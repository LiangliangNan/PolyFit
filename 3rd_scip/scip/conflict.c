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

/**@file   conflict.c
 * @brief  methods and datastructures for conflict analysis
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Jakob Witzig
 *
 * This file implements a conflict analysis method like the one used in modern
 * SAT solvers like zchaff. The algorithm works as follows:
 *
 * Given is a set of bound changes that are not allowed being applied simultaneously, because they
 * render the current node infeasible (e.g. because a single constraint is infeasible in the these
 * bounds, or because the LP relaxation is infeasible).  The goal is to deduce a clause on variables
 * -- a conflict clause -- representing the "reason" for this conflict, i.e., the branching decisions
 * or the deductions (applied e.g. in domain propagation) that lead to the conflict. This clause can
 * then be added to the constraint set to help cutting off similar parts of the branch and bound
 * tree, that would lead to the same conflict.  A conflict clause can also be generated, if the
 * conflict was detected by a locally valid constraint. In this case, the resulting conflict clause
 * is also locally valid in the same depth as the conflict detecting constraint. If all involved
 * variables are binary, a linear (set covering) constraint can be generated, otherwise a bound
 * disjunction constraint is generated. Details are given in
 *
 * Tobias Achterberg, Conflict Analysis in Mixed Integer Programming@n
 * Discrete Optimization, 4, 4-20 (2007)
 *
 * See also @ref CONF. Here is an outline of the algorithm:
 *
 * -#  Put all the given bound changes to a priority queue, which is ordered,
 *     such that the bound change that was applied last due to branching or deduction
 *     is at the top of the queue. The variables in the queue are always active
 *     problem variables. Because binary variables are preferred over general integer
 *     variables, integer variables are put on the priority queue prior to the binary
 *     variables. Create an empty conflict set.
 * -#  Remove the top bound change b from the priority queue.
 * -#  Perform the following case distinction:
 *     -#  If the remaining queue is non-empty, and bound change b' (the one that is now
 *         on the top of the queue) was applied at the same depth level as b, and if
 *         b was a deduction with known inference reason, and if the inference constraint's
 *         valid depth is smaller or equal to the conflict detecting constraint's valid
 *         depth:
 *          - Resolve bound change b by asking the constraint that inferred the
 *            bound change to put all the bound changes on the priority queue, that
 *            lead to the deduction of b.
 *            Note that these bound changes have at most the same inference depth
 *            level as b, and were deduced earlier than b.
 *     -#  Otherwise, the bound change b was a branching decision or a deduction with
 *         missing inference reason, or the inference constraint's validity is more local
 *         than the one of the conflict detecting constraint.
 *          - If a the bound changed corresponds to a binary variable, add it or its
 *            negation to the conflict set, depending on which of them is currently fixed to
 *            FALSE (i.e., the conflict set consists of literals that cannot be FALSE
 *            altogether at the same time).
 *          - Otherwise put the bound change into the conflict set.
 *         Note that if the bound change was a branching, all deduced bound changes
 *         remaining in the priority queue have smaller inference depth level than b,
 *         since deductions are always applied after the branching decisions. However,
 *         there is the possibility, that b was a deduction, where the inference
 *         reason was not given or the inference constraint was too local.
 *         With this lack of information, we must treat the deduced bound change like
 *         a branching, and there may exist other deduced bound changes of the same
 *         inference depth level in the priority queue.
 * -#  If priority queue is non-empty, goto step 2.
 * -#  The conflict set represents the conflict clause saying that at least one
 *     of the conflict variables must take a different value. The conflict set is then passed
 *     to the conflict handlers, that may create a corresponding constraint (e.g. a logicor
 *     constraint or bound disjunction constraint) out of these conflict variables and
 *     add it to the problem.
 *
 * If all deduced bound changes come with (global) inference information, depending on
 * the conflict analyzing strategy, the resulting conflict set has the following property:
 *  - 1-FirstUIP: In the depth level where the conflict was found, at most one variable
 *    assigned at that level is member of the conflict set. This conflict variable is the
 *    first unique implication point of its depth level (FUIP).
 *  - All-FirstUIP: For each depth level, at most one variable assigned at that level is
 *    member of the conflict set. This conflict variable is the first unique implication
 *    point of its depth level (FUIP).
 *
 * The user has to do the following to get the conflict analysis running in its
 * current implementation:
 *  - A constraint handler or propagator supporting the conflict analysis must implement
 *    the CONSRESPROP/PROPRESPROP call, that processes a bound change inference b and puts all
 *    the reason bounds leading to the application of b with calls to
 *    SCIPaddConflictBound() on the conflict queue (algorithm step 3.(a)).
 *  - If the current bounds lead to a deduction of a bound change (e.g. in domain
 *    propagation), a constraint handler should call SCIPinferVarLbCons() or
 *    SCIPinferVarUbCons(), thus providing the constraint that infered the bound change.
 *    A propagator should call SCIPinferVarLbProp() or SCIPinferVarUbProp() instead,
 *    thus providing a pointer to itself.
 *  - If (in the current bounds) an infeasibility is detected, the constraint handler or
 *    propagator should
 *     1. call SCIPinitConflictAnalysis() to initialize the conflict queue,
 *     2. call SCIPaddConflictBound() for each bound that lead to the conflict,
 *     3. call SCIPanalyzeConflictCons() or SCIPanalyzeConflict() to analyze the conflict
 *        and add an appropriate conflict constraint.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/visual.h"
#include "scip/history.h"
#include "scip/paramset.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/tree.h"
#include "scip/sol.h"
#include "scip/scip.h"
#include "scip/conflict.h"
#include "scip/cons.h"
#include "scip/prop.h"
#include "scip/presolve.h"
#include "scip/debug.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/cuts.h"
#include "lpi/lpi.h"

#include "scip/struct_conflict.h"
#include "scip/cons_linear.h"

#define BOUNDSWITCH                0.51 /**< threshold for bound switching - see cuts.c */
#define POSTPROCESS               FALSE /**< apply postprocessing to the cut - see cuts.c */
#define USEVBDS                   FALSE /**< use variable bounds - see cuts.c */
#define ALLOWLOCAL                FALSE /**< allow to generate local cuts - see cuts. */
#define MINFRAC                   0.05  /**< minimal fractionality of floor(rhs) - see cuts.c */
#define MAXFRAC                   0.999 /**< maximal fractionality of floor(rhs) - see cuts.c */

/*#define SCIP_CONFGRAPH*/


#ifdef SCIP_CONFGRAPH
/*
 * Output of Conflict Graph
 */

#include <stdio.h>

static FILE*             confgraphfile = NULL;              /**< output file for current conflict graph */
static SCIP_BDCHGINFO*   confgraphcurrentbdchginfo = NULL;  /**< currently resolved bound change */
static int               confgraphnconflictsets = 0;        /**< number of conflict sets marked in the graph */

/** writes a node section to the conflict graph file */
static
void confgraphWriteNode(
   void*                 idptr,              /**< id of the node */
   const char*           label,              /**< label of the node */
   const char*           nodetype,           /**< type of the node */
   const char*           fillcolor,          /**< color of the node's interior */
   const char*           bordercolor         /**< color of the node's border */
   )
{
   assert(confgraphfile != NULL);

   SCIPgmlWriteNode(confgraphfile, (unsigned int)(size_t)idptr, label, nodetype, fillcolor, bordercolor);
}

/** writes an edge section to the conflict graph file */
static
void confgraphWriteEdge(
   void*                 source,             /**< source node of the edge */
   void*                 target,             /**< target node of the edge */
   const char*           color               /**< color of the edge */
   )
{
   assert(confgraphfile != NULL);

#ifndef SCIP_CONFGRAPH_EDGE
   SCIPgmlWriteArc(confgraphfile, (unsigned int)(size_t)source, (unsigned int)(size_t)target, NULL, color);
#else
   SCIPgmlWriteEdge(confgraphfile, (unsigned int)(size_t)source, (unsigned int)(size_t)target, NULL, color);
#endif
}

/** creates a file to output the current conflict graph into; adds the conflict vertex to the graph */
static
SCIP_RETCODE confgraphCreate(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   char fname[SCIP_MAXSTRLEN];

   assert(conflict != NULL);
   assert(confgraphfile == NULL);

   (void) SCIPsnprintf(fname, SCIP_MAXSTRLEN, "conf%p%d.gml", conflict, conflict->count);
   SCIPinfoMessage(set->scip, NULL, "storing conflict graph in file <%s>\n", fname);

   confgraphfile = fopen(fname, "w");

   if( confgraphfile == NULL )
   {
      SCIPerrorMessage("cannot open graph file <%s>\n", fname);
      SCIPABORT(); /*lint !e527*/
      return SCIP_WRITEERROR;
   }

   SCIPgmlWriteOpening(confgraphfile, TRUE);

   confgraphWriteNode(NULL, "conflict", "ellipse", "#ff0000", "#000000");

   confgraphcurrentbdchginfo = NULL;

   return SCIP_OKAY;
}

/** closes conflict graph file */
static
void confgraphFree(
   void
   )
{
   if( confgraphfile != NULL )
   {
      SCIPgmlWriteClosing(confgraphfile);

      fclose(confgraphfile);

      confgraphfile = NULL;
      confgraphnconflictsets = 0;
   }
}

/** adds a bound change node to the conflict graph and links it to the currently resolved bound change */
static
void confgraphAddBdchg(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to add to the conflict graph */
   )
{
   const char* colors[] = {
      "#8888ff", /* blue for constraint resolving */
      "#ffff00", /* yellow for propagator resolving */
      "#55ff55"  /* green branching decision */
   };
   char label[SCIP_MAXSTRLEN];
   char depth[SCIP_MAXSTRLEN];
   int col;


   switch( SCIPbdchginfoGetChgtype(bdchginfo) )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
      col = 2;
      break;
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      col = 0;
      break;
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      col = (SCIPbdchginfoGetInferProp(bdchginfo) == NULL ? 1 : 0);
      break;
   default:
      SCIPerrorMessage("invalid bound change type\n");
      col = 0;
      SCIPABORT();
      break;
   }

   if( SCIPbdchginfoGetDepth(bdchginfo) == INT_MAX )
      (void) SCIPsnprintf(depth, SCIP_MAXSTRLEN, "dive");
   else
      (void) SCIPsnprintf(depth, SCIP_MAXSTRLEN, "%d", SCIPbdchginfoGetDepth(bdchginfo));
   (void) SCIPsnprintf(label, SCIP_MAXSTRLEN, "%s %s %g\n[%s:%d]", SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
      SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
      SCIPbdchginfoGetNewbound(bdchginfo), depth, SCIPbdchginfoGetPos(bdchginfo));
   confgraphWriteNode(bdchginfo, label, "ellipse", colors[col], "#000000");
   confgraphWriteEdge(bdchginfo, confgraphcurrentbdchginfo, "#000000");
}

/** links the already existing bound change node to the currently resolved bound change */
static
void confgraphLinkBdchg(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to add to the conflict graph */
   )
{
   confgraphWriteEdge(bdchginfo, confgraphcurrentbdchginfo, "#000000");
}

/** marks the given bound change to be the currently resolved bound change */
static
void confgraphSetCurrentBdchg(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to add to the conflict graph */
   )
{
   confgraphcurrentbdchginfo = bdchginfo;
}

/** marks given conflict set in the conflict graph */
static
void confgraphMarkConflictset(
   SCIP_CONFLICTSET*     conflictset         /**< conflict set */
   )
{
   char label[SCIP_MAXSTRLEN];
   int i;

   assert(conflictset != NULL);

   confgraphnconflictsets++;
   (void) SCIPsnprintf(label, SCIP_MAXSTRLEN, "conf %d (%d)", confgraphnconflictsets, conflictset->validdepth);
   confgraphWriteNode((void*)(size_t)confgraphnconflictsets, label, "rectangle", "#ff00ff", "#000000");
   for( i = 0; i < conflictset->nbdchginfos; ++i )
      confgraphWriteEdge((void*)(size_t)confgraphnconflictsets, conflictset->bdchginfos[i], "#ff00ff");
}

#endif

/*
 * Conflict Handler
 */

/** compares two conflict handlers w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPconflicthdlrComp)
{  /*lint --e{715}*/
   return ((SCIP_CONFLICTHDLR*)elem2)->priority - ((SCIP_CONFLICTHDLR*)elem1)->priority;
}

/** comparison method for sorting conflict handler w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPconflicthdlrCompName)
{
   return strcmp(SCIPconflicthdlrGetName((SCIP_CONFLICTHDLR*)elem1), SCIPconflicthdlrGetName((SCIP_CONFLICTHDLR*)elem2));
}

/** method to call, when the priority of a conflict handler was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdConflicthdlrPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetConflicthdlrPriority() to mark the conflicthdlrs unsorted */
   SCIP_CALL( SCIPsetConflicthdlrPriority(scip, (SCIP_CONFLICTHDLR*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given conflict handler to a new scip */
SCIP_RETCODE SCIPconflicthdlrCopyInclude(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( conflicthdlr->conflictcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including conflict handler %s in subscip %p\n", SCIPconflicthdlrGetName(conflicthdlr), (void*)set->scip);
      SCIP_CALL( conflicthdlr->conflictcopy(set->scip, conflicthdlr) );
   }

   return SCIP_OKAY;
}

/** creates a conflict handler */
SCIP_RETCODE SCIPconflicthdlrCreate(
   SCIP_CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of conflict handler */
   const char*           desc,               /**< description of conflict handler */
   int                   priority,           /**< priority of the conflict handler */
   SCIP_DECL_CONFLICTCOPY((*conflictcopy)),  /**< copy method of conflict handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONFLICTFREE((*conflictfree)),  /**< destructor of conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit)),  /**< initialize conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit)),  /**< deinitialize conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol)),/**< solving process initialization method of conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol)),/**< solving process deinitialization method of conflict handler */
   SCIP_DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(conflicthdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   SCIP_ALLOC( BMSallocMemory(conflicthdlr) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*conflicthdlr)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*conflicthdlr)->desc, desc, strlen(desc)+1) );
   (*conflicthdlr)->priority = priority;
   (*conflicthdlr)->conflictcopy = conflictcopy;
   (*conflicthdlr)->conflictfree = conflictfree;
   (*conflicthdlr)->conflictinit = conflictinit;
   (*conflicthdlr)->conflictexit = conflictexit;
   (*conflicthdlr)->conflictinitsol = conflictinitsol;
   (*conflicthdlr)->conflictexitsol = conflictexitsol;
   (*conflicthdlr)->conflictexec = conflictexec;
   (*conflicthdlr)->conflicthdlrdata = conflicthdlrdata;
   (*conflicthdlr)->initialized = FALSE;

   SCIP_CALL( SCIPclockCreate(&(*conflicthdlr)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflicthdlr)->conflicttime, SCIP_CLOCKTYPE_DEFAULT) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "conflict/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of conflict handler <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc, &(*conflicthdlr)->priority, TRUE, \
         priority, INT_MIN, INT_MAX, paramChgdConflicthdlrPriority, (SCIP_PARAMDATA*)(*conflicthdlr)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of conflict handler */
SCIP_RETCODE SCIPconflicthdlrFree(
   SCIP_CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(*conflicthdlr != NULL);
   assert(!(*conflicthdlr)->initialized);
   assert(set != NULL);

   /* call destructor of conflict handler */
   if( (*conflicthdlr)->conflictfree != NULL )
   {
      SCIP_CALL( (*conflicthdlr)->conflictfree(set->scip, *conflicthdlr) );
   }

   SCIPclockFree(&(*conflicthdlr)->conflicttime);
   SCIPclockFree(&(*conflicthdlr)->setuptime);

   BMSfreeMemoryArray(&(*conflicthdlr)->name);
   BMSfreeMemoryArray(&(*conflicthdlr)->desc);
   BMSfreeMemory(conflicthdlr);

   return SCIP_OKAY;
}

/** calls initialization method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrInit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   if( conflicthdlr->initialized )
   {
      SCIPerrorMessage("conflict handler <%s> already initialized\n", conflicthdlr->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(conflicthdlr->setuptime);
      SCIPclockReset(conflicthdlr->conflicttime);
   }

   /* call initialization method of conflict handler */
   if( conflicthdlr->conflictinit != NULL )
   {
      /* start timing */
      SCIPclockStart(conflicthdlr->setuptime, set);

      SCIP_CALL( conflicthdlr->conflictinit(set->scip, conflicthdlr) );

      /* stop timing */
      SCIPclockStop(conflicthdlr->setuptime, set);
   }
   conflicthdlr->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrExit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   if( !conflicthdlr->initialized )
   {
      SCIPerrorMessage("conflict handler <%s> not initialized\n", conflicthdlr->name);
      return SCIP_INVALIDCALL;
   }

   /* call deinitialization method of conflict handler */
   if( conflicthdlr->conflictexit != NULL )
   {
      /* start timing */
      SCIPclockStart(conflicthdlr->setuptime, set);

      SCIP_CALL( conflicthdlr->conflictexit(set->scip, conflicthdlr) );

      /* stop timing */
      SCIPclockStop(conflicthdlr->setuptime, set);
   }
   conflicthdlr->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs conflict handler that the branch and bound process is being started */
SCIP_RETCODE SCIPconflicthdlrInitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   /* call solving process initialization method of conflict handler */
   if( conflicthdlr->conflictinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(conflicthdlr->setuptime, set);

      SCIP_CALL( conflicthdlr->conflictinitsol(set->scip, conflicthdlr) );

      /* stop timing */
      SCIPclockStop(conflicthdlr->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs conflict handler that the branch and bound process data is being freed */
SCIP_RETCODE SCIPconflicthdlrExitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of conflict handler */
   if( conflicthdlr->conflictexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(conflicthdlr->setuptime, set);

      SCIP_CALL( conflicthdlr->conflictexitsol(set->scip, conflicthdlr) );

      /* stop timing */
      SCIPclockStop(conflicthdlr->setuptime, set);
   }

   return SCIP_OKAY;
}

/** calls execution method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrExec(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node to add conflict constraint to */
   SCIP_NODE*            validnode,          /**< node at which the constraint is valid */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound change resembling the conflict set */
   SCIP_Real*            relaxedbds,         /**< array with relaxed bounds which are efficient to create a valid conflict */
   int                   nbdchginfos,        /**< number of bound changes in the conflict set */
   SCIP_CONFTYPE         conftype,           /**< type of the conflict */
   SCIP_Bool             usescutoffbound,    /**< depends the conflict on the cutoff bound? */
   SCIP_Bool             resolved,           /**< was the conflict set already used to create a constraint? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{

   assert(conflicthdlr != NULL);
   assert(set != NULL);
   assert(bdchginfos != NULL || nbdchginfos == 0);
   assert(result != NULL);

   /* call solution start method of conflict handler */
   *result = SCIP_DIDNOTRUN;
   if( conflicthdlr->conflictexec != NULL )
   {
      /* start timing */
      SCIPclockStart(conflicthdlr->conflicttime, set);

      SCIP_CALL( conflicthdlr->conflictexec(set->scip, conflicthdlr, node, validnode, bdchginfos, relaxedbds, nbdchginfos,
            conftype, usescutoffbound, set->conf_separate, (SCIPnodeGetDepth(validnode) > 0), set->conf_dynamic,
            set->conf_removable, resolved, result) );

      /* stop timing */
      SCIPclockStop(conflicthdlr->conflicttime, set);

      if( *result != SCIP_CONSADDED
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         SCIPerrorMessage("execution method of conflict handler <%s> returned invalid result <%d>\n",
            conflicthdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** gets user data of conflict handler */
SCIP_CONFLICTHDLRDATA* SCIPconflicthdlrGetData(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->conflicthdlrdata;
}

/** sets user data of conflict handler; user has to free old data in advance! */
void SCIPconflicthdlrSetData(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< new conflict handler user data */
   )
{
   assert(conflicthdlr != NULL);

   conflicthdlr->conflicthdlrdata = conflicthdlrdata;
}

/** set copy method of conflict handler */
void SCIPconflicthdlrSetCopy(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTCOPY((*conflictcopy))   /**< copy method of the conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   conflicthdlr->conflictcopy = conflictcopy;
}

/** set destructor of conflict handler */
void SCIPconflicthdlrSetFree(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTFREE((*conflictfree))   /**< destructor of conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   conflicthdlr->conflictfree = conflictfree;
}

/** set initialization method of conflict handler */
void SCIPconflicthdlrSetInit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit))   /**< initialization method conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   conflicthdlr->conflictinit = conflictinit;
}

/** set deinitialization method of conflict handler */
void SCIPconflicthdlrSetExit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit))   /**< deinitialization method conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   conflicthdlr->conflictexit = conflictexit;
}

/** set solving process initialization method of conflict handler */
void SCIPconflicthdlrSetInitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol))/**< solving process initialization method of conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   conflicthdlr->conflictinitsol = conflictinitsol;
}

/** set solving process deinitialization method of conflict handler */
void SCIPconflicthdlrSetExitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol))/**< solving process deinitialization method of conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   conflicthdlr->conflictexitsol = conflictexitsol;
}

/** gets name of conflict handler */
const char* SCIPconflicthdlrGetName(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->name;
}

/** gets description of conflict handler */
const char* SCIPconflicthdlrGetDesc(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->desc;
}

/** gets priority of conflict handler */
int SCIPconflicthdlrGetPriority(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->priority;
}

/** sets priority of conflict handler */
void SCIPconflicthdlrSetPriority(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the conflict handler */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   conflicthdlr->priority = priority;
   set->conflicthdlrssorted = FALSE;
}

/** is conflict handler initialized? */
SCIP_Bool SCIPconflicthdlrIsInitialized(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->initialized;
}

/** enables or disables all clocks of \p conflicthdlr, depending on the value of the flag */
void SCIPconflicthdlrEnableOrDisableClocks(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< the conflict handler for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the conflict handler be enabled? */
   )
{
   assert(conflicthdlr != NULL);

   SCIPclockEnableOrDisable(conflicthdlr->setuptime, enable);
   SCIPclockEnableOrDisable(conflicthdlr->conflicttime, enable);
}

/** gets time in seconds used in this conflict handler for setting up for next stages */
SCIP_Real SCIPconflicthdlrGetSetupTime(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return SCIPclockGetTime(conflicthdlr->setuptime);
}

/** gets time in seconds used in this conflict handler */
SCIP_Real SCIPconflicthdlrGetTime(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return SCIPclockGetTime(conflicthdlr->conflicttime);
}

/*
 * Conflict LP Bound Changes
 */


/** create conflict LP bound change data structure */
static
SCIP_RETCODE lpbdchgsCreate(
   SCIP_LPBDCHGS**       lpbdchgs,           /**< pointer to store the conflict LP bound change data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   ncols               /**< number of columns */
   )
{
   SCIP_CALL( SCIPsetAllocBuffer(set, lpbdchgs) );

   SCIP_CALL( SCIPsetAllocBufferArray(set, &(*lpbdchgs)->bdchginds, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &(*lpbdchgs)->bdchglbs, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &(*lpbdchgs)->bdchgubs, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &(*lpbdchgs)->bdchgcolinds, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &(*lpbdchgs)->usedcols, ncols) );
   BMSclearMemoryArray((*lpbdchgs)->usedcols, ncols);

   (*lpbdchgs)->nbdchgs = 0;

   return SCIP_OKAY;
}

/** reset conflict LP bound change data structure */
static
void lpbdchgsReset(
   SCIP_LPBDCHGS*        lpbdchgs,           /**< conflict LP bound change data structure */
   int                   ncols               /**< number of columns */
   )
{
   assert(lpbdchgs != NULL);

   BMSclearMemoryArray(lpbdchgs->usedcols, ncols);
   lpbdchgs->nbdchgs = 0;
}

/** free conflict LP bound change data structure */
static
void lpbdchgsFree(
   SCIP_LPBDCHGS**       lpbdchgs,           /**< pointer to store the conflict LP bound change data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIPsetFreeBufferArray(set, &(*lpbdchgs)->usedcols);
   SCIPsetFreeBufferArray(set, &(*lpbdchgs)->bdchgcolinds);
   SCIPsetFreeBufferArray(set, &(*lpbdchgs)->bdchgubs);
   SCIPsetFreeBufferArray(set, &(*lpbdchgs)->bdchglbs);
   SCIPsetFreeBufferArray(set, &(*lpbdchgs)->bdchginds);

   SCIPsetFreeBuffer(set, lpbdchgs);
}

/*
 * Proof Sets
 */

/** resets the data structure of a proofset */
static
void proofsetClear(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   proofset->nnz = 0;
   proofset->rhs = 0.0;
   proofset->conflicttype = SCIP_CONFTYPE_UNKNOWN;
}

/** creates a proofset */
static
SCIP_RETCODE proofsetCreate(
   SCIP_PROOFSET**       proofset,           /**< proof set */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(proofset != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, proofset) );
   (*proofset)->vals = NULL;
   (*proofset)->inds = NULL;
   (*proofset)->rhs = 0.0;
   (*proofset)->nnz = 0;
   (*proofset)->size = 0;
   (*proofset)->conflicttype = SCIP_CONFTYPE_UNKNOWN;

   return SCIP_OKAY;
}

/** creates and clears the proofset */
static
SCIP_RETCODE conflictInitProofset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(conflict != NULL);
   assert(blkmem != NULL);

   SCIP_CALL( proofsetCreate(&conflict->proofset, blkmem) );

   return SCIP_OKAY;
}

/** frees a proofset */
static
void proofsetFree(
   SCIP_PROOFSET**       proofset,           /**< proof set */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(proofset != NULL);
   assert(*proofset != NULL);
   assert(blkmem != NULL);

   BMSfreeBlockMemoryArrayNull(blkmem, &(*proofset)->vals, (*proofset)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*proofset)->inds, (*proofset)->size);
   BMSfreeBlockMemory(blkmem, proofset);
   (*proofset) = NULL;
}

#ifdef SCIP_DEBUG
static
void proofsetPrint(
   SCIP_PROOFSET*        proofset,
   SCIP_SET*             set,
   SCIP_PROB*            transprob
   )
{
   SCIP_VAR** vars;
   int i;

   assert(proofset != NULL);

   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   printf("proofset: ");
   for( i = 0; i < proofset->nnz; i++ )
      printf("%+.15g <%s> ", proofset->vals[i], SCIPvarGetName(vars[proofset->inds[i]]));
   printf(" <= %.15g\n", proofset->rhs);
}
#endif

/** return the indices of variables in the proofset */
static
int* proofsetGetInds(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   return proofset->inds;
}

/** return coefficient of variable in the proofset with given probindex */
static
SCIP_Real* proofsetGetVals(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   return proofset->vals;
}

/** return the right-hand side if a proofset */
static
SCIP_Real proofsetGetRhs(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   return proofset->rhs;
}

/** returns the number of variables in the proofset */
static
int proofsetGetNVars(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   return proofset->nnz;
}

/** returns the number of variables in the proofset */
static
SCIP_CONFTYPE proofsetGetConftype(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   return proofset->conflicttype;
}

/** adds given data as aggregation row to the proofset */
static
SCIP_RETCODE proofsetAddSparseData(
   SCIP_PROOFSET*        proofset,           /**< proof set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Real*            vals,               /**< variable coefficients */
   int*                  inds,               /**< variable array */
   int                   nnz,                /**< size of variable and coefficient array */
   SCIP_Real             rhs                 /**< right-hand side of the aggregation row */
   )
{
   assert(proofset != NULL);
   assert(blkmem != NULL);

   if( proofset->size == 0 )
   {
      assert(proofset->vals == NULL);
      assert(proofset->inds == NULL);

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &proofset->vals, vals, nnz) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &proofset->inds, inds, nnz) );

      proofset->size = nnz;
   }
   else
   {
      int i;

      assert(proofset->vals != NULL);
      assert(proofset->inds != NULL);

      if( proofset->size < nnz )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &proofset->vals, proofset->size, nnz) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &proofset->inds, proofset->size, nnz) );
         proofset->size = nnz;
      }

      for( i = 0; i < nnz; i++ )
      {
         proofset->vals[i] = vals[i];
         proofset->inds[i] = inds[i];
      }
   }

   proofset->rhs = rhs;
   proofset->nnz = nnz;

   return SCIP_OKAY;
}

/** adds an aggregation row to the proofset */
static
SCIP_RETCODE proofsetAddAggrrow(
   SCIP_PROOFSET*        proofset,           /**< proof set */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_AGGRROW*         aggrrow             /**< aggregation row to add */
   )
{
   SCIP_Real* vals;
   int* inds;
   int nnz;
   int i;

   assert(proofset != NULL);
   assert(set != NULL);

   inds = SCIPaggrRowGetInds(aggrrow);
   assert(inds != NULL);

   nnz = SCIPaggrRowGetNNz(aggrrow);
   assert(nnz > 0);

   SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, nnz) );

   for( i = 0; i < nnz; i++ )
   {
      vals[i] = SCIPaggrRowGetProbvarValue(aggrrow, inds[i]);
   }

   SCIP_CALL( proofsetAddSparseData(proofset, blkmem, vals, inds, nnz, SCIPaggrRowGetRhs(aggrrow)) );

   SCIPsetFreeBufferArray(set, &vals);

   return SCIP_OKAY;
}

/** Removes a given variable @p var from position @p pos from the proofset and updates the right-hand side according
 *  to sign of the coefficient, i.e., rhs -= coef * bound, where bound = lb if coef >= 0 and bound = ub, otherwise.
 *
 *  @note: The list of non-zero indices and coefficients will be updated by swapping the last non-zero index to @p pos.
 */
static
void proofsetCancelVarWithBound(
   SCIP_PROOFSET*        proofset,
   SCIP_SET*             set,
   SCIP_VAR*             var,
   int                   pos,
   SCIP_Bool*            valid
   )
{
   assert(proofset != NULL);
   assert(var != NULL);
   assert(pos >= 0 && pos < proofset->nnz);
   assert(valid != NULL);

   *valid = TRUE;

   /* cancel with lower bound */
   if( proofset->vals[pos] > 0.0 )
   {
      proofset->rhs -= proofset->vals[pos] * SCIPvarGetLbGlobal(var);
   }
   /* cancel with upper bound */
   else
   {
      assert(proofset->vals[pos] < 0.0);
      proofset->rhs -= proofset->vals[pos] * SCIPvarGetUbGlobal(var);
   }

   --proofset->nnz;

   proofset->vals[pos] = proofset->vals[proofset->nnz];
   proofset->inds[pos] = proofset->inds[proofset->nnz];
   proofset->vals[proofset->nnz] = 0.0;
   proofset->inds[proofset->nnz] = 0;

   if( SCIPsetIsInfinity(set, proofset->rhs) )
      *valid = FALSE;
}

/*
 * Conflict Sets
 */

/** resizes the array of the temporary bound change informations to be able to store at least num bound change entries */
static
SCIP_RETCODE conflictEnsureTmpbdchginfosMem(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in arrays */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);

   if( num > conflict->tmpbdchginfossize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conflict->tmpbdchginfos, newsize) );
      conflict->tmpbdchginfossize = newsize;
   }
   assert(num <= conflict->tmpbdchginfossize);

   return SCIP_OKAY;
}

/** creates a temporary bound change information object that is destroyed after the conflict sets are flushed */
static
SCIP_RETCODE conflictCreateTmpBdchginfo(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< active variable that changed the bounds */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound for var: lower or upper bound */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BDCHGINFO**      bdchginfo           /**< pointer to store bound change information */
   )
{
   assert(conflict != NULL);

   SCIP_CALL( conflictEnsureTmpbdchginfosMem(conflict, set, conflict->ntmpbdchginfos+1) );
   SCIP_CALL( SCIPbdchginfoCreate(&conflict->tmpbdchginfos[conflict->ntmpbdchginfos], blkmem,
         var, boundtype, oldbound, newbound) );
   *bdchginfo = conflict->tmpbdchginfos[conflict->ntmpbdchginfos];
   conflict->ntmpbdchginfos++;

   return SCIP_OKAY;
}

/** frees all temporarily created bound change information data */
static
void conflictFreeTmpBdchginfos(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int i;

   assert(conflict != NULL);

   for( i = 0; i < conflict->ntmpbdchginfos; ++i )
      SCIPbdchginfoFree(&conflict->tmpbdchginfos[i], blkmem);
   conflict->ntmpbdchginfos = 0;
}

/** clears the given conflict set */
static
void conflictsetClear(
   SCIP_CONFLICTSET*     conflictset         /**< conflict set */
   )
{
   assert(conflictset != NULL);

   conflictset->nbdchginfos = 0;
   conflictset->validdepth = 0;
   conflictset->insertdepth = 0;
   conflictset->conflictdepth = 0;
   conflictset->repropdepth = 0;
   conflictset->repropagate = TRUE;
   conflictset->usescutoffbound = FALSE;
   conflictset->conflicttype = SCIP_CONFTYPE_UNKNOWN;
}

/** creates an empty conflict set */
static
SCIP_RETCODE conflictsetCreate(
   SCIP_CONFLICTSET**    conflictset,        /**< pointer to store the conflict set */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(conflictset != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, conflictset) );
   (*conflictset)->bdchginfos = NULL;
   (*conflictset)->relaxedbds = NULL;
   (*conflictset)->sortvals = NULL;
   (*conflictset)->bdchginfossize = 0;

   conflictsetClear(*conflictset);

   return SCIP_OKAY;
}

/** creates a copy of the given conflict set, allocating an additional amount of memory */
static
SCIP_RETCODE conflictsetCopy(
   SCIP_CONFLICTSET**    targetconflictset,  /**< pointer to store the conflict set */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_CONFLICTSET*     sourceconflictset,  /**< source conflict set */
   int                   nadditionalelems    /**< number of additional elements to allocate memory for */
   )
{
   int targetsize;

   assert(targetconflictset != NULL);
   assert(sourceconflictset != NULL);

   targetsize = sourceconflictset->nbdchginfos + nadditionalelems;
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, targetconflictset) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*targetconflictset)->bdchginfos, targetsize) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*targetconflictset)->relaxedbds, targetsize) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*targetconflictset)->sortvals, targetsize) );
   (*targetconflictset)->bdchginfossize = targetsize;

   BMScopyMemoryArray((*targetconflictset)->bdchginfos, sourceconflictset->bdchginfos, sourceconflictset->nbdchginfos);
   BMScopyMemoryArray((*targetconflictset)->relaxedbds, sourceconflictset->relaxedbds, sourceconflictset->nbdchginfos);
   BMScopyMemoryArray((*targetconflictset)->sortvals, sourceconflictset->sortvals, sourceconflictset->nbdchginfos);

   (*targetconflictset)->nbdchginfos = sourceconflictset->nbdchginfos;
   (*targetconflictset)->validdepth = sourceconflictset->validdepth;
   (*targetconflictset)->insertdepth = sourceconflictset->insertdepth;
   (*targetconflictset)->conflictdepth = sourceconflictset->conflictdepth;
   (*targetconflictset)->repropdepth = sourceconflictset->repropdepth;
   (*targetconflictset)->usescutoffbound = sourceconflictset->usescutoffbound;
   (*targetconflictset)->conflicttype = sourceconflictset->conflicttype;

   return SCIP_OKAY;
}

/** frees a conflict set */
static
void conflictsetFree(
   SCIP_CONFLICTSET**    conflictset,        /**< pointer to the conflict set */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(conflictset != NULL);
   assert(*conflictset != NULL);

   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictset)->bdchginfos, (*conflictset)->bdchginfossize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictset)->relaxedbds, (*conflictset)->bdchginfossize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictset)->sortvals, (*conflictset)->bdchginfossize);
   BMSfreeBlockMemory(blkmem, conflictset);
}

/** resizes the arrays of the conflict set to be able to store at least num bound change entries */
static
SCIP_RETCODE conflictsetEnsureBdchginfosMem(
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in arrays */
   )
{
   assert(conflictset != NULL);
   assert(set != NULL);

   if( num > conflictset->bdchginfossize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictset->bdchginfos, conflictset->bdchginfossize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictset->relaxedbds, conflictset->bdchginfossize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictset->sortvals, conflictset->bdchginfossize, newsize) );
      conflictset->bdchginfossize = newsize;
   }
   assert(num <= conflictset->bdchginfossize);

   return SCIP_OKAY;
}

/** calculates the score of the conflict set
 *
 *  the score is weighted sum of number of bound changes, repropagation depth, and valid depth
 */
static
SCIP_Real conflictsetCalcScore(
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflictset != NULL);

   return -(set->conf_weightsize * conflictset->nbdchginfos
         + set->conf_weightrepropdepth * conflictset->repropdepth
         + set->conf_weightvaliddepth * conflictset->validdepth);
}

/** calculates the score of a bound change within a conflict */
static
SCIP_Real calcBdchgScore(
   SCIP_Real             prooflhs,           /**< lhs of proof constraint */
   SCIP_Real             proofact,           /**< activity of the proof constraint */
   SCIP_Real             proofactdelta,      /**< activity change */
   SCIP_Real             proofcoef,          /**< coefficient in proof constraint */
   int                   depth,              /**< bound change depth */
   int                   currentdepth,       /**< current depth */
   SCIP_VAR*             var,                /**< variable corresponding to bound change */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_COL* col;
   SCIP_Real score;

   score = set->conf_proofscorefac * (1.0 - proofactdelta/(prooflhs - proofact));
   score = MAX(score, 0.0);
   score += set->conf_depthscorefac * (SCIP_Real)(depth+1)/(SCIP_Real)(currentdepth+1);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      col = SCIPvarGetCol(var);
   else
      col = NULL;

   if( proofcoef > 0.0 )
   {
      if( col != NULL && SCIPcolGetNNonz(col) > 0 )
         score += set->conf_uplockscorefac * (SCIP_Real)(SCIPvarGetNLocksUp(var))/(SCIP_Real)(SCIPcolGetNNonz(col));
      else
         score += set->conf_uplockscorefac * SCIPvarGetNLocksUp(var);
   }
   else
   {
      if( col != NULL && SCIPcolGetNNonz(col) > 0 )
         score += set->conf_downlockscorefac * (SCIP_Real)(SCIPvarGetNLocksDown(var))/(SCIP_Real)(SCIPcolGetNNonz(col));
      else
         score += set->conf_downlockscorefac * SCIPvarGetNLocksDown(var);
   }

   return score;
}

/** check if the bound change info (which is the potential next candidate which is queued) is valid for the current
 *  conflict analysis; a bound change info can get invalid if after this one was added to the queue, a weaker bound
 *  change was added to the queue (due the bound widening idea) which immediately makes this bound change redundant; due
 *  to the priority we did not removed that bound change info since that cost O(log(n)); hence we have to skip/ignore it
 *  now
 *
 *  The following situations can occur before for example the bound change info (x >= 3) is potentially popped from the
 *  queue.
 *
 *  Postcondition: the reason why (x >= 3) was queued is that at this time point no lower bound of x was involved yet in
 *                 the current conflict or the lower bound which was involved until then was stronger, e.g., (x >= 2).
 *
 *  1) during the time until (x >= 3) gets potentially popped no weaker lower bound was added to the queue, in that case
 *     the conflictlbcount is valid and conflictlb is 3; that is (var->conflictlbcount == conflict->count &&
 *     var->conflictlb == 3)
 *
 *  2) a weaker bound change info gets queued (e.g., x >= 4); this bound change is popped before (x >= 3) since it has
 *     higher priority (which is the time stamp of the bound change info and (x >= 4) has to be done after (x >= 3)
 *     during propagation or branching)
 *
 *    a) if (x >= 4) is popped and added to the conflict set the conflictlbcount is still valid and conflictlb is at
 *      most 4; that is (var->conflictlbcount == conflict->count && var->conflictlb >= 4); it follows that any bound
 *      change info which is stronger than (x >= 4) gets ignored (for example x >= 2)
 *
 *    b) if (x >= 4) is popped and resolved without introducing a new lower bound on x until (x >= 3) is a potentially
 *       candidate the conflictlbcount indicates that bound change is currently not present; that is
 *       (var->conflictlbcount != conflict->count)
 *
 *    c) if (x >= 4) is popped and resolved and a new lower bound on x (e.g., x >= 2) is introduced until (x >= 3) is
 *       pooped, the conflictlbcount indicates that bound change is currently present; that is (var->conflictlbcount ==
 *       conflict->count); however the (x >= 3) only has be explained if conflictlb matches that one; that is
 *       (var->conflictlb == bdchginfo->newbound); otherwise it redundant/invalid.
 */
static
SCIP_Bool bdchginfoIsInvalid(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   SCIP_VAR* var;

   assert(bdchginfo != NULL);

   var = SCIPbdchginfoGetVar(bdchginfo);
   assert(var != NULL);

   /* the bound change info of a binary (domained) variable can never be invalid since the concepts of relaxed bounds
    * and bound widening do not make sense for these type of variables
    */
   if( SCIPvarIsBinary(var) )
      return FALSE;

   /* check if the bdchginfo is invaild since a tight/weaker bound change was already explained */
   if( SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER )
   {
      if( var->conflictlbcount != conflict->count || var->conflictlb != SCIPbdchginfoGetNewbound(bdchginfo) ) /*lint !e777*/
      {
         assert(!SCIPvarIsBinary(var));
         return TRUE;
      }
   }
   else
   {
      assert(SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER);

      if( var->conflictubcount != conflict->count || var->conflictub != SCIPbdchginfoGetNewbound(bdchginfo) ) /*lint !e777*/
      {
         assert(!SCIPvarIsBinary(var));
         return TRUE;
      }
   }

   return FALSE;
}

/** adds a bound change to a conflict set */
static
SCIP_RETCODE conflictsetAddBound(
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change to add to the conflict set */
   SCIP_Real             relaxedbd           /**< relaxed bound */
   )
{
   SCIP_BDCHGINFO** bdchginfos;
   SCIP_Real* relaxedbds;
   int* sortvals;
   SCIP_VAR* var;
   SCIP_BOUNDTYPE boundtype;
   int idx;
   int sortval;
   int pos;

   assert(conflictset != NULL);
   assert(bdchginfo != NULL);

   /* allocate memory for additional element */
   SCIP_CALL( conflictsetEnsureBdchginfosMem(conflictset, blkmem, set, conflictset->nbdchginfos+1) );

   /* insert the new bound change in the arrays sorted by increasing variable index and by bound type */
   bdchginfos = conflictset->bdchginfos;
   relaxedbds = conflictset->relaxedbds;
   sortvals = conflictset->sortvals;
   var = SCIPbdchginfoGetVar(bdchginfo);
   boundtype = SCIPbdchginfoGetBoundtype(bdchginfo);
   idx = SCIPvarGetIndex(var);
   assert(idx < INT_MAX/2);
   assert((int)boundtype == 0 || (int)boundtype == 1);
   sortval = 2*idx + (int)boundtype; /* first sorting criteria: variable index, second criteria: boundtype */

   /* insert new element into the sorted arrays; if an element exits with the same value insert the new element afterwards
    *
    * @todo check if it better (faster) to first search for the position O(log n) and compare the sort values and if
    *       they are equal just replace the element and if not run the insert method O(n)
    */

   SCIPsortedvecInsertIntPtrReal(sortvals, (void**)bdchginfos, relaxedbds, sortval, (void*)bdchginfo, relaxedbd, &conflictset->nbdchginfos, &pos);
   assert(pos == conflictset->nbdchginfos - 1 || sortval < sortvals[pos+1]);

   /* merge multiple bound changes */
   if( pos > 0 && sortval == sortvals[pos-1] )
   {
      /* this is a multiple bound change */
      if( SCIPbdchginfoIsTighter(bdchginfo, bdchginfos[pos-1]) )
      {
         /* remove the "old" bound change since the "new" one in tighter */
         SCIPsortedvecDelPosIntPtrReal(sortvals, (void**)bdchginfos, relaxedbds, pos-1, &conflictset->nbdchginfos);
      }
      else if( SCIPbdchginfoIsTighter(bdchginfos[pos-1], bdchginfo) )
      {
         /* remove the "new"  bound change since the "old" one is tighter */
         SCIPsortedvecDelPosIntPtrReal(sortvals, (void**)bdchginfos, relaxedbds, pos, &conflictset->nbdchginfos);
      }
      else
      {
         /* both bound change are equivalent; hence, keep the worse relaxed bound and remove one of them */
         relaxedbds[pos-1] = boundtype == SCIP_BOUNDTYPE_LOWER ? MAX(relaxedbds[pos-1], relaxedbd) : MIN(relaxedbds[pos-1], relaxedbd);
         SCIPsortedvecDelPosIntPtrReal(sortvals, (void**)bdchginfos, relaxedbds, pos, &conflictset->nbdchginfos);
      }
   }

   return SCIP_OKAY;
}

/** adds given bound changes to a conflict set */
static
SCIP_RETCODE conflictsetAddBounds(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound changes to add to the conflict set */
   int                   nbdchginfos         /**< number of bound changes to add */
   )
{
   SCIP_BDCHGINFO** confbdchginfos;
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_Real* confrelaxedbds;
   int* confsortvals;
   int confnbdchginfos;
   int idx;
   int sortval;
   int i;
   SCIP_BOUNDTYPE boundtype;

   assert(conflict != NULL);
   assert(conflictset != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(bdchginfos != NULL || nbdchginfos == 0);

   /* nothing to add */
   if( nbdchginfos == 0 )
      return SCIP_OKAY;

   assert(bdchginfos != NULL);

   /* only one element to add, use the single insertion method */
   if( nbdchginfos == 1 )
   {
      bdchginfo = bdchginfos[0];
      assert(bdchginfo != NULL);

      if( !bdchginfoIsInvalid(conflict, bdchginfo) )
      {
         SCIP_CALL( conflictsetAddBound(conflictset, blkmem, set, bdchginfo, SCIPbdchginfoGetRelaxedBound(bdchginfo)) );
      }
      else
      {
         SCIPsetDebugMsg(set, "-> bound change info [%d:<%s> %s %g] is invaild -> ignore it\n", SCIPbdchginfoGetDepth(bdchginfo),
            SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo));
      }

      return SCIP_OKAY;
   }

   confnbdchginfos = conflictset->nbdchginfos;

   /* allocate memory for additional element */
   SCIP_CALL( conflictsetEnsureBdchginfosMem(conflictset, blkmem, set, confnbdchginfos + nbdchginfos) );

   confbdchginfos = conflictset->bdchginfos;
   confrelaxedbds = conflictset->relaxedbds;
   confsortvals = conflictset->sortvals;

   assert(SCIP_BOUNDTYPE_LOWER == FALSE);/*lint !e641*/
   assert(SCIP_BOUNDTYPE_UPPER == TRUE);/*lint !e641*/

   for( i = 0; i < nbdchginfos; ++i )
   {
      bdchginfo = bdchginfos[i];
      assert(bdchginfo != NULL);

      /* add only valid bound change infos */
      if( !bdchginfoIsInvalid(conflict, bdchginfo) )
      {
         /* calculate sorting value */
         boundtype = SCIPbdchginfoGetBoundtype(bdchginfo);
         assert(SCIPbdchginfoGetVar(bdchginfo) != NULL);

         idx = SCIPvarGetIndex(SCIPbdchginfoGetVar(bdchginfo));
         assert(idx < INT_MAX/2);

         assert((int)boundtype == 0 || (int)boundtype == 1);
         sortval = 2*idx + (int)boundtype; /* first sorting criteria: variable index, second criteria: boundtype */

         /* add new element */
         confbdchginfos[confnbdchginfos] = bdchginfo;
         confrelaxedbds[confnbdchginfos] = SCIPbdchginfoGetRelaxedBound(bdchginfo);
         confsortvals[confnbdchginfos] = sortval;
         ++confnbdchginfos;
      }
      else
      {
         SCIPsetDebugMsg(set, "-> bound change info [%d:<%s> %s %g] is invaild -> ignore it\n", SCIPbdchginfoGetDepth(bdchginfo),
            SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo));
      }
   }
   assert(confnbdchginfos <= conflictset->nbdchginfos + nbdchginfos);

   /* sort and merge the new conflict set */
   if( confnbdchginfos > conflictset->nbdchginfos )
   {
      int k = 0;

      /* sort array */
      SCIPsortIntPtrReal(confsortvals, (void**)confbdchginfos, confrelaxedbds, confnbdchginfos);

      i = 1;
      /* merge multiple bound changes */
      while( i < confnbdchginfos )
      {
         assert(i > k);

         /* is this a multiple bound change */
         if( confsortvals[k] == confsortvals[i] )
         {
            if( SCIPbdchginfoIsTighter(confbdchginfos[k], confbdchginfos[i]) )
               ++i;
            else if( SCIPbdchginfoIsTighter(confbdchginfos[i], confbdchginfos[k]) )
            {
               /* replace worse bound change info by tighter bound change info */
               confbdchginfos[k] = confbdchginfos[i];
               confrelaxedbds[k] = confrelaxedbds[i];
               confsortvals[k] = confsortvals[i];
               ++i;
            }
            else
            {
               assert(confsortvals[k] == confsortvals[i]);

               /* both bound change are equivalent; hence, keep the worse relaxed bound and remove one of them */
               confrelaxedbds[k] = (confsortvals[k] % 2 == 0) ? MAX(confrelaxedbds[k], confrelaxedbds[i]) : MIN(confrelaxedbds[k], confrelaxedbds[i]);
               ++i;
            }
         }
         else
         {
            /* all bound change infos must be valid */
            assert(!bdchginfoIsInvalid(conflict, confbdchginfos[k]));

            ++k;
            /* move next comparison element to the correct position */
            if( k != i )
            {
               confbdchginfos[k] = confbdchginfos[i];
               confrelaxedbds[k] = confrelaxedbds[i];
               confsortvals[k] = confsortvals[i];
            }
            ++i;
         }
      }
      /* last bound change infos must also be valid */
      assert(!bdchginfoIsInvalid(conflict, confbdchginfos[k]));
      /* the number of bound change infos cannot be decreased, it would mean that the conflict set was not merged
       * before
       */
      assert(conflictset->nbdchginfos <= k + 1 );
      assert(k + 1 <= confnbdchginfos);

      conflictset->nbdchginfos = k + 1;
   }

   return SCIP_OKAY;
}

/** calculates the conflict and the repropagation depths of the conflict set */
static
void conflictsetCalcConflictDepth(
   SCIP_CONFLICTSET*     conflictset         /**< conflict set */
   )
{
   int maxdepth[2];
   int i;

   assert(conflictset != NULL);
   assert(conflictset->validdepth <= conflictset->insertdepth);

   /* get the depth of the last and last but one bound change */
   maxdepth[0] = conflictset->validdepth;
   maxdepth[1] = conflictset->validdepth;
   for( i = 0; i < conflictset->nbdchginfos; ++i )
   {
      int depth;

      depth = SCIPbdchginfoGetDepth(conflictset->bdchginfos[i]);
      assert(depth >= 0);
      if( depth > maxdepth[0] )
      {
         maxdepth[1] = maxdepth[0];
         maxdepth[0] = depth;
      }
      else if( depth > maxdepth[1] )
         maxdepth[1] = depth;
   }
   assert(maxdepth[0] >= maxdepth[1]);

   conflictset->conflictdepth = maxdepth[0];
   conflictset->repropdepth = maxdepth[1];
}

/** identifies the depth, at which the conflict set should be added:
 *  - if the branching rule operates on variables only, and if all branching variables up to a certain
 *    depth level are member of the conflict, the conflict constraint can only be violated in the subtree
 *    of the node at that depth, because in all other nodes, at least one of these branching variables
 *    violates its conflicting bound, such that the conflict constraint is feasible
 *  - if there is at least one branching variable in a node, we assume, that this branching was performed
 *    on variables, and that the siblings of this node are disjunct w.r.t. the branching variables' fixings
 *  - we have to add the conflict set at least in the valid depth of the initial conflict set,
 *    so we start searching at the first branching after this depth level, i.e. validdepth+1
 */
static
SCIP_RETCODE conflictsetCalcInsertDepth(
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   SCIP_Bool* branchingincluded;
   int currentdepth;
   int i;

   assert(conflictset != NULL);
   assert(set != NULL);
   assert(tree != NULL);

   /* the conflict set must not be inserted prior to its valid depth */
   conflictset->insertdepth = conflictset->validdepth;
   assert(conflictset->insertdepth >= 0);

   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);

   /* mark the levels for which a branching variable is included in the conflict set */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &branchingincluded, currentdepth+2) );
   BMSclearMemoryArray(branchingincluded, currentdepth+2);
   for( i = 0; i < conflictset->nbdchginfos; ++i )
   {
      int depth;

      depth = SCIPbdchginfoGetDepth(conflictset->bdchginfos[i]);
      depth = MIN(depth, currentdepth+1); /* put diving/probing/strong branching changes in this depth level */
      branchingincluded[depth] = TRUE;
   }

   /* skip additional depth levels where branching on the conflict variables was applied */
   while( conflictset->insertdepth < currentdepth && branchingincluded[conflictset->insertdepth+1] )
      conflictset->insertdepth++;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &branchingincluded);

   assert(conflictset->validdepth <= conflictset->insertdepth && conflictset->insertdepth <= currentdepth);

   return SCIP_OKAY;
}

/** checks whether the first conflict set is redundant to the second one */
static
SCIP_Bool conflictsetIsRedundant(
   SCIP_CONFLICTSET*     conflictset1,       /**< first conflict conflict set */
   SCIP_CONFLICTSET*     conflictset2        /**< second conflict conflict set */
   )
{
   int i1;
   int i2;

   assert(conflictset1 != NULL);
   assert(conflictset2 != NULL);

   /* if conflictset1 has smaller validdepth, it is definitely not redundant to conflictset2 */
   if( conflictset1->validdepth < conflictset2->validdepth )
      return FALSE;

   /* check, if all bound changes in conflictset2 are also present at least as tight in conflictset1;
    * we can stop immediately, if more bound changes are remaining in conflictset2 than in conflictset1
    */
   for( i1 = 0, i2 = 0; i2 < conflictset2->nbdchginfos && conflictset1->nbdchginfos - i1 >= conflictset2->nbdchginfos - i2;
        ++i1, ++i2 )
   {
      int sortval;

      assert(i2 == 0 || conflictset2->sortvals[i2-1] < conflictset2->sortvals[i2]);

      sortval = conflictset2->sortvals[i2];
      for( ; i1 < conflictset1->nbdchginfos && conflictset1->sortvals[i1] < sortval; ++i1 )
      {
         /* while scanning conflictset1, check consistency */
         assert(i1 == 0 || conflictset1->sortvals[i1-1] < conflictset1->sortvals[i1]);
      }
      if( i1 >= conflictset1->nbdchginfos || conflictset1->sortvals[i1] > sortval
         || SCIPbdchginfoIsTighter(conflictset2->bdchginfos[i2], conflictset1->bdchginfos[i1]) )
         return FALSE;
   }

   return (i2 == conflictset2->nbdchginfos);
}

#ifdef SCIP_DEBUG
/** prints a conflict set to the screen */
static
void conflictsetPrint(
   SCIP_CONFLICTSET*     conflictset         /**< conflict set */
   )
{
   int i;

   assert(conflictset != NULL);
   for( i = 0; i < conflictset->nbdchginfos; ++i )
   {
      SCIPdebugPrintf(" [%d:<%s> %s %g(%g)]", SCIPbdchginfoGetDepth(conflictset->bdchginfos[i]),
         SCIPvarGetName(SCIPbdchginfoGetVar(conflictset->bdchginfos[i])),
         SCIPbdchginfoGetBoundtype(conflictset->bdchginfos[i]) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
         SCIPbdchginfoGetNewbound(conflictset->bdchginfos[i]), conflictset->relaxedbds[i]);
   }
   SCIPdebugPrintf("\n");
}
#endif

/** resizes proofsets array to be able to store at least num entries */
static
SCIP_RETCODE conflictEnsureProofsetsMem(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);

   if( num > conflict->proofsetssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conflict->proofsets, newsize) );
      conflict->proofsetssize = newsize;
   }
   assert(num <= conflict->proofsetssize);

   return SCIP_OKAY;
}

/** resizes conflictsets array to be able to store at least num entries */
static
SCIP_RETCODE conflictEnsureConflictsetsMem(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);

   if( num > conflict->conflictsetssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conflict->conflictsets, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&conflict->conflictsetscores, newsize) );
      conflict->conflictsetssize = newsize;
   }
   assert(num <= conflict->conflictsetssize);

   return SCIP_OKAY;
}

/** add a proofset to the list of all proofsets */
static
SCIP_RETCODE conflictInsertProofset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROOFSET*        proofset            /**< proof set to add */
   )
{
   assert(conflict != NULL);
   assert(proofset != NULL);

   /* insert proofset into the sorted proofsets array */
   SCIP_CALL( conflictEnsureProofsetsMem(conflict, set, conflict->nproofsets + 1) );

   conflict->proofsets[conflict->nproofsets] = proofset;
   ++conflict->nproofsets;

   return SCIP_OKAY;
}

/** inserts conflict set into sorted conflictsets array and deletes the conflict set pointer */
static
SCIP_RETCODE conflictInsertConflictset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONFLICTSET**    conflictset         /**< pointer to conflict set to insert */
   )
{
   SCIP_Real score;
   int pos;
   int i;
   int j;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(conflictset != NULL);
   assert(*conflictset != NULL);
   assert((*conflictset)->validdepth <= (*conflictset)->insertdepth);
   assert(set->conf_allowlocal || (*conflictset)->validdepth == 0);

   /* calculate conflict and repropagation depth */
   conflictsetCalcConflictDepth(*conflictset);

   /* if we apply repropagations, the conflict set should be inserted at most at its repropdepth */
   if( set->conf_repropagate )
      (*conflictset)->insertdepth = MIN((*conflictset)->insertdepth, (*conflictset)->repropdepth);
   else
      (*conflictset)->repropdepth = INT_MAX;
   assert((*conflictset)->insertdepth <= (*conflictset)->repropdepth);

   SCIPsetDebugMsg(set, "inserting conflict set (valid: %d, insert: %d, conf: %d, reprop: %d):\n",
      (*conflictset)->validdepth, (*conflictset)->insertdepth, (*conflictset)->conflictdepth, (*conflictset)->repropdepth);
   SCIPdebug(conflictsetPrint(*conflictset));

   /* get the score of the conflict set */
   score = conflictsetCalcScore(*conflictset, set);

   /* check, if conflict set is redundant to a better conflict set */
   for( pos = 0; pos < conflict->nconflictsets && score < conflict->conflictsetscores[pos]; ++pos )
   {
      /* check if conflict set is redundant with respect to conflictsets[pos] */
      if( conflictsetIsRedundant(*conflictset, conflict->conflictsets[pos]) )
      {
         SCIPsetDebugMsg(set, " -> conflict set is redundant to: ");
         SCIPdebug(conflictsetPrint(conflict->conflictsets[pos]));
         conflictsetFree(conflictset, blkmem);
         return SCIP_OKAY;
      }

      /**@todo like in sepastore.c: calculate overlap between conflictsets -> large overlap reduces score */

   }

   /* insert conflictset into the sorted conflictsets array */
   SCIP_CALL( conflictEnsureConflictsetsMem(conflict, set, conflict->nconflictsets + 1) );
   for( i = conflict->nconflictsets; i > pos; --i )
   {
      assert(score >= conflict->conflictsetscores[i-1]);
      conflict->conflictsets[i] = conflict->conflictsets[i-1];
      conflict->conflictsetscores[i] = conflict->conflictsetscores[i-1];
   }
   conflict->conflictsets[pos] = *conflictset;
   conflict->conflictsetscores[pos] = score;
   conflict->nconflictsets++;

   /* remove worse conflictsets that are redundant to the new conflictset */
   for( i = pos+1, j = pos+1; i < conflict->nconflictsets; ++i )
   {
      if( conflictsetIsRedundant(conflict->conflictsets[i], *conflictset) )
      {
         SCIPsetDebugMsg(set, " -> conflict set dominates: ");
         SCIPdebug(conflictsetPrint(conflict->conflictsets[i]));
         conflictsetFree(&conflict->conflictsets[i], blkmem);
      }
      else
      {
         assert(j <= i);
         conflict->conflictsets[j] = conflict->conflictsets[i];
         conflict->conflictsetscores[j] = conflict->conflictsetscores[i];
         j++;
      }
   }
   assert(j <= conflict->nconflictsets);
   conflict->nconflictsets = j;

#ifdef SCIP_CONFGRAPH
   confgraphMarkConflictset(*conflictset);
#endif

   *conflictset = NULL; /* ownership of pointer is now in the conflictsets array */

   return SCIP_OKAY;
}

/** calculates the maximal size of conflict sets to be used */
static
int conflictCalcMaxsize(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   int maxsize;

   assert(set != NULL);
   assert(prob != NULL);

   maxsize = (int)(set->conf_maxvarsfac * (prob->nvars - prob->ncontvars));
   maxsize = MAX(maxsize, set->conf_minmaxvars);

   return maxsize;
}

/** increases the conflict score of the variable in the given direction */
static
SCIP_RETCODE incVSIDS(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound for which the score should be increased */
   SCIP_Real             value,              /**< value of the bound */
   SCIP_Real             weight              /**< weight of this VSIDS updates */
   )
{
   SCIP_BRANCHDIR branchdir;

   assert(var != NULL);
   assert(stat != NULL);

   /* weight the VSIDS by the given weight */
   weight *= stat->vsidsweight;

   if( SCIPsetIsZero(set, weight) )
      return SCIP_OKAY;

   branchdir = (boundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS); /*lint !e641*/
   SCIP_CALL( SCIPvarIncVSIDS(var, blkmem, set, stat, branchdir, value, weight) );
   SCIPhistoryIncVSIDS(stat->glbhistory, branchdir,  weight);
   SCIPhistoryIncVSIDS(stat->glbhistorycrun, branchdir,  weight);

   return SCIP_OKAY;
}

/** update conflict statistics */
static
SCIP_RETCODE updateStatistics(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set to add to the tree */
   int                   insertdepth         /**< depth level at which the conflict set should be added */
   )
{
   if( insertdepth > 0 )
   {
      conflict->nappliedlocconss++;
      conflict->nappliedlocliterals += conflictset->nbdchginfos;
   }
   else
   {
      int i;
      int conflictlength;
      conflictlength = conflictset->nbdchginfos;

      for( i = 0; i < conflictlength; i++ )
      {
         SCIP_VAR* var;
         SCIP_BRANCHDIR branchdir;
         SCIP_BOUNDTYPE boundtype;
         SCIP_Real bound;

         assert(stat != NULL);

         var = conflictset->bdchginfos[i]->var;
         boundtype = SCIPbdchginfoGetBoundtype(conflictset->bdchginfos[i]);
         bound = conflictset->relaxedbds[i];

         branchdir = (boundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS); /*lint !e641*/

         SCIP_CALL( SCIPvarIncNActiveConflicts(var, blkmem, set, stat,  branchdir, bound, (SCIP_Real)conflictlength) );
         SCIPhistoryIncNActiveConflicts(stat->glbhistory, branchdir, (SCIP_Real)conflictlength);
         SCIPhistoryIncNActiveConflicts(stat->glbhistorycrun, branchdir, (SCIP_Real)conflictlength);

         /* each variable which is part of the conflict gets an increase in the VSIDS */
         SCIP_CALL( incVSIDS(var, blkmem, set, stat, boundtype, bound, set->conf_conflictweight) );
      }
      conflict->nappliedglbconss++;
      conflict->nappliedglbliterals += conflictset->nbdchginfos;
   }

   return SCIP_OKAY;
}


/** check conflict set for redundancy, other conflicts in the same conflict analysis could have led to global reductions
 *  an made this conflict set redundant
 */
static
SCIP_Bool checkRedundancy(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONFLICTSET*     conflictset         /**< conflict set */
   )
{
   SCIP_BDCHGINFO** bdchginfos;
   SCIP_VAR* var;
   SCIP_Real* relaxedbds;
   SCIP_Real bound;
   int v;

   assert(set != NULL);
   assert(conflictset != NULL);

   bdchginfos = conflictset->bdchginfos;
   relaxedbds = conflictset->relaxedbds;
   assert(bdchginfos != NULL);
   assert(relaxedbds != NULL);

   /* check all boundtypes and bounds for redundancy */
   for( v = conflictset->nbdchginfos - 1; v >= 0; --v )
   {
      var = SCIPbdchginfoGetVar(bdchginfos[v]);
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) >= 0);

      /* check if the relaxed bound is really a relaxed bound */
      assert(SCIPbdchginfoGetBoundtype(bdchginfos[v]) == SCIP_BOUNDTYPE_LOWER || SCIPsetIsGE(set, relaxedbds[v], SCIPbdchginfoGetNewbound(bdchginfos[v])));
      assert(SCIPbdchginfoGetBoundtype(bdchginfos[v]) == SCIP_BOUNDTYPE_UPPER || SCIPsetIsLE(set, relaxedbds[v], SCIPbdchginfoGetNewbound(bdchginfos[v])));

      bound = relaxedbds[v];

      if( SCIPbdchginfoGetBoundtype(bdchginfos[v]) == SCIP_BOUNDTYPE_UPPER )
      {
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            assert(SCIPsetIsIntegral(set, bound));
            bound += 1.0;
         }

         /* check if the bound is already fulfilled globally */
         if( SCIPsetIsFeasGE(set, SCIPvarGetLbGlobal(var), bound) )
            return TRUE;
      }
      else
      {
         assert(SCIPbdchginfoGetBoundtype(bdchginfos[v]) == SCIP_BOUNDTYPE_LOWER);

         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            assert(SCIPsetIsIntegral(set, bound));
            bound -= 1.0;
         }

         /* check if the bound is already fulfilled globally */
         if( SCIPsetIsFeasLE(set, SCIPvarGetUbGlobal(var), bound) )
            return TRUE;
      }
   }

   return FALSE;
}

/** find global fixings which can be derived from the new conflict set */
static
SCIP_RETCODE detectImpliedBounds(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set to add to the tree */
   int*                  nbdchgs,            /**< number of global deducted bound changes due to the conflict set */
   int*                  nredvars,           /**< number of redundant and removed variables from conflict set */
   SCIP_Bool*            redundant           /**< did we found a global reduction on a conflict set variable, which makes this conflict redundant */
   )
{
   SCIP_BDCHGINFO** bdchginfos;
   SCIP_Real* relaxedbds;
   SCIP_VAR* var;
   SCIP_Bool* boundtypes;
   SCIP_Real* bounds;
   SCIP_Longint* nbinimpls;
   int* sortvals;
   SCIP_Real bound;
   SCIP_Bool isupper;
   int ntrivialredvars;
   int nbdchginfos;
   int nzeroimpls;
   int v;

   assert(set != NULL);
   assert(prob != NULL);
   assert(SCIPprobIsTransformed(prob));
   assert(conflictset != NULL);
   assert(nbdchgs != NULL);
   assert(nredvars != NULL);
   /* only check conflict sets with more than one variable */
   assert(conflictset->nbdchginfos > 1);

   *nbdchgs = 0;
   *nredvars = 0;

   /* due to other conflict in the same conflict analysis, this conflict set might have become redundant */
   *redundant = checkRedundancy(set, conflictset);

   if( *redundant )
      return SCIP_OKAY;

   bdchginfos = conflictset->bdchginfos;
   relaxedbds = conflictset->relaxedbds;
   nbdchginfos = conflictset->nbdchginfos;
   sortvals = conflictset->sortvals;

   assert(bdchginfos != NULL);
   assert(relaxedbds != NULL);
   assert(sortvals != NULL);

   /* check if the boolean representation of boundtypes matches the 'standard' definition */
   assert(SCIP_BOUNDTYPE_LOWER == FALSE); /*lint !e641*/
   assert(SCIP_BOUNDTYPE_UPPER == TRUE); /*lint !e641*/

   ntrivialredvars = 0;

   /* due to multiple conflict sets for one conflict, it can happen, that we already have redundant information in the
    * conflict set
    */
   for( v = nbdchginfos - 1; v >= 0; --v )
   {
      var = SCIPbdchginfoGetVar(bdchginfos[v]);
      bound = relaxedbds[v];
      isupper = (SCIP_Bool) SCIPboundtypeOpposite(SCIPbdchginfoGetBoundtype(bdchginfos[v]));

      /* for integral variable we can increase/decrease the conflicting bound */
      if( SCIPvarIsIntegral(var) )
         bound += (isupper ? -1.0 : +1.0);

      /* if conflict variable cannot fulfill the conflict we can remove it */
      if( (isupper && SCIPsetIsFeasLT(set, bound, SCIPvarGetLbGlobal(var))) ||
         (!isupper && SCIPsetIsFeasGT(set, bound, SCIPvarGetUbGlobal(var))) )
      {
         SCIPsetDebugMsg(set, "remove redundant variable <%s> from conflict set\n", SCIPvarGetName(var));

         bdchginfos[v] = bdchginfos[nbdchginfos - 1];
         relaxedbds[v] = relaxedbds[nbdchginfos - 1];
         sortvals[v] = sortvals[nbdchginfos - 1];

         --nbdchginfos;
         ++ntrivialredvars;
      }
   }
   assert(ntrivialredvars + nbdchginfos == conflictset->nbdchginfos);

   SCIPsetDebugMsg(set, "trivially removed %d redundant of %d variables from conflictset (%p)\n", ntrivialredvars, conflictset->nbdchginfos, (void*)conflictset);
   conflictset->nbdchginfos = nbdchginfos;

   /* all variables where removed, the conflict cannot be fulfilled, i.e., we have an infeasibility proof */
   if( conflictset->nbdchginfos == 0 )
      return SCIP_OKAY;

   /* do not check to big or trivial conflicts */
   if( conflictset->nbdchginfos > set->conf_maxvarsdetectimpliedbounds || conflictset->nbdchginfos == 1 )
   {
      *nredvars = ntrivialredvars;
      return SCIP_OKAY;
   }

   /* create array of boundtypes, and bound values in conflict set */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &boundtypes, nbdchginfos) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bounds, nbdchginfos) );
   /* memory for the estimates for binary implications used for sorting */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &nbinimpls, nbdchginfos) );

   nzeroimpls = 0;

   /* collect estimates and initialize variables, boundtypes, and bounds array */
   for( v = 0; v < nbdchginfos; ++v )
   {
      var = SCIPbdchginfoGetVar(bdchginfos[v]);
      boundtypes[v] = (SCIP_Bool) SCIPboundtypeOpposite(SCIPbdchginfoGetBoundtype(bdchginfos[v]));
      bounds[v] = relaxedbds[v];

      assert(SCIPvarGetProbindex(var) >= 0);

      /* check if the relaxed bound is really a relaxed bound */
      assert(SCIPbdchginfoGetBoundtype(bdchginfos[v]) == SCIP_BOUNDTYPE_LOWER || SCIPsetIsGE(set, relaxedbds[v], SCIPbdchginfoGetNewbound(bdchginfos[v])));
      assert(SCIPbdchginfoGetBoundtype(bdchginfos[v]) == SCIP_BOUNDTYPE_UPPER || SCIPsetIsLE(set, relaxedbds[v], SCIPbdchginfoGetNewbound(bdchginfos[v])));

      /* for continuous variables, we can only use the relaxed version of the bounds negation: !(x <= u) -> x >= u */
      if( SCIPvarIsBinary(var) )
      {
         if( !boundtypes[v] )
         {
            assert(SCIPsetIsZero(set, bounds[v]));
            bounds[v] = 1.0;
            nbinimpls[v] = (SCIP_Longint)SCIPvarGetNCliques(var, TRUE) * 2;
         }
         else
         {
            assert(SCIPsetIsEQ(set, bounds[v], 1.0));
            bounds[v] = 0.0;
            nbinimpls[v] = (SCIP_Longint)SCIPvarGetNCliques(var, FALSE) * 2;
         }
      }
      else if( SCIPvarIsIntegral(var) )
      {
         assert(SCIPsetIsIntegral(set, bounds[v]));

         bounds[v] += ((!boundtypes[v]) ? +1.0 : -1.0);
         nbinimpls[v] = (boundtypes[v] ? SCIPvarGetNVlbs(var) : SCIPvarGetNVubs(var));
      }
      else if( ((!boundtypes[v]) && SCIPsetIsFeasEQ(set, SCIPvarGetLbGlobal(var), bounds[v]))
         || ((boundtypes[v]) && SCIPsetIsFeasEQ(set, SCIPvarGetUbGlobal(var), bounds[v])) )
      {
         /* the literal is satisfied in global bounds (may happen due to weak "negation" of continuous variables)
          * -> discard the conflict constraint
          */
         break;
      }
      else
      {
         nbinimpls[v] = (boundtypes[v] ? SCIPvarGetNVlbs(var) : SCIPvarGetNVubs(var));
      }

      if( nbinimpls[v] == 0 )
         ++nzeroimpls;
   }

   /* starting to derive global bound changes */
   if( v == nbdchginfos && ((!set->conf_fullshortenconflict && nzeroimpls < 2) || (set->conf_fullshortenconflict && nzeroimpls < nbdchginfos)) )
   {
      SCIP_VAR** vars;
      SCIP_Bool* redundants;
      SCIP_Bool glbinfeas;

      /* sort variables in increasing order of binary implications to gain speed later on */
      SCIPsortLongPtrRealRealBool(nbinimpls, (void**)bdchginfos, relaxedbds, bounds, boundtypes, v);

      SCIPsetDebugMsg(set, "checking for global reductions and redundant conflict variables(in %s) on conflict:\n", SCIPprobGetName(prob));
      SCIPsetDebugMsg(set, "[");
      for( v = 0; v < nbdchginfos; ++v )
      {
         SCIPsetDebugMsgPrint(set, "%s %s %g", SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfos[v])), (!boundtypes[v]) ? ">=" : "<=", bounds[v]);
         if( v < nbdchginfos - 1 )
            SCIPsetDebugMsgPrint(set, ", ");
      }
      SCIPsetDebugMsgPrint(set, "]\n");

      SCIP_CALL( SCIPsetAllocBufferArray(set, &vars, v) );
      SCIP_CALL( SCIPsetAllocCleanBufferArray(set, &redundants, v) );

      /* initialize conflict variable data */
      for( v = 0; v < nbdchginfos; ++v )
         vars[v] = SCIPbdchginfoGetVar(bdchginfos[v]);

      SCIP_CALL( SCIPshrinkDisjunctiveVarSet(set->scip, vars, bounds, boundtypes, redundants, nbdchginfos, nredvars, \
            nbdchgs, redundant, &glbinfeas, set->conf_fullshortenconflict) );

      if( glbinfeas )
      {
         SCIPsetDebugMsg(set, "conflict set (%p) led to global infeasibility\n", (void*) conflictset);
         goto TERMINATE;
      }

#ifdef SCIP_DEBUG
      if( *nbdchgs > 0 )
      {
         SCIPsetDebugMsg(set, "conflict set (%p) led to %d global bound reductions\n", (void*) conflictset, *nbdchgs);
      }
#endif

      /* remove as redundant marked variables */
      if( *redundant )
      {
         SCIPsetDebugMsg(set, "conflict set (%p) is redundant because at least one global reduction, fulfills the conflict constraint\n", (void*)conflictset);

         BMSclearMemoryArray(redundants, nbdchginfos);
      }
      else if( *nredvars > 0 )
      {
         assert(bdchginfos == conflictset->bdchginfos);
         assert(relaxedbds == conflictset->relaxedbds);
         assert(sortvals == conflictset->sortvals);

         for( v = nbdchginfos - 1; v >= 0; --v )
         {
            /* if conflict variable was marked to be redundant remove it */
            if( redundants[v] )
            {
               SCIPsetDebugMsg(set, "remove redundant variable <%s> from conflict set\n", SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfos[v])));

               bdchginfos[v] = bdchginfos[nbdchginfos - 1];
               relaxedbds[v] = relaxedbds[nbdchginfos - 1];
               sortvals[v] = sortvals[nbdchginfos - 1];

               /* reset redundants[v] to 0 */
               redundants[v] = 0;

               --nbdchginfos;
            }
         }
         assert((*nredvars) + nbdchginfos == conflictset->nbdchginfos);

         SCIPsetDebugMsg(set, "removed %d redundant of %d variables from conflictset (%p)\n", (*nredvars), conflictset->nbdchginfos, (void*)conflictset);
         conflictset->nbdchginfos = nbdchginfos;
      }

     TERMINATE:
      SCIPsetFreeCleanBufferArray(set, &redundants);
      SCIPsetFreeBufferArray(set, &vars);
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &nbinimpls);
   SCIPsetFreeBufferArray(set, &bounds);
   SCIPsetFreeBufferArray(set, &boundtypes);

   *nredvars += ntrivialredvars;

   return SCIP_OKAY;
}

/** tighten the bound of a singleton variable in a constraint
 *
 *  if the bound is contradicting with a global bound we cannot tighten the bound directly.
 *  in this case we need to create and add a constraint of size one such that propagating this constraint will
 *  enforce the infeasibility.
 */
static
SCIP_RETCODE tightenSingleVar(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_TREE*            tree,               /**< tree data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             val,                /**< coefficient of the variable */
   SCIP_Real             rhs,                /**< rhs of the constraint */
   SCIP_CONFTYPE         prooftype            /**< type of the proof */
   )
{
   SCIP_Real newbound;
   SCIP_BOUNDTYPE boundtype;

   assert(tree != NULL);

   /* if variable and coefficient are integral the rhs can be rounded down */
   if( SCIPvarIsIntegral(var) && SCIPsetIsIntegral(set, val) )
   {
      newbound = SCIPsetFeasFloor(set, rhs)/val;
      boundtype = (val > 0.0 ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
      SCIPvarAdjustBd(var, set, boundtype, &newbound);
   }
   else
   {
      newbound = rhs/val;
      boundtype = (val > 0.0 ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
      SCIPvarAdjustBd(var, set, boundtype, &newbound);
   }

   /* skip numerical unstable bound changes */
   if( (boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsLE(set, newbound, SCIPvarGetLbGlobal(var)))
      || (boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsGE(set, newbound, SCIPvarGetUbGlobal(var))) )
   {
      return SCIP_OKAY;
   }

   /* the new bound contradicts a global bound, we can cutoff the root node immediately */
   if( (boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsGT(set, newbound, SCIPvarGetUbGlobal(var)))
      || (boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsLT(set, newbound, SCIPvarGetLbGlobal(var))) )
   {
      SCIPsetDebugMsg(set, "detect global infeasibility at var <%s>: locdom=[%g,%g] glbdom=[%g,%g] new %s bound=%g\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetLbGlobal(var),
            SCIPvarGetUbGlobal(var), (boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper"), newbound);

      SCIP_CALL( SCIPnodeCutoff(tree->path[0], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
   }
   else
   {
      SCIPsetDebugMsg(set, "change %s bound of <%s>: %g -> %g\n", (boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper"),
            SCIPvarGetName(var), (boundtype == SCIP_BOUNDTYPE_LOWER ? SCIPvarGetLbGlobal(var) : SCIPvarGetUbGlobal(var)),
            newbound);

      if( lp->strongbranching )
      {
         SCIP_CONS* cons;
         SCIP_Real conslhs;
         SCIP_Real consrhs;
         char name[SCIP_MAXSTRLEN];

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "pc_fix_%s", SCIPvarGetName(var));

         if( boundtype == SCIP_BOUNDTYPE_UPPER )
         {
	    conslhs = -SCIPsetInfinity(set);
            consrhs = newbound;
         }
	 else
         {
	    conslhs = newbound;
            consrhs = SCIPsetInfinity(set);
	 }

         SCIP_CALL( SCIPcreateConsLinear(set->scip, &cons, name, 0, NULL, NULL, conslhs, consrhs,
               FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

         SCIP_CALL( SCIPaddCoefLinear(set->scip, cons, var, 1.0) );

         SCIP_CALL( SCIPprobAddCons(transprob, set, stat, cons) );
         SCIP_CALL( SCIPconsRelease(&cons, blkmem, set) );
      }
      else
      {
         SCIP_CALL( SCIPnodeAddBoundchg(tree->root, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, \
               eventqueue, cliquetable, var, newbound, boundtype, FALSE) );

         /* mark the node in the repropdepth to be propagated again */
         SCIPnodePropagateAgain(tree->path[0], set, stat, tree);
      }
   }

   ++conflict->nglbchgbds;

   if( prooftype == SCIP_CONFTYPE_INFEASLP || prooftype == SCIP_CONFTYPE_ALTINFPROOF )
   {
      ++conflict->dualrayinfnnonzeros; /* we count a global bound reduction as size 1 */
      ++conflict->ndualrayinfsuccess;
      ++conflict->ndualrayinfglobal;
      ++conflict->ninflpsuccess;
   }
   else
   {
      ++conflict->dualraybndnnonzeros; /* we count a global bound reduction as size 1 */
      ++conflict->ndualraybndsuccess;
      ++conflict->ndualraybndglobal;
      ++conflict->nboundlpsuccess;
   }

   return SCIP_OKAY;
}

/** calculates the minimal activity of a given aggregation row */
static
SCIP_Real aggrRowGetMinActivity(
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables (or NULL for global bounds) */
   SCIP_Real*            curvarubs           /**< current upper bounds of active problem variables (or NULL for global bounds) */
   )
{
   SCIP_VAR** vars;
   SCIP_Real minact = 0.0;
   int* inds;
   int nnz;
   int i;

   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   nnz = SCIPaggrRowGetNNz(aggrrow);
   inds = SCIPaggrRowGetInds(aggrrow);

   for( i = 0; i < nnz; i++ )
   {
      SCIP_Real val;
      int v = inds[i];

      assert(SCIPvarGetProbindex(vars[v]) == v);

      val = SCIPaggrRowGetProbvarValue(aggrrow, v);

      /* calculate the minimal activity */
      if( val > 0.0 )
         minact += val * (curvarlbs == NULL ? SCIPvarGetLbGlobal(vars[v]) : curvarlbs[v]);
      else
         minact += val * (curvarubs == NULL ? SCIPvarGetUbGlobal(vars[v]) : curvarubs[v]);
   }

   return minact;
}

/** calculates the minimal activity of a given set of bounds and coefficients */
static
SCIP_Real getMinActivity(
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_Real*            coefs,              /**< coefficients in sparse representation */
   int*                  inds,               /**< non-zero indices */
   int                   nnz,                /**< number of non-zero indices */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables (or NULL for global bounds) */
   SCIP_Real*            curvarubs           /**< current upper bounds of active problem variables (or NULL for global bounds) */
   )
{
   SCIP_VAR** vars;
   SCIP_Real minact = 0.0;
   int i;

   assert(coefs != NULL);
   assert(inds != NULL);

   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   for( i = 0; i < nnz; i++ )
   {
      SCIP_Real val;
      int v = inds[i];

      assert(SCIPvarGetProbindex(vars[v]) == v);

      val = coefs[i];

      /* calculate the minimal activity */
      if( val > 0.0 )
         minact += val * (curvarlbs == NULL ? SCIPvarGetLbGlobal(vars[v]) : curvarlbs[v]);
      else
         minact += val * (curvarubs == NULL ? SCIPvarGetUbGlobal(vars[v]) : curvarubs[v]);
   }

   return minact;
}

/** calculates the minimal activity of a given set of bounds and coefficients */
static
SCIP_Real getMaxActivity(
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_Real*            coefs,              /**< coefficients in sparse representation */
   int*                  inds,               /**< non-zero indices */
   int                   nnz,                /**< number of non-zero indices */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables (or NULL for global bounds) */
   SCIP_Real*            curvarubs           /**< current upper bounds of active problem variables (or NULL for global bounds) */
   )
{
   SCIP_VAR** vars;
   SCIP_Real maxact = 0.0;
   int i;

   assert(coefs != NULL);
   assert(inds != NULL);

   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   for( i = 0; i < nnz; i++ )
   {
      SCIP_Real val;
      int v = inds[i];

      assert(SCIPvarGetProbindex(vars[v]) == v);

      val = coefs[i];

      /* calculate the minimal activity */
      if( val < 0.0 )
         maxact += val * (curvarlbs == NULL ? SCIPvarGetLbGlobal(vars[v]) : curvarlbs[v]);
      else
         maxact += val * (curvarubs == NULL ? SCIPvarGetUbGlobal(vars[v]) : curvarubs[v]);
   }

   return maxact;
}

static
SCIP_RETCODE propagateLongProof(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_TREE*            tree,               /**< tree data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real*            coefs,              /**< coefficients in sparse representation */
   int*                  inds,               /**< non-zero indices */
   int                   nnz,                /**< number of non-zero indices */
   SCIP_Real             rhs,                /**< right-hand side */
   SCIP_CONFTYPE         conflicttype        /**< type of the conflict */
   )
{
   SCIP_VAR** vars;
   SCIP_Real minact;
   int i;

   assert(coefs != NULL);
   assert(inds != NULL);
   assert(nnz >= 0);

   vars = SCIPprobGetVars(transprob);
   minact = getMinActivity(transprob, coefs, inds, nnz, NULL, NULL);

   for( i = 0; i < nnz; i++ )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_Real resminact;
      SCIP_Real lb;
      SCIP_Real ub;
      int pos;

      pos = inds[i];
      val = coefs[i];
      var = vars[pos];
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      assert(!SCIPsetIsZero(set, val));

      resminact = minact;

      /* we got a potential new upper bound */
      if( val > 0.0 )
      {
         SCIP_Real newub;

         resminact -= (val * lb);
         newub = (rhs - resminact)/val;

         if( SCIPsetIsInfinity(set, newub) )
            continue;

         /* we cannot tighten the upper bound */
         if( SCIPsetIsGE(set, newub, ub) )
            continue;

         SCIP_CALL( tightenSingleVar(conflict, set, stat, tree, blkmem, origprob, transprob, reopt, lp, branchcand, \
               eventqueue, cliquetable, var, val, rhs-resminact, conflicttype) );
      }
      /* we got a potential new lower bound */
      else
      {
         SCIP_Real newlb;

         resminact -= (val * ub);
         newlb = (rhs - resminact)/val;

         if( SCIPsetIsInfinity(set, -newlb) )
            continue;

         /* we cannot tighten the lower bound */
         if( SCIPsetIsLE(set, newlb, lb) )
            continue;

         SCIP_CALL( tightenSingleVar(conflict, set, stat, tree, blkmem, origprob, transprob, reopt, lp, branchcand, \
               eventqueue, cliquetable, var, val, rhs-resminact, conflicttype) );
      }

      /* the minimal activity should stay unchanged because we tightened the bound that doesn't contribute to the
       * minimal activity
       */
      assert(SCIPsetIsEQ(set, minact, getMinActivity(transprob, coefs, inds, nnz, NULL, NULL)));
   }

   return SCIP_OKAY;
}


/** creates a proof constraint and tries to add it to the storage */
static
SCIP_RETCODE createAndAddProofcons(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict pool data */
   SCIP_PROOFSET*        proofset,           /**< proof set */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< tree data */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_CONS* cons;
   SCIP_CONS* upgdcons;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int* inds;
   SCIP_Real rhs;
   SCIP_Real fillin;
   SCIP_Real globalminactivity;
   SCIP_Bool toolong;
   SCIP_CONFTYPE conflicttype;
   char name[SCIP_MAXSTRLEN];
   int nnz;
   int i;

   assert(conflict != NULL);
   assert(conflictstore != NULL);
   assert(proofset != NULL);

   nnz = proofsetGetNVars(proofset);

   if( nnz == 0 )
      return SCIP_OKAY;

   vars = SCIPprobGetVars(transprob);

   rhs = proofsetGetRhs(proofset);
   assert(!SCIPsetIsInfinity(set, rhs));

   coefs = proofsetGetVals(proofset);
   assert(coefs != NULL);

   inds = proofsetGetInds(proofset);
   assert(inds != NULL);

   conflicttype = proofsetGetConftype(proofset);

   if( conflicttype == SCIP_CONFTYPE_ALTINFPROOF || conflicttype == SCIP_CONFTYPE_ALTBNDPROOF )
   {
      SCIP_Real globalmaxactivity = getMaxActivity(transprob, coefs, inds, nnz, NULL, NULL);

      /* check whether the alternative proof is redundant */
      if( SCIPsetIsLE(set, globalmaxactivity, rhs) )
         return SCIP_OKAY;
   }

   /* check whether the constraint proves global infeasibility */
   globalminactivity = getMinActivity(transprob, coefs, inds, nnz, NULL, NULL);
   if( SCIPsetIsGT(set, globalminactivity, rhs) )
   {
      SCIPsetDebugMsg(set, "detect global infeasibility: minactivity=%g, rhs=%g\n", globalminactivity, rhs);

      SCIP_CALL( SCIPnodeCutoff(tree->path[0], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );

      goto UPDATESTATISTICS;
   }

   if( set->conf_minmaxvars >= nnz  )
      toolong = FALSE;
   else
   {
      SCIP_Real maxnnz;

      if( transprob->startnconss < 100 )
         maxnnz = 0.85 * transprob->nvars;
      else
         maxnnz = (SCIP_Real)transprob->nvars;

      fillin = nnz;
      if( conflicttype == SCIP_CONFTYPE_INFEASLP || conflicttype == SCIP_CONFTYPE_ALTINFPROOF )
      {
         fillin += SCIPconflictstoreGetNDualInfProofs(conflictstore) * SCIPconflictstoreGetAvgNnzDualInfProofs(conflictstore);
         fillin /= (SCIPconflictstoreGetNDualInfProofs(conflictstore) + 1.0);
         toolong = (fillin > MIN(2.0 * stat->avgnnz, maxnnz));
      }
      else
      {
         assert(conflicttype == SCIP_CONFTYPE_BNDEXCEEDING || conflicttype == SCIP_CONFTYPE_ALTBNDPROOF);

         fillin += SCIPconflictstoreGetNDualBndProofs(conflictstore) * SCIPconflictstoreGetAvgNnzDualBndProofs(conflictstore);
         fillin /= (SCIPconflictstoreGetNDualBndProofs(conflictstore) + 1.0);
         toolong = (fillin > MIN(1.5 * stat->avgnnz, maxnnz));
      }

      toolong = (toolong && (nnz > set->conf_maxvarsfac * transprob->nvars));
   }

   /* don't store global dual proofs that are to long / have to much non-zeros */
   if( toolong )
   {
      SCIP_CALL( propagateLongProof(conflict, set, stat, reopt, tree, blkmem, origprob, transprob, lp, branchcand,
            eventqueue, cliquetable, coefs, inds, nnz, rhs, conflicttype) );
      return SCIP_OKAY;
   }

   if( conflicttype == SCIP_CONFTYPE_INFEASLP || conflicttype == SCIP_CONFTYPE_ALTINFPROOF )
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "dualproof_inf_%d", conflict->ndualrayinfsuccess);
   else if( conflicttype == SCIP_CONFTYPE_BNDEXCEEDING || conflicttype == SCIP_CONFTYPE_ALTBNDPROOF )
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "dualproof_bnd_%d", conflict->ndualraybndsuccess);
   else
      return SCIP_INVALIDCALL;

   SCIP_CALL( SCIPcreateConsLinear(set->scip, &cons, name, 0, NULL, NULL, -SCIPsetInfinity(set), rhs,
         FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   for( i = 0; i < nnz; i++ )
   {
      int v = inds[i];
      SCIP_CALL( SCIPaddCoefLinear(set->scip, cons, vars[v], coefs[i]) );
   }

   /* do not upgrade linear constraints of size 1 */
   if( nnz > 1 )
   {
      upgdcons = NULL;
      /* try to automatically convert a linear constraint into a more specific and more specialized constraint */
      SCIP_CALL( SCIPupgradeConsLinear(set->scip, cons, &upgdcons) );
      if( upgdcons != NULL )
      {
         SCIP_CALL( SCIPreleaseCons(set->scip, &cons) );
         cons = upgdcons;

         if( conflicttype == SCIP_CONFTYPE_INFEASLP )
            conflicttype = SCIP_CONFTYPE_ALTINFPROOF;
         else if( conflicttype == SCIP_CONFTYPE_BNDEXCEEDING )
            conflicttype = SCIP_CONFTYPE_ALTBNDPROOF;
      }
   }

   /* mark constraint to be a conflict */
   SCIPconsMarkConflict(cons);

   /* add constraint to storage */
   if( conflicttype == SCIP_CONFTYPE_INFEASLP || conflicttype == SCIP_CONFTYPE_ALTINFPROOF )
   {
      /* add constraint based on dual ray to storage */
      SCIP_CALL( SCIPconflictstoreAddDualraycons(conflictstore, cons, blkmem, set, stat, transprob, reopt) );
   }
   else
   {
      SCIP_Real scale = 1.0;
      SCIP_Bool updateside = FALSE;

      /* In some cases the constraint could not be updated to a more special type. However, it is possible that
       * constraint got scaled. Therefore, we need to be very careful when updating the lhs/rhs after the incumbent
       * solution has improved.
       */
      if( conflicttype == SCIP_CONFTYPE_BNDEXCEEDING )
      {
         SCIP_Real side;

#ifndef NDEBUG
         SCIP_CONSHDLR* conshdlr = SCIPconsGetHdlr(cons);

         assert(conshdlr != NULL);
         assert(strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0);
#endif
         side = SCIPgetLhsLinear(set->scip, cons);

         if( !SCIPsetIsInfinity(set, -side) )
         {
            if( SCIPsetIsZero(set, side) )
            {
               scale = 1.0;
            }
            else
            {
               scale = proofsetGetRhs(proofset) / side;
               assert(SCIPsetIsNegative(set, scale));
            }
         }
         else
         {
            side = SCIPgetRhsLinear(set->scip, cons);
            assert(!SCIPsetIsInfinity(set, side));

            if( SCIPsetIsZero(set, side) )
            {
               scale = 1.0;
            }
            else
            {
               scale = proofsetGetRhs(proofset) / side;
               assert(SCIPsetIsPositive(set, scale));
            }
         }
         updateside = TRUE;
      }

      /* add constraint based on dual solution to storage */
      SCIP_CALL( SCIPconflictstoreAddDualsolcons(conflictstore, cons, blkmem, set, stat, transprob, reopt, scale, updateside) );
   }

   /* add the constraint to the global problem */
   SCIP_CALL( SCIPprobAddCons(transprob, set, stat, cons) );

   SCIPsetDebugMsg(set, "added proof-constraint to node %p in depth 0 (nproofconss %d)\n", (void*)tree->path[0],
         SCIPconflictstoreGetNDualInfProofs(conflictstore));

   /* release the constraint */
   SCIP_CALL( SCIPreleaseCons(set->scip, &cons) );

  UPDATESTATISTICS:
   /* update statistics */
   if( conflicttype == SCIP_CONFTYPE_INFEASLP || conflicttype == SCIP_CONFTYPE_ALTINFPROOF )
   {
      conflict->dualrayinfnnonzeros += nnz;
      ++conflict->ndualrayinfglobal;
      ++conflict->ndualrayinfsuccess;
   }
   else
   {
      assert(conflicttype == SCIP_CONFTYPE_BNDEXCEEDING || conflicttype == SCIP_CONFTYPE_ALTBNDPROOF);
      conflict->dualraybndnnonzeros += nnz;
      ++conflict->ndualraybndglobal;
      ++conflict->ndualraybndsuccess;
   }
   return SCIP_OKAY;
}

/* create proof constraints out of proof sets */
static
SCIP_RETCODE conflictFlushProofset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(conflict != NULL);

   if( proofsetGetConftype(conflict->proofset) != SCIP_CONFTYPE_UNKNOWN )
   {
      /* only one variable has a coefficient different to zero, we add this bound change instead of a constraint */
      if( proofsetGetNVars(conflict->proofset) == 1 )
      {
         SCIP_VAR** vars;
         SCIP_Real* coefs;
         int* inds;
         SCIP_Real rhs;

         vars = SCIPprobGetVars(transprob);

         coefs = proofsetGetVals(conflict->proofset);
         inds = proofsetGetInds(conflict->proofset);
         rhs = proofsetGetRhs(conflict->proofset);

         SCIP_CALL( tightenSingleVar(conflict, set, stat, tree, blkmem, origprob, transprob, reopt, lp, \
               branchcand, eventqueue, cliquetable, vars[inds[0]], coefs[0], rhs, conflict->proofset->conflicttype) );
      }
      else
      {
         SCIP_Bool skipinitialproof = FALSE;

         /* prefer an infeasibility proof
          *
          * todo: check whether this is really what we want
          */
         if( set->conf_prefinfproof && proofsetGetConftype(conflict->proofset) == SCIP_CONFTYPE_BNDEXCEEDING )
         {
            int i;

            for( i = 0; i < conflict->nproofsets; i++ )
            {
               if( proofsetGetConftype(conflict->proofsets[i]) == SCIP_CONFTYPE_INFEASLP )
               {
                  skipinitialproof = TRUE;
                  break;
               }
            }
         }

         if( !skipinitialproof )
         {
            /* create and add the original proof */
            SCIP_CALL( createAndAddProofcons(conflict, conflictstore, conflict->proofset, set, stat, origprob, transprob, \
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, blkmem) );
         }
      }

      /* clear the proof set anyway */
      proofsetClear(conflict->proofset);
   }

   if( conflict->nproofsets > 0 )
   {
      int i;

      for( i = 0; i < conflict->nproofsets; i++ )
      {
         assert(conflict->proofsets[i] != NULL);
         assert(proofsetGetConftype(conflict->proofsets[i]) != SCIP_CONFTYPE_UNKNOWN);

         /* only one variable has a coefficient different to zero, we add this bound change instead of a constraint */
         if( proofsetGetNVars(conflict->proofsets[i]) == 1 )
         {
            SCIP_VAR** vars;
            SCIP_Real* coefs;
            int* inds;
            SCIP_Real rhs;

            vars = SCIPprobGetVars(transprob);

            coefs = proofsetGetVals(conflict->proofsets[i]);
            inds = proofsetGetInds(conflict->proofsets[i]);
            rhs = proofsetGetRhs(conflict->proofsets[i]);

            SCIP_CALL( tightenSingleVar(conflict, set, stat, tree, blkmem, origprob, transprob, reopt, lp,
                  branchcand, eventqueue, cliquetable, vars[inds[0]], coefs[0], rhs,
                  conflict->proofsets[i]->conflicttype) );
         }
         else
         {
            /* create and add proof constraint */
            SCIP_CALL( createAndAddProofcons(conflict, conflictstore, conflict->proofsets[i], set, stat, origprob, \
                  transprob, tree, reopt, lp, branchcand, eventqueue, cliquetable, blkmem) );
         }
      }

      /* free all proofsets */
      for( i = 0; i < conflict->nproofsets; i++ )
         proofsetFree(&conflict->proofsets[i], blkmem);

      conflict->nproofsets = 0;
   }

   return SCIP_OKAY;
}

/** adds the given conflict set as conflict constraint to the problem */
static
SCIP_RETCODE conflictAddConflictCons(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set to add to the tree */
   int                   insertdepth,        /**< depth level at which the conflict set should be added */
   SCIP_Bool*            success             /**< pointer to store whether the addition was successful */
   )
{
   SCIP_Bool redundant;
   int h;

   assert(conflict != NULL);
   assert(tree != NULL);
   assert(tree->path != NULL);
   assert(conflictset != NULL);
   assert(conflictset->validdepth <= insertdepth);
   assert(success != NULL);

   *success = FALSE;
   redundant = FALSE;

   /* try to derive global bound changes and shorten the conflictset by using implication and clique and variable bound
    * information
    */
   if( conflictset->nbdchginfos > 1 && insertdepth == 0 && !lp->strongbranching )
   {
      int nbdchgs;
      int nredvars;
#ifdef SCIP_DEBUG
      int oldnbdchginfos = conflictset->nbdchginfos;
#endif
      assert(conflictset->validdepth == 0);

      SCIPclockStart(conflict->dIBclock, set);

      /* find global bound changes which can be derived from the new conflict set */
      SCIP_CALL( detectImpliedBounds(set, transprob, conflictset, &nbdchgs, &nredvars, &redundant) );

      /* all variables where removed, we have an infeasibility proof */
      if( conflictset->nbdchginfos == 0 )
         return SCIP_OKAY;

      /* debug check for reduced conflict set */
      if( nredvars > 0 )
      {
         /* check conflict set on debugging solution */
         SCIP_CALL( SCIPdebugCheckConflict(blkmem, set, tree->root, conflictset->bdchginfos, conflictset->relaxedbds, conflictset->nbdchginfos) ); /*lint !e506 !e774*/
      }

#ifdef SCIP_DEBUG
      SCIPsetDebugMsg(set, " -> conflict set removed %d redundant variables (old nvars %d, new nvars = %d)\n", nredvars, oldnbdchginfos, conflictset->nbdchginfos);
      SCIPsetDebugMsg(set, " -> conflict set led to %d global bound changes %s(cdpt:%d, fdpt:%d, confdpt:%d, len:%d):\n",
         nbdchgs, redundant ? "(conflict became redundant) " : "", SCIPtreeGetCurrentDepth(tree), SCIPtreeGetFocusDepth(tree),
         conflictset->conflictdepth, conflictset->nbdchginfos);
      conflictsetPrint(conflictset);
#endif

      SCIPclockStop(conflict->dIBclock, set);

      if( redundant )
      {
         if( nbdchgs > 0 )
            *success = TRUE;

         return SCIP_OKAY;
      }
   }

   /* in case the conflict set contains only one bound change which is globally valid we apply that bound change
    * directly (except if we are in strong branching or diving - in this case a bound change would yield an unflushed LP
    * and is not handled when restoring the information)
    *
    * @note A bound change can only be applied if it is are related to the active node or if is a global bound
    *       change. Bound changes which are related to any other node cannot be handled at point due to the internal
    *       data structure
    */
   if( conflictset->nbdchginfos == 1 && insertdepth == 0 && !lp->strongbranching && !lp->diving )
   {
      SCIP_VAR* var;
      SCIP_Real bound;
      SCIP_BOUNDTYPE boundtype;

      var = conflictset->bdchginfos[0]->var;
      assert(var != NULL);

      boundtype = SCIPboundtypeOpposite((SCIP_BOUNDTYPE) conflictset->bdchginfos[0]->boundtype);
      bound = conflictset->relaxedbds[0];

      /* for continuous variables, we can only use the relaxed version of the bounds negation: !(x <= u) -> x >= u */
      if( SCIPvarIsIntegral(var) )
      {
         assert(SCIPsetIsIntegral(set, bound));
         bound += (boundtype == SCIP_BOUNDTYPE_LOWER ? +1.0 : -1.0);
      }

      SCIPsetDebugMsg(set, " -> apply global bound change: <%s> %s %g\n",
         SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", bound);

      SCIP_CALL( SCIPnodeAddBoundchg(tree->path[conflictset->validdepth], blkmem, set, stat, transprob, origprob, tree,
            reopt, lp, branchcand, eventqueue, cliquetable, var, bound, boundtype, FALSE) );

      *success = TRUE;
      SCIP_CALL( updateStatistics(conflict, blkmem, set, stat, conflictset, insertdepth) );
   }
   else
   {
      /* sort conflict handlers by priority */
      SCIPsetSortConflicthdlrs(set);

      /* call conflict handlers to create a conflict constraint */
      for( h = 0; h < set->nconflicthdlrs; ++h )
      {
         SCIP_RESULT result;

         assert(conflictset->conflicttype != SCIP_CONFTYPE_UNKNOWN);

         SCIP_CALL( SCIPconflicthdlrExec(set->conflicthdlrs[h], set, tree->path[insertdepth],
               tree->path[conflictset->validdepth], conflictset->bdchginfos, conflictset->relaxedbds,
               conflictset->nbdchginfos, conflictset->conflicttype, conflictset->usescutoffbound, *success, &result) );
         if( result == SCIP_CONSADDED )
         {
            *success = TRUE;
            SCIP_CALL( updateStatistics(conflict, blkmem, set, stat, conflictset, insertdepth) );
         }

         SCIPsetDebugMsg(set, " -> call conflict handler <%s> (prio=%d) to create conflict set with %d bounds returned result %d\n",
            SCIPconflicthdlrGetName(set->conflicthdlrs[h]), SCIPconflicthdlrGetPriority(set->conflicthdlrs[h]),
            conflictset->nbdchginfos, result);
      }
   }

   return SCIP_OKAY;
}

/** adds the collected conflict constraints to the corresponding nodes; the best set->conf_maxconss conflict constraints
 *  are added to the node of their validdepth; additionally (if not yet added, and if repropagation is activated), the
 *  conflict constraint that triggers the earliest repropagation is added to the node of its validdepth
 */
SCIP_RETCODE SCIPconflictFlushConss(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);

   /* is there anything to do? */
   if( conflict->nconflictsets > 0 )
   {
      SCIP_CONFLICTSET* repropconflictset;
      int nconflictsetsused;
      int focusdepth;
#ifndef NDEBUG
      int currentdepth;
#endif
      int cutoffdepth;
      int repropdepth;
      int maxconflictsets;
      int maxsize;
      int i;

      /* calculate the maximal number of conflict sets to accept, and the maximal size of each accepted conflict set */
      maxconflictsets = (set->conf_maxconss == -1 ? INT_MAX : set->conf_maxconss);
      maxsize = conflictCalcMaxsize(set, transprob);

      focusdepth = SCIPtreeGetFocusDepth(tree);
#ifndef NDEBUG
      currentdepth = SCIPtreeGetCurrentDepth(tree);
      assert(focusdepth <= currentdepth);
      assert(currentdepth == tree->pathlen-1);
#endif

      SCIPsetDebugMsg(set, "flushing %d conflict sets at focus depth %d (maxconflictsets: %d, maxsize: %d)\n",
         conflict->nconflictsets, focusdepth, maxconflictsets, maxsize);

      /* mark the focus node to have produced conflict sets in the visualization output */
      SCIPvisualFoundConflict(stat->visual, stat, tree->path[focusdepth]);

      /* insert the conflict sets at the corresponding nodes */
      nconflictsetsused = 0;
      cutoffdepth = INT_MAX;
      repropdepth = INT_MAX;
      repropconflictset = NULL;
      for( i = 0; i < conflict->nconflictsets && nconflictsetsused < maxconflictsets; ++i )
      {
         SCIP_CONFLICTSET* conflictset;

         conflictset = conflict->conflictsets[i];
         assert(conflictset != NULL);
         assert(0 <= conflictset->validdepth);
         assert(conflictset->validdepth <= conflictset->insertdepth);
         assert(conflictset->insertdepth <= focusdepth);
         assert(conflictset->insertdepth <= conflictset->repropdepth);
         assert(conflictset->repropdepth <= currentdepth || conflictset->repropdepth == INT_MAX); /* INT_MAX for dive/probing/strong */
         assert(conflictset->conflictdepth <= currentdepth || conflictset->conflictdepth == INT_MAX); /* INT_MAX for dive/probing/strong */

         /* ignore conflict sets that are only valid at a node that was already cut off */
         if( conflictset->insertdepth >= cutoffdepth )
         {
            SCIPsetDebugMsg(set, " -> ignoring conflict set with insertdepth %d >= cutoffdepth %d\n",
               conflictset->validdepth, cutoffdepth);
            continue;
         }

         /* if no conflict bounds exist, the node and its sub tree in the conflict set's valid depth can be
          * cut off completely
          */
         if( conflictset->nbdchginfos == 0 )
         {
            SCIPsetDebugMsg(set, " -> empty conflict set in depth %d cuts off sub tree at depth %d\n",
               focusdepth, conflictset->validdepth);

            SCIP_CALL( SCIPnodeCutoff(tree->path[conflictset->validdepth], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
            cutoffdepth = conflictset->validdepth;
            continue;
         }

         /* if the conflict set is too long, use the conflict set only if it decreases the repropagation depth */
         if( conflictset->nbdchginfos > maxsize )
         {
            SCIPsetDebugMsg(set, " -> conflict set is too long: %d > %d literals\n", conflictset->nbdchginfos, maxsize);
            if( set->conf_keepreprop && conflictset->repropagate && conflictset->repropdepth < repropdepth )
            {
               repropdepth = conflictset->repropdepth;
               repropconflictset = conflictset;
            }
         }
         else
         {
            SCIP_Bool success;

            /* call conflict handlers to create a conflict constraint */
            SCIP_CALL( conflictAddConflictCons(conflict, blkmem, set, stat, transprob, origprob, tree, reopt, lp, \
                  branchcand, eventqueue, cliquetable, conflictset, conflictset->insertdepth, &success) );

            /* if no conflict bounds exist, the node and its sub tree in the conflict set's valid depth can be
             * cut off completely
             */
            if( conflictset->nbdchginfos == 0 )
            {
               assert(!success);

               SCIPsetDebugMsg(set, " -> empty conflict set in depth %d cuts off sub tree at depth %d\n",
                  focusdepth, conflictset->validdepth);

               SCIP_CALL( SCIPnodeCutoff(tree->path[conflictset->validdepth], set, stat, tree, transprob, origprob, \
                     reopt, lp, blkmem) );
               cutoffdepth = conflictset->validdepth;
               continue;
            }

            if( success )
            {
               SCIPsetDebugMsg(set, " -> conflict set %d/%d added (cdpt:%d, fdpt:%d, insert:%d, valid:%d, conf:%d, reprop:%d, len:%d):\n",
                  nconflictsetsused+1, maxconflictsets, SCIPtreeGetCurrentDepth(tree), SCIPtreeGetFocusDepth(tree),
                  conflictset->insertdepth, conflictset->validdepth, conflictset->conflictdepth, conflictset->repropdepth,
                  conflictset->nbdchginfos);
               SCIPdebug(conflictsetPrint(conflictset));

               if( conflictset->repropagate && conflictset->repropdepth <= repropdepth )
               {
                  repropdepth = conflictset->repropdepth;
                  repropconflictset = NULL;
               }
               nconflictsetsused++;
            }
         }
      }

      /* reactivate propagation on the first node where one of the new conflict sets trigger a deduction */
      if( set->conf_repropagate && repropdepth < cutoffdepth && repropdepth < tree->pathlen )
      {
         assert(0 <= repropdepth && repropdepth < tree->pathlen);
         assert((int) tree->path[repropdepth]->depth == repropdepth);

         /* if the conflict constraint of smallest repropagation depth was not yet added, insert it now */
         if( repropconflictset != NULL )
         {
            SCIP_Bool success;

            assert(repropconflictset->repropagate);
            assert(repropconflictset->repropdepth == repropdepth);

            SCIP_CALL( conflictAddConflictCons(conflict, blkmem, set, stat, transprob, origprob, tree, reopt, lp, \
                  branchcand, eventqueue, cliquetable, repropconflictset, repropdepth, &success) );

            /* if no conflict bounds exist, the node and its sub tree in the conflict set's valid depth can be
             * cut off completely
             */
            if( repropconflictset->nbdchginfos == 0 )
            {
               assert(!success);

               SCIPsetDebugMsg(set, " -> empty reprop conflict set in depth %d cuts off sub tree at depth %d\n",
                  focusdepth, repropconflictset->validdepth);

               SCIP_CALL( SCIPnodeCutoff(tree->path[repropconflictset->validdepth], set, stat, tree, transprob, \
                     origprob, reopt, lp, blkmem) );
            }

#ifdef SCIP_DEBUG
            if( success )
            {
               SCIPsetDebugMsg(set, " -> additional reprop conflict set added (cdpt:%d, fdpt:%d, insert:%d, valid:%d, conf:%d, reprop:%d, len:%d):\n",
                  SCIPtreeGetCurrentDepth(tree), SCIPtreeGetFocusDepth(tree),
                  repropconflictset->insertdepth, repropconflictset->validdepth, repropconflictset->conflictdepth,
                  repropconflictset->repropdepth, repropconflictset->nbdchginfos);
               SCIPdebug(conflictsetPrint(repropconflictset));
            }
#endif
         }

         /* mark the node in the repropdepth to be propagated again */
         SCIPnodePropagateAgain(tree->path[repropdepth], set, stat, tree);

         SCIPsetDebugMsg(set, "marked node %p in depth %d to be repropagated due to conflicts found in depth %d\n",
            (void*)tree->path[repropdepth], repropdepth, focusdepth);
      }

      /* free the conflict store */
      for( i = 0; i < conflict->nconflictsets; ++i )
      {
         conflictsetFree(&conflict->conflictsets[i], blkmem);
      }
      conflict->nconflictsets = 0;
   }

   /* free all temporarily created bound change information data */
   conflictFreeTmpBdchginfos(conflict, blkmem);

   return SCIP_OKAY;
}

/** returns the current number of conflict sets in the conflict set storage */
int SCIPconflictGetNConflicts(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nconflictsets;
}

/** returns the total number of conflict constraints that were added to the problem */
SCIP_Longint SCIPconflictGetNAppliedConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbconss + conflict->nappliedlocconss;
}

/** returns the total number of literals in conflict constraints that were added to the problem */
SCIP_Longint SCIPconflictGetNAppliedLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbliterals + conflict->nappliedlocliterals;
}

/** returns the total number of global bound changes applied by the conflict analysis */
SCIP_Longint SCIPconflictGetNGlobalChgBds(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nglbchgbds;
}

/** returns the total number of conflict constraints that were added globally to the problem */
SCIP_Longint SCIPconflictGetNAppliedGlobalConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbconss;
}

/** returns the total number of literals in conflict constraints that were added globally to the problem */
SCIP_Longint SCIPconflictGetNAppliedGlobalLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbliterals;
}

/** returns the total number of local bound changes applied by the conflict analysis */
SCIP_Longint SCIPconflictGetNLocalChgBds(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlocchgbds;
}

/** returns the total number of conflict constraints that were added locally to the problem */
SCIP_Longint SCIPconflictGetNAppliedLocalConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedlocconss;
}

/** returns the total number of literals in conflict constraints that were added locally to the problem */
SCIP_Longint SCIPconflictGetNAppliedLocalLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedlocliterals;
}




/*
 * Propagation Conflict Analysis
 */

/** returns whether bound change has a valid reason that can be resolved in conflict analysis */
static
SCIP_Bool bdchginfoIsResolvable(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   return (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
      || (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_PROPINFER
         && SCIPbdchginfoGetInferProp(bdchginfo) != NULL));
}

/** compares two conflict set entries, such that bound changes infered later are
 *  ordered prior to ones that were infered earlier
 */
static
SCIP_DECL_SORTPTRCOMP(conflictBdchginfoComp)
{  /*lint --e{715}*/
   SCIP_BDCHGINFO* bdchginfo1;
   SCIP_BDCHGINFO* bdchginfo2;

   bdchginfo1 = (SCIP_BDCHGINFO*)elem1;
   bdchginfo2 = (SCIP_BDCHGINFO*)elem2;
   assert(bdchginfo1 != NULL);
   assert(bdchginfo2 != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo1));
   assert(!SCIPbdchginfoIsRedundant(bdchginfo2));

   if( !SCIPbdchgidxIsEarlierNonNull(SCIPbdchginfoGetIdx(bdchginfo1), SCIPbdchginfoGetIdx(bdchginfo2)) )
      return -1;
   else
      return +1;
}

/** return TRUE if conflict analysis is applicable; In case the function return FALSE there is no need to initialize the
 *  conflict analysis since it will not be applied
 */
SCIP_Bool SCIPconflictApplicable(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   /* check, if propagation conflict analysis is enabled */
   if( !set->conf_enable || !set->conf_useprop )
      return FALSE;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return FALSE;

   return TRUE;
}

/** creates conflict analysis data for propagation conflicts */
SCIP_RETCODE SCIPconflictCreate(
   SCIP_CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflict != NULL);

   SCIP_ALLOC( BMSallocMemory(conflict) );

   SCIP_CALL( SCIPclockCreate(&(*conflict)->dIBclock, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->propanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->inflpanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->boundlpanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->sbanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->pseudoanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );

   /* enable or disable timing depending on the parameter statistic timing */
   SCIPconflictEnableOrDisableClocks((*conflict), set->time_statistictiming);

   SCIP_CALL( SCIPpqueueCreate(&(*conflict)->bdchgqueue, set->mem_arraygrowinit, set->mem_arraygrowfac,
         conflictBdchginfoComp) );
   SCIP_CALL( SCIPpqueueCreate(&(*conflict)->forcedbdchgqueue, set->mem_arraygrowinit, set->mem_arraygrowfac,
         conflictBdchginfoComp) );
   SCIP_CALL( conflictsetCreate(&(*conflict)->conflictset, blkmem) );
   (*conflict)->conflictsets = NULL;
   (*conflict)->conflictsetscores = NULL;
   (*conflict)->tmpbdchginfos = NULL;
   (*conflict)->conflictsetssize = 0;
   (*conflict)->nconflictsets = 0;
   (*conflict)->proofsets = NULL;
   (*conflict)->proofsetssize = 0;
   (*conflict)->nproofsets = 0;
   (*conflict)->tmpbdchginfossize = 0;
   (*conflict)->ntmpbdchginfos = 0;
   (*conflict)->count = 0;
   (*conflict)->nglbchgbds = 0;
   (*conflict)->nappliedglbconss = 0;
   (*conflict)->nappliedglbliterals = 0;
   (*conflict)->nlocchgbds = 0;
   (*conflict)->nappliedlocconss = 0;
   (*conflict)->nappliedlocliterals = 0;
   (*conflict)->npropcalls = 0;
   (*conflict)->npropsuccess = 0;
   (*conflict)->npropconfconss = 0;
   (*conflict)->npropconfliterals = 0;
   (*conflict)->npropreconvconss = 0;
   (*conflict)->npropreconvliterals = 0;
   (*conflict)->ninflpcalls = 0;
   (*conflict)->ninflpsuccess = 0;
   (*conflict)->ninflpconfconss = 0;
   (*conflict)->ninflpconfliterals = 0;
   (*conflict)->ninflpreconvconss = 0;
   (*conflict)->ninflpreconvliterals = 0;
   (*conflict)->ninflpiterations = 0;
   (*conflict)->nboundlpcalls = 0;
   (*conflict)->nboundlpsuccess = 0;
   (*conflict)->nboundlpconfconss = 0;
   (*conflict)->nboundlpconfliterals = 0;
   (*conflict)->nboundlpreconvconss = 0;
   (*conflict)->nboundlpreconvliterals = 0;
   (*conflict)->nboundlpiterations = 0;
   (*conflict)->nsbcalls = 0;
   (*conflict)->nsbsuccess = 0;
   (*conflict)->nsbconfconss = 0;
   (*conflict)->nsbconfliterals = 0;
   (*conflict)->nsbreconvconss = 0;
   (*conflict)->nsbreconvliterals = 0;
   (*conflict)->nsbiterations = 0;
   (*conflict)->npseudocalls = 0;
   (*conflict)->npseudosuccess = 0;
   (*conflict)->npseudoconfconss = 0;
   (*conflict)->npseudoconfliterals = 0;
   (*conflict)->npseudoreconvconss = 0;
   (*conflict)->npseudoreconvliterals = 0;
   (*conflict)->ndualrayinfglobal = 0;
   (*conflict)->ndualrayinfsuccess = 0;
   (*conflict)->dualrayinfnnonzeros = 0;
   (*conflict)->ndualraybndglobal = 0;
   (*conflict)->ndualraybndsuccess = 0;
   (*conflict)->dualraybndnnonzeros = 0;

   SCIP_CALL( conflictInitProofset((*conflict), blkmem) );

   return SCIP_OKAY;
}

/** frees conflict analysis data for propagation conflicts */
SCIP_RETCODE SCIPconflictFree(
   SCIP_CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(conflict != NULL);
   assert(*conflict != NULL);
   assert((*conflict)->nconflictsets == 0);
   assert((*conflict)->ntmpbdchginfos == 0);

#ifdef SCIP_CONFGRAPH
   confgraphFree();
#endif

   SCIPclockFree(&(*conflict)->dIBclock);
   SCIPclockFree(&(*conflict)->propanalyzetime);
   SCIPclockFree(&(*conflict)->inflpanalyzetime);
   SCIPclockFree(&(*conflict)->boundlpanalyzetime);
   SCIPclockFree(&(*conflict)->sbanalyzetime);
   SCIPclockFree(&(*conflict)->pseudoanalyzetime);
   SCIPpqueueFree(&(*conflict)->bdchgqueue);
   SCIPpqueueFree(&(*conflict)->forcedbdchgqueue);
   conflictsetFree(&(*conflict)->conflictset, blkmem);
   proofsetFree(&(*conflict)->proofset, blkmem);

   BMSfreeMemoryArrayNull(&(*conflict)->conflictsets);
   BMSfreeMemoryArrayNull(&(*conflict)->conflictsetscores);
   BMSfreeMemoryArrayNull(&(*conflict)->proofsets);
   BMSfreeMemoryArrayNull(&(*conflict)->tmpbdchginfos);
   BMSfreeMemory(conflict);

   return SCIP_OKAY;
}

/** clears the conflict queue and the current conflict set */
static
void conflictClear(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   SCIPpqueueClear(conflict->bdchgqueue);
   SCIPpqueueClear(conflict->forcedbdchgqueue);
   conflictsetClear(conflict->conflictset);
}

/** initializes the propagation conflict analysis by clearing the conflict candidate queue */
SCIP_RETCODE SCIPconflictInit(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_CONFTYPE         conftype,           /**< type of the conflict */
   SCIP_Bool             usescutoffbound     /**< depends the conflict on a cutoff bound? */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);

   SCIPsetDebugMsg(set, "initializing conflict analysis\n");

   /* clear the conflict candidate queue and the conflict set */
   conflictClear(conflict);

   /* set conflict type */
   assert(conftype == SCIP_CONFTYPE_BNDEXCEEDING || conftype == SCIP_CONFTYPE_INFEASLP
       || conftype == SCIP_CONFTYPE_PROPAGATION);
   conflict->conflictset->conflicttype = conftype;

   /* set whether a cutoff bound is involved */
   conflict->conflictset->usescutoffbound = usescutoffbound;

   /* increase the conflict counter, such that binary variables of new conflict set and new conflict queue are labeled
    * with this new counter
    */
   conflict->count++;
   if( conflict->count == 0 ) /* make sure, 0 is not a valid conflict counter (may happen due to integer overflow) */
      conflict->count = 1;

   /* increase the conflict score weight for history updates of future conflict reasons */
   if( stat->nnodes > stat->lastconflictnode )
   {
      assert(0.0 < set->conf_scorefac && set->conf_scorefac <= 1.0);
      stat->vsidsweight /= set->conf_scorefac;
      assert(stat->vsidsweight > 0.0);

      /* if the conflict score for the next conflict exceeds 1000.0, rescale all history conflict scores */
      if( stat->vsidsweight >= 1000.0 )
      {
         int v;

         for( v = 0; v < prob->nvars; ++v )
         {
            SCIP_CALL( SCIPvarScaleVSIDS(prob->vars[v], 1.0/stat->vsidsweight) );
         }
         SCIPhistoryScaleVSIDS(stat->glbhistory, 1.0/stat->vsidsweight);
         SCIPhistoryScaleVSIDS(stat->glbhistorycrun, 1.0/stat->vsidsweight);
         stat->vsidsweight = 1.0;
      }
      stat->lastconflictnode = stat->nnodes;
   }

#ifdef SCIP_CONFGRAPH
   confgraphFree();
   SCIP_CALL( confgraphCreate(set, conflict) );
#endif

   return SCIP_OKAY;
}

/** marks bound to be present in the current conflict and returns whether a bound which is at least as tight was already
 *  member of the current conflict (i.e., the given bound change does not need to be added)
 */
static
SCIP_Bool conflictMarkBoundCheckPresence(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change to add to the conflict set */
   SCIP_Real             relaxedbd           /**< relaxed bound */
   )
{
   SCIP_VAR* var;
   SCIP_Real newbound;

   assert(conflict != NULL);

   var = SCIPbdchginfoGetVar(bdchginfo);
   newbound = SCIPbdchginfoGetNewbound(bdchginfo);
   assert(var != NULL);

   switch( SCIPbdchginfoGetBoundtype(bdchginfo) )
   {
   case SCIP_BOUNDTYPE_LOWER:
      /* check if the variables lower bound is already member of the conflict */
      if( var->conflictlbcount == conflict->count )
      {
         /* the variable is already member of the conflict; hence check if the new bound is redundant */
         if( var->conflictlb > newbound )
         {
            SCIPsetDebugMsg(set, "ignoring redundant bound change <%s> >= %g since a stronger lower bound exist <%s> >= %g\n",
               SCIPvarGetName(var), newbound, SCIPvarGetName(var), var->conflictlb);
            return TRUE;
         }
         else if( var->conflictlb == newbound ) /*lint !e777*/
         {
            SCIPsetDebugMsg(set, "ignoring redundant bound change <%s> >= %g since this lower bound is already present\n", SCIPvarGetName(var), newbound);
            SCIPsetDebugMsg(set, "adjust relaxed lower bound <%g> -> <%g>\n", var->conflictlb, relaxedbd);
            var->conflictrelaxedlb = MAX(var->conflictrelaxedlb, relaxedbd);
            return TRUE;
         }
      }

      /* add the variable lower bound to the current conflict */
      var->conflictlbcount = conflict->count;

      /* remember the lower bound and relaxed bound to allow only better/tighter lower bounds for that variables
       * w.r.t. this conflict
       */
      var->conflictlb = newbound;
      var->conflictrelaxedlb = relaxedbd;

      return FALSE;

   case SCIP_BOUNDTYPE_UPPER:
      /* check if the variables upper bound is already member of the conflict */
      if( var->conflictubcount == conflict->count )
      {
         /* the variable is already member of the conflict; hence check if the new bound is redundant */
         if( var->conflictub < newbound )
         {
            SCIPsetDebugMsg(set, "ignoring redundant bound change <%s> <= %g since a stronger upper bound exist <%s> <= %g\n",
               SCIPvarGetName(var), newbound, SCIPvarGetName(var), var->conflictub);
            return TRUE;
         }
         else if( var->conflictub == newbound ) /*lint !e777*/
         {
            SCIPsetDebugMsg(set, "ignoring redundant bound change <%s> <= %g since this upper bound is already present\n", SCIPvarGetName(var), newbound);
            SCIPsetDebugMsg(set, "adjust relaxed upper bound <%g> -> <%g>\n", var->conflictub, relaxedbd);
            var->conflictrelaxedub = MIN(var->conflictrelaxedub, relaxedbd);
            return TRUE;
         }
      }

      /* add the variable upper bound to the current conflict */
      var->conflictubcount = conflict->count;

      /* remember the upper bound and relaxed bound to allow only better/tighter upper bounds for that variables
       * w.r.t. this conflict
       */
      var->conflictub = newbound;
      var->conflictrelaxedub = relaxedbd;

      return FALSE;

   default:
      SCIPerrorMessage("invalid bound type %d\n", SCIPbdchginfoGetBoundtype(bdchginfo));
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }
}

/** puts bound change into the current conflict set */
static
SCIP_RETCODE conflictAddConflictBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change to add to the conflict set */
   SCIP_Real             relaxedbd           /**< relaxed bound */
   )
{
   assert(conflict != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   /* check if the relaxed bound is really a relaxed bound */
   assert(SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER || SCIPsetIsGE(set, relaxedbd, SCIPbdchginfoGetNewbound(bdchginfo)));
   assert(SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER || SCIPsetIsLE(set, relaxedbd, SCIPbdchginfoGetNewbound(bdchginfo)));

   SCIPsetDebugMsg(set, "putting bound change <%s> %s %g(%g) at depth %d to current conflict set\n",
      SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
      SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", SCIPbdchginfoGetNewbound(bdchginfo),
      relaxedbd, SCIPbdchginfoGetDepth(bdchginfo));

   /* mark the bound to be member of the conflict and check if a bound which is at least as tight is already member of
    * the conflict
    */
   if( !conflictMarkBoundCheckPresence(conflict, set, bdchginfo, relaxedbd) )
   {
      /* add the bound change to the current conflict set */
      SCIP_CALL( conflictsetAddBound(conflict->conflictset, blkmem, set, bdchginfo, relaxedbd) );

#ifdef SCIP_CONFGRAPH
      if( bdchginfo != confgraphcurrentbdchginfo )
         confgraphAddBdchg(bdchginfo);
#endif
   }
#ifdef SCIP_CONFGRAPH
   else
      confgraphLinkBdchg(bdchginfo);
#endif

   return SCIP_OKAY;
}

/** returns whether the negation of the given bound change would lead to a globally valid literal */
static
SCIP_Bool isBoundchgUseless(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   SCIP_VAR* var;
   SCIP_BOUNDTYPE boundtype;
   SCIP_Real bound;

   var = SCIPbdchginfoGetVar(bdchginfo);
   boundtype = SCIPbdchginfoGetBoundtype(bdchginfo);
   bound = SCIPbdchginfoGetNewbound(bdchginfo);

   return (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
      && ((boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasGE(set, bound, SCIPvarGetUbGlobal(var)))
         || (boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasLE(set, bound, SCIPvarGetLbGlobal(var)))));
}

/** adds given bound change information to the conflict candidate queue */
static
SCIP_RETCODE conflictQueueBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change information */
   SCIP_Real             relaxedbd           /**< relaxed bound */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);
   assert(bdchginfo != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   /* check if the relaxed bound is really a relaxed bound */
   assert(SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER || SCIPsetIsGE(set, relaxedbd, SCIPbdchginfoGetNewbound(bdchginfo)));
   assert(SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER || SCIPsetIsLE(set, relaxedbd, SCIPbdchginfoGetNewbound(bdchginfo)));

   /* mark the bound to be member of the conflict and check if a bound which is at least as tight is already member of
    * the conflict
    */
   if( !conflictMarkBoundCheckPresence(conflict, set, bdchginfo, relaxedbd) )
   {
      /* insert the bound change into the conflict queue */
      if( (!set->conf_preferbinary || SCIPvarIsBinary(SCIPbdchginfoGetVar(bdchginfo)))
         && !isBoundchgUseless(set, bdchginfo) )
      {
         SCIP_CALL( SCIPpqueueInsert(conflict->bdchgqueue, (void*)bdchginfo) );
      }
      else
      {
         SCIP_CALL( SCIPpqueueInsert(conflict->forcedbdchgqueue, (void*)bdchginfo) );
      }

#ifdef SCIP_CONFGRAPH
      confgraphAddBdchg(bdchginfo);
#endif
   }
#ifdef SCIP_CONFGRAPH
   else
      confgraphLinkBdchg(bdchginfo);
#endif

   return SCIP_OKAY;
}

/** convert variable and bound change to active variable */
static
SCIP_RETCODE convertToActiveVar(
   SCIP_VAR**            var,                /**< pointer to variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BOUNDTYPE*       boundtype,          /**< pointer to type of bound that was changed: lower or upper bound */
   SCIP_Real*            bound               /**< pointer to bound to convert, or NULL */
   )
{
   SCIP_Real scalar;
   SCIP_Real constant;

   scalar = 1.0;
   constant = 0.0;

   /* transform given varibale to active varibale */
   SCIP_CALL( SCIPvarGetProbvarSum(var, set, &scalar, &constant) );
   assert(SCIPvarGetStatus(*var) == SCIP_VARSTATUS_FIXED || scalar != 0.0); /*lint !e777*/

   if( SCIPvarGetStatus(*var) == SCIP_VARSTATUS_FIXED )
      return SCIP_OKAY;

   /* if the scalar of the aggregation is negative, we have to switch the bound type */
   if( scalar < 0.0 )
      (*boundtype) = SCIPboundtypeOpposite(*boundtype);

   if( bound != NULL )
   {
      (*bound) -= constant;
      (*bound) /= scalar;
   }

   return SCIP_OKAY;
}

/** adds variable's bound to conflict candidate queue */
static
SCIP_RETCODE conflictAddBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound that was changed: lower or upper bound */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change info, or NULL */
   SCIP_Real             relaxedbd           /**< relaxed bound */
   )
{
   assert(SCIPvarIsActive(var));
   assert(bdchginfo != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   SCIPsetDebugMsg(set, " -> adding bound <%s> %s %.15g(%.15g) [status:%d, type:%d, depth:%d, pos:%d, reason:<%s>, info:%d] to candidates\n",
      SCIPvarGetName(var),
      boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
      SCIPbdchginfoGetNewbound(bdchginfo), relaxedbd,
      SCIPvarGetStatus(var), SCIPvarGetType(var),
      SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
      SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
      : (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
         ? SCIPconsGetName(SCIPbdchginfoGetInferCons(bdchginfo))
         : (SCIPbdchginfoGetInferProp(bdchginfo) != NULL ? SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo))
            : "none")),
      SCIPbdchginfoGetChgtype(bdchginfo) != SCIP_BOUNDCHGTYPE_BRANCHING ? SCIPbdchginfoGetInferInfo(bdchginfo) : -1);

   /* the local bound change may be resolved and has to be put on the candidate queue;
    * we even put bound changes without inference information on the queue in order to automatically
    * eliminate multiple insertions of the same bound change
    */
   assert(SCIPbdchginfoGetVar(bdchginfo) == var);
   assert(SCIPbdchginfoGetBoundtype(bdchginfo) == boundtype);
   assert(SCIPbdchginfoGetDepth(bdchginfo) >= 0);
   assert(SCIPbdchginfoGetPos(bdchginfo) >= 0);

   /* the relaxed bound should be a relaxation */
   assert(boundtype == SCIP_BOUNDTYPE_LOWER ? SCIPsetIsLE(set, relaxedbd, SCIPbdchginfoGetNewbound(bdchginfo)) : SCIPsetIsGE(set, relaxedbd, SCIPbdchginfoGetNewbound(bdchginfo)));

   /* the relaxed bound should be worse then the old bound of the bound change info */
   assert(boundtype == SCIP_BOUNDTYPE_LOWER ? SCIPsetIsGT(set, relaxedbd, SCIPbdchginfoGetOldbound(bdchginfo)) : SCIPsetIsLT(set, relaxedbd, SCIPbdchginfoGetOldbound(bdchginfo)));

   /* put bound change information into priority queue */
   SCIP_CALL( conflictQueueBound(conflict, set, bdchginfo, relaxedbd) );

   /* each variable which is add to the conflict graph gets an increase in the VSIDS
    *
    * @note That is different to the VSIDS preseted in the literature
    */
   SCIP_CALL( incVSIDS(var, blkmem, set, stat, boundtype, relaxedbd, set->conf_conflictgraphweight) );

   return SCIP_OKAY;
}

/** adds variable's bound to conflict candidate queue */
SCIP_RETCODE SCIPconflictAddBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound that was changed: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{
   SCIP_BDCHGINFO* bdchginfo;

   assert(conflict != NULL);
   assert(stat != NULL);
   assert(var != NULL);

   /* convert bound to active problem variable */
   SCIP_CALL( convertToActiveVar(&var, set, &boundtype, NULL) );

   /* we can ignore fixed variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
      return SCIP_OKAY;

   /* if the variable is multi-aggregated, add the bounds of all aggregation variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_VAR** vars;
      SCIP_Real* scalars;
      int nvars;
      int i;

      vars = SCIPvarGetMultaggrVars(var);
      scalars = SCIPvarGetMultaggrScalars(var);
      nvars = SCIPvarGetMultaggrNVars(var);
      for( i = 0; i < nvars; ++i )
      {
         SCIP_CALL( SCIPconflictAddBound(conflict, blkmem, set, stat, vars[i],
               (scalars[i] < 0.0 ? SCIPboundtypeOpposite(boundtype) : boundtype), bdchgidx) );
      }

      return SCIP_OKAY;
   }
   assert(SCIPvarIsActive(var));

   /* get bound change information */
   bdchginfo = SCIPvarGetBdchgInfo(var, boundtype, bdchgidx, FALSE);

   /* if bound of variable was not changed (this means it is still the global bound), we can ignore the conflicting
    * bound
    */
   if( bdchginfo == NULL )
      return SCIP_OKAY;

   assert(SCIPbdchgidxIsEarlier(SCIPbdchginfoGetIdx(bdchginfo), bdchgidx));

   SCIP_CALL( conflictAddBound(conflict, blkmem, set, stat, var, boundtype, bdchginfo, SCIPbdchginfoGetNewbound(bdchginfo)) );

   return SCIP_OKAY;
}

/** adds variable's bound to conflict candidate queue */
SCIP_RETCODE SCIPconflictAddRelaxedBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound that was changed: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_Real             relaxedbd           /**< the relaxed bound */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   int nbdchgs;

   assert(conflict != NULL);
   assert(stat != NULL);
   assert(var != NULL);

   if( !SCIPvarIsActive(var) )
   {
      /* convert bound to active problem variable */
      SCIP_CALL( convertToActiveVar(&var, set, &boundtype, &relaxedbd) );

      /* we can ignore fixed variables */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
         return SCIP_OKAY;

      /* if the variable is multi-aggregated, add the bounds of all aggregation variables */
      if(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      {
         SCIPsetDebugMsg(set, "ignoring relaxed bound information since variable <%s> is multi-aggregated active\n", SCIPvarGetName(var));

         SCIP_CALL( SCIPconflictAddBound(conflict, blkmem, set, stat, var, boundtype, bdchgidx) );

         return SCIP_OKAY;
      }
   }
   assert(SCIPvarIsActive(var));

   /* get bound change information */
   bdchginfo = SCIPvarGetBdchgInfo(var, boundtype, bdchgidx, FALSE);

   /* if bound of variable was not changed (this means it is still the global bound), we can ignore the conflicting
    * bound
    */
   if( bdchginfo == NULL )
      return SCIP_OKAY;

   /* check that the bound change info is not a temporary one */
   assert(SCIPbdchgidxGetPos(&bdchginfo->bdchgidx) >= 0);

   /* get the position of the bound change information within the bound change array of the variable */
   nbdchgs = (int) bdchginfo->pos;
   assert(nbdchgs >= 0);

   /* if the relaxed bound should be ignored, set the relaxed bound to the bound given by the bdchgidx; that ensures
    * that the loop(s) below will be skipped
    */
   if( set->conf_ignorerelaxedbd )
      relaxedbd = SCIPbdchginfoGetNewbound(bdchginfo);

   /* search for the bound change information which includes the relaxed bound */
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_Real newbound;

      /* adjust relaxed lower bound w.r.t. variable type */
      SCIPvarAdjustLb(var, set, &relaxedbd);

      /* due to numericis we compare the relaxed lower bound to the one present at the particular time point and take
       * the better one
       */
      newbound = SCIPbdchginfoGetNewbound(bdchginfo);
      relaxedbd = MIN(relaxedbd, newbound);

      /* check if relaxed lower bound is smaller or equal to global lower bound; if so we can ignore the conflicting
       * bound
       */
      if( SCIPsetIsLE(set, relaxedbd, SCIPvarGetLbGlobal(var)) )
         return SCIP_OKAY;

      while( nbdchgs > 0 )
      {
         assert(SCIPsetIsLE(set, relaxedbd, SCIPbdchginfoGetNewbound(bdchginfo)));

         /* check if the old lower bound is greater than or equal to relaxed lower bound; if not we found the bound
          * change info which we need to report
          */
         if( SCIPsetIsGT(set, relaxedbd, SCIPbdchginfoGetOldbound(bdchginfo)) )
            break;

         bdchginfo = SCIPvarGetBdchgInfoLb(var, nbdchgs-1);

         SCIPsetDebugMsg(set, "lower bound change %d oldbd=%.15g, newbd=%.15g, depth=%d, pos=%d, redundant=%u\n",
            nbdchgs, SCIPbdchginfoGetOldbound(bdchginfo), SCIPbdchginfoGetNewbound(bdchginfo),
            SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPbdchginfoIsRedundant(bdchginfo));

         /* if bound change is redundant (this means it now a global bound), we can ignore the conflicting bound */
         if( SCIPbdchginfoIsRedundant(bdchginfo) )
            return SCIP_OKAY;

         nbdchgs--;
      }
      assert(SCIPsetIsGT(set, relaxedbd, SCIPbdchginfoGetOldbound(bdchginfo)));
   }
   else
   {
      SCIP_Real newbound;

      assert(boundtype == SCIP_BOUNDTYPE_UPPER);

      /* adjust relaxed upper bound w.r.t. variable type */
      SCIPvarAdjustUb(var, set, &relaxedbd);

      /* due to numericis we compare the relaxed upper bound to the one present at the particular time point and take
       * the better one
       */
      newbound = SCIPbdchginfoGetNewbound(bdchginfo);
      relaxedbd = MAX(relaxedbd, newbound);

      /* check if relaxed upper bound is greater or equal to global upper bound; if so we can ignore the conflicting
       * bound
       */
      if( SCIPsetIsGE(set, relaxedbd, SCIPvarGetUbGlobal(var)) )
         return SCIP_OKAY;

      while( nbdchgs > 0 )
      {
         assert(SCIPsetIsGE(set, relaxedbd, SCIPbdchginfoGetNewbound(bdchginfo)));

         /* check if the old upper bound is smaller than or equal to the relaxed upper bound; if not we found the
          * bound change info which we need to report
          */
         if( SCIPsetIsLT(set, relaxedbd, SCIPbdchginfoGetOldbound(bdchginfo)) )
            break;

         bdchginfo = SCIPvarGetBdchgInfoUb(var, nbdchgs-1);

         SCIPsetDebugMsg(set, "upper bound change %d oldbd=%.15g, newbd=%.15g, depth=%d, pos=%d, redundant=%u\n",
            nbdchgs, SCIPbdchginfoGetOldbound(bdchginfo), SCIPbdchginfoGetNewbound(bdchginfo),
            SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPbdchginfoIsRedundant(bdchginfo));

         /* if bound change is redundant (this means it now a global bound), we can ignore the conflicting bound */
         if( SCIPbdchginfoIsRedundant(bdchginfo) )
            return SCIP_OKAY;

         nbdchgs--;
      }
      assert(SCIPsetIsLT(set, relaxedbd, SCIPbdchginfoGetOldbound(bdchginfo)));
   }

   assert(SCIPbdchgidxIsEarlier(SCIPbdchginfoGetIdx(bdchginfo), bdchgidx));

   /* put bound change information into priority queue */
   SCIP_CALL( conflictAddBound(conflict, blkmem, set, stat, var, boundtype, bdchginfo, relaxedbd) );

   return SCIP_OKAY;
}

/** checks if the given variable is already part of the current conflict set or queued for resolving with the same or
 *  even stronger bound
 */
SCIP_RETCODE SCIPconflictIsVarUsed(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound for which the score should be increased */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_Bool*            used                /**< pointer to store if the variable is already used */
   )
{
   SCIP_Real newbound;

   /* convert bound to active problem variable */
   SCIP_CALL( convertToActiveVar(&var, set, &boundtype, NULL) );

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      *used = FALSE;
   else
   {
      assert(SCIPvarIsActive(var));
      assert(var != NULL);

      switch( boundtype )
      {
      case SCIP_BOUNDTYPE_LOWER:

         newbound = SCIPgetVarLbAtIndex(set->scip, var, bdchgidx, FALSE);

         if( var->conflictlbcount == conflict->count && var->conflictlb >= newbound )
         {
            SCIPsetDebugMsg(set, "already queued bound change <%s> >= %g\n", SCIPvarGetName(var), newbound);
            *used = TRUE;
         }
         else
            *used = FALSE;
         break;
      case SCIP_BOUNDTYPE_UPPER:

         newbound = SCIPgetVarUbAtIndex(set->scip, var, bdchgidx, FALSE);

         if( var->conflictubcount == conflict->count && var->conflictub <= newbound )
         {
            SCIPsetDebugMsg(set, "already queued bound change <%s> <= %g\n", SCIPvarGetName(var), newbound);
            *used = TRUE;
         }
         else
            *used = FALSE;
         break;
      default:
         SCIPerrorMessage("invalid bound type %d\n", boundtype);
         SCIPABORT();
         *used = FALSE; /*lint !e527*/
      }
   }

   return SCIP_OKAY;
}

/** returns the conflict lower bound if the variable is present in the current conflict set; otherwise the global lower
 *  bound
 */
SCIP_Real SCIPconflictGetVarLb(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   if( var->conflictlbcount == conflict->count )
   {
      assert(EPSGE(var->conflictlb, var->conflictrelaxedlb, 1e-09));
      return var->conflictrelaxedlb;
   }

   return SCIPvarGetLbGlobal(var);
}

/** returns the conflict upper bound if the variable is present in the current conflict set; otherwise the global upper
 *  bound
 */
SCIP_Real SCIPconflictGetVarUb(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   if( var->conflictubcount == conflict->count )
   {
      assert(EPSLE(var->conflictub, var->conflictrelaxedub, 1e-09));
      return var->conflictrelaxedub;
   }

   return SCIPvarGetUbGlobal(var);
}

/** removes and returns next conflict analysis candidate from the candidate queue */
static
SCIP_BDCHGINFO* conflictRemoveCand(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_VAR* var;

   assert(conflict != NULL);

   if( SCIPpqueueNElems(conflict->forcedbdchgqueue) > 0 )
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueRemove(conflict->forcedbdchgqueue));
   else
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueRemove(conflict->bdchgqueue));

   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   /* if we have a candidate this one should be valid for the current conflict analysis */
   assert(!bdchginfoIsInvalid(conflict, bdchginfo));

   /* mark the bound change to be no longer in the conflict (it will be either added again to the conflict set or
    * replaced by resolving, which might add a weaker change on the same bound to the queue)
    */
   var = SCIPbdchginfoGetVar(bdchginfo);
   if( SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER )
   {
      var->conflictlbcount = 0;
      var->conflictrelaxedlb = SCIP_REAL_MIN;
   }
   else
   {
      assert(SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER);
      var->conflictubcount = 0;
      var->conflictrelaxedub = SCIP_REAL_MAX;
   }

#ifdef SCIP_CONFGRAPH
   confgraphSetCurrentBdchg(bdchginfo);
#endif

   return bdchginfo;
}

/** returns next conflict analysis candidate from the candidate queue without removing it */
static
SCIP_BDCHGINFO* conflictFirstCand(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   SCIP_BDCHGINFO* bdchginfo;

   assert(conflict != NULL);

   if( SCIPpqueueNElems(conflict->forcedbdchgqueue) > 0 )
   {
      /* get next potetioal candidate */
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueFirst(conflict->forcedbdchgqueue));

      /* check if this candidate is valid */
      if( bdchginfoIsInvalid(conflict, bdchginfo) )
      {
         SCIPdebugMessage("bound change info [%d:<%s> %s %g] is invaild -> pop it from the force queue\n", SCIPbdchginfoGetDepth(bdchginfo),
            SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo));

         /* pop the invalid bound change info from the queue */
         (void)(SCIPpqueueRemove(conflict->forcedbdchgqueue));

         /* call method recursively to get next conflict analysis candidate */
         bdchginfo = conflictFirstCand(conflict);
      }
   }
   else
   {
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueFirst(conflict->bdchgqueue));

      /* check if this candidate is valid */
      if( bdchginfo != NULL && bdchginfoIsInvalid(conflict, bdchginfo) )
      {
         SCIPdebugMessage("bound change info [%d:<%s> %s %g] is invaild -> pop it from the queue\n", SCIPbdchginfoGetDepth(bdchginfo),
            SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo));

         /* pop the invalid bound change info from the queue */
         (void)(SCIPpqueueRemove(conflict->bdchgqueue));

         /* call method recursively to get next conflict analysis candidate */
         bdchginfo = conflictFirstCand(conflict);
      }
   }
   assert(bdchginfo == NULL || !SCIPbdchginfoIsRedundant(bdchginfo));

   return bdchginfo;
}

/** adds the current conflict set (extended by all remaining bound changes in the queue) to the pool of conflict sets */
static
SCIP_RETCODE conflictAddConflictset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   validdepth,         /**< minimal depth level at which the conflict set is valid */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   SCIP_Bool             repropagate,        /**< should the constraint trigger a repropagation? */
   SCIP_Bool*            success,            /**< pointer to store whether the conflict set is valid */
   int*                  nliterals           /**< pointer to store the number of literals in the generated conflictset */
   )
{
   SCIP_CONFLICTSET* conflictset;
   SCIP_BDCHGINFO** bdchginfos;
   int nbdchginfos;
   int currentdepth;
   int focusdepth;

   assert(conflict != NULL);
   assert(conflict->conflictset != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(success != NULL);
   assert(nliterals != NULL);
   assert(SCIPpqueueNElems(conflict->forcedbdchgqueue) == 0);

   *success = FALSE;
   *nliterals = 0;

   /* check, whether local conflicts are allowed */
   validdepth = MAX(validdepth, conflict->conflictset->validdepth);
   if( !set->conf_allowlocal && validdepth > 0 )
      return SCIP_OKAY;

   focusdepth = SCIPtreeGetFocusDepth(tree);
   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(focusdepth <= currentdepth);
   assert(0 <= conflict->conflictset->validdepth && conflict->conflictset->validdepth <= currentdepth);
   assert(0 <= validdepth && validdepth <= currentdepth);

   /* get the elements of the bound change queue */
   bdchginfos = (SCIP_BDCHGINFO**)SCIPpqueueElems(conflict->bdchgqueue);
   nbdchginfos = SCIPpqueueNElems(conflict->bdchgqueue);

   /* create a copy of the current conflict set, allocating memory for the additional elements of the queue */
   SCIP_CALL( conflictsetCopy(&conflictset, blkmem, conflict->conflictset, nbdchginfos) );
   conflictset->validdepth = validdepth;
   conflictset->repropagate = repropagate;

   /* add the valid queue elements to the conflict set  */
   SCIPsetDebugMsg(set, "adding %d variables from the queue as temporary conflict variables\n", nbdchginfos);
   SCIP_CALL( conflictsetAddBounds(conflict, conflictset, blkmem, set, bdchginfos, nbdchginfos) );

   /* calculate the depth, at which the conflictset should be inserted */
   SCIP_CALL( conflictsetCalcInsertDepth(conflictset, set, tree) );
   assert(conflictset->validdepth <= conflictset->insertdepth && conflictset->insertdepth <= currentdepth);
   SCIPsetDebugMsg(set, " -> conflict with %d literals found at depth %d is active in depth %d and valid in depth %d\n",
      conflictset->nbdchginfos, currentdepth, conflictset->insertdepth, conflictset->validdepth);

   /* if all branching variables are in the conflict set, the conflict set is of no use;
    * don't use conflict sets that are only valid in the probing path but not in the problem tree
    */
   if( (diving || conflictset->insertdepth < currentdepth) && conflictset->insertdepth <= focusdepth )
   {
      /* if the conflict should not be located only in the subtree where it is useful, put it to its valid depth level */
      if( !set->conf_settlelocal )
         conflictset->insertdepth = conflictset->validdepth;

      *nliterals = conflictset->nbdchginfos;
      SCIPsetDebugMsg(set, " -> final conflict set has %d literals\n", *nliterals);

      /* check conflict set on debugging solution */
      SCIP_CALL( SCIPdebugCheckConflict(blkmem, set, tree->path[validdepth], \
            conflictset->bdchginfos, conflictset->relaxedbds, conflictset->nbdchginfos) ); /*lint !e506 !e774*/

      /* move conflictset to the conflictset storage */
      SCIP_CALL( conflictInsertConflictset(conflict, blkmem, set, &conflictset) );
      *success = TRUE;
   }
   else
   {
      /* free the temporary conflict set */
      conflictsetFree(&conflictset, blkmem);
   }

   return SCIP_OKAY;
}

/** tries to resolve given bound change
 *   - resolutions on local constraints are only applied, if the constraint is valid at the
 *     current minimal valid depth level, because this depth level is the topmost level to add the conflict
 *     constraint to anyways
 *
 *  @note it is sufficient to explain the relaxed bound change
 */
static
SCIP_RETCODE conflictResolveBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change to resolve */
   SCIP_Real             relaxedbd,          /**< the relaxed bound */
   int                   validdepth,         /**< minimal depth level at which the conflict is valid */
   SCIP_Bool*            resolved            /**< pointer to store whether the bound change was resolved */
   )
{
   SCIP_VAR* actvar;
   SCIP_CONS* infercons;
   SCIP_PROP* inferprop;
   SCIP_RESULT result;

#ifndef NDEBUG
   int nforcedbdchgqueue;
   int nbdchgqueue;

   /* store the current size of the conflict queues */
   assert(conflict != NULL);
   nforcedbdchgqueue = SCIPpqueueNElems(conflict->forcedbdchgqueue);
   nbdchgqueue = SCIPpqueueNElems(conflict->bdchgqueue);
#else
   assert(conflict != NULL);
#endif

   assert(resolved != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   *resolved = FALSE;

   actvar = SCIPbdchginfoGetVar(bdchginfo);
   assert(actvar != NULL);
   assert(SCIPvarIsActive(actvar));

#ifdef SCIP_DEBUG
   {
      int i;
      SCIPsetDebugMsg(set, "processing next conflicting bound (depth: %d, valid depth: %d, bdchgtype: %s [%s], vartype: %d): [<%s> %s %g(%g)]\n",
         SCIPbdchginfoGetDepth(bdchginfo), validdepth,
         SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
            : SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER ? "cons" : "prop",
               SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "-"
                        : SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
                          ? SCIPconsGetName(SCIPbdchginfoGetInferCons(bdchginfo))
                                : SCIPbdchginfoGetInferProp(bdchginfo) == NULL ? "-"
                                      : SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo)),
                                        SCIPvarGetType(actvar), SCIPvarGetName(actvar),
                                        SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
                                              SCIPbdchginfoGetNewbound(bdchginfo), relaxedbd);
      SCIPsetDebugMsg(set, " - conflict set       :");

      for( i = 0; i < conflict->conflictset->nbdchginfos; ++i )
      {
         SCIPsetDebugMsgPrint(set, " [%d:<%s> %s %g(%g)]", SCIPbdchginfoGetDepth(conflict->conflictset->bdchginfos[i]),
               SCIPvarGetName(SCIPbdchginfoGetVar(conflict->conflictset->bdchginfos[i])),
               SCIPbdchginfoGetBoundtype(conflict->conflictset->bdchginfos[i]) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
                     SCIPbdchginfoGetNewbound(conflict->conflictset->bdchginfos[i]), conflict->conflictset->relaxedbds[i]);
      }
      SCIPsetDebugMsgPrint(set, "\n");
      SCIPsetDebugMsg(set, " - forced candidates  :");

      for( i = 0; i < SCIPpqueueNElems(conflict->forcedbdchgqueue); ++i )
      {
         SCIP_BDCHGINFO* info = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->forcedbdchgqueue)[i]);
         SCIPsetDebugMsgPrint(set, " [%d:<%s> %s %g(%g)]", SCIPbdchginfoGetDepth(info), SCIPvarGetName(SCIPbdchginfoGetVar(info)),
               bdchginfoIsInvalid(conflict, info) ? "<!>" : SCIPbdchginfoGetBoundtype(info) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
                     SCIPbdchginfoGetNewbound(info), SCIPbdchginfoGetRelaxedBound(info));
      }
      SCIPsetDebugMsgPrint(set, "\n");
      SCIPsetDebugMsg(set, " - optional candidates:");

      for( i = 0; i < SCIPpqueueNElems(conflict->bdchgqueue); ++i )
      {
         SCIP_BDCHGINFO* info = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->bdchgqueue)[i]);
         SCIPsetDebugMsgPrint(set, " [%d:<%s> %s %g(%g)]", SCIPbdchginfoGetDepth(info), SCIPvarGetName(SCIPbdchginfoGetVar(info)),
               bdchginfoIsInvalid(conflict, info) ? "<!>" : SCIPbdchginfoGetBoundtype(info) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
                     SCIPbdchginfoGetNewbound(info), SCIPbdchginfoGetRelaxedBound(info));
      }
      SCIPsetDebugMsgPrint(set, "\n");
   }
#endif

   /* check, if the bound change can and should be resolved:
    *  - resolutions on local constraints should only be applied, if the constraint is valid at the
    *    current minimal valid depth level (which is initialized with the valid depth level of the initial
    *    conflict set), because this depth level is the topmost level to add the conflict constraint to anyways
    */
   switch( SCIPbdchginfoGetChgtype(bdchginfo) )
   {
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      infercons = SCIPbdchginfoGetInferCons(bdchginfo);
      assert(infercons != NULL);

      if( SCIPconsIsGlobal(infercons) || SCIPconsGetValidDepth(infercons) <= validdepth )
      {
         SCIP_VAR* infervar;
         int inferinfo;
         SCIP_BOUNDTYPE inferboundtype;
         SCIP_BDCHGIDX* bdchgidx;

         /* resolve bound change by asking the constraint that infered the bound to put all bounds that were
          * the reasons for the conflicting bound change on the priority queue
          */
         infervar = SCIPbdchginfoGetInferVar(bdchginfo);
         inferinfo = SCIPbdchginfoGetInferInfo(bdchginfo);
         inferboundtype = SCIPbdchginfoGetInferBoundtype(bdchginfo);
         bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
         assert(infervar != NULL);

         SCIPsetDebugMsg(set, "resolving bound <%s> %s %g(%g) [status:%d, type:%d, depth:%d, pos:%d]: <%s> %s %g [cons:<%s>(%s), info:%d]\n",
            SCIPvarGetName(actvar),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo), relaxedbd,
            SCIPvarGetStatus(actvar), SCIPvarGetType(actvar),
            SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPvarGetName(infervar),
            inferboundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPgetVarBdAtIndex(set->scip, infervar, inferboundtype, bdchgidx, TRUE),
            SCIPconsGetName(infercons),
            SCIPconsIsGlobal(infercons) ? "global" : "local",
            inferinfo);

         /* in case the inference variables is not an active variables, we need to transform the relaxed bound */
         if( actvar != infervar )
         {
            SCIP_VAR* var;
            SCIP_Real scalar;
            SCIP_Real constant;

            assert(SCIPvarGetStatus(infervar) == SCIP_VARSTATUS_AGGREGATED
               || SCIPvarGetStatus(infervar) == SCIP_VARSTATUS_NEGATED
               || (SCIPvarGetStatus(infervar) == SCIP_VARSTATUS_MULTAGGR && SCIPvarGetMultaggrNVars(infervar) == 1));

            scalar = 1.0;
            constant = 0.0;

            var = infervar;

            /* transform given varibale to active varibale */
            SCIP_CALL( SCIPvarGetProbvarSum(&var, set, &scalar, &constant) );
            assert(var == actvar);

            relaxedbd *= scalar;
            relaxedbd += constant;
         }

         SCIP_CALL( SCIPconsResolvePropagation(infercons, set, infervar, inferinfo, inferboundtype, bdchgidx, relaxedbd, &result) );
         *resolved = (result == SCIP_SUCCESS);
      }
      break;

   case SCIP_BOUNDCHGTYPE_PROPINFER:
      inferprop = SCIPbdchginfoGetInferProp(bdchginfo);
      if( inferprop != NULL )
      {
         SCIP_VAR* infervar;
         int inferinfo;
         SCIP_BOUNDTYPE inferboundtype;
         SCIP_BDCHGIDX* bdchgidx;

         /* resolve bound change by asking the propagator that infered the bound to put all bounds that were
          * the reasons for the conflicting bound change on the priority queue
          */
         infervar = SCIPbdchginfoGetInferVar(bdchginfo);
         inferinfo = SCIPbdchginfoGetInferInfo(bdchginfo);
         inferboundtype = SCIPbdchginfoGetInferBoundtype(bdchginfo);
         bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
         assert(infervar != NULL);

         SCIPsetDebugMsg(set, "resolving bound <%s> %s %g(%g) [status:%d, depth:%d, pos:%d]: <%s> %s %g [prop:<%s>, info:%d]\n",
            SCIPvarGetName(actvar),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo), relaxedbd,
            SCIPvarGetStatus(actvar), SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPvarGetName(infervar),
            inferboundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPgetVarBdAtIndex(set->scip, infervar, inferboundtype, bdchgidx, TRUE),
            SCIPpropGetName(inferprop), inferinfo);

         SCIP_CALL( SCIPpropResolvePropagation(inferprop, set, infervar, inferinfo, inferboundtype, bdchgidx, relaxedbd, &result) );
         *resolved = (result == SCIP_SUCCESS);
      }
      break;

   case SCIP_BOUNDCHGTYPE_BRANCHING:
      assert(!(*resolved));
      break;

   default:
      SCIPerrorMessage("invalid bound change type <%d>\n", SCIPbdchginfoGetChgtype(bdchginfo));
      return SCIP_INVALIDDATA;
   }

   SCIPsetDebugMsg(set, "resolving status: %u\n", *resolved);

#ifndef NDEBUG
   /* subtract the size of the conflicq queues */
   nforcedbdchgqueue -= SCIPpqueueNElems(conflict->forcedbdchgqueue);
   nbdchgqueue -= SCIPpqueueNElems(conflict->bdchgqueue);

   /* in case the bound change was not resolved, the conflict queues should have the same size (contents) */
   assert((*resolved) || (nforcedbdchgqueue == 0 && nbdchgqueue == 0));
#endif

   return SCIP_OKAY;
}

/** if only one conflicting bound change of the last depth level was used, and if this can be resolved,
 *  creates GRASP-like reconvergence conflict constraints in the conflict graph up to the branching variable of this
 *  depth level
 */
static
SCIP_RETCODE conflictCreateReconvergenceConss(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_BDCHGINFO*       firstuip,           /**< first UIP of conflict graph */
   int*                  nreconvconss,       /**< pointer to store the number of generated reconvergence constraints */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence constraints */
   )
{
   SCIP_BDCHGINFO* uip;
   SCIP_CONFTYPE conftype;
   SCIP_Bool usescutoffbound;
   int firstuipdepth;
   int focusdepth;
   int currentdepth;
   int maxvaliddepth;

   assert(conflict != NULL);
   assert(firstuip != NULL);
   assert(nreconvconss != NULL);
   assert(nreconvliterals != NULL);
   assert(!SCIPbdchginfoIsRedundant(firstuip));

   focusdepth = SCIPtreeGetFocusDepth(tree);
   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(focusdepth <= currentdepth);

   /* check, whether local constraints are allowed; however, don't generate reconvergence constraints that are only valid
    * in the probing path and not in the problem tree (i.e. that exceed the focusdepth)
    */
   maxvaliddepth = (set->conf_allowlocal ? MIN(currentdepth-1, focusdepth) : 0);
   if( validdepth > maxvaliddepth )
      return SCIP_OKAY;

   firstuipdepth = SCIPbdchginfoGetDepth(firstuip);

   conftype = conflict->conflictset->conflicttype;
   usescutoffbound = conflict->conflictset->usescutoffbound;

   /* for each succeeding UIP pair of the last depth level, create one reconvergence constraint */
   uip = firstuip;
   while( uip != NULL && SCIPbdchginfoGetDepth(uip) == SCIPbdchginfoGetDepth(firstuip) && bdchginfoIsResolvable(uip) )
   {
      SCIP_BDCHGINFO* oppositeuip;
      SCIP_BDCHGINFO* bdchginfo;
      SCIP_BDCHGINFO* nextuip;
      SCIP_VAR* uipvar;
      SCIP_Real oppositeuipbound;
      SCIP_BOUNDTYPE oppositeuipboundtype;
      int nresolutions;

      assert(!SCIPbdchginfoIsRedundant(uip));

      SCIPsetDebugMsg(set, "creating reconvergence constraint for UIP <%s> %s %g in depth %d pos %d\n",
         SCIPvarGetName(SCIPbdchginfoGetVar(uip)), SCIPbdchginfoGetBoundtype(uip) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
         SCIPbdchginfoGetNewbound(uip), SCIPbdchginfoGetDepth(uip), SCIPbdchginfoGetPos(uip));

      /* initialize conflict data */
      SCIP_CALL( SCIPconflictInit(conflict, set, stat, prob, conftype, usescutoffbound) );

      conflict->conflictset->conflicttype = conftype;
      conflict->conflictset->usescutoffbound = usescutoffbound;

      /* create a temporary bound change information for the negation of the UIP's bound change;
       * this bound change information is freed in the SCIPconflictFlushConss() call;
       * for reconvergence constraints for continuous variables we can only use the "negation" !(x <= u) == (x >= u);
       * during conflict analysis, we treat a continuous bound "x >= u" in the conflict set as "x > u", and in the
       * generated constraint this is negated again to "x <= u" which is correct.
       */
      uipvar = SCIPbdchginfoGetVar(uip);
      oppositeuipboundtype = SCIPboundtypeOpposite(SCIPbdchginfoGetBoundtype(uip));
      oppositeuipbound = SCIPbdchginfoGetNewbound(uip);
      if( SCIPvarIsIntegral(uipvar) )
      {
         assert(SCIPsetIsIntegral(set, oppositeuipbound));
         oppositeuipbound += (oppositeuipboundtype == SCIP_BOUNDTYPE_LOWER ? +1.0 : -1.0);
      }
      SCIP_CALL( conflictCreateTmpBdchginfo(conflict, blkmem, set, uipvar, oppositeuipboundtype, \
            oppositeuipboundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_REAL_MIN : SCIP_REAL_MAX, oppositeuipbound, &oppositeuip) );

      /* put the negated UIP into the conflict set */
      SCIP_CALL( conflictAddConflictBound(conflict, blkmem, set, oppositeuip, oppositeuipbound) );

      /* put positive UIP into priority queue */
      SCIP_CALL( conflictQueueBound(conflict, set, uip, SCIPbdchginfoGetNewbound(uip) ) );

      /* resolve the queue until the next UIP is reached */
      bdchginfo = conflictFirstCand(conflict);
      nextuip = NULL;
      nresolutions = 0;
      while( bdchginfo != NULL && validdepth <= maxvaliddepth )
      {
         SCIP_BDCHGINFO* nextbdchginfo;
         SCIP_Real relaxedbd;
         SCIP_Bool forceresolve;
         int bdchgdepth;

         /* check if the next bound change must be resolved in every case */
         forceresolve = (SCIPpqueueNElems(conflict->forcedbdchgqueue) > 0);

         /* remove currently processed candidate and get next conflicting bound from the conflict candidate queue before
          * we remove the candidate we have to collect the relaxed bound since removing the candidate from the queue
          * invalidates the relaxed bound
          */
         assert(bdchginfo == conflictFirstCand(conflict));
         relaxedbd = SCIPbdchginfoGetRelaxedBound(bdchginfo);
         bdchginfo = conflictRemoveCand(conflict);
         nextbdchginfo = conflictFirstCand(conflict);
         bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);
         assert(bdchginfo != NULL);
         assert(!SCIPbdchginfoIsRedundant(bdchginfo));
         assert(nextbdchginfo == NULL || SCIPbdchginfoGetDepth(bdchginfo) >= SCIPbdchginfoGetDepth(nextbdchginfo)
            || forceresolve);
         assert(bdchgdepth <= firstuipdepth);

         /* bound changes that are higher in the tree than the valid depth of the conflict can be ignored;
          * multiple insertions of the same bound change can be ignored
          */
         if( bdchgdepth > validdepth && bdchginfo != nextbdchginfo )
         {
            SCIP_VAR* actvar;
            SCIP_Bool resolved;

            actvar = SCIPbdchginfoGetVar(bdchginfo);
            assert(actvar != NULL);
            assert(SCIPvarIsActive(actvar));

            /* check if we have to resolve the bound change in this depth level
             *  - the starting uip has to be resolved
             *  - a bound change should be resolved, if it is in the fuip's depth level and not the
             *    next uip (i.e., if it is not the last bound change in the fuip's depth level)
             *  - a forced bound change must be resolved in any case
             */
            resolved = FALSE;
            if( bdchginfo == uip
               || (bdchgdepth == firstuipdepth
                  && nextbdchginfo != NULL
                  && SCIPbdchginfoGetDepth(nextbdchginfo) == bdchgdepth)
               || forceresolve )
            {
               SCIP_CALL( conflictResolveBound(conflict, set, bdchginfo, relaxedbd, validdepth, &resolved) );
            }

            if( resolved )
               nresolutions++;
            else if( forceresolve )
            {
               /* variable cannot enter the conflict clause: we have to make the conflict clause local, s.t.
                * the unresolved bound change is active in the whole sub tree of the conflict clause
                */
               assert(bdchgdepth >= validdepth);
               validdepth = bdchgdepth;

               SCIPsetDebugMsg(set, "couldn't resolve forced bound change on <%s> -> new valid depth: %d\n",
                  SCIPvarGetName(actvar), validdepth);
            }
            else if( bdchginfo != uip )
            {
               assert(conflict->conflictset != NULL);
               assert(conflict->conflictset->nbdchginfos >= 1); /* starting UIP is already member of the conflict set */

               /* if this is the first variable of the conflict set besides the current starting UIP, it is the next
                * UIP (or the first unresolvable bound change)
                */
               if( bdchgdepth == firstuipdepth && conflict->conflictset->nbdchginfos == 1 )
               {
                  assert(nextuip == NULL);
                  nextuip = bdchginfo;
               }

               /* put bound change into the conflict set */
               SCIP_CALL( conflictAddConflictBound(conflict, blkmem, set, bdchginfo, relaxedbd) );
               assert(conflict->conflictset->nbdchginfos >= 2);
            }
            else
               assert(conflictFirstCand(conflict) == NULL); /* the starting UIP was not resolved */
         }

         /* get next conflicting bound from the conflict candidate queue (this does not need to be nextbdchginfo, because
          * due to resolving the bound changes, a variable could be added to the queue which must be
          * resolved before nextbdchginfo)
          */
         bdchginfo = conflictFirstCand(conflict);
      }
      assert(nextuip != uip);

      /* if only one propagation was resolved, the reconvergence constraint is already member of the constraint set
       * (it is exactly the constraint that produced the propagation)
       */
      if( nextuip != NULL && nresolutions >= 2 && bdchginfo == NULL && validdepth <= maxvaliddepth )
      {
         int nlits;
         SCIP_Bool success;

         assert(SCIPbdchginfoGetDepth(nextuip) == SCIPbdchginfoGetDepth(uip));

         /* check conflict graph frontier on debugging solution */
         SCIP_CALL( SCIPdebugCheckConflictFrontier(blkmem, set, tree->path[validdepth], \
               bdchginfo, conflict->conflictset->bdchginfos, conflict->conflictset->relaxedbds, \
               conflict->conflictset->nbdchginfos, conflict->bdchgqueue, conflict->forcedbdchgqueue) ); /*lint !e506 !e774*/

         SCIPsetDebugMsg(set, "creating reconvergence constraint from UIP <%s> to UIP <%s> in depth %d with %d literals after %d resolutions\n",
            SCIPvarGetName(SCIPbdchginfoGetVar(uip)), SCIPvarGetName(SCIPbdchginfoGetVar(nextuip)),
            SCIPbdchginfoGetDepth(uip), conflict->conflictset->nbdchginfos, nresolutions);

         /* call the conflict handlers to create a conflict set */
         SCIP_CALL( conflictAddConflictset(conflict, blkmem, set, stat, tree, validdepth, diving, FALSE, &success, &nlits) );
         if( success )
         {
            (*nreconvconss)++;
            (*nreconvliterals) += nlits;
         }
      }

      /* clear the conflict candidate queue and the conflict set (to make sure, oppositeuip is not referenced anymore) */
      conflictClear(conflict);

      uip = nextuip;
   }

   conflict->conflictset->conflicttype = conftype;
   conflict->conflictset->usescutoffbound = usescutoffbound;


   return SCIP_OKAY;
}

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound() and
 *  SCIPconflictAddRelaxedBound(), and on success, calls the conflict handlers to create a conflict constraint out of
 *  the resulting conflict set; afterwards the conflict queue and the conflict set is cleared
 */
static
SCIP_RETCODE conflictAnalyze(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool             mustresolve,        /**< should the conflict set only be used, if a resolution was applied? */
   int*                  nconss,             /**< pointer to store the number of generated conflict constraints */
   int*                  nliterals,          /**< pointer to store the number of literals in generated conflict constraints */
   int*                  nreconvconss,       /**< pointer to store the number of generated reconvergence constraints */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence constraints */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BDCHGINFO** firstuips;
   int nfirstuips;
   int focusdepth;
   int currentdepth;
   int maxvaliddepth;
   int resolvedepth;
   int nresolutions;
   int lastconsnresolutions;
   int lastconsresoldepth;

   assert(conflict != NULL);
   assert(conflict->conflictset != NULL);
   assert(conflict->conflictset->nbdchginfos >= 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(0 <= validdepth && validdepth <= SCIPtreeGetCurrentDepth(tree));
   assert(nconss != NULL);
   assert(nliterals != NULL);
   assert(nreconvconss != NULL);
   assert(nreconvliterals != NULL);

   focusdepth = SCIPtreeGetFocusDepth(tree);
   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(focusdepth <= currentdepth);

   resolvedepth = ((set->conf_fuiplevels >= 0 && set->conf_fuiplevels <= currentdepth)
      ? currentdepth - set->conf_fuiplevels + 1 : 0);
   assert(0 <= resolvedepth && resolvedepth <= currentdepth + 1);

   /* if we must resolve at least one bound change, find the first UIP at least in the last depth level */
   if( mustresolve )
      resolvedepth = MIN(resolvedepth, currentdepth);

   SCIPsetDebugMsg(set, "analyzing conflict with %d+%d conflict candidates and starting conflict set of size %d in depth %d (resolvedepth=%d)\n",
      SCIPpqueueNElems(conflict->forcedbdchgqueue), SCIPpqueueNElems(conflict->bdchgqueue),
      conflict->conflictset->nbdchginfos, currentdepth, resolvedepth);

   *nconss = 0;
   *nliterals = 0;
   *nreconvconss = 0;
   *nreconvliterals = 0;

   /* check, whether local conflicts are allowed; however, don't generate conflict constraints that are only valid in the
    * probing path and not in the problem tree (i.e. that exceed the focusdepth)
    */
   maxvaliddepth = (set->conf_allowlocal ? MIN(currentdepth-1, focusdepth) : 0);
   if( validdepth > maxvaliddepth )
      return SCIP_OKAY;

   /* allocate temporary memory for storing first UIPs (in each depth level, at most two bound changes can be flagged
    * as UIP, namely a binary and a non-binary bound change)
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &firstuips, 2*(currentdepth+1)) ); /*lint !e647*/

   /* process all bound changes in the conflict candidate queue */
   nresolutions = 0;
   lastconsnresolutions = (mustresolve ? 0 : -1);
   lastconsresoldepth = (mustresolve ? currentdepth : INT_MAX);
   bdchginfo = conflictFirstCand(conflict);
   nfirstuips = 0;

   /* check if the initial reason on debugging solution */
   SCIP_CALL( SCIPdebugCheckConflictFrontier(blkmem, set, tree->path[validdepth], \
         NULL, conflict->conflictset->bdchginfos, conflict->conflictset->relaxedbds, conflict->conflictset->nbdchginfos, \
         conflict->bdchgqueue, conflict->forcedbdchgqueue) ); /*lint !e506 !e774*/

   while( bdchginfo != NULL && validdepth <= maxvaliddepth )
   {
      SCIP_BDCHGINFO* nextbdchginfo;
      SCIP_Real relaxedbd;
      SCIP_Bool forceresolve;
      int bdchgdepth;

      assert(!SCIPbdchginfoIsRedundant(bdchginfo));

      /* check if the next bound change must be resolved in every case */
      forceresolve = (SCIPpqueueNElems(conflict->forcedbdchgqueue) > 0);

      /* resolve next bound change in queue */
      bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);
      assert(0 <= bdchgdepth && bdchgdepth <= currentdepth);
      assert(SCIPvarIsActive(SCIPbdchginfoGetVar(bdchginfo)));
      assert(bdchgdepth < tree->pathlen);
      assert(tree->path[bdchgdepth] != NULL);
      assert(tree->path[bdchgdepth]->domchg != NULL);
      assert(SCIPbdchginfoGetPos(bdchginfo) < (int)tree->path[bdchgdepth]->domchg->domchgbound.nboundchgs);
      assert(tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].var
         == SCIPbdchginfoGetVar(bdchginfo));
      assert(tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].newbound
         == SCIPbdchginfoGetNewbound(bdchginfo)
         || (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER
            ? SCIPvarGetLbGlobal(SCIPbdchginfoGetVar(bdchginfo)) : SCIPvarGetUbGlobal(SCIPbdchginfoGetVar(bdchginfo)))
         == SCIPbdchginfoGetNewbound(bdchginfo)); /*lint !e777*/
      assert((SCIP_BOUNDTYPE)tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].boundtype
         == SCIPbdchginfoGetBoundtype(bdchginfo));

      /* create intermediate conflict constraint */
      assert(nresolutions >= lastconsnresolutions);
      if( !forceresolve )
      {
         if( nresolutions == lastconsnresolutions )
            lastconsresoldepth = bdchgdepth; /* all intermediate depth levels consisted of only unresolved bound changes */
         else if( bdchgdepth < lastconsresoldepth && (set->conf_interconss == -1 || *nconss < set->conf_interconss) )
         {
            int nlits;
            SCIP_Bool success;

            /* call the conflict handlers to create a conflict set */
            SCIPsetDebugMsg(set, "creating intermediate conflictset after %d resolutions up to depth %d (valid at depth %d): %d conflict bounds, %d bounds in queue\n",
               nresolutions, bdchgdepth, validdepth, conflict->conflictset->nbdchginfos,
               SCIPpqueueNElems(conflict->bdchgqueue));

            SCIP_CALL( conflictAddConflictset(conflict, blkmem, set, stat, tree, validdepth, diving, TRUE, &success, &nlits) );
            lastconsnresolutions = nresolutions;
            lastconsresoldepth = bdchgdepth;
            if( success )
            {
               (*nconss)++;
               (*nliterals) += nlits;
            }
         }
      }

      /* remove currently processed candidate and get next conflicting bound from the conflict candidate queue before
       * we remove the candidate we have to collect the relaxed bound since removing the candidate from the queue
       * invalidates the relaxed bound
       */
      assert(bdchginfo == conflictFirstCand(conflict));
      relaxedbd = SCIPbdchginfoGetRelaxedBound(bdchginfo);
      bdchginfo = conflictRemoveCand(conflict);
      nextbdchginfo = conflictFirstCand(conflict);
      assert(bdchginfo != NULL);
      assert(!SCIPbdchginfoIsRedundant(bdchginfo));
      assert(nextbdchginfo == NULL || SCIPbdchginfoGetDepth(bdchginfo) >= SCIPbdchginfoGetDepth(nextbdchginfo)
         || forceresolve);

      /* we don't need to resolve bound changes that are already active in the valid depth of the current conflict set,
       * because the conflict set can only be added locally at the valid depth, and all bound changes applied in this
       * depth or earlier can be removed from the conflict constraint, since they are already applied in the constraint's
       * subtree;
       * if the next bound change on the remaining queue is equal to the current bound change,
       * this is a multiple insertion in the conflict candidate queue and we can ignore the current
       * bound change
       */
      if( bdchgdepth > validdepth && bdchginfo != nextbdchginfo )
      {
         SCIP_VAR* actvar;
         SCIP_Bool resolved;

         actvar = SCIPbdchginfoGetVar(bdchginfo);
         assert(actvar != NULL);
         assert(SCIPvarIsActive(actvar));

         /* check if we want to resolve the bound change in this depth level
          *  - bound changes should be resolved, if
          *     (i)   we must apply at least one resolution and didn't resolve a bound change yet, or
          *     (ii)  their depth level is at least equal to the minimal resolving depth, and
          *           they are not the last remaining conflicting bound change in their depth level
          *     (iii) the bound change resolving is forced (i.e., the forced queue was non-empty)
          */
         resolved = FALSE;
         if( (mustresolve && nresolutions == 0)
            || (bdchgdepth >= resolvedepth
               && nextbdchginfo != NULL
               && SCIPbdchginfoGetDepth(nextbdchginfo) == bdchgdepth)
            || forceresolve )
         {
            SCIP_CALL( conflictResolveBound(conflict, set, bdchginfo, relaxedbd, validdepth, &resolved) );
         }

         if( resolved )
            nresolutions++;
         else if( forceresolve )
         {
            /* variable cannot enter the conflict clause: we have to make the conflict clause local, s.t.
             * the unresolved bound change is active in the whole sub tree of the conflict clause
             */
            assert(bdchgdepth >= validdepth);
            validdepth = bdchgdepth;

            SCIPsetDebugMsg(set, "couldn't resolve forced bound change on <%s> -> new valid depth: %d\n",
               SCIPvarGetName(actvar), validdepth);
         }
         else
         {
            /* if this is a UIP (the last bound change in its depth level), it can be used to generate a
             * UIP reconvergence constraint
             */
            if( nextbdchginfo == NULL || SCIPbdchginfoGetDepth(nextbdchginfo) != bdchgdepth )
            {
               assert(nfirstuips < 2*(currentdepth+1));
               firstuips[nfirstuips] = bdchginfo;
               nfirstuips++;
            }

            /* put variable into the conflict set, using the literal that is currently fixed to FALSE */
            SCIP_CALL( conflictAddConflictBound(conflict, blkmem, set, bdchginfo, relaxedbd) );
         }
      }

      /* check conflict graph frontier on debugging solution */
      SCIP_CALL( SCIPdebugCheckConflictFrontier(blkmem, set, tree->path[validdepth], \
            bdchginfo, conflict->conflictset->bdchginfos, conflict->conflictset->relaxedbds, conflict->conflictset->nbdchginfos, \
            conflict->bdchgqueue, conflict->forcedbdchgqueue) ); /*lint !e506 !e774*/

      /* get next conflicting bound from the conflict candidate queue (this needs not to be nextbdchginfo, because
       * due to resolving the bound changes, a bound change could be added to the queue which must be
       * resolved before nextbdchginfo)
       */
      bdchginfo = conflictFirstCand(conflict);
   }

   /* check, if a valid conflict set was found */
   if( bdchginfo == NULL
      && nresolutions > lastconsnresolutions
      && validdepth <= maxvaliddepth
      && (!mustresolve || nresolutions > 0 || conflict->conflictset->nbdchginfos == 0)
      && SCIPpqueueNElems(conflict->forcedbdchgqueue) == 0 )
   {
      int nlits;
      SCIP_Bool success;

      /* call the conflict handlers to create a conflict set */
      SCIP_CALL( conflictAddConflictset(conflict, blkmem, set, stat, tree, validdepth, diving, TRUE, &success, &nlits) );
      if( success )
      {
         (*nconss)++;
         (*nliterals) += nlits;
      }
   }

   /* produce reconvergence constraints defined by succeeding UIP's of the last depth level */
   if( set->conf_reconvlevels != 0 && validdepth <= maxvaliddepth )
   {
      int reconvlevels;
      int i;

      reconvlevels = (set->conf_reconvlevels == -1 ? INT_MAX : set->conf_reconvlevels);
      for( i = 0; i < nfirstuips; ++i )
      {
         if( SCIPbdchginfoHasInferenceReason(firstuips[i])
            && currentdepth - SCIPbdchginfoGetDepth(firstuips[i]) < reconvlevels )
         {
            SCIP_CALL( conflictCreateReconvergenceConss(conflict, blkmem, set, stat, prob, tree, diving, \
                  validdepth, firstuips[i], nreconvconss, nreconvliterals) );
         }
      }
   }

   /* free the temporary memory */
   SCIPsetFreeBufferArray(set, &firstuips);

   /* clear the conflict candidate queue and the conflict set */
   conflictClear(conflict);

   return SCIP_OKAY;
}

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set;
 *  updates statistics for propagation conflict analysis
 */
SCIP_RETCODE SCIPconflictAnalyze(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   int nconss;
   int nliterals;
   int nreconvconss;
   int nreconvliterals;

   assert(conflict != NULL);
   assert(conflict->conflictset != NULL);
   assert(set != NULL);
   assert(prob != NULL);

   if( success != NULL )
      *success = FALSE;

   /* check if the conflict analysis is applicable */
   if( !SCIPconflictApplicable(set) )
      return SCIP_OKAY;

   /* check, if the conflict set will get too large with high probability */
   if( conflict->conflictset->nbdchginfos + SCIPpqueueNElems(conflict->bdchgqueue)
      + SCIPpqueueNElems(conflict->forcedbdchgqueue) >= 2*conflictCalcMaxsize(set, prob) )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "analyzing conflict after infeasible propagation in depth %d\n", SCIPtreeGetCurrentDepth(tree));

   /* start timing */
   SCIPclockStart(conflict->propanalyzetime, set);

   conflict->npropcalls++;

   /* analyze the conflict set, and create a conflict constraint on success */
   SCIP_CALL( conflictAnalyze(conflict, blkmem, set, stat, prob, tree, FALSE, validdepth, TRUE, &nconss, &nliterals, \
         &nreconvconss, &nreconvliterals) );
   conflict->npropsuccess += (nconss > 0 ? 1 : 0);
   conflict->npropconfconss += nconss;
   conflict->npropconfliterals += nliterals;
   conflict->npropreconvconss += nreconvconss;
   conflict->npropreconvliterals += nreconvliterals;
   if( success != NULL )
      *success = (nconss > 0);

   /* stop timing */
   SCIPclockStop(conflict->propanalyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for preprocessing global conflict constraint before appliance */
SCIP_Real SCIPconflictGetGlobalApplTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->dIBclock);
}

/** gets time in seconds used for analyzing propagation conflicts */
SCIP_Real SCIPconflictGetPropTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->propanalyzetime);
}

/** gets number of calls to propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropcalls;
}

/** gets number of calls to propagation conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNPropSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropsuccess;
}

/** gets number of conflict constraints detected in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropconfconss;
}

/** gets total number of literals in conflict constraints created in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropconfliterals;
}

/** gets number of reconvergence constraints detected in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropreconvconss;
}

/** gets total number of literals in reconvergence constraints created in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropreconvliterals;
}




/*
 * Infeasible LP Conflict Analysis
 */

/** ensures, that side change arrays can store at least num entries */
static
SCIP_RETCODE ensureSidechgsSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   int**                 sidechginds,        /**< pointer to side change index array */
   SCIP_Real**           sidechgoldlhss,     /**< pointer to side change old left hand sides array */
   SCIP_Real**           sidechgoldrhss,     /**< pointer to side change old right hand sides array */
   SCIP_Real**           sidechgnewlhss,     /**< pointer to side change new left hand sides array */
   SCIP_Real**           sidechgnewrhss,     /**< pointer to side change new right hand sides array */
   int*                  sidechgssize,       /**< pointer to size of side change arrays */
   int                   num                 /**< minimal number of entries to be able to store in side change arrays */
   )
{
   assert(sidechginds != NULL);
   assert(sidechgoldlhss != NULL);
   assert(sidechgoldrhss != NULL);
   assert(sidechgnewlhss != NULL);
   assert(sidechgnewrhss != NULL);
   assert(sidechgssize != NULL);

   if( num > *sidechgssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechginds, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgoldlhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgoldrhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgnewlhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgnewrhss, newsize) );
      *sidechgssize = newsize;
   }
   assert(num <= *sidechgssize);

   return SCIP_OKAY;
}

/** adds removal of row's side to side change arrays; finite sides are only replaced by near infinite sides, such
 *  that the row's sense in the LP solver is not changed
 */
static
SCIP_RETCODE addSideRemoval(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row,                /**< LP row to change the sides for */
   SCIP_Real             lpiinfinity,        /**< value treated as infinity in LP solver */
   int**                 sidechginds,        /**< pointer to side change index array */
   SCIP_Real**           sidechgoldlhss,     /**< pointer to side change old left hand sides array */
   SCIP_Real**           sidechgoldrhss,     /**< pointer to side change old right hand sides array */
   SCIP_Real**           sidechgnewlhss,     /**< pointer to side change new left hand sides array */
   SCIP_Real**           sidechgnewrhss,     /**< pointer to side change new right hand sides array */
   int*                  sidechgssize,       /**< pointer to size of side change arrays */
   int*                  nsidechgs           /**< pointer to number of used slots in side change arrays */
   )
{
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real constant;

   assert(sidechginds != NULL);
   assert(sidechgoldlhss != NULL);
   assert(sidechgoldrhss != NULL);
   assert(sidechgnewlhss != NULL);
   assert(sidechgnewrhss != NULL);
   assert(sidechgssize != NULL);
   assert(nsidechgs != NULL);

   lhs = SCIProwGetLhs(row);
   rhs = SCIProwGetRhs(row);
   constant = SCIProwGetConstant(row);
   assert(!SCIPsetIsInfinity(set, -lhs) || !SCIPsetIsInfinity(set, rhs));

   /* get memory to store additional side change */
   SCIP_CALL( ensureSidechgsSize(set, sidechginds, sidechgoldlhss, sidechgoldrhss, sidechgnewlhss, sidechgnewrhss, \
         sidechgssize, (*nsidechgs)+1) );
   assert(*nsidechgs < *sidechgssize);
   assert(*sidechginds != NULL);
   assert(*sidechgoldlhss != NULL);
   assert(*sidechgoldrhss != NULL);
   assert(*sidechgnewlhss != NULL);
   assert(*sidechgnewrhss != NULL);

   /* store side change */
   (*sidechginds)[*nsidechgs] = SCIProwGetLPPos(row);
   if( SCIPsetIsInfinity(set, -lhs) )
   {
      (*sidechgoldlhss)[*nsidechgs] = -lpiinfinity;
      (*sidechgnewlhss)[*nsidechgs] = -lpiinfinity;
   }
   else
   {
      (*sidechgoldlhss)[*nsidechgs] = lhs - constant;
      (*sidechgnewlhss)[*nsidechgs] = -lpiinfinity;
   }
   if( SCIPsetIsInfinity(set, rhs) )
   {
      (*sidechgoldrhss)[*nsidechgs] = lpiinfinity;
      (*sidechgnewrhss)[*nsidechgs] = lpiinfinity;
   }
   else
   {
      (*sidechgoldrhss)[*nsidechgs] = rhs - constant;
      (*sidechgnewrhss)[*nsidechgs] = lpiinfinity;
   }
   (*nsidechgs)++;

   return SCIP_OKAY;
}

/** inserts variable's new bounds into bound change arrays */
static
SCIP_RETCODE addBdchg(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to change the LP bounds for */
   SCIP_Real             newlb,              /**< new lower bound */
   SCIP_Real             newub,              /**< new upper bound */
   SCIP_LPBDCHGS*        oldlpbdchgs,        /**< old LP bound changes used for reset the LP bound change */
   SCIP_LPBDCHGS*        relaxedlpbdchgs,    /**< relaxed LP bound changes used for reset the LP bound change */
   SCIP_LPI*             lpi                 /**< pointer to LPi to access infinity of LP solver; necessary to set correct value */
   )
{
   assert(newlb <= newub);
   assert(oldlpbdchgs != NULL);
   assert(relaxedlpbdchgs != NULL);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
   {
      SCIP_COL* col;
      int idx;
      int c;

      col = SCIPvarGetCol(var);
      c = SCIPcolGetLPPos(col);

      if( c >= 0 )
      {
         /* store old bound change for resetting the LP later */
         if( !oldlpbdchgs->usedcols[c] )
         {
            idx = oldlpbdchgs->nbdchgs;
            oldlpbdchgs->usedcols[c] = TRUE;
            oldlpbdchgs->bdchgcolinds[c] = idx;
            oldlpbdchgs->nbdchgs++;

            oldlpbdchgs->bdchginds[idx] = c;
            oldlpbdchgs->bdchglbs[idx] = SCIPvarGetLbLP(var, set);
            oldlpbdchgs->bdchgubs[idx] = SCIPvarGetUbLP(var, set);
         }
         assert(oldlpbdchgs->bdchginds[oldlpbdchgs->bdchgcolinds[c]] == c);
         assert((SCIPlpiIsInfinity(lpi, -oldlpbdchgs->bdchglbs[oldlpbdchgs->bdchgcolinds[c]]) && SCIPsetIsInfinity(set, -SCIPvarGetLbLP(var, set))) ||
            SCIPsetIsEQ(set, oldlpbdchgs->bdchglbs[oldlpbdchgs->bdchgcolinds[c]], SCIPvarGetLbLP(var, set)));
         assert((SCIPlpiIsInfinity(lpi, oldlpbdchgs->bdchgubs[oldlpbdchgs->bdchgcolinds[c]]) && SCIPsetIsInfinity(set, SCIPvarGetUbLP(var, set))) ||
            SCIPsetIsEQ(set, oldlpbdchgs->bdchgubs[oldlpbdchgs->bdchgcolinds[c]], SCIPvarGetUbLP(var, set)));

         /* store bound change for conflict analysis */
         if( !relaxedlpbdchgs->usedcols[c] )
         {
            idx = relaxedlpbdchgs->nbdchgs;
            relaxedlpbdchgs->usedcols[c] = TRUE;
            relaxedlpbdchgs->bdchgcolinds[c] = idx;
            relaxedlpbdchgs->nbdchgs++;

            /* remember the positive for later further bound widenings */
            relaxedlpbdchgs->bdchginds[idx] = c;
         }
         else
         {
            idx = relaxedlpbdchgs->bdchgcolinds[c];
            assert(relaxedlpbdchgs->bdchginds[idx] == c);

            /* the new bound should be the same or more relaxed */
            assert(relaxedlpbdchgs->bdchglbs[idx] >= newlb ||
               (SCIPlpiIsInfinity(lpi, -relaxedlpbdchgs->bdchglbs[idx]) && SCIPsetIsInfinity(set, -newlb)));
            assert(relaxedlpbdchgs->bdchgubs[idx] <= newub ||
               (SCIPlpiIsInfinity(lpi, relaxedlpbdchgs->bdchgubs[idx]) && SCIPsetIsInfinity(set, newub)));
         }

         /* set the new bounds for the LP with the correct infinity value */
         relaxedlpbdchgs->bdchglbs[idx] = SCIPsetIsInfinity(set, -newlb) ? -SCIPlpiInfinity(lpi) : newlb;
         relaxedlpbdchgs->bdchgubs[idx] = SCIPsetIsInfinity(set, newub) ? SCIPlpiInfinity(lpi) : newub;
         if( SCIPsetIsInfinity(set, -oldlpbdchgs->bdchglbs[idx]) )
            oldlpbdchgs->bdchglbs[idx] = -SCIPlpiInfinity(lpi);
         if( SCIPsetIsInfinity(set, oldlpbdchgs->bdchgubs[idx]) )
            oldlpbdchgs->bdchgubs[idx] = SCIPlpiInfinity(lpi);
      }
   }

   return SCIP_OKAY;
}

/** ensures, that candidate array can store at least num entries */
static
SCIP_RETCODE ensureCandsSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR***           cands,              /**< pointer to candidate array */
   SCIP_Real**           candscores,         /**< pointer to candidate score array */
   SCIP_Real**           newbounds,          /**< pointer to candidate new bounds array */
   SCIP_Real**           proofactdeltas,     /**< pointer to candidate proof delta array */
   int*                  candssize,          /**< pointer to size of array */
   int                   num                 /**< minimal number of candidates to store in array */
   )
{
   assert(cands != NULL);
   assert(candssize != NULL);

   if( num > *candssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_CALL( SCIPsetReallocBufferArray(set, cands, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, candscores, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, newbounds, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, proofactdeltas, newsize) );
      *candssize = newsize;
   }
   assert(num <= *candssize);

   return SCIP_OKAY;
}

/** adds variable to candidate list, if the current best bound corresponding to the proof coefficient is local;
 *  returns the array position in the candidate list, where the new candidate was inserted, or -1 if the
 *  variable can relaxed to global bounds immediately without increasing the proof's activity;
 *  the candidates are sorted with respect to the following two criteria:
 *  - prefer bound changes that have been applied deeper in the tree, to get a more global conflict
 *  - prefer variables with small Farkas coefficient to get rid of as many bound changes as possible
 */
static
SCIP_RETCODE addCand(
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_VAR*             var,                /**< variable to add to candidate array */
   int                   lbchginfopos,       /**< positions of currently active lower bound change information in variable's array */
   int                   ubchginfopos,       /**< positions of currently active upper bound change information in variable's array */
   SCIP_Real             proofcoef,          /**< coefficient of variable in infeasibility/bound proof */
   SCIP_Real             prooflhs,           /**< left hand side of infeasibility/bound proof */
   SCIP_Real             proofact,           /**< activity of infeasibility/bound proof row */
   SCIP_VAR***           cands,              /**< pointer to candidate array for undoing bound changes */
   SCIP_Real**           candscores,         /**< pointer to candidate score array for undoing bound changes */
   SCIP_Real**           newbounds,          /**< pointer to candidate new bounds array for undoing bound changes */
   SCIP_Real**           proofactdeltas,     /**< pointer to proof activity increase array for undoing bound changes */
   int*                  candssize,          /**< pointer to size of cands arrays */
   int*                  ncands,             /**< pointer to count number of candidates in bound change list */
   int                   firstcand           /**< position of first unprocessed bound change candidate */
   )
{
   SCIP_Real oldbound;
   SCIP_Real newbound;
   SCIP_Real proofactdelta;
   SCIP_Real score;
   int depth;
   int i;
   SCIP_Bool resolvable;

   assert(set != NULL);
   assert(var != NULL);
   assert(-1 <= lbchginfopos && lbchginfopos <= var->nlbchginfos);
   assert(-1 <= ubchginfopos && ubchginfopos <= var->nubchginfos);
   assert(!SCIPsetIsZero(set, proofcoef));
   assert(SCIPsetIsGT(set, prooflhs, proofact));
   assert(cands != NULL);
   assert(candscores != NULL);
   assert(newbounds != NULL);
   assert(proofactdeltas != NULL);
   assert(candssize != NULL);
   assert(ncands != NULL);
   assert(*ncands <= *candssize);
   assert(0 <= firstcand && firstcand <= *ncands);

   /* in the infeasibility or dual bound proof, the variable's bound is chosen to maximize the proof's activity */
   if( proofcoef > 0.0 )
   {
      assert(ubchginfopos >= 0); /* otherwise, undoBdchgsProof() should already have relaxed the local bound */

      /* calculate the difference of current bound to the previous bound the variable was set to */
      if( ubchginfopos == var->nubchginfos )
      {
         /* current bound is the strong branching or diving bound */
         oldbound = SCIPvarGetUbLP(var, set);
         newbound = SCIPvarGetUbLocal(var);
         depth = currentdepth+1;
         resolvable = FALSE;
      }
      else
      {
         /* current bound is the result of a local bound change */
         resolvable = bdchginfoIsResolvable(&var->ubchginfos[ubchginfopos]);
         depth = var->ubchginfos[ubchginfopos].bdchgidx.depth;
         oldbound = var->ubchginfos[ubchginfopos].newbound;
         newbound = var->ubchginfos[ubchginfopos].oldbound;
      }
   }
   else
   {
      assert(lbchginfopos >= 0); /* otherwise, undoBdchgsProof() should already have relaxed the local bound */

      /* calculate the difference of current bound to the previous bound the variable was set to */
      if( lbchginfopos == var->nlbchginfos )
      {
         /* current bound is the strong branching or diving bound */
         oldbound = SCIPvarGetLbLP(var, set);
         newbound = SCIPvarGetLbLocal(var);
         depth = currentdepth+1;
         resolvable = FALSE;
      }
      else
      {
         /* current bound is the result of a local bound change */
         resolvable = bdchginfoIsResolvable(&var->lbchginfos[lbchginfopos]);
         depth = var->lbchginfos[lbchginfopos].bdchgidx.depth;
         oldbound = var->lbchginfos[lbchginfopos].newbound;
         newbound = var->lbchginfos[lbchginfopos].oldbound;
      }
   }

   /* calculate the increase in the proof's activity */
   proofactdelta = (newbound - oldbound)*proofcoef;
   assert(proofactdelta > 0.0);


   /* calculate score for undoing the bound change */
   score = calcBdchgScore(prooflhs, proofact, proofactdelta, proofcoef, depth, currentdepth, var, set);

   if( !resolvable )
   {
      score += 10.0;
      if( !SCIPvarIsBinary(var) )
         score += 10.0;
   }

   /* get enough memory to store new candidate */
   SCIP_CALL( ensureCandsSize(set, cands, candscores, newbounds, proofactdeltas, candssize, (*ncands)+1) );
   assert(*cands != NULL);
   assert(*candscores != NULL);
   assert(*newbounds != NULL);
   assert(*proofactdeltas != NULL);

   SCIPsetDebugMsg(set, " -> local <%s> %s %g, relax <%s> %s %g, proofcoef=%g, dpt=%d, resolve=%u, delta=%g, score=%g\n",
      SCIPvarGetName(var), proofcoef > 0.0 ? "<=" : ">=", oldbound,
      SCIPvarGetName(var), proofcoef > 0.0 ? "<=" : ">=", newbound,
      proofcoef, depth, resolvable, proofactdelta, score);

   /* insert variable in candidate list without touching the already processed candidates */
   for( i = *ncands; i > firstcand && score > (*candscores)[i-1]; --i )
   {
      (*cands)[i] = (*cands)[i-1];
      (*candscores)[i] = (*candscores)[i-1];
      (*newbounds)[i] = (*newbounds)[i-1];
      (*proofactdeltas)[i] = (*proofactdeltas)[i-1];
   }
   (*cands)[i] = var;
   (*candscores)[i] = score;
   (*newbounds)[i] = newbound;
   (*proofactdeltas)[i] = proofactdelta;
   (*ncands)++;

   return SCIP_OKAY;
}

/** after changing the global bound of a variable, the bdchginfos that are now redundant are replaced with
 *  oldbound = newbound = global bound; if the current bdchginfo is of such kind, the bound is equal to the
 *  global bound and we can ignore it by installing a -1 as the corresponding bound change info position
 */
static
void skipRedundantBdchginfos(
   SCIP_VAR*             var,                /**< problem variable */
   int*                  lbchginfopos,       /**< pointer to lower bound change information position */
   int*                  ubchginfopos        /**< pointer to upper bound change information position */
   )
{
   assert(var != NULL);
   assert(lbchginfopos != NULL);
   assert(ubchginfopos != NULL);
   assert(-1 <= *lbchginfopos && *lbchginfopos <= var->nlbchginfos);
   assert(-1 <= *ubchginfopos && *ubchginfopos <= var->nubchginfos);
   assert(*lbchginfopos == -1 || *lbchginfopos == var->nlbchginfos
      || var->lbchginfos[*lbchginfopos].redundant
      == (var->lbchginfos[*lbchginfopos].oldbound == var->lbchginfos[*lbchginfopos].newbound)); /*lint !e777*/
   assert(*ubchginfopos == -1 || *ubchginfopos == var->nubchginfos
      || var->ubchginfos[*ubchginfopos].redundant
      == (var->ubchginfos[*ubchginfopos].oldbound == var->ubchginfos[*ubchginfopos].newbound)); /*lint !e777*/

   if( *lbchginfopos >= 0 && *lbchginfopos < var->nlbchginfos && var->lbchginfos[*lbchginfopos].redundant )
   {
      assert(SCIPvarGetLbGlobal(var) == var->lbchginfos[*lbchginfopos].oldbound); /*lint !e777*/
      *lbchginfopos = -1;
   }
   if( *ubchginfopos >= 0 && *ubchginfopos < var->nubchginfos && var->ubchginfos[*ubchginfopos].redundant )
   {
      assert(SCIPvarGetUbGlobal(var) == var->ubchginfos[*ubchginfopos].oldbound); /*lint !e777*/
      *ubchginfopos = -1;
   }
}

/** undoes bound changes on variables, still leaving the given infeasibility proof valid */
static
SCIP_RETCODE undoBdchgsProof(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_Real*            proofcoefs,         /**< coefficients in infeasibility proof */
   SCIP_Real             prooflhs,           /**< left hand side of proof */
   SCIP_Real*            proofact,           /**< current activity of proof */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   SCIP_LPBDCHGS*        oldlpbdchgs,        /**< old LP bound changes used for reset the LP bound change, or NULL */
   SCIP_LPBDCHGS*        relaxedlpbdchgs,    /**< relaxed LP bound changes used for reset the LP bound change, or NULL */
   SCIP_Bool*            resolve,            /**< pointer to store whether the changed LP should be resolved again, or NULL */
   SCIP_LPI*             lpi                 /**< pointer to LPi to access infinity of LP solver; necessary to set correct values */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** cands;
   SCIP_Real* candscores;
   SCIP_Real* newbounds;
   SCIP_Real* proofactdeltas;
   int nvars;
   int ncands;
   int candssize;
   int v;
   int i;

   assert(prob != NULL);
   assert(proofcoefs != NULL);
   assert(SCIPsetIsFeasGT(set, prooflhs, (*proofact)));
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);

   if( resolve != NULL )
      *resolve = FALSE;

   vars = prob->vars;
   nvars = prob->nvars;
   assert(nvars == 0 || vars != NULL);

   /* calculate the order in which the bound changes are tried to be undone, and relax all bounds if this doesn't
    * increase the proof's activity
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cands, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &candscores, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &newbounds, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &proofactdeltas, nvars) );
   ncands = 0;
   candssize = nvars;
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Bool relaxed;

      var = vars[v];

      /* after changing the global bound of a variable, the bdchginfos that are now redundant are replaced with
       * oldbound = newbound = global bound; if the current bdchginfo is of such kind, the bound is equal to the
       * global bound and we can ignore it
       */
      skipRedundantBdchginfos(var, &lbchginfoposs[v], &ubchginfoposs[v]);

      /* ignore variables already relaxed to global bounds */
      if( (lbchginfoposs[v] == -1 && ubchginfoposs[v] == -1) )
      {
         proofcoefs[v] = 0.0;
         continue;
      }

      /* relax bounds that are not used in the proof to the global bounds */
      relaxed = FALSE;
      if( !SCIPsetIsNegative(set, proofcoefs[v]) )
      {
         /* the lower bound is not used */
         if( lbchginfoposs[v] >= 0 )
         {
            SCIPsetDebugMsg(set, " -> relaxing variable <%s>[%g,%g] to [%g,%g]: proofcoef=%g, %g <= %g\n",
               SCIPvarGetName(var), curvarlbs[v], curvarubs[v], SCIPvarGetLbGlobal(var), curvarubs[v],
               proofcoefs[v], prooflhs, (*proofact));
            curvarlbs[v] = SCIPvarGetLbGlobal(var);
            lbchginfoposs[v] = -1;
            relaxed = TRUE;
         }
      }
      if( !SCIPsetIsPositive(set, proofcoefs[v]) )
      {
         /* the upper bound is not used */
         if( ubchginfoposs[v] >= 0 )
         {
            SCIPsetDebugMsg(set, " -> relaxing variable <%s>[%g,%g] to [%g,%g]: proofcoef=%g, %g <= %g\n",
               SCIPvarGetName(var), curvarlbs[v], curvarubs[v], curvarlbs[v], SCIPvarGetUbGlobal(var),
               proofcoefs[v], prooflhs, (*proofact));
            curvarubs[v] = SCIPvarGetUbGlobal(var);
            ubchginfoposs[v] = -1;
            relaxed = TRUE;
         }
      }
      if( relaxed && oldlpbdchgs != NULL )
      {
         SCIP_CALL( addBdchg(set, var, curvarlbs[v], curvarubs[v], oldlpbdchgs, relaxedlpbdchgs, lpi) );
      }

      /* add bound to candidate list */
      if( lbchginfoposs[v] >= 0 || ubchginfoposs[v] >= 0 )
      {
         SCIP_CALL( addCand(set, currentdepth, var, lbchginfoposs[v], ubchginfoposs[v], proofcoefs[v],
               prooflhs, (*proofact), &cands, &candscores, &newbounds, &proofactdeltas, &candssize, &ncands, 0) );
      }
      /* we can set the proof coefficient to zero, because the variable is not needed */
      else
         proofcoefs[v] = 0.0;
   }

   /* try to undo remaining local bound changes while still keeping the proof row violated:
    * bound changes can be undone, if prooflhs > proofact + proofactdelta;
    * afterwards, the current proof activity has to be updated
    */
   for( i = 0; i < ncands; ++i )
   {
      assert(proofactdeltas[i] > 0.0);
      assert((lbchginfoposs[SCIPvarGetProbindex(cands[i])] >= 0) != (ubchginfoposs[SCIPvarGetProbindex(cands[i])] >= 0));

      /* when relaxing a constraint we still need to stay infeasible; therefore we need to do the comparison in
       * feasibility tolerance because if 'prooflhs' is (feas-))equal to 'proofact + proofactdeltas[i]' it would mean
       * that there is no violation
       */
      if( SCIPsetIsFeasGT(set, prooflhs, (*proofact) + proofactdeltas[i]) )
      {
         v = SCIPvarGetProbindex(cands[i]);
         assert(0 <= v && v < nvars);
         assert((lbchginfoposs[v] >= 0) != (ubchginfoposs[v] >= 0));

         SCIPsetDebugMsg(set, " -> relaxing variable <%s>[%g,%g] to [%g,%g]: proofcoef=%g, %g <= %g + %g\n",
            SCIPvarGetName(cands[i]), curvarlbs[v], curvarubs[v],
            proofcoefs[v] > 0.0 ? curvarlbs[v] : newbounds[i],
            proofcoefs[v] > 0.0 ? newbounds[i] : curvarubs[v],
            proofcoefs[v], prooflhs, (*proofact), proofactdeltas[i]);

         assert((SCIPsetIsPositive(set, proofcoefs[v]) && SCIPsetIsGT(set, newbounds[i], curvarubs[v]))
            || (SCIPsetIsNegative(set, proofcoefs[v]) && SCIPsetIsLT(set, newbounds[i], curvarlbs[v])));
         assert((SCIPsetIsPositive(set, proofcoefs[v])
               && SCIPsetIsEQ(set, proofactdeltas[i], (newbounds[i] - curvarubs[v])*proofcoefs[v]))
            || (SCIPsetIsNegative(set, proofcoefs[v])
               && SCIPsetIsEQ(set, proofactdeltas[i], (newbounds[i] - curvarlbs[v])*proofcoefs[v])));
         assert(!SCIPsetIsZero(set, proofcoefs[v]));

         if( proofcoefs[v] > 0.0 )
         {
            assert(ubchginfoposs[v] >= 0);
            assert(lbchginfoposs[v] == -1);
            curvarubs[v] = newbounds[i];
            ubchginfoposs[v]--;
         }
         else
         {
            assert(lbchginfoposs[v] >= 0);
            assert(ubchginfoposs[v] == -1);
            curvarlbs[v] = newbounds[i];
            lbchginfoposs[v]--;
         }
         if( oldlpbdchgs != NULL )
         {
            SCIP_CALL( addBdchg(set, cands[i], curvarlbs[v], curvarubs[v], oldlpbdchgs, relaxedlpbdchgs, lpi) );
         }
         (*proofact) += proofactdeltas[i];
         if( resolve != NULL && SCIPvarIsInLP(cands[i]) )
            *resolve = TRUE;

         /* after changing the global bound of a variable, the bdchginfos that are now redundant are replaced with
          * oldbound = newbound = global bound; if the current bdchginfo is of such kind, the bound is equal to the
          * global bound and we can ignore it
          */
         skipRedundantBdchginfos(cands[i], &lbchginfoposs[v], &ubchginfoposs[v]);

         /* insert the new local bound of the variable into the candidate list */
         if( lbchginfoposs[v] >= 0 || ubchginfoposs[v] >= 0 )
         {
            SCIP_CALL( addCand(set, currentdepth, cands[i], lbchginfoposs[v], ubchginfoposs[v], proofcoefs[v],
                  prooflhs, (*proofact), &cands, &candscores, &newbounds, &proofactdeltas, &candssize, &ncands, i+1) );
         }
         else
            proofcoefs[v] = 0.0;
      }
   }

   /* free the buffer for the sorted bound change candidates */
   SCIPsetFreeBufferArray(set, &proofactdeltas);
   SCIPsetFreeBufferArray(set, &newbounds);
   SCIPsetFreeBufferArray(set, &candscores);
   SCIPsetFreeBufferArray(set, &cands);

   return SCIP_OKAY;
}

/* because calculations might cancel out some values, we stop the infeasibility analysis if a value is bigger than
 * 2^53 = 9007199254740992
 */
#define NUMSTOP 9007199254740992.0

/** analyzes an infeasible LP and undoes additional bound changes while staying infeasible */
static
SCIP_RETCODE undoBdchgsDualfarkas(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   SCIP_LPBDCHGS*        oldlpbdchgs,        /**< old LP bound changes used for reset the LP bound change, or NULL */
   SCIP_LPBDCHGS*        relaxedlpbdchgs,    /**< relaxed LP bound changes used for reset the LP bound change, or NULL */
   SCIP_Bool*            valid,              /**< pointer to store whether the unfixings are valid */
   SCIP_Bool*            resolve,            /**< pointer to store whether the changed LP should be resolved again */
   SCIP_Real*            farkascoefs,        /**< coefficients in the proof constraint */
   SCIP_Real             farkaslhs,          /**< lhs of the proof constraint */
   SCIP_Real*            farkasactivity      /**< maximal activity of the proof constraint */
   )
{
   SCIP_LPI* lpi;

   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);
   assert(valid != NULL);
   assert(resolve != NULL);

   SCIPsetDebugMsg(set, "undoing bound changes in infeasible LP: cutoff=%g\n", lp->cutoffbound);

   *valid = FALSE;
   *resolve = FALSE;

   lpi = SCIPlpGetLPI(lp);

   /* check, if the Farkas row is still violated (using current bounds and ignoring local rows) */
   if( SCIPsetIsFeasGT(set, farkaslhs, *farkasactivity) )
   {
      /* undo bound changes while keeping the infeasibility proof valid */
      SCIP_CALL( undoBdchgsProof(set, prob, currentdepth, farkascoefs, farkaslhs, farkasactivity, \
            curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs, oldlpbdchgs, relaxedlpbdchgs, resolve, lpi) );

      *valid = TRUE;

      /* resolving does not make sense: the old dual ray is still valid -> resolving will not change the solution */
      *resolve = FALSE;
   }

   return SCIP_OKAY;
}

/** analyzes an LP exceeding the objective limit and undoes additional bound changes while staying beyond the
 *  objective limit
 */
static
SCIP_RETCODE undoBdchgsDualsol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   SCIP_LPBDCHGS*        oldlpbdchgs,        /**< old LP bound changes used for reset the LP bound change, or NULL */
   SCIP_LPBDCHGS*        relaxedlpbdchgs,    /**< relaxed LP bound changes used for reset the LP bound change, or NULL */
   SCIP_Bool*            valid,              /**< pointer to store whether the unfixings are valid */
   SCIP_Bool*            resolve,            /**< pointer to store whether the changed LP should be resolved again */
   SCIP_Real*            dualcoefs,          /**< coefficients in the proof constraint */
   SCIP_Real             duallhs,            /**< lhs of the proof constraint */
   SCIP_Real*            dualactivity        /**< maximal activity of the proof constraint */
   )
{
   SCIP_LPI* lpi;

   assert(set != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);
   assert(valid != NULL);
   assert(resolve != NULL);

   *valid = FALSE;
   *resolve = FALSE;

   SCIPsetDebugMsg(set, "undoing bound changes in LP exceeding cutoff: cutoff=%g\n", lp->cutoffbound);

   /* get LP solver interface */
   lpi = SCIPlpGetLPI(lp);

   /* check, if the dual row is still violated (using current bounds and ignoring local rows) */
   if( SCIPsetIsFeasGT(set, duallhs, *dualactivity) )
   {
      /* undo bound changes while keeping the infeasibility proof valid */
      SCIP_CALL( undoBdchgsProof(set, prob, currentdepth, dualcoefs, duallhs, dualactivity, curvarlbs, curvarubs, \
            lbchginfoposs, ubchginfoposs, oldlpbdchgs, relaxedlpbdchgs, resolve, lpi) );

      *valid = TRUE;
   }

   return SCIP_OKAY;
}

/** applies conflict analysis starting with given bound changes, that could not be undone during previous
 *  infeasibility analysis
 */
static
SCIP_RETCODE conflictAnalyzeRemainingBdchgs(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   int*                  nconss,             /**< pointer to store the number of generated conflict constraints */
   int*                  nliterals,          /**< pointer to store the number of literals in generated conflict constraints */
   int*                  nreconvconss,       /**< pointer to store the number of generated reconvergence constraints */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence constraints */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_CONFTYPE conftype;
   SCIP_Bool usescutoffbound;
   int nvars;
   int v;
   int nbdchgs;
   int maxsize;

   assert(prob != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);
   assert(nconss != NULL);
   assert(nliterals != NULL);
   assert(nreconvconss != NULL);
   assert(nreconvliterals != NULL);

   *nconss = 0;
   *nliterals = 0;
   *nreconvconss = 0;
   *nreconvliterals = 0;

   vars = prob->vars;
   nvars = prob->nvars;
   assert(nvars == 0 || vars != NULL);

   maxsize = 2*conflictCalcMaxsize(set, prob);

   /* initialize conflict data */
   conftype = conflict->conflictset->conflicttype;
   usescutoffbound = conflict->conflictset->usescutoffbound;

   SCIP_CALL( SCIPconflictInit(conflict, set, stat, prob, conftype, usescutoffbound) );

   conflict->conflictset->conflicttype = conftype;
   conflict->conflictset->usescutoffbound = usescutoffbound;

   /* add remaining bound changes to conflict queue */
   SCIPsetDebugMsg(set, "initial conflict set after undoing bound changes:\n");

   nbdchgs = 0;
   for( v = 0; v < nvars && nbdchgs < maxsize; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      assert(var->nlbchginfos >= 0);
      assert(var->nubchginfos >= 0);
      assert(-1 <= lbchginfoposs[v] && lbchginfoposs[v] <= var->nlbchginfos);
      assert(-1 <= ubchginfoposs[v] && ubchginfoposs[v] <= var->nubchginfos);

      if( lbchginfoposs[v] == var->nlbchginfos || ubchginfoposs[v] == var->nubchginfos )
      {
         SCIP_BDCHGINFO* bdchginfo;
         SCIP_Real relaxedbd;

         /* the strong branching or diving bound stored in the column is responsible for the conflict:
          * it cannot be resolved and therefore has to be directly put into the conflict set
          */
         assert((lbchginfoposs[v] == var->nlbchginfos) != (ubchginfoposs[v] == var->nubchginfos)); /* only one can be tight in the dual! */
         assert(lbchginfoposs[v] < var->nlbchginfos || SCIPvarGetLbLP(var, set) > SCIPvarGetLbLocal(var));
         assert(ubchginfoposs[v] < var->nubchginfos || SCIPvarGetUbLP(var, set) < SCIPvarGetUbLocal(var));

         /* create an artificial bound change information for the diving/strong branching bound change;
          * they are freed in the SCIPconflictFlushConss() call
          */
         if( lbchginfoposs[v] == var->nlbchginfos )
         {
            SCIP_CALL( conflictCreateTmpBdchginfo(conflict, blkmem, set, var, SCIP_BOUNDTYPE_LOWER,
                  SCIPvarGetLbLocal(var), SCIPvarGetLbLP(var, set), &bdchginfo) );
            relaxedbd = SCIPvarGetLbLP(var, set);
         }
         else
         {
            SCIP_CALL( conflictCreateTmpBdchginfo(conflict, blkmem, set, var, SCIP_BOUNDTYPE_UPPER,
                  SCIPvarGetUbLocal(var), SCIPvarGetUbLP(var, set), &bdchginfo) );
            relaxedbd = SCIPvarGetUbLP(var, set);
         }

         /* put variable into the conflict set */
         SCIPsetDebugMsg(set, "   force: <%s> %s %g [status: %d, type: %d, dive/strong]\n",
            SCIPvarGetName(var), lbchginfoposs[v] == var->nlbchginfos ? ">=" : "<=",
            lbchginfoposs[v] == var->nlbchginfos ? SCIPvarGetLbLP(var, set) : SCIPvarGetUbLP(var, set),
            SCIPvarGetStatus(var), SCIPvarGetType(var));
         SCIP_CALL( conflictAddConflictBound(conflict, blkmem, set, bdchginfo, relaxedbd) );

         /* each variable which is add to the conflict graph gets an increase in the VSIDS
          *
          * @note That is different to the VSIDS preseted in the literature
          */
         SCIP_CALL( incVSIDS(var, blkmem, set, stat, SCIPbdchginfoGetBoundtype(bdchginfo), relaxedbd, set->conf_conflictgraphweight) );
         nbdchgs++;
      }
      else
      {
         /* put remaining bound changes into conflict candidate queue */
         if( lbchginfoposs[v] >= 0 )
         {
            SCIP_CALL( conflictAddBound(conflict, blkmem, set, stat, var, SCIP_BOUNDTYPE_LOWER, \
                  &var->lbchginfos[lbchginfoposs[v]], SCIPbdchginfoGetNewbound(&var->lbchginfos[lbchginfoposs[v]])) );
            nbdchgs++;
         }
         if( ubchginfoposs[v] >= 0 )
         {
            assert(!SCIPbdchginfoIsRedundant(&var->ubchginfos[ubchginfoposs[v]]));
            SCIP_CALL( conflictAddBound(conflict, blkmem, set, stat, var, SCIP_BOUNDTYPE_UPPER, \
                  &var->ubchginfos[ubchginfoposs[v]], SCIPbdchginfoGetNewbound(&var->ubchginfos[ubchginfoposs[v]])) );
            nbdchgs++;
         }
      }
   }

   if( v == nvars )
   {
      /* analyze the conflict set, and create conflict constraints on success */
      SCIP_CALL( conflictAnalyze(conflict, blkmem, set, stat, prob, tree, diving, 0, FALSE, nconss, nliterals, \
            nreconvconss, nreconvliterals) );
   }

   return SCIP_OKAY;
}

/** calculates a Farkas proof from the current dual LP solution */
static
SCIP_RETCODE getFarkasProof(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPI*             lpi,                /**< LPI data */
   SCIP_AGGRROW*         farkasrow,          /**< aggregated row representing the proof */
   SCIP_Real*            farkasact,          /**< maximal activity of the proof constraint */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_Bool*            valid               /**< pointer store whether the proof constraint is valid */
   )
{
   SCIP_ROW** rows;
   SCIP_Real* dualfarkas;
   SCIP_ROW* row;
   int nrows;
   int r;

   assert(SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsDualFeasible(lpi));
   assert(SCIPlpiIsPrimalInfeasible(lpi) || !SCIPlpDivingObjChanged(lp));

   /* get LP rows and problem variables */
   rows = SCIPlpGetRows(lp);
   nrows = SCIPlpGetNRows(lp);
   assert(nrows == 0 || rows != NULL);
   assert(nrows == lp->nlpirows);

   /* it can happen that infeasibility is detetected within LP presolve. in that case, the LP solver may not be able to
    * to return the dual ray.
    */
   if( !SCIPlpiHasDualRay(lpi) )
   {
      *valid = FALSE;
      return SCIP_OKAY;
   }

   assert(farkasrow != NULL);

   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualfarkas, nrows) );
   BMSclearMemoryArray(dualfarkas, nrows);

   /* get dual Farkas values of rows */
   SCIP_CALL( SCIPlpiGetDualfarkas(lpi, dualfarkas) );

   /* calculate the Farkas row */
   (*valid) = TRUE;
   for( r = 0; r < nrows; ++r )
   {
      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->vals != NULL);
      assert(row == lp->lpirows[r]);

      /* ignore local rows and rows with Farkas value 0.0 */
      if( !row->local && !SCIPsetIsDualfeasZero(set, dualfarkas[r]) )
      {
#ifndef NDEBUG
         {
            SCIP_Real lpilhs;
            SCIP_Real lpirhs;

            SCIP_CALL( SCIPlpiGetSides(lpi, r, r, &lpilhs, &lpirhs) );
            assert((SCIPsetIsInfinity(set, -lpilhs) && SCIPsetIsInfinity(set, -row->lhs))
               || SCIPsetIsRelEQ(set, lpilhs, row->lhs - row->constant));
            assert((SCIPsetIsInfinity(set, lpirhs) && SCIPsetIsInfinity(set, row->rhs))
               || SCIPsetIsRelEQ(set, lpirhs, row->rhs - row->constant));
         }
#endif

         /* add row side to Farkas row lhs: dualfarkas > 0 -> lhs, dualfarkas < 0 -> rhs */
         if( dualfarkas[r] > 0.0 )
         {
            /* check if sign of dual Farkas value is valid */
            if( SCIPsetIsInfinity(set, -row->lhs) )
               continue;

            /* due to numerical reasons we want to stop */
            if( REALABS(dualfarkas[r] * (row->lhs - row->constant)) > NUMSTOP )
            {
               (*valid) = FALSE;
               goto TERMINATE;
            }

            SCIP_CALL( SCIPaggrRowAddRow(set->scip, farkasrow, row, -dualfarkas[r], -1) );
         }
         else
         {
            /* check if sign of dual Farkas value is valid */
            if( SCIPsetIsInfinity(set, row->rhs) )
               continue;

            /* due to numerical reasons we want to stop */
            if( REALABS(dualfarkas[r] * (row->rhs - row->constant)) > NUMSTOP )
            {
               (*valid) = FALSE;
               goto TERMINATE;
            }

            SCIP_CALL( SCIPaggrRowAddRow(set->scip, farkasrow, row, -dualfarkas[r], +1) );
         }
         SCIPsetDebugMsg(set, " -> farkasrhs: %g<%s>[%g,%g] -> %g\n", dualfarkas[r], SCIProwGetName(row),
            row->lhs - row->constant, row->rhs - row->constant, SCIPaggrRowGetRhs(farkasrow));

         /* due to numerical reasons we want to stop */
         if( REALABS(SCIPaggrRowGetRhs(farkasrow)) > NUMSTOP )
         {
            (*valid) = FALSE;
            goto TERMINATE;
         }
      }
#ifdef SCIP_DEBUG
      else if( !SCIPsetIsZero(set, dualfarkas[r]) )
      {
         SCIPsetDebugMsg(set, " -> ignoring %s row <%s> with dual Farkas value %.10f (lhs=%g, rhs=%g)\n",
            row->local ? "local" : "global", SCIProwGetName(row), dualfarkas[r],
            row->lhs - row->constant, row->rhs - row->constant);
      }
#endif
   }

   /* remove all coefficients that are too close to zero */
   SCIPaggrRowRemoveZeros(set->scip, farkasrow, valid);

   if( !(*valid) )
      goto TERMINATE;

   /* calculate the current Farkas activity, always using the best bound w.r.t. the Farkas coefficient */
   *farkasact = aggrRowGetMinActivity(prob, farkasrow, curvarlbs, curvarubs);

   SCIPsetDebugMsg(set, " -> farkasact=%g farkasrhs=%g, \n", (*farkasact), SCIPaggrRowGetRhs(farkasrow));

   /* the constructed proof is not valid, this can happen due to numerical reasons,
    * e.g., we only consider rows r with !SCIPsetIsZero(set, dualfarkas[r])
    */
   if( SCIPsetIsFeasLE(set, *farkasact, SCIPaggrRowGetRhs(farkasrow)) )
   {
      (*valid) = FALSE;
      SCIPsetDebugMsg(set, " -> proof is not valid: %g <= %g\n", *farkasact, SCIPaggrRowGetRhs(farkasrow));
   }

TERMINATE:
   SCIPsetFreeBufferArray(set, &dualfarkas);

   return SCIP_OKAY;
}

/** calculates a Farkas proof from the current dual LP solution */
static
SCIP_RETCODE getDualProof(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPI*             lpi,                /**< LPI data */
   SCIP_AGGRROW*         farkasrow,          /**< aggregated row representing the proof */
   SCIP_Real*            farkasact,          /**< maximal activity of the proof constraint */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_Bool*            valid               /**< pointer store whether the proof constraint is valid */
   )
{
   SCIP_RETCODE retcode;
   SCIP_ROW** rows;
   SCIP_ROW* row;
   SCIP_Real* primsols;
   SCIP_Real* dualsols;
   SCIP_Real* redcosts;
   SCIP_Real maxabsdualsol;
   int nrows;
   int ncols;
   int r;

   assert(set != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(valid != NULL);

   *valid = TRUE;

   /* get LP rows and problem variables */
   rows = SCIPlpGetRows(lp);
   nrows = SCIPlpGetNRows(lp);
   ncols = SCIPlpGetNCols(lp);
   assert(nrows == 0 || rows != NULL);
   assert(nrows == lp->nlpirows);

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &primsols, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualsols, nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &redcosts, ncols) );

   /* get solution from LPI */
   retcode = SCIPlpiGetSol(lpi, NULL, primsols, dualsols, NULL, redcosts);
   if( retcode == SCIP_LPERROR ) /* on an error in the LP solver, just abort the conflict analysis */
   {
      (*valid) = FALSE;
      goto TERMINATE;
   }
   SCIP_CALL( retcode );
#ifdef SCIP_DEBUG
   {
      SCIP_Real objval;
      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      SCIPsetDebugMsg(set, " -> LP objval: %g\n", objval);
   }
#endif

   /* check whether the dual solution is numerically stable */
   maxabsdualsol = 0;
   for( r = 0; r < nrows; r++ )
   {
      SCIP_Real absdualsol = REALABS(dualsols[r]);

      if( absdualsol > maxabsdualsol )
         maxabsdualsol = absdualsol;
   }

   /* don't consider dual solution with maxabsdualsol > 1e+07, this would almost cancel out the objective constraint */
   if( maxabsdualsol > 1e+07 )
   {
      (*valid) = FALSE;
      goto TERMINATE;
   }

   /* clear the proof */
   SCIPaggrRowClear(farkasrow);

   /* Let y be the dual solution and r be the reduced cost vector. Let z be defined as
    *    z_i := y_i if i is a global row,
    *    z_i := 0   if i is a local row.
    * Define the set X := {x | lhs <= Ax <= rhs, lb <= x <= ub, c^Tx <= c*}, with c* being the current primal bound.
    * Then the following inequalities are valid for all x \in X:
    *                                 - c* <= -c^Tx
    *   <=>                     z^TAx - c* <= (z^TA - c^T) x
    *   <=>                     z^TAx - c* <= (y^TA - c^T - (y-z)^TA) x
    *   <=>                     z^TAx - c* <= (-r^T - (y-z)^TA) x         (dual feasibility of (y,r): y^TA + r^T == c^T)
    * Because lhs <= Ax <= rhs and lb <= x <= ub, the inequality can be relaxed to give
    *     min{z^Tq | lhs <= q <= rhs} - c* <= max{(-r^T - (y-z)^TA) x | lb <= x <= ub}, or X = {}.
    *
    * The resulting dual row is:  z^T{lhs,rhs} - c* <= (-r^T - (y-z)^TA){lb,ub},
    * where lhs, rhs, lb, and ub are selected in order to maximize the feasibility of the row.
    */

   /* add the objective function to the aggregation row with current cutoff bound as right-hand side
    *
    * use a slightly tighter cutoff bound, because solutions with equal objective value should also be declared
    * infeasible
    */
   SCIP_CALL( SCIPaggrRowAddObjectiveFunction(set->scip, farkasrow, lp->cutoffbound - SCIPsetSumepsilon(set), 1.0) );

   /* dual row: z^T{lhs,rhs} - c* <= (-r^T - (y-z)^TA){lb,ub}
    * process rows: add z^T{lhs,rhs} to the dual row's left hand side, and -(y-z)^TA to the dual row's coefficients
    */
   for( r = 0; r < nrows; ++r )
   {
      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->vals != NULL);
      assert(row == lp->lpirows[r]);

      /* ignore dual solution values of 0.0 (in this case: y_i == z_i == 0) */
      if( SCIPsetIsDualfeasZero(set, dualsols[r]) )
         continue;

      /* check dual feasibility */
      if( (SCIPsetIsInfinity(set, -row->lhs) && dualsols[r] > 0.0) || (SCIPsetIsInfinity(set, row->rhs) && dualsols[r] < 0.0) )
      {
         SCIPsetDebugMsg(set, " -> infeasible dual solution %g in row <%s>: lhs=%g, rhs=%g\n",
            dualsols[r], SCIProwGetName(row), row->lhs, row->rhs);
         (*valid) = FALSE;
         goto TERMINATE;
      }

      /* local rows add up to the dual row's coefficients (because z_i == 0 => -(y_i - z_i) == -y_i),
       * global rows add up to the dual row's left hand side (because z_i == y_i != 0)
       */
      if( row->local )
      {
#if 0
         /* add -y_i A_i to coefficients of dual row */
         for( i = 0; i < row->len; ++i )
         {
            v = SCIPvarGetProbindex(SCIPcolGetVar(row->cols[i]));
            assert(0 <= v && v < prob->nvars);
            dualcoefs[v] -= dualsols[r] * row->vals[i]; // TODO use row
         }
#endif
         SCIPsetDebugMsg(set, " -> local row <%s>: dual=%g\n", SCIProwGetName(row), dualsols[r]);
      }
      else
      {
         /* add minimal value to dual row's left hand side: z_i == y_i > 0 -> lhs, z_i == y_i < 0 -> rhs */
         if( dualsols[r] > 0.0 )
         {
            assert(!SCIPsetIsInfinity(set, -row->lhs));
            SCIP_CALL( SCIPaggrRowAddRow(set->scip, farkasrow, row, -dualsols[r], -1) );
         }
         else
         {
            assert(!SCIPsetIsInfinity(set, row->rhs));
            SCIP_CALL( SCIPaggrRowAddRow(set->scip, farkasrow, row, -dualsols[r], +1) );
         }
         SCIPsetDebugMsg(set, " -> global row <%s>[%g,%g]: dual=%g -> dualrhs=%g\n",
            SCIProwGetName(row), row->lhs - row->constant, row->rhs - row->constant, -dualsols[r], SCIPaggrRowGetRhs(farkasrow));
      }
   }

   /* remove all nearly zero coefficients */
   SCIPaggrRowRemoveZeros(set->scip, farkasrow, valid);

   if( !(*valid) )
      goto TERMINATE;

   /* check validity of the proof */
   *farkasact = aggrRowGetMinActivity(prob, farkasrow, curvarlbs, curvarubs);

   if( SCIPsetIsLE(set, *farkasact, SCIPaggrRowGetRhs(farkasrow)) )
   {
      *valid = FALSE;
      SCIPsetDebugMsg(set, " -> proof is not valid: %g <= %g\n", *farkasact, SCIPaggrRowGetRhs(farkasrow));
   }

  TERMINATE:
   SCIPsetFreeBufferArray(set, &redcosts);
   SCIPsetFreeBufferArray(set, &dualsols);
   SCIPsetFreeBufferArray(set, &primsols);

   return SCIP_OKAY;
}

#ifdef SCIP_DEBUG
static
void debugPrintViolationInfo(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             minact,             /**< min activity */
   SCIP_Real             rhs,                /**< right hand side */
   const char*           infostr             /**< additional info for this debug message, or NULL */
   )
{
   SCIPsetDebugMsg(set, "-> %sminact=%.15g rhs=%.15g violation=%.15g\n",infostr != NULL ? infostr : "" , minact, rhs, minact - rhs);
}
#else
#define debugPrintViolationInfo(...) /**/
#endif

/** apply coefficient tightening */
static
void tightenCoefficients(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROOFSET*        proofset,           /**< proof set */
   int*                  nchgcoefs,          /**< pointer to store number of changed coefficients */
   SCIP_Bool*            redundant           /**< pointer to store whether the proof set is redundant */
   )
{
#ifdef SCIP_DEBUG
   SCIP_Real absmax = 0.0;
   SCIP_Real absmin = SCIPsetInfinity(set);
   int i;

   for( i = 0; i < proofset->nnz; i++ )
   {
      absmax = MAX(absmax, REALABS(proofset->vals[i]));
      absmin = MIN(absmin, REALABS(proofset->vals[i]));
   }
#endif

   (*redundant) = SCIPcutsTightenCoefficients(set->scip, FALSE, proofset->vals, &proofset->rhs, proofset->inds, &proofset->nnz, nchgcoefs);

#ifdef SCIP_DEBUG
   {
      SCIP_Real newabsmax = 0.0;
      SCIP_Real newabsmin = SCIPsetInfinity(set);

      for( i = 0; i < proofset->nnz; i++ )
      {
         newabsmax = MAX(newabsmax, REALABS(proofset->vals[i]));
         newabsmin = MIN(newabsmin, REALABS(proofset->vals[i]));
      }

      SCIPsetDebugMsg(set, "coefficient tightening: [%.15g,%.15g] -> [%.15g,%.15g] (nnz: %d, nchg: %d rhs: %.15g)\n",
            absmin, absmax, newabsmin, newabsmax, proofsetGetNVars(proofset), *nchgcoefs, proofsetGetRhs(proofset));
      printf("coefficient tightening: [%.15g,%.15g] -> [%.15g,%.15g] (nnz: %d, nchg: %d rhs: %.15g)\n",
            absmin, absmax, newabsmin, newabsmax, proofsetGetNVars(proofset), *nchgcoefs, proofsetGetRhs(proofset));
   }
#endif
}

/** try to generate alternative proofs by applying subadditive functions */
static
SCIP_RETCODE separateAlternativeProofs(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< tree data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_AGGRROW*         proofrow,           /**< proof rows data */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_CONFTYPE         conflicttype        /**< type of the conflict */
   )
{
   SCIP_VAR** vars;
   SCIP_SOL* refsol;
   SCIP_Real* cutcoefs;
   SCIP_Real cutefficacy;
   SCIP_Real cutrhs;
   SCIP_Real proofefficiacy;
   SCIP_Real efficiacynorm;
   SCIP_Bool islocal;
   SCIP_Bool cutsuccess;
   SCIP_Bool success;
   int* cutinds;
   int* inds;
   int cutnnz;
   int nnz;
   int nvars;
   int i;

   vars = SCIPprobGetVars(transprob);
   nvars = SCIPprobGetNVars(transprob);

   inds = SCIPaggrRowGetInds(proofrow);
   nnz = SCIPaggrRowGetNNz(proofrow);

   proofefficiacy = aggrRowGetMinActivity(transprob, proofrow, curvarlbs, curvarubs) - SCIPaggrRowGetRhs(proofrow);

   efficiacynorm = SCIPaggrRowCalcEfficacyNorm(set->scip, proofrow);
   proofefficiacy /= MAX(1e-6, efficiacynorm);

   /* create reference solution */
   SCIP_CALL( SCIPcreateSol(set->scip, &refsol, NULL) );

   /* initialize with average solution */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPsolSetVal(refsol, set, stat, tree, vars[i], SCIPvarGetAvgSol(vars[i])) );
   }

   /* set all variables that are part of the proof to its active local bound */
   for( i = 0; i < nnz; i++ )
   {
      SCIP_Real val = SCIPaggrRowGetProbvarValue(proofrow, inds[i]);

      if( val > 0.0 )
      {
         SCIP_CALL( SCIPsolSetVal(refsol, set, stat, tree, vars[inds[i]], curvarubs[inds[i]]) );
      }
      else
      {
         SCIP_CALL( SCIPsolSetVal(refsol, set, stat, tree, vars[inds[i]], curvarlbs[inds[i]]) );
      }
   }

   SCIP_CALL( SCIPsetAllocBufferArray(set, &cutcoefs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cutinds, nvars) );

   cutnnz = 0;
   cutefficacy = -SCIPsetInfinity(set);

   /* apply flow cover */
   SCIP_CALL( SCIPcalcFlowCover(set->scip, refsol, POSTPROCESS, BOUNDSWITCH, ALLOWLOCAL, proofrow, \
         cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, NULL, &islocal, &cutsuccess) );
   success = cutsuccess;

   /* apply MIR */
   SCIP_CALL( SCIPcutGenerationHeuristicCMIR(set->scip, refsol, POSTPROCESS, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, INT_MAX, \
         NULL, NULL, MINFRAC, MAXFRAC, proofrow, cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, NULL, \
         &islocal, &cutsuccess) );
   success = (success || cutsuccess);

   /* replace the current proof */
   if( success && !islocal && SCIPsetIsPositive(set, cutefficacy) && cutefficacy * nnz > proofefficiacy * cutnnz )
   {
      SCIP_PROOFSET* alternativeproofset;
      SCIP_Bool redundant;
      int nchgcoefs;

      SCIP_CALL( proofsetCreate(&alternativeproofset, blkmem) );
      alternativeproofset->conflicttype = (conflicttype == SCIP_CONFTYPE_INFEASLP ? SCIP_CONFTYPE_ALTINFPROOF : SCIP_CONFTYPE_ALTBNDPROOF);

      SCIP_CALL( proofsetAddSparseData(alternativeproofset, blkmem, cutcoefs, cutinds, cutnnz, cutrhs) );

      /* apply coefficient tightening */
      tightenCoefficients(set, alternativeproofset, &nchgcoefs, &redundant);

      if( !redundant )
      {
         SCIP_CALL( conflictInsertProofset(conflict, set, alternativeproofset) );
      }
      else
      {
         proofsetFree(&alternativeproofset, blkmem);
      }
   }

   SCIPsetFreeBufferArray(set, &cutinds);
   SCIPsetFreeBufferArray(set, &cutcoefs);

   SCIP_CALL( SCIPfreeSol(set->scip, &refsol) );

   return SCIP_OKAY;
}

/** tighten a given infeasibility proof a^Tx <= b with minact > b w.r.t. local bounds
 *
 *  1) Apply cut generating functions
 *    - c-MIR
 *    - Flow-cover
 *    - TODO: implement other subadditive functions
 *  2) Remove continuous variables contributing with its global bound
 *    - TODO: implement a variant of non-zero-cancellation
 */
static
SCIP_RETCODE tightenDualproof(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< tree data */
   SCIP_AGGRROW*         proofrow,           /**< aggregated row representing the proof */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_Bool             initialproof        /**< do we analyze the initial reason of infeasibility? */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int* inds;
   SCIP_PROOFSET* proofset;
   SCIP_Bool valid;
   SCIP_Bool redundant;
   int nnz;
   int nchgcoefs;
   int nbinvars;
   int ncontvars;
   int nintvars;
   int i;

   assert(conflict->proofset != NULL);

   vars = SCIPprobGetVars(transprob);
   nbinvars = 0;
   nintvars = 0;
   ncontvars = 0;

   inds = SCIPaggrRowGetInds(proofrow);
   nnz = SCIPaggrRowGetNNz(proofrow);

   /* count number of binary, integer, and continuous variables */
   for( i = 0; i < nnz; i++ )
   {
      assert(SCIPvarGetProbindex(vars[inds[i]]) == inds[i]);

      if( SCIPvarIsBinary(vars[inds[i]]) )
         ++nbinvars;
      else if( SCIPvarIsIntegral(vars[inds[i]]) )
         ++nintvars;
      else
         ++ncontvars;
   }

   SCIPsetDebugMsg(set, "start dualray tightening:\n");
   SCIPsetDebugMsg(set, "-> tighten dual ray: nvars=%d (bin=%d, int=%d, cont=%d)\n",
         nnz, nbinvars, nintvars, ncontvars);
   debugPrintViolationInfo(set, aggrRowGetMinActivity(transprob, proofrow, curvarlbs, curvarubs), SCIPaggrRowGetRhs(proofrow), NULL);

   /* try to find an alternative proof of local infeasibility that is stronger */
   if( set->conf_sepaaltproofs )
   {
      SCIP_CALL( separateAlternativeProofs(conflict, set, stat, transprob, tree, blkmem, proofrow, curvarlbs, curvarubs,
            conflict->conflictset->conflicttype) );
   }

   if( initialproof )
      proofset = conflict->proofset;
   else
   {
      SCIP_CALL( proofsetCreate(&proofset, blkmem) );
      SCIP_CALL( conflictInsertProofset(conflict, set, proofset) );
   }

   /* start with a proofset containing all variables with a non-zero coefficient in the dual proof */
   SCIP_CALL( proofsetAddAggrrow(proofset, set, blkmem, proofrow) );
   proofset->conflicttype = conflict->conflictset->conflicttype;

   /* get proof data */
   vals = proofsetGetVals(proofset);
   inds = proofsetGetInds(proofset);
   nnz = proofsetGetNVars(proofset);

   /* remove continuous variable contributing with their global bound
    *
    * todo: check whether we also want to do that for bound exceeding proofs, but then we cannot update the
    *       conflict anymore
    */
   if( proofset->conflicttype == SCIP_CONFTYPE_INFEASLP )
   {
      /* remove all continuous variables that have equal global and local bounds (ub or lb depend on the sign)
       * from the proof
       */

      for( i = 0; i < nnz && nnz > 1; )
      {
         SCIP_Real val;
         int idx = inds[i];

         assert(vars[idx] != NULL);

         val = vals[i];
         assert(!SCIPsetIsZero(set, val));

         /* skip integral variables */
         if( SCIPvarGetType(vars[idx]) != SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(vars[idx]) != SCIP_VARTYPE_IMPLINT )
         {
            i++;
            continue;
         }
         else
         {
            SCIP_Real glbbd;
            SCIP_Real locbd;

            /* get appropriate global and local bounds */
            glbbd = (val < 0.0 ? SCIPvarGetUbGlobal(vars[idx]) : SCIPvarGetLbGlobal(vars[idx]));
            locbd = (val < 0.0 ? curvarubs[idx] : curvarlbs[idx]);

            if( !SCIPsetIsEQ(set, glbbd, locbd) )
            {
               i++;
               continue;
            }

            SCIPsetDebugMsg(set, "-> remove continuous variable <%s>: glb=[%g,%g], loc=[%g,%g], val=%g\n",
                  SCIPvarGetName(vars[idx]), SCIPvarGetLbGlobal(vars[idx]), SCIPvarGetUbGlobal(vars[idx]),
                  SCIPvarGetLbLocal(vars[idx]), SCIPvarGetUbLocal(vars[idx]), val);

            proofsetCancelVarWithBound(proofset, set, vars[idx], i, &valid);
            assert(valid); /* this should be always fulfilled at this place */

            --nnz;
         }
      }
   }

   /* apply coefficient tightening to initial proof */
   tightenCoefficients(set, proofset, &nchgcoefs, &redundant);

   /* it can happen that the constraints is almost globally redundant w.r.t to the maximal activity,
    * e.g., due to numerics. in this case, we want to discard the proof
    */
   if( redundant )
   {
#ifndef NDEBUG
      SCIP_Real eps = MIN(0.01, 10.0*set->num_feastol);
      assert(proofset->rhs - getMaxActivity(transprob, proofset->vals, proofset->inds, proofset->nnz, NULL, NULL) < eps);
#endif
      proofsetClear(proofset);
   }
   else if( nchgcoefs > 0 )
   {
      if( proofset->conflicttype == SCIP_CONFTYPE_INFEASLP )
         proofset->conflicttype = SCIP_CONFTYPE_ALTINFPROOF;
      else if( proofset->conflicttype == SCIP_CONFTYPE_BNDEXCEEDING )
         proofset->conflicttype = SCIP_CONFTYPE_ALTBNDPROOF;
   }

   return SCIP_OKAY;
}

/** perform conflict analysis based on a dual unbounded ray
 *
 *  given an aggregation of rows lhs <= a^Tx such that lhs > maxactivity. if the constraint has size one we add a
 *  bound change instead of the constraint.
 */
static
SCIP_RETCODE conflictAnalyzeDualProof(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< tree data */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_AGGRROW*         proofrow,           /**< aggregated row representing the proof */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_Bool             initialproof,       /**< do we analyze the initial reason of infeasibility? */
   SCIP_Bool*            globalinfeasible,   /**< pointer to store whether global infeasibility could be proven */
   SCIP_Bool*            success             /**< pointer to store success result */
   )
{
   SCIP_Real rhs;
   SCIP_Real minact;
   int nnz;

   assert(set != NULL);
   assert(transprob != NULL);

   /* get sparse data */
   nnz = SCIPaggrRowGetNNz(proofrow);
   rhs = SCIPaggrRowGetRhs(proofrow);

   *globalinfeasible = FALSE;
   *success = FALSE;

   /* get minimal activity w.r.t. local bounds */
   minact = aggrRowGetMinActivity(transprob, proofrow, curvarlbs, curvarubs);

   /* only run is the proof proves local infeasibility */
   if( SCIPsetIsFeasLE(set, minact, rhs) )
      return SCIP_OKAY;

   /* if the farkas-proof is empty, the node and its sub tree can be cut off completely */
   if( nnz == 0 )
   {
      SCIPsetDebugMsg(set, " -> empty farkas-proof in depth %d cuts off sub tree at depth 0\n", SCIPtreeGetFocusDepth(tree));

      SCIP_CALL( SCIPnodeCutoff(tree->path[0], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );

      *globalinfeasible = TRUE;
      *success = TRUE;

      ++conflict->ndualrayinfsuccess;

      return SCIP_OKAY;
   }

   /* try to enforce the constraint based on a dual ray */
   SCIP_CALL( tightenDualproof(conflict, set, stat, blkmem, transprob, tree, proofrow, curvarlbs, curvarubs, initialproof) );

   if( *globalinfeasible )
   {
      SCIPsetDebugMsg(set, "detect global: cutoff root node\n");
      SCIP_CALL( SCIPnodeCutoff(tree->path[0], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
      *success = TRUE;

      ++conflict->ndualrayinfsuccess;
   }

   return SCIP_OKAY;
}

/** try to find a subset of changed bounds leading to an infeasible LP
 *
 *  1. call undoBdchgsDualfarkas() or undoBdchgsDualsol()
 *     -> update lb/ubchginfoposs arrays
 *     -> store additional changes in bdchg and curvarlbs/ubs arrays
 *     -> apply additional changes to the LPI
 *  2. (optional) if additional bound changes were undone:
 *     -> resolve LP
 *     -> goto 1.
 *  3. redo all bound changes in the LPI to restore the LPI to its original state
 *  4. analyze conflict
 *     -> put remaining changed bounds (see lb/ubchginfoposs arrays) into starting conflict set
 */
static
SCIP_RETCODE runBoundHeuristic(
   SCIP_CONFLICT*        conflict,           /**< conflict data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPI*             lpi,                /**< LPI data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Real*            proofcoefs,         /**< coefficients in the proof constraint */
   SCIP_Real*            prooflhs,           /**< lhs of the proof constraint */
   SCIP_Real*            proofactivity,      /**< maximal activity of the proof constraint */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   int*                  iterations,         /**< pointer to store the total number of LP iterations used */
   SCIP_Bool             marklpunsolved,     /**< whether LP should be marked unsolved after analysis (needed for strong branching) */
   SCIP_Bool*            dualraysuccess,     /**< pointer to store success result of dualray analysis */
   SCIP_Bool*            valid               /**< pointer to store whether the result is still a valid proof */
   )
{
   SCIP_LPBDCHGS* oldlpbdchgs;
   SCIP_LPBDCHGS* relaxedlpbdchgs;
   SCIP_Bool solvelp;
   SCIP_Bool resolve;
   int ncols;

   assert(set != NULL);

   /* get number of columns in the LP */
   ncols = SCIPlpGetNCols(lp);

   /* get temporary memory for remembering bound changes on LPI columns */
   SCIP_CALL( lpbdchgsCreate(&oldlpbdchgs, set, ncols) );
   SCIP_CALL( lpbdchgsCreate(&relaxedlpbdchgs, set, ncols) );

   /* undo as many bound changes as possible with the current LP solution */
   resolve = FALSE;
   if( (*valid) )
   {
      int currentdepth;
      currentdepth = SCIPtreeGetCurrentDepth(tree);

      if( SCIPlpiIsPrimalInfeasible(lpi) )
      {
         SCIP_CALL( undoBdchgsDualfarkas(set, transprob, lp, currentdepth, curvarlbs, curvarubs, lbchginfoposs, \
               ubchginfoposs, oldlpbdchgs, relaxedlpbdchgs, valid, &resolve, proofcoefs, *prooflhs, proofactivity) );
      }
      else
      {
         assert(SCIPlpiIsDualFeasible(lpi) || SCIPlpiIsObjlimExc(lpi));
         SCIP_CALL( undoBdchgsDualsol(set, transprob, lp, currentdepth, curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs, \
               oldlpbdchgs, relaxedlpbdchgs, valid, &resolve, proofcoefs, *prooflhs, proofactivity) );
      }
   }

   /* check if we want to solve the LP */
   assert(SCIPprobAllColsInLP(transprob, set, lp));
   solvelp = (set->conf_maxlploops != 0 && set->conf_lpiterations != 0);

   if( (*valid) && resolve && solvelp )
   {
      SCIP_RETCODE retcode;
      SCIP_ROW** rows;
      int* sidechginds;
      SCIP_Real* sidechgoldlhss;
      SCIP_Real* sidechgoldrhss;
      SCIP_Real* sidechgnewlhss;
      SCIP_Real* sidechgnewrhss;
      SCIP_Real lpiinfinity;
      SCIP_Bool globalinfeasible;
      int maxlploops;
      int lpiterations;
      int sidechgssize;
      int nsidechgs;
      int nrows;
      int nloops;
      int r;

      /* get infinity value of LP solver */
      lpiinfinity = SCIPlpiInfinity(lpi);

      /* temporarily disable objective limit and install an iteration limit */
      maxlploops = (set->conf_maxlploops >= 0 ? set->conf_maxlploops : INT_MAX);
      lpiterations = (set->conf_lpiterations >= 0 ? set->conf_lpiterations : INT_MAX);
      SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_OBJLIM, lpiinfinity) );
      SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, lpiterations) );

      /* get LP rows */
      rows = SCIPlpGetRows(lp);
      nrows = SCIPlpGetNRows(lp);
      assert(nrows == 0 || rows != NULL);

      /* get temporary memory for remembering side changes on LPI rows */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechginds, nrows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgoldlhss, nrows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgoldrhss, nrows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgnewlhss, nrows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgnewrhss, nrows) );
      sidechgssize = nrows;
      nsidechgs = 0;

      /* remove all local rows by setting their sides to infinity;
       * finite sides are only changed to near infinity, such that the row's sense in the LP solver
       * is not affected (e.g. CPLEX cannot handle free rows)
       */
      for( r = 0; r < nrows; ++r )
      {
         assert(SCIProwGetLPPos(rows[r]) == r);

         if( SCIProwIsLocal(rows[r]) )
         {
            SCIPsetDebugMsg(set, " -> removing local row <%s> [%g,%g]\n",
               SCIProwGetName(rows[r]), SCIProwGetLhs(rows[r]), SCIProwGetRhs(rows[r]));
            SCIP_CALL( addSideRemoval(set, rows[r], lpiinfinity, &sidechginds, &sidechgoldlhss, &sidechgoldrhss,
                  &sidechgnewlhss, &sidechgnewrhss, &sidechgssize, &nsidechgs) );
         }
      }

      /* apply changes of local rows to the LP solver */
      if( nsidechgs > 0 )
      {
         SCIP_CALL( SCIPlpiChgSides(lpi, nsidechgs, sidechginds, sidechgnewlhss, sidechgnewrhss) );
      }

      /* undo as many additional bound changes as possible by resolving the LP */
      assert((*valid));
      assert(resolve);
      nloops = 0;
      globalinfeasible = FALSE;
      while( (*valid) && resolve && nloops < maxlploops )
      {
         int iter;

         assert(!globalinfeasible);

         nloops++;
         resolve = FALSE;

         SCIPsetDebugMsg(set, "infeasible LP conflict analysis loop %d (changed col bounds: %d)\n", nloops, relaxedlpbdchgs->nbdchgs);

         /* apply bound changes to the LP solver */
         assert(relaxedlpbdchgs->nbdchgs >= 0);
         if( relaxedlpbdchgs->nbdchgs > 0 )
         {
            SCIPsetDebugMsg(set, " -> applying %d bound changes to the LP solver\n", relaxedlpbdchgs->nbdchgs);
            SCIP_CALL( SCIPlpiChgBounds(lpi, relaxedlpbdchgs->nbdchgs, relaxedlpbdchgs->bdchginds, \
                  relaxedlpbdchgs->bdchglbs, relaxedlpbdchgs->bdchgubs) );

            /* reset conflict LP bound change data structure */
            lpbdchgsReset(relaxedlpbdchgs, ncols);
         }

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve LP */
         retcode = SCIPlpiSolveDual(lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         if( retcode == SCIP_LPERROR )
         {
            (*valid) = FALSE;
            break;
         }
         SCIP_CALL( retcode );

         /* count number of LP iterations */
         SCIP_CALL( SCIPlpiGetIterations(lpi, &iter) );
         (*iterations) += iter;
         stat->nconflictlps++;
         stat->nconflictlpiterations += iter;
         SCIPsetDebugMsg(set, " -> resolved LP in %d iterations (total: %" SCIP_LONGINT_FORMAT ") (infeasible:%u)\n",
            iter, stat->nconflictlpiterations, SCIPlpiIsPrimalInfeasible(lpi));

         /* evaluate result */
         if( SCIPlpiIsDualFeasible(lpi) || SCIPlpiIsObjlimExc(lpi) )
         {
            SCIP_Real objval;

            SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
            (*valid) = (objval >= lp->lpiobjlim && !SCIPlpDivingObjChanged(lp));
         }
         else
            (*valid) = SCIPlpiIsPrimalInfeasible(lpi);

         if( (*valid) )
         {
            int currentdepth;
            currentdepth = SCIPtreeGetCurrentDepth(tree);

            /* undo additional bound changes */
            if( SCIPlpiIsPrimalInfeasible(lpi) )
            {
               SCIP_AGGRROW* farkasrow;
               int* inds;
               int nnz;
               int v;

#ifndef NDEBUG
               SCIP_VAR** vars = SCIPprobGetVars(transprob);
#endif

               SCIP_CALL( SCIPaggrRowCreate(set->scip, &farkasrow) );

               /* the original LP exceeds the current cutoff bound, thus, we have not constructed the Farkas proof */
               SCIP_CALL( getFarkasProof(set, transprob, lp, lpi, farkasrow, proofactivity, curvarlbs, curvarubs, valid) );

               /* the constructed Farkas proof is not valid, we need to break here */
               if( !(*valid) )
               {
                  SCIPaggrRowFree(set->scip, &farkasrow);
                  break;
               }

               /* start dual ray analysis */
               if( set->conf_useinflp == 'd' || set->conf_useinflp == 'b' )
               {
                  /* change the conflict type */
                  SCIP_CONFTYPE oldconftype = conflict->conflictset->conflicttype;
                  conflict->conflictset->conflicttype = SCIP_CONFTYPE_INFEASLP;

                  /* start dual ray analysis */
                  SCIP_CALL( conflictAnalyzeDualProof(conflict, set, stat, blkmem, origprob, transprob, tree, reopt, lp, \
                        farkasrow, curvarlbs, curvarubs, FALSE, &globalinfeasible, dualraysuccess) );

                  conflict->conflictset->conflicttype = oldconftype;
               }

               if( globalinfeasible )
               {
                  SCIPaggrRowFree(set->scip, &farkasrow);
                  goto TERMINATE;
               }

               BMSclearMemoryArray(proofcoefs, SCIPprobGetNVars(transprob));
               (*prooflhs) = -SCIPaggrRowGetRhs(farkasrow);
               (*proofactivity) = -(*proofactivity);

               inds = SCIPaggrRowGetInds(farkasrow);
               nnz = SCIPaggrRowGetNNz(farkasrow);

               for( v = 0; v < nnz; v++ )
               {
                  int i = inds[v];

                  assert(SCIPvarGetProbindex(vars[i]) == inds[v]);

                  proofcoefs[i] = -SCIPaggrRowGetProbvarValue(farkasrow, i);
               }

               /* free aggregation rows */
               SCIPaggrRowFree(set->scip, &farkasrow);

               SCIP_CALL( undoBdchgsDualfarkas(set, transprob, lp, currentdepth, curvarlbs, curvarubs, lbchginfoposs, \
                     ubchginfoposs,  oldlpbdchgs, relaxedlpbdchgs, valid, &resolve, proofcoefs, (*prooflhs), proofactivity) );
            }
            else
            {
               SCIP_AGGRROW* proofrow;
               int* inds;
               int nnz;
               int v;

#ifndef NDEBUG
               SCIP_VAR** vars = SCIPprobGetVars(transprob);
#endif

               assert(SCIPlpiIsDualFeasible(lpi) || SCIPlpiIsObjlimExc(lpi));

               SCIP_CALL( SCIPaggrRowCreate(set->scip, &proofrow) );

               SCIP_CALL( getDualProof(set, transprob, lp, lpi, proofrow, proofactivity, curvarlbs, curvarubs, valid) );

               /* the constructed dual proof is not valid, we need to break here */
               if( !(*valid) )
               {
                  SCIPaggrRowFree(set->scip, &proofrow);
                  break;
               }
               /* in contrast to the infeasible case we don't want to analyze the (probably identical) proof again. */

               BMSclearMemoryArray(proofcoefs, SCIPprobGetNVars(transprob));
               (*prooflhs) = -SCIPaggrRowGetRhs(proofrow);
               (*proofactivity) = -(*proofactivity);

               inds = SCIPaggrRowGetInds(proofrow);
               nnz = SCIPaggrRowGetNNz(proofrow);

               for( v = 0; v < nnz; v++ )
               {
                  int i = inds[v];

                  assert(SCIPvarGetProbindex(vars[i]) == inds[v]);

                  proofcoefs[i] = -SCIPaggrRowGetProbvarValue(proofrow, i);
               }

               /* free aggregation rows */
               SCIPaggrRowFree(set->scip, &proofrow);

               SCIP_CALL( undoBdchgsDualsol(set, transprob, lp, currentdepth, curvarlbs, curvarubs, lbchginfoposs, \
                     ubchginfoposs, oldlpbdchgs, relaxedlpbdchgs, valid, &resolve, proofcoefs, *prooflhs, proofactivity) );
            }
         }
         assert(!resolve || (*valid));
         assert(!resolve || relaxedlpbdchgs->nbdchgs > 0);
         SCIPsetDebugMsg(set, " -> finished infeasible LP conflict analysis loop %d (iter: %d, nbdchgs: %d)\n",
            nloops, iter, relaxedlpbdchgs->nbdchgs);
      }

      SCIPsetDebugMsg(set, "finished undoing bound changes after %d loops (valid=%u, nbdchgs: %d)\n",
         nloops, (*valid), oldlpbdchgs->nbdchgs);

   TERMINATE:
      /* reset variables to local bounds */
      if( oldlpbdchgs->nbdchgs > 0 )
      {
         SCIP_CALL( SCIPlpiChgBounds(lpi, oldlpbdchgs->nbdchgs, oldlpbdchgs->bdchginds, oldlpbdchgs->bdchglbs, oldlpbdchgs->bdchgubs) );
      }

      /* reset changes of local rows */
      if( nsidechgs > 0 )
      {
         SCIP_CALL( SCIPlpiChgSides(lpi, nsidechgs, sidechginds, sidechgoldlhss, sidechgoldrhss) );
      }

      /* mark the LP unsolved */
      if( oldlpbdchgs->nbdchgs > 0 || nsidechgs > 0 )
      {
         /* The LPI data are out of sync with LP data. Thus, the LP should be marked
          * unsolved. However, for strong branching calls, the LP has to have status 'solved'; in
          * this case, marklpunsolved is FALSE and synchronization is performed later. */
         if ( marklpunsolved )
         {
            lp->solved = FALSE;
            lp->primalfeasible = FALSE;
            lp->primalchecked = FALSE;
            lp->dualfeasible = FALSE;
            lp->dualchecked = FALSE;
            lp->lpobjval = SCIP_INVALID;
            lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         }
      }

      /* reinstall old objective and iteration limits in LP solver */
      SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_OBJLIM, lp->lpiobjlim) );
      SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, lp->lpiitlim) );

      /* free temporary memory */
      SCIPsetFreeBufferArray(set, &sidechgnewrhss);
      SCIPsetFreeBufferArray(set, &sidechgnewlhss);
      SCIPsetFreeBufferArray(set, &sidechgoldrhss);
      SCIPsetFreeBufferArray(set, &sidechgoldlhss);
      SCIPsetFreeBufferArray(set, &sidechginds);
   }

   /* free temporary memory */
   lpbdchgsFree(&relaxedlpbdchgs, set);
   lpbdchgsFree(&oldlpbdchgs, set);

   return SCIP_OKAY;
}

/** actually performs analysis of infeasible LP */
static
SCIP_RETCODE conflictAnalyzeLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
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
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   SCIP_Bool*            dualraysuccess,     /**< pointer to store success result of dualray analysis */
   int*                  iterations,         /**< pointer to store the total number of LP iterations used */
   int*                  nconss,             /**< pointer to store the number of generated conflict constraints */
   int*                  nliterals,          /**< pointer to store the number of literals in generated conflict constraints */
   int*                  nreconvconss,       /**< pointer to store the number of generated reconvergence constraints */
   int*                  nreconvliterals,    /**< pointer to store the number of literals generated reconvergence constraints */
   SCIP_Bool             marklpunsolved      /**< whether LP should be marked unsolved after analysis (needed for strong branching) */
   )
{
   SCIP_VAR** vars;
   SCIP_AGGRROW* farkasrow;
   SCIP_LPI* lpi;
   SCIP_Bool valid;
   SCIP_Bool globalinfeasible;
   int* lbchginfoposs;
   int* ubchginfoposs;
   int nvars;
   int v;
   SCIP_Real* curvarlbs;
   SCIP_Real* curvarubs;
   SCIP_Real farkasactivity;

   assert(conflict != NULL);
   assert(conflict->nconflictsets == 0);
   assert(set != NULL);
   assert(SCIPprobAllColsInLP(transprob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(iterations != NULL);
   assert(nconss != NULL);
   assert(nliterals != NULL);
   assert(nreconvconss != NULL);
   assert(nreconvliterals != NULL);

   *iterations = 0;
   *nconss = 0;
   *nliterals = 0;
   *nreconvconss = 0;
   *nreconvliterals = 0;

   vars = transprob->vars;
   nvars = transprob->nvars;

   valid = TRUE;

   /* get LP solver interface */
   lpi = SCIPlpGetLPI(lp);
   assert(SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsDualFeasible(lpi));
   assert(SCIPlpiIsPrimalInfeasible(lpi) || !SCIPlpDivingObjChanged(lp));

   if( !SCIPlpiIsPrimalInfeasible(lpi) )
   {
      SCIP_Real objval;

      assert(!SCIPlpDivingObjChanged(lp));

      /* make sure, a dual feasible solution exists, that exceeds the objective limit;
       * With FASTMIP setting, CPLEX does not apply the final pivot to reach the dual solution exceeding the objective
       * limit. Therefore, we have to either turn off FASTMIP and resolve the problem or continue solving it without
       * objective limit for at least one iteration. It seems that the strategy to continue with FASTMIP for one
       * additional simplex iteration yields better results.
       */
      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      if( objval < lp->lpiobjlim )
      {
         SCIP_RETCODE retcode;

         /* temporarily disable objective limit and install an iteration limit */
         SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_OBJLIM, SCIPlpiInfinity(lpi)) );
         SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, 1) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve LP */
         retcode = SCIPlpiSolveDual(lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         valid = (retcode != SCIP_LPERROR);
         if( valid )
         {
            int iter;

            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lpi, &iter) );
            (*iterations) += iter;
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            SCIPsetDebugMsg(set, " -> resolved objlim exceeding LP in %d iterations (total: %" SCIP_LONGINT_FORMAT ") (infeasible:%u, objlim: %u, optimal:%u)\n",
               iter, stat->nconflictlpiterations, SCIPlpiIsPrimalInfeasible(lpi), SCIPlpiIsObjlimExc(lpi), SCIPlpiIsOptimal(lpi));
            valid = (SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsDualFeasible(lpi));
         }

         /* reinstall old objective and iteration limits in LP solver */
         SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_OBJLIM, lp->lpiobjlim) );
         SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, lp->lpiitlim) );

         /* abort, if the LP produced an error */
         if( !valid )
            return SCIP_OKAY;
      }
   }
   assert(SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsDualFeasible(lpi));

   if( !SCIPlpiIsPrimalInfeasible(lpi) )
   {
      SCIP_Real objval;

      assert(!SCIPlpDivingObjChanged(lp));

      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      if( objval < lp->lpiobjlim )
      {
         SCIPsetDebugMsg(set, " -> LP does not exceed the cutoff bound: obj=%g, cutoff=%g\n", objval, lp->lpiobjlim);
         return SCIP_OKAY;
      }
      else
      {
         SCIPsetDebugMsg(set, " -> LP exceeds the cutoff bound: obj=%g, cutoff=%g\n", objval, lp->lpiobjlim);
      }
   }

   assert(valid);

   SCIP_CALL( SCIPaggrRowCreate(set->scip, &farkasrow) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lbchginfoposs, transprob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ubchginfoposs, transprob->nvars) );

   farkasactivity = 0.0;

   /* get temporary memory for remembering variables' current bounds and corresponding bound change information
    * positions in variable's bound change information arrays
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarlbs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarubs, nvars) );


   /* get current bounds and current positions in lb/ubchginfos arrays of variables */
   valid = TRUE;
   for( v = 0; v < nvars && valid; ++v )
   {
      SCIP_VAR* var;

      var = vars[v];

      curvarlbs[v] = SCIPvarGetLbLP(var, set);
      curvarubs[v] = SCIPvarGetUbLP(var, set);
      lbchginfoposs[v] = var->nlbchginfos-1;
      ubchginfoposs[v] = var->nubchginfos-1;
      assert(diving || SCIPsetIsEQ(set, curvarlbs[v], SCIPvarGetLbLocal(var)));
      assert(diving || SCIPsetIsEQ(set, curvarubs[v], SCIPvarGetUbLocal(var)));

      /* check, if last bound changes were due to strong branching or diving */
      if( diving )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
         if( SCIPsetIsGT(set, curvarlbs[v], lb) )
            lbchginfoposs[v] = var->nlbchginfos;
         else if( SCIPsetIsLT(set, curvarlbs[v], lb) )
         {
            /* the bound in the diving LP was relaxed -> the LP is not a subproblem of the current node -> abort! */
            /**@todo we could still analyze such a conflict, but we would have to take care with our data structures */
            valid = FALSE;
         }
         if( SCIPsetIsLT(set, curvarubs[v], ub) )
            ubchginfoposs[v] = var->nubchginfos;
         else if( SCIPsetIsGT(set, curvarubs[v], ub) )
         {
            /* the bound in the diving LP was relaxed -> the LP is not a subproblem of the current node -> abort! */
            /**@todo we could still analyze such a conflict, but we would have to take care with our data structures */
            valid = FALSE;
         }
      }
   }

   if( !valid )
      goto TERMINATE;

   /* the LP is prooven to be infeasible */
   if( SCIPlpiIsPrimalInfeasible(lpi) )
   {
      SCIP_CALL( getFarkasProof(set, transprob, lp, lpi, farkasrow, &farkasactivity, curvarlbs, curvarubs, &valid) );
   }
   /* the LP is dual feasible and/or exceeds the current incumbant solution */
   else
   {
      assert(SCIPlpiIsDualFeasible(lpi) || SCIPlpiIsObjlimExc(lpi));
      SCIP_CALL( getDualProof(set, transprob, lp, lpi, farkasrow, &farkasactivity, curvarlbs, curvarubs, &valid) );
   }

   if( !valid )
      goto TERMINATE;

   globalinfeasible = FALSE;

   /* start dual proof analysis */
   if( ((set->conf_useinflp == 'b' || set->conf_useinflp == 'd') && conflict->conflictset->conflicttype == SCIP_CONFTYPE_INFEASLP)
      || ((set->conf_useboundlp == 'b' || set->conf_useboundlp == 'd') && conflict->conflictset->conflicttype == SCIP_CONFTYPE_BNDEXCEEDING) )
   {
      /* start dual ray analysis */
      SCIP_CALL( conflictAnalyzeDualProof(conflict, set, stat, blkmem, origprob, transprob, tree, reopt, lp, farkasrow, \
            curvarlbs, curvarubs, TRUE, &globalinfeasible, dualraysuccess) );
   }

   assert(valid);

   if( !globalinfeasible && (((set->conf_useinflp == 'b' || set->conf_useinflp == 'c') && conflict->conflictset->conflicttype == SCIP_CONFTYPE_INFEASLP)
      || ((set->conf_useboundlp == 'b' || set->conf_useboundlp == 'c') && conflict->conflictset->conflicttype == SCIP_CONFTYPE_BNDEXCEEDING)) )
   {
      SCIP_Real* farkascoefs;
      SCIP_Real farkaslhs;
      int* inds;
      int nnz;

#ifdef SCIP_DEBUG
      {
         SCIP_Real objlim;
         SCIPsetDebugMsg(set, "analyzing conflict on infeasible LP (infeasible: %u, objlimexc: %u, optimal:%u) in depth %d (diving: %u)\n",
               SCIPlpiIsPrimalInfeasible(lpi), SCIPlpiIsObjlimExc(lpi), SCIPlpiIsOptimal(lpi), SCIPtreeGetCurrentDepth(tree), diving);

         SCIP_CALL( SCIPlpiGetRealpar(lpi, SCIP_LPPAR_OBJLIM, &objlim) );
         SCIPsetDebugMsg(set, " -> objective limit in LP solver: %g (in LP: %g)\n", objlim, lp->lpiobjlim);
      }
#endif

      SCIP_CALL( SCIPsetAllocBufferArray(set, &farkascoefs, SCIPprobGetNVars(transprob)) );
      BMSclearMemoryArray(farkascoefs, SCIPprobGetNVars(transprob));

      farkaslhs = -SCIPaggrRowGetRhs(farkasrow);
      farkasactivity = -farkasactivity;

      inds = SCIPaggrRowGetInds(farkasrow);
      nnz = SCIPaggrRowGetNNz(farkasrow);

      for( v = 0; v < nnz; v++ )
      {
         int i = inds[v];

         assert(SCIPvarGetProbindex(vars[i]) == inds[v]);

         farkascoefs[i] = -SCIPaggrRowGetProbvarValue(farkasrow, i);
      }

      SCIP_CALL( runBoundHeuristic(conflict, set, stat, origprob, transprob, tree, reopt, lp, lpi, blkmem, farkascoefs,
            &farkaslhs, &farkasactivity, curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs, iterations, marklpunsolved,
            dualraysuccess, &valid) );

      SCIPsetFreeBufferArray(set, &farkascoefs);

      if( !valid )
         goto FLUSHPROOFSETS;

      /* analyze the conflict starting with remaining bound changes */
      SCIP_CALL( conflictAnalyzeRemainingBdchgs(conflict, blkmem, set, stat, transprob, tree, diving, \
            lbchginfoposs, ubchginfoposs, nconss, nliterals, nreconvconss, nreconvliterals) );

      /* flush conflict set storage */
      SCIP_CALL( SCIPconflictFlushConss(conflict, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, \
            eventqueue, cliquetable) );
   }

  FLUSHPROOFSETS:
   /* flush proof set */
   if( proofsetGetNVars(conflict->proofset) > 0 || conflict->nproofsets > 0 )
   {
      SCIP_CALL( conflictFlushProofset(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, \
            branchcand, eventqueue, cliquetable) );
   }

  TERMINATE:
   SCIPsetFreeBufferArray(set, &curvarubs);
   SCIPsetFreeBufferArray(set, &curvarlbs);
   SCIPsetFreeBufferArray(set, &ubchginfoposs);
   SCIPsetFreeBufferArray(set, &lbchginfoposs);
   SCIPaggrRowFree(set->scip, &farkasrow);

   return SCIP_OKAY;
}

/** analyzes an infeasible LP to find out the bound changes on variables that were responsible for the infeasibility;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for infeasible LP conflict analysis
 */
static
SCIP_RETCODE conflictAnalyzeInfeasibleLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
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
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_Bool dualraysuccess = FALSE;
   SCIP_Longint olddualproofsuccess;
   int iterations;
   int nconss;
   int nliterals;
   int nreconvconss;
   int nreconvliterals;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(SCIPprobAllColsInLP(transprob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */

   assert(success == NULL || *success == FALSE);

   /* check, if infeasible LP conflict analysis is enabled */
   if( !set->conf_enable || set->conf_useinflp == 'o' )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "analyzing conflict on infeasible LP in depth %d (solstat: %d, objchanged: %u)\n",
      SCIPtreeGetCurrentDepth(tree), SCIPlpGetSolstat(lp), SCIPlpDivingObjChanged(lp));

   /* start timing */
   SCIPclockStart(conflict->inflpanalyzetime, set);
   conflict->ninflpcalls++;

   conflict->conflictset->conflicttype = SCIP_CONFTYPE_INFEASLP;

   olddualproofsuccess = conflict->ndualrayinfsuccess;

   /* perform conflict analysis */
   SCIP_CALL( conflictAnalyzeLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, \
         cliquetable, SCIPlpDiving(lp), &dualraysuccess, &iterations, &nconss, &nliterals, &nreconvconss, &nreconvliterals, TRUE) );
   conflict->ninflpsuccess += ((nconss > 0 || conflict->ndualrayinfsuccess > olddualproofsuccess) ? 1 : 0);
   conflict->ninflpiterations += iterations;
   conflict->ninflpconfconss += nconss;
   conflict->ninflpconfliterals += nliterals;
   conflict->ninflpreconvconss += nreconvconss;
   conflict->ninflpreconvliterals += nreconvliterals;
   if( success != NULL )
      *success = (nconss > 0 || conflict->ndualrayinfsuccess > olddualproofsuccess);

   /* stop timing */
   SCIPclockStop(conflict->inflpanalyzetime, set);

   return SCIP_OKAY;
}

/** analyzes a bound exceeding LP to find out the bound changes on variables that were responsible for exceeding the
 *  primal bound;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for bound exceeding LP conflict analysis
 */
static
SCIP_RETCODE conflictAnalyzeBoundexceedingLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
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
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_Bool dualraysuccess;
   SCIP_Longint oldnsuccess;
   int iterations;
   int nconss;
   int nliterals;
   int nreconvconss;
   int nreconvliterals;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(!SCIPlpDivingObjChanged(lp));
   assert(SCIPprobAllColsInLP(transprob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */

   assert(success == NULL || *success == FALSE);

   /* check, if bound exceeding LP conflict analysis is enabled */
   if( !set->conf_enable || set->conf_useboundlp == 'o')
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "analyzing conflict on bound exceeding LP in depth %d (solstat: %d)\n",
      SCIPtreeGetCurrentDepth(tree), SCIPlpGetSolstat(lp));

   /* start timing */
   SCIPclockStart(conflict->boundlpanalyzetime, set);
   conflict->nboundlpcalls++;

   /* mark the conflict to depend on the cutoff bound */
   conflict->conflictset->conflicttype = SCIP_CONFTYPE_BNDEXCEEDING;
   conflict->conflictset->usescutoffbound = TRUE;

   oldnsuccess = conflict->ndualraybndsuccess + conflict->ndualrayinfsuccess;

   /* perform conflict analysis */
   SCIP_CALL( conflictAnalyzeLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, \
         cliquetable, SCIPlpDiving(lp), &dualraysuccess, &iterations, &nconss, &nliterals, &nreconvconss, &nreconvliterals, TRUE) );
   conflict->nboundlpsuccess += ((nconss > 0 || conflict->ndualraybndsuccess + conflict->ndualrayinfsuccess > oldnsuccess) ? 1 : 0);
   conflict->nboundlpiterations += iterations;
   conflict->nboundlpconfconss += nconss;
   conflict->nboundlpconfliterals += nliterals;
   conflict->nboundlpreconvconss += nreconvconss;
   conflict->nboundlpreconvliterals += nreconvliterals;
   if( success != NULL )
      *success = (nconss > 0 || conflict->ndualraybndsuccess + conflict->ndualrayinfsuccess > oldnsuccess);

   /* stop timing */
   SCIPclockStop(conflict->boundlpanalyzetime, set);

   return SCIP_OKAY;
}

/** analyzes an infeasible or bound exceeding LP to find out the bound changes on variables that were responsible for the
 *  infeasibility or for exceeding the primal bound;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for infeasible or bound exceeding LP conflict analysis;
 *  may only be called if SCIPprobAllColsInLP()
 */
SCIP_RETCODE SCIPconflictAnalyzeLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
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
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_LPSOLVALS storedsolvals;
   SCIP_COLSOLVALS* storedcolsolvals;
   SCIP_ROWSOLVALS* storedrowsolvals;
   int c;
   int r;

   if( success != NULL )
      *success = FALSE;

   /* check if the conflict analysis is applicable */
   if( !set->conf_enable || (set->conf_useinflp == 'o' && set->conf_useboundlp == 'o') )
      return SCIP_OKAY;

   /* in rare cases, it might happen that the solution stati of the LP and the LPI are out of sync; in particular this
    * happens when a new incumbent which cuts off the current node is found during the LP solving loop; in this case the
    * LP has status objlimit, but if diving has been used, the LPI only has the basis information, but is not solved
    *
    * @todo: alternatively, solve the LPI
    */
   if( !SCIPlpiWasSolved(SCIPlpGetLPI(lp)) )
      return SCIP_OKAY;

   /* LP conflict analysis is only valid, if all variables are known */
   assert( SCIPprobAllColsInLP(transprob, set, lp) );
   assert( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT
      || (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && set->lp_disablecutoff == 1) );

   /* save status */
   storedsolvals.lpsolstat = lp->lpsolstat;
   storedsolvals.lpobjval = lp->lpobjval;
   storedsolvals.primalfeasible = lp->primalfeasible;
   storedsolvals.primalchecked = lp->primalchecked;
   storedsolvals.dualfeasible = lp->dualfeasible;
   storedsolvals.dualchecked = lp->dualchecked;
   storedsolvals.solisbasic = lp->solisbasic;
   storedsolvals.lpissolved = lp->solved;

   /* store solution values */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &storedcolsolvals, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &storedrowsolvals, lp->nrows) );
   for (c = 0; c < lp->ncols; ++c)
   {
      SCIP_COL* col;

      col = lp->cols[c];
      assert( col != NULL );

      storedcolsolvals[c].primsol = col->primsol;
      storedcolsolvals[c].redcost = col->redcost;
      storedcolsolvals[c].basisstatus = col->basisstatus; /*lint !e641 !e732*/
   }
   for (r = 0; r < lp->nrows; ++r)
   {
      SCIP_ROW* row;

      row = lp->rows[r];
      assert( row != NULL );

      if ( lp->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE )
         storedrowsolvals[r].dualsol = row->dualfarkas;
      else
      {
         assert( lp->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT ||
            (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && set->lp_disablecutoff == 1) );
         storedrowsolvals[r].dualsol = row->dualsol;
      }
      storedrowsolvals[r].activity = row->activity;
      storedrowsolvals[r].basisstatus = row->basisstatus; /*lint !e641 !e732*/
   }

   /* check, if the LP was infeasible or bound exceeding */
   if( SCIPlpiIsPrimalInfeasible(SCIPlpGetLPI(lp)) )
   {
      SCIP_CALL( conflictAnalyzeInfeasibleLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, \
            reopt, lp, branchcand, eventqueue, cliquetable, success) );
   }
   else
   {
      SCIP_CALL( conflictAnalyzeBoundexceedingLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, \
            reopt, lp, branchcand, eventqueue, cliquetable, success) );
   }

   /* possibly restore solution values */
   if( lp->flushed && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_NOTSOLVED )
   {
      /* restore status */
      lp->lpsolstat = storedsolvals.lpsolstat;
      lp->lpobjval = storedsolvals.lpobjval;
      lp->primalfeasible = storedsolvals.primalfeasible;
      lp->primalchecked = storedsolvals.primalchecked;
      lp->dualfeasible = storedsolvals.dualfeasible;
      lp->dualchecked = storedsolvals.dualchecked;
      lp->solisbasic = storedsolvals.solisbasic;
      lp->solved = storedsolvals.lpissolved;

      for (c = 0; c < lp->ncols; ++c)
      {
         SCIP_COL* col;

         col = lp->cols[c];
         assert( col != NULL );
         col->primsol = storedcolsolvals[c].primsol;
         col->redcost = storedcolsolvals[c].redcost;
         col->basisstatus = storedcolsolvals[c].basisstatus; /*lint !e641 !e732*/
      }
      for (r = 0; r < lp->nrows; ++r)
      {
         SCIP_ROW* row;

         row = lp->rows[r];
         assert( row != NULL );

         if ( lp->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE )
            row->dualfarkas = storedrowsolvals[r].dualsol;
         else
         {
            assert( lp->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT );
            row->dualsol = storedrowsolvals[r].dualsol;
         }
         row->activity = storedrowsolvals[r].activity;
         row->basisstatus = storedrowsolvals[r].basisstatus; /*lint !e641 !e732*/
      }
   }
   SCIPsetFreeBufferArray(set, &storedrowsolvals);
   SCIPsetFreeBufferArray(set, &storedcolsolvals);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing infeasible LP conflicts */
SCIP_Real SCIPconflictGetInfeasibleLPTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->inflpanalyzetime);
}

/** gets number of calls to infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpcalls;
}

/** gets number of calls to infeasible LP conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNInfeasibleLPSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpsuccess;
}

/** gets number of conflict constraints detected in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpconfconss;
}

/** gets total number of literals in conflict constraints created in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpconfliterals;
}

/** gets number of reconvergence constraints detected in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpreconvconss;
}

/** gets total number of literals in reconvergence constraints created in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpreconvliterals;
}

/** gets number of LP iterations in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpiterations;
}

/** gets time in seconds used for analyzing bound exceeding LP conflicts */
SCIP_Real SCIPconflictGetBoundexceedingLPTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->boundlpanalyzetime);
}

/** gets number of calls to bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpcalls;
}

/** gets number of calls to bound exceeding LP conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNBoundexceedingLPSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpsuccess;
}

/** gets number of conflict constraints detected in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpconfconss;
}

/** gets total number of literals in conflict constraints created in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpconfliterals;
}

/** gets number of reconvergence constraints detected in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpreconvconss;
}

/** gets total number of literals in reconvergence constraints created in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpreconvliterals;
}

/** gets number of LP iterations in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpiterations;
}




/*
 * infeasible strong branching conflict analysis
 */

/** analyses infeasible strong branching sub problems for conflicts */
SCIP_RETCODE SCIPconflictAnalyzeStrongbranch(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_COL*             col,                /**< LP column with at least one infeasible strong branching subproblem */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict          /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   )
{
   int* cstat;
   int* rstat;
   SCIP_RETCODE retcode;
   SCIP_Bool resolve;
   SCIP_Real oldlb;
   SCIP_Real oldub;
   SCIP_Real newlb;
   SCIP_Real newub;
   SCIP_Bool dualraysuccess;
   int iter;
   int nconss;
   int nliterals;
   int nreconvconss;
   int nreconvliterals;

   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPprobAllColsInLP(transprob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */
   assert(col != NULL);
   assert((col->sbdownvalid && SCIPsetIsGE(set, col->sbdown, lp->cutoffbound)
         && SCIPsetFeasCeil(set, col->primsol-1.0) >= col->lb - 0.5)
      || (col->sbupvalid && SCIPsetIsGE(set, col->sbup, lp->cutoffbound)
         && SCIPsetFeasFloor(set, col->primsol+1.0) <= col->ub + 0.5));
   assert(SCIPtreeGetCurrentDepth(tree) > 0);

   if( downconflict != NULL )
      *downconflict = FALSE;
   if( upconflict != NULL )
      *upconflict = FALSE;

   /* check, if infeasible LP conflict analysis is enabled */
   if( !set->conf_enable || !set->conf_usesb )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   /* inform the LPI that strong branch is (temporarily) finished */
   SCIP_CALL( SCIPlpiEndStrongbranch(lp->lpi) );

   /* start timing */
   SCIPclockStart(conflict->sbanalyzetime, set);

   /* get temporary memory for storing current LP basis */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cstat, lp->nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rstat, lp->nlpirows) );

   /* get current LP basis */
   SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );

   /* remember old bounds */
   oldlb = col->lb;
   oldub = col->ub;

   resolve = FALSE;

   /* is down branch infeasible? */
   if( col->sbdownvalid && SCIPsetIsGE(set, col->sbdown, lp->cutoffbound) )
   {
      newub = SCIPsetFeasCeil(set, col->primsol-1.0);
      if( newub >= col->lb - 0.5 )
      {
         SCIPsetDebugMsg(set, "analyzing conflict on infeasible downwards strongbranch for variable <%s>[%g,%g] in depth %d\n",
            SCIPvarGetName(SCIPcolGetVar(col)), SCIPvarGetLbLocal(SCIPcolGetVar(col)), SCIPvarGetUbLocal(SCIPcolGetVar(col)),
            SCIPtreeGetCurrentDepth(tree));

         conflict->conflictset->conflicttype = SCIP_CONFTYPE_INFEASLP;
         conflict->nsbcalls++;

         /* change the upper bound */
         col->ub = newub;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve the LP */
         retcode = SCIPlpiSolveDual(lp->lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         if( retcode != SCIP_LPERROR )
         {
            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lp->lpi, &iter) );
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            conflict->nsbiterations += iter;
            SCIPsetDebugMsg(set, " -> resolved downwards strong branching LP in %d iterations\n", iter);

            /* perform conflict analysis on infeasible LP; last parameter guarantees status 'solved' on return */
            SCIP_CALL( conflictAnalyzeLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, reopt, \
                  lp, branchcand, eventqueue, cliquetable, TRUE, &dualraysuccess, &iter, &nconss, &nliterals, \
                  &nreconvconss, &nreconvliterals, FALSE) );
            conflict->nsbsuccess += ((nconss > 0 || dualraysuccess) ? 1 : 0);
            conflict->nsbiterations += iter;
            conflict->nsbconfconss += nconss;
            conflict->nsbconfliterals += nliterals;
            conflict->nsbreconvconss += nreconvconss;
            conflict->nsbreconvliterals += nreconvliterals;
            if( downconflict != NULL )
               *downconflict = (nconss > 0);
         }

         /* reset the upper bound */
         col->ub = oldub;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* reset LP basis */
         SCIP_CALL( SCIPlpiSetBase(lp->lpi, cstat, rstat) );

         /* mark the LP to be resolved at the end */
         resolve = TRUE;
      }
   }

   /* is up branch infeasible? */
   if( col->sbupvalid && SCIPsetIsGE(set, col->sbup, lp->cutoffbound) )
   {
      newlb = SCIPsetFeasFloor(set, col->primsol+1.0);
      if( newlb <= col->ub + 0.5 )
      {
         SCIPsetDebugMsg(set, "analyzing conflict on infeasible upwards strongbranch for variable <%s>[%g,%g] in depth %d\n",
            SCIPvarGetName(SCIPcolGetVar(col)), SCIPvarGetLbLocal(SCIPcolGetVar(col)), SCIPvarGetUbLocal(SCIPcolGetVar(col)),
            SCIPtreeGetCurrentDepth(tree));

         conflict->conflictset->conflicttype = SCIP_CONFTYPE_INFEASLP;
         conflict->nsbcalls++;

         /* change the lower bound */
         col->lb = newlb;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve the LP */
         retcode = SCIPlpiSolveDual(lp->lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         if( retcode != SCIP_LPERROR )
         {
            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lp->lpi, &iter) );
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            conflict->nsbiterations += iter;
            SCIPsetDebugMsg(set, " -> resolved upwards strong branching LP in %d iterations\n", iter);

            /* perform conflict analysis on infeasible LP; last parameter guarantees status 'solved' on return */
            SCIP_CALL( conflictAnalyzeLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, reopt, \
                  lp, branchcand, eventqueue, cliquetable, TRUE, &dualraysuccess, &iter, &nconss, &nliterals, \
                  &nreconvconss, &nreconvliterals, FALSE) );
            conflict->nsbsuccess += ((nconss > 0 || dualraysuccess) ? 1 : 0);
            conflict->nsbiterations += iter;
            conflict->nsbconfconss += nconss;
            conflict->nsbconfliterals += nliterals;
            conflict->nsbreconvconss += nreconvconss;
            conflict->nsbreconvliterals += nreconvliterals;
            if( upconflict != NULL )
               *upconflict = (nconss > 0);
         }

         /* reset the lower bound */
         col->lb = oldlb;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* reset LP basis */
         SCIP_CALL( SCIPlpiSetBase(lp->lpi, cstat, rstat) );

         /* mark the LP to be resolved at the end */
         resolve = TRUE;
      }
   }

   /* free temporary memory for storing current LP basis */
   SCIPsetFreeBufferArray(set, &rstat);
   SCIPsetFreeBufferArray(set, &cstat);

   assert(lp->flushed);

   /* resolve LP if something has changed in order to synchronize LPI and LP */
   if ( resolve )
   {
      /* start LP timer */
      SCIPclockStart(stat->conflictlptime, set);

      /* resolve the LP */
      SCIP_CALL( SCIPlpiSolveDual(lp->lpi) );

      /* stop LP timer */
      SCIPclockStop(stat->conflictlptime, set);
   }

   /* stop timing */
   SCIPclockStop(conflict->sbanalyzetime, set);

   /* inform the LPI that strong branch starts (again) */
   SCIP_CALL( SCIPlpiStartStrongbranch(lp->lpi) );

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing infeasible strong branching conflicts */
SCIP_Real SCIPconflictGetStrongbranchTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->sbanalyzetime);
}

/** gets number of successful calls to infeasible dualray analysis */
SCIP_Longint SCIPconflictGetNDualrayInfSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ndualrayinfsuccess;
}

/** gets number of globally valid dualray constraints */
SCIP_Longint SCIPconflictGetNDualrayInfGlobal(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ndualrayinfglobal;
}

/** gets average length of infeasible dualrays */
SCIP_Longint SCIPconflictGetNDualrayInfNonzeros(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->dualrayinfnnonzeros;
}

/** gets number of successfully analyzed dual proofs of boundexceeding LPs */
SCIP_Longint SCIPconflictGetNDualrayBndSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ndualraybndsuccess;
}

/** gets number of globally applied dual proofs of boundexceeding LPs */
SCIP_Longint SCIPconflictGetNDualrayBndGlobal(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ndualraybndglobal;
}

/** gets average length of dual proofs of boundexceeding LPs */
SCIP_Longint SCIPconflictGetNDualrayBndNonzeros(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->dualraybndnnonzeros;
}

/** gets number of calls to infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbcalls;
}

/** gets number of calls to infeasible strong branching conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNStrongbranchSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbsuccess;
}

/** gets number of conflict constraints detected in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbconfconss;
}

/** gets total number of literals in conflict constraints created in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbconfliterals;
}

/** gets number of reconvergence constraints detected in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbreconvconss;
}

/** gets total number of literals in reconvergence constraints created in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbreconvliterals;
}

/** gets number of LP iterations in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbiterations;
}




/*
 * pseudo solution conflict analysis
 */

/** analyzes a pseudo solution with objective value exceeding the current cutoff to find out the bound changes on
 *  variables that were responsible for the objective value degradation;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for pseudo solution conflict analysis
 */
SCIP_RETCODE SCIPconflictAnalyzePseudo(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
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
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real* curvarlbs;
   SCIP_Real* curvarubs;
   int* lbchginfoposs;
   int* ubchginfoposs;
   SCIP_Real* pseudocoefs;
   SCIP_Real pseudolhs;
   SCIP_Real pseudoact;
   int nvars;
   int v;

   assert(conflict != NULL);
   assert(conflict->nconflictsets == 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(lp != NULL);
   assert(!SCIPsetIsInfinity(set, -SCIPlpGetPseudoObjval(lp, set, transprob)));
   assert(!SCIPsetIsInfinity(set, lp->cutoffbound));

   if( success != NULL )
      *success = FALSE;

   /* check, if pseudo solution conflict analysis is enabled */
   if( !set->conf_enable || !set->conf_usepseudo )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "analyzing pseudo solution (obj: %g) that exceeds objective limit (%g)\n",
      SCIPlpGetPseudoObjval(lp, set, transprob), lp->cutoffbound);

   conflict->conflictset->conflicttype = SCIP_CONFTYPE_BNDEXCEEDING;
   conflict->conflictset->usescutoffbound = TRUE;

   /* start timing */
   SCIPclockStart(conflict->pseudoanalyzetime, set);
   conflict->npseudocalls++;

   vars = transprob->vars;
   nvars = transprob->nvars;
   assert(nvars == 0 || vars != NULL);

   /* The current primal bound c* gives an upper bound for the current pseudo objective value:
    *   min{c^T x | lb <= x <= ub} <= c*.
    * We have to transform this row into a >= inequality in order to use methods above:
    *                          -c* <= max{-c^T x | lb <= x <= ub}.
    * In the local subproblem, this row is violated. We want to undo bound changes while still keeping the
    * row violated.
    */

   /* get temporary memory for remembering variables' current bounds and corresponding bound change information
    * positions in variable's bound change information arrays
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarlbs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarubs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lbchginfoposs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ubchginfoposs, nvars) );

   /* get temporary memory for infeasibility proof coefficients */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &pseudocoefs, nvars) );

   /* use a slightly tighter cutoff bound, because solutions with equal objective value should also be declared
    * infeasible
    */
   pseudolhs = -(lp->cutoffbound - SCIPsetSumepsilon(set));

   /* store the objective values as infeasibility proof coefficients, and recalculate the pseudo activity */
   pseudoact = 0.0;
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      pseudocoefs[v] = -SCIPvarGetObj(var);
      curvarlbs[v] = SCIPvarGetLbLocal(var);
      curvarubs[v] = SCIPvarGetUbLocal(var);
      if( pseudocoefs[v] > 0.0 )
         pseudoact += pseudocoefs[v] * curvarubs[v];
      else
         pseudoact += pseudocoefs[v] * curvarlbs[v];
      lbchginfoposs[v] = var->nlbchginfos-1;
      ubchginfoposs[v] = var->nubchginfos-1;
   }
   assert(SCIPsetIsFeasEQ(set, pseudoact, -SCIPlpGetPseudoObjval(lp, set, transprob)));
   SCIPsetDebugMsg(set, "  -> recalculated pseudo infeasibility proof:  %g <= %g\n", pseudolhs, pseudoact);

   /* check, if the pseudo row is still violated (after recalculation of pseudo activity) */
   if( SCIPsetIsFeasGT(set, pseudolhs, pseudoact) )
   {
      int nconss;
      int nliterals;
      int nreconvconss;
      int nreconvliterals;

      /* undo bound changes without destroying the infeasibility proof */
      SCIP_CALL( undoBdchgsProof(set, transprob, SCIPtreeGetCurrentDepth(tree), pseudocoefs, pseudolhs, &pseudoact,
            curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs, NULL, NULL, NULL, lp->lpi) );

      /* analyze conflict on remaining bound changes */
      SCIP_CALL( conflictAnalyzeRemainingBdchgs(conflict, blkmem, set, stat, transprob, tree, FALSE, \
            lbchginfoposs, ubchginfoposs, &nconss, &nliterals, &nreconvconss, &nreconvliterals) );
      conflict->npseudosuccess += (nconss > 0 ? 1 : 0);
      conflict->npseudoconfconss += nconss;
      conflict->npseudoconfliterals += nliterals;
      conflict->npseudoreconvconss += nreconvconss;
      conflict->npseudoreconvliterals += nreconvliterals;
      if( success != NULL )
         *success = (nconss > 0);
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &pseudocoefs);
   SCIPsetFreeBufferArray(set, &ubchginfoposs);
   SCIPsetFreeBufferArray(set, &lbchginfoposs);
   SCIPsetFreeBufferArray(set, &curvarubs);
   SCIPsetFreeBufferArray(set, &curvarlbs);

   /* flush conflict set storage */
   SCIP_CALL( SCIPconflictFlushConss(conflict, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cliquetable) );

   /* stop timing */
   SCIPclockStop(conflict->pseudoanalyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing pseudo solution conflicts */
SCIP_Real SCIPconflictGetPseudoTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->pseudoanalyzetime);
}

/** gets number of calls to pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudocalls;
}

/** gets number of calls to pseudo solution conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNPseudoSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudosuccess;
}

/** gets number of conflict constraints detected in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoconfconss;
}

/** gets total number of literals in conflict constraints created in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoconfliterals;
}

/** gets number of reconvergence constraints detected in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoreconvconss;
}

/** gets total number of literals in reconvergence constraints created in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoreconvliterals;
}


/** enables or disables all clocks of \p conflict, depending on the value of the flag */
void SCIPconflictEnableOrDisableClocks(
   SCIP_CONFLICT*        conflict,           /**< the conflict analysis data for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the conflict analysis data be enabled? */
   )
{
   assert(conflict != NULL);

   SCIPclockEnableOrDisable(conflict->boundlpanalyzetime, enable);
   SCIPclockEnableOrDisable(conflict->dIBclock, enable);
   SCIPclockEnableOrDisable(conflict->inflpanalyzetime, enable);
   SCIPclockEnableOrDisable(conflict->propanalyzetime, enable);
   SCIPclockEnableOrDisable(conflict->pseudoanalyzetime, enable);
   SCIPclockEnableOrDisable(conflict->sbanalyzetime, enable);
}

