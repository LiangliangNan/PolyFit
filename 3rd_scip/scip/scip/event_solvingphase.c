/* * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * * * * * * */
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

/**@file   event_solvingphase.c
 * @brief  event handler for solving phase dependent parameter adjustment
 * @author Gregor Hendel
 *
 * this event handler provides methods to support parameter adjustment at every new of the three solving phases:
 *   - Feasibility phase - before the first solution is found
 *   - Improvement phase - after the first solution was found until an optimal solution is found or believed to be found
 *   - Proof phase - the remaining time of the solution process after an optimal or believed-to-be optimal incumbent has been found.
 *
 * Of course, this event handler cannot detect by itself whether a given incumbent is optimal prior to termination of the
 * solution process. It rather uses heuristic transitions based on properties of the search tree in order to
 * determine the appropriate stage. Settings files can be passed to this event handler for each of the three phases.
 *
 * This approach of phase-based parameter adjustment was first presented in
 *
 * Gregor Hendel
 * Empirical Analysis of Solving Phases in Mixed-Integer Programming
 * Master thesis, Technical University Berlin (2014)
 *
 * with the main results also available from
 *
 * Gregor Hendel
 * Exploiting solving phases in mixed-integer programs (2015)
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_solvingphase.h"
#include "string.h"
#include "scip/pub_misc.h"

#define EVENTHDLR_NAME         "solvingphase"
#define EVENTHDLR_DESC         "event handler to adjust settings depending on current stage"

#define EVENTHDLR_EVENT SCIP_EVENTTYPE_BESTSOLFOUND | SCIP_EVENTTYPE_NODEBRANCHED | SCIP_EVENTTYPE_NODEFOCUSED /**< the actual event to be caught */
#define TRANSITIONMETHODS          "elor" /**< which heuristic transition method: (e)stimate based, (l)ogarithmic regression based, (o)ptimal value based (cheat!),
                                            * (r)ank-1 node based? */
#define DEFAULT_SETNAME               "-" /**< default settings file name for solving phase setting files */
#define DEFAULT_TRANSITIONMETHOD      'r' /**< the default transition method */
#define DEFAULT_NODEOFFSET            50L /**< default node offset before transition to proof phase is active */
#define DEFAULT_FALLBACK            FALSE /**< should the phase transition fall back to suboptimal phase? */
#define DEFAULT_INTERRUPTOPTIMAL    FALSE /**< should solving process be interrupted if optimal solution was found? */

#define DEFAULT_ENABLED             FALSE /**< should the event handler be executed? */
#define DEFAULT_TESTMODE            FALSE /**< should the event handler test the criteria? */

#define DEFAULT_USERESTART1TO2      FALSE /**< should a restart be applied between the feasibility and improvement phase? */
#define DEFAULT_USERESTART2TO3      FALSE /**< should a restart be applied between the improvement and the proof phase? */
#define DEFAULT_USEEMPHSETTINGS     TRUE  /**< should emphasis settings be used for the different solving phases, or settings files? */

/* logarithmic regression settings */
#define DEFAULT_LOGREGRESSION_XTYPE   'n' /**< default type to use for log regression - (t)ime, (n)odes, (l)p iterations */
#define LOGREGRESSION_XTYPES        "lnt" /**< available types for log regression - (t)ime, (n)odes, (l)p iterations */
/*
 * Data structures
 */

/** enumerator to represent the current solving phase */
enum SolvingPhase
{
   SOLVINGPHASE_UNINITIALIZED = -1,          /**< solving phase has not been initialized yet */
   SOLVINGPHASE_FEASIBILITY   = 0,           /**< no solution was found until now */
   SOLVINGPHASE_IMPROVEMENT   = 1,           /**< current incumbent solution is suboptimal */
   SOLVINGPHASE_PROOF         = 2            /**< current incumbent is optimal */
};
typedef enum SolvingPhase SOLVINGPHASE;

/** depth information structure */
struct DepthInfo
{
   int                   nsolvednodes;       /**< number of nodes that were solved so far at this depth */
   SCIP_Real             minestimate;        /**< the minimum estimate of a solved node */
   SCIP_NODE**           minnodes;           /**< points to the rank-1 nodes at this depth (open nodes whose estimate is lower than current
                                                  minimum estimate over solved nodes) */
   int                   nminnodes;          /**< the number of minimum nodes */
   int                   minnodescapacity;   /**< the capacity of the min nodes array */
};

typedef struct DepthInfo DEPTHINFO;

/** event handler data */
struct SCIP_EventhdlrData
{
   char                  logregression_xtype;/**< type to use for log regression - (t)ime, (n)odes, (l)p iterations */
   SCIP_Bool             enabled;            /**< should the event handler be executed? */
   char*                 feassetname;        /**< settings file parameter for the feasibility phase -- precedence over emphasis settings */
   char*                 improvesetname;     /**< settings file parameter for the improvement phase -- precedence over emphasis settings  */
   char*                 proofsetname;       /**< settings file parameter for the proof phase -- precedence over emphasis settings */
   SCIP_Real             optimalvalue;       /**< value of optimal solution of the problem */
   SCIP_Longint          nnodesleft;         /**< store the number of open nodes that are considered internally to update data */
   SOLVINGPHASE          solvingphase;       /**< the current solving phase */
   char                  transitionmethod;   /**< transition method from improvement phase -> proof phase?
                                               *  (e)stimate based, (l)ogarithmic regression based, (o)ptimal value based (cheat!),
                                               *  (r)ank-1 node based */
   SCIP_Longint          nodeoffset;         /**< node offset for triggering rank-1 node based phased transition */
   SCIP_Longint          lastndelayedcutoffs;/**< the number of delayed cutoffs since the last update of a focus node */
   SCIP_Bool             fallback;           /**< should the phase transition fall back to improvement phase? */
   SCIP_Bool             interruptoptimal;   /**< interrupt after optimal solution was found */
   SCIP_Bool             userestart1to2;     /**< should a restart be applied between the feasibility and improvement phase? */
   SCIP_Bool             userestart2to3;     /**< should a restart be applied between the improvement and the proof phase? */
   SCIP_Bool             useemphsettings;    /**< should emphasis settings for the solving phases be used, or settings files? */

   SCIP_Bool             testmode;           /**< should transitions be tested only, but not triggered? */
   SCIP_Bool             rank1reached;       /**< has the rank-1 transition into proof phase been reached? */
   SCIP_Bool             estimatereached;    /**< has the best-estimate transition been reached? */
   SCIP_Bool             optimalreached;     /**< is the incumbent already optimal? */
   SCIP_Bool             logreached;         /**< has a logarithmic phase transition been reached? */
   SCIP_Bool             newbestsol;         /**< has a new incumbent been found since the last node was solved? */

   SCIP_REGRESSION*      regression;         /**< regression data for log linear regression of the incumbent solutions */
   SCIP_Real             lastx;              /**< X-value of last observation */
   SCIP_Real             lasty;              /**< Y-value of last observation */
   SCIP_PARAM**          nondefaultparams;   /**< parameters with non-default values during problem initialization */
   int                   nnondefaultparams;  /**< number of parameters with non-default values during problem initialization  */
   int                   nondefaultparamssize;/**< capacity of the array of non-default parameters */
   int                   eventfilterpos;     /**< the event filter position, or -1, if event has not (yet) been caught */
   DEPTHINFO**           depthinfos;         /**< array of depth infos for every depth of the search tree */
   int                   maxdepth;           /**< maximum depth so far */
   int                   nrank1nodes;        /**< number of rank-1 nodes */
   int                   nnodesbelowincumbent;/**< number of open nodes with an estimate lower than the current incumbent */
};


/*
 * methods for rank-1 and active estimate transition
 */

/** nodes are sorted first by their estimates, and if estimates are equal, by their number */
static
SCIP_DECL_SORTPTRCOMP(sortCompTreeinfo)
{
   SCIP_NODE* node1;
   SCIP_NODE* node2;
   SCIP_Real estim1;
   SCIP_Real estim2;
   node1 = (SCIP_NODE*)elem1;
   node2 = (SCIP_NODE*)elem2;

   estim1 = SCIPnodeGetEstimate(node1);
   estim2 = SCIPnodeGetEstimate(node2);

   /* compare estimates */
   if( estim1 < estim2 )
      return -1;
   else if( estim1 > estim2 )
      return 1;
   else
   {
      SCIP_Longint number1;
      SCIP_Longint number2;

      number1 = SCIPnodeGetNumber(node1);
      number2 = SCIPnodeGetNumber(node2);

      /* compare numbers */
      if( number1 < number2 )
         return -1;
      else if( number1 > number2 )
         return 1;
   }

   return 0;
}

/** insert an array of open nodes (leaves/siblings/children) into the event handler data structures and update the transition information */
static
SCIP_RETCODE addNodesInformation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   SCIP_NODE**           nodes,              /**< array of nodes */
   int                   nnodes              /**< number of nodes */
   )
{
   int n;

   assert(nnodes == 0 || nodes != NULL);
   assert(scip != NULL);
   assert(eventhdlrdata->depthinfos != NULL);

   /* store every relevant node in the data structure for its depth */
   for( n = 0; n < nnodes; ++n )
   {
      SCIP_NODE* node = nodes[n];
      DEPTHINFO* depthinfo = eventhdlrdata->depthinfos[SCIPnodeGetDepth(node)];
      SCIP_Real estim = SCIPnodeGetEstimate(node);

      assert(SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD || SCIPnodeGetType(node) == SCIP_NODETYPE_LEAF
            || SCIPnodeGetType(node) == SCIP_NODETYPE_SIBLING);

      /* an open node has rank 1 if it has an estimate at least as small as the best solved node at this depth */
      if( depthinfo->nsolvednodes == 0 || SCIPisGE(scip, depthinfo->minestimate, SCIPnodeGetEstimate(node)) )
      {
         int pos;

         /* allocate additional memory to hold new node */
         if( depthinfo->nminnodes == depthinfo->minnodescapacity )
         {
            int oldcapacity = depthinfo->minnodescapacity;
            depthinfo->minnodescapacity *= 2;
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &depthinfo->minnodes, oldcapacity, depthinfo->minnodescapacity) );
         }

         /* find correct insert position */
         SCIPsortedvecInsertPtr((void **)depthinfo->minnodes, sortCompTreeinfo, (void*)node, &depthinfo->nminnodes, &pos);
         assert(pos >= 0 && pos < depthinfo->nminnodes);
         assert(depthinfo->minnodes[pos] == node);

         /* update rank 1 node information */
         ++eventhdlrdata->nrank1nodes;
      }

      /* update active estimate information by bookkeeping nodes with an estimate smaller than the current incumbent */
      if( SCIPisLT(scip, estim, SCIPgetUpperbound(scip) ) )
         ++eventhdlrdata->nnodesbelowincumbent;
   }

   /* update the number of open search nodes */
   eventhdlrdata->nnodesleft += nnodes;

   return SCIP_OKAY;
}

/** remove a node from the data structures of the event handler */
static
void removeNode(
   SCIP_NODE*            node,               /**< node that should be removed */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   DEPTHINFO* depthinfo;
   int pos;
   SCIP_Bool contained;

   assert(node != NULL);

   /* get depth information for the depth of this node */
   depthinfo = eventhdlrdata->depthinfos[SCIPnodeGetDepth(node)];

   /* no node is saved at this depth */
   if( depthinfo->nminnodes == 0 )
      return;

   /* search for the node by using binary search */
   contained = SCIPsortedvecFindPtr((void **)depthinfo->minnodes, sortCompTreeinfo, (void *)node, depthinfo->nminnodes, &pos);

   /* remove the node if it is contained */
   if( contained )
   {
      SCIPsortedvecDelPosPtr((void **)depthinfo->minnodes, sortCompTreeinfo, pos, &(depthinfo->nminnodes));
      --eventhdlrdata->nrank1nodes;
   }
}

/** returns the current number of rank 1 nodes in the tree */
static
int getNRank1Nodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(SCIPfindEventhdlr(scip, EVENTHDLR_NAME));

   /* return the stored number of rank 1 nodes only during solving stage */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      return eventhdlrdata->nrank1nodes;
   else
      return -1;
}

/** returns the current number of open nodes which have an estimate lower than the incumbent solution */
static
int getNNodesBelowIncumbent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(SCIPfindEventhdlr(scip, EVENTHDLR_NAME));

   /* return the stored number of nodes only during solving stage */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      return eventhdlrdata->nnodesbelowincumbent;
   else
      return -1;
}

/** discards all previous node information and renews it */
static
SCIP_RETCODE recomputeNodeInformation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_NODE** leaves;
   SCIP_NODE** children;
   SCIP_NODE** siblings;

   int nleaves;
   int nchildren;
   int nsiblings;
   int d;

   /* the required node information is only available after solving started */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   assert(eventhdlrdata != NULL);

   /* reset depth information */
   for( d = 0; d < eventhdlrdata->maxdepth; ++d )
      eventhdlrdata->depthinfos[d]->nminnodes = 0;

   eventhdlrdata->nrank1nodes = 0;
   eventhdlrdata->nnodesbelowincumbent = 0;
   eventhdlrdata->nnodesleft = 0;

   nleaves = nchildren = nsiblings = 0;

   /* get leaves, children, and sibling arrays and update the event handler data structures */
   SCIP_CALL( SCIPgetOpenNodesData(scip, &leaves, &children, &siblings, &nleaves, &nchildren, &nsiblings) );

   SCIP_CALL ( addNodesInformation(scip, eventhdlrdata, children, nchildren) );

   SCIP_CALL ( addNodesInformation(scip, eventhdlrdata, siblings, nsiblings) );

   SCIP_CALL ( addNodesInformation(scip, eventhdlrdata, leaves, nleaves) );

   /* information needs to be recomputed from scratch if a new incumbent is found */
   eventhdlrdata->newbestsol = FALSE;

   return SCIP_OKAY;
}

/** allocates memory for a depth info */
static
SCIP_RETCODE createDepthinfo(
   SCIP*                 scip,               /**< SCIP data structure */
   DEPTHINFO**           depthinfo           /**< pointer to depth information structure */
   )
{
   assert(scip != NULL);
   assert(depthinfo != NULL);

   /* allocate the necessary memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, depthinfo) );

   /* reset the depth information */
   (*depthinfo)->minestimate = SCIPinfinity(scip);
   (*depthinfo)->nsolvednodes = 0;
   (*depthinfo)->nminnodes = 0;
   (*depthinfo)->minnodescapacity = 2;

   /* allocate array to store nodes */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*depthinfo)->minnodes, (*depthinfo)->minnodescapacity) );

   return SCIP_OKAY;
}

/** frees depth information data structure */
static
SCIP_RETCODE freeDepthinfo(
   SCIP*                 scip,               /**< SCIP data structure */
   DEPTHINFO**           depthinfo           /**< pointer to depth information structure */
   )
{
   assert(scip != NULL);
   assert(depthinfo != NULL);
   assert(*depthinfo != NULL);
   assert((*depthinfo)->minnodes != NULL);

   /* free nodes data structure and then the structure itself */
   SCIPfreeBlockMemoryArray(scip, &(*depthinfo)->minnodes, (*depthinfo)->minnodescapacity);
   SCIPfreeBlockMemory(scip, depthinfo);

   return SCIP_OKAY;
}

/** removes the node itself and updates the data if this node defined an active estimate globally or locally at its depth level */
static
void releaseNodeFromDepthInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   SCIP_NODE*            node                /**< node to be removed from the data structures of the event handler */
   )
{
   DEPTHINFO* depthinfo;

   assert(scip != NULL);
   assert(node != NULL);
   assert(eventhdlrdata != NULL);

   /* get the correct depth info at the node depth */
   depthinfo = eventhdlrdata->depthinfos[SCIPnodeGetDepth(node)];
   assert(depthinfo != NULL);

   /* remove the node from the data structures */
   removeNode(node, eventhdlrdata);

   /* compare the node estimate to the minimum estimate of the particular depth */
   if( SCIPisLT(scip, SCIPnodeGetEstimate(node), depthinfo->minestimate) )
      depthinfo->minestimate = SCIPnodeGetEstimate(node);

   /* decrease counter of active estimate nodes if node has an estimate that is below the current incumbent */
   if( SCIPisLT(scip, SCIPnodeGetEstimate(node), SCIPgetUpperbound(scip)) && SCIPnodeGetDepth(node) > 0 )
      eventhdlrdata->nnodesbelowincumbent--;

   /* loop over remaining, unsolved nodes and decide whether they are still rank-1 nodes */
   while( depthinfo->nminnodes > 0 && SCIPisGT(scip, SCIPnodeGetEstimate(depthinfo->minnodes[depthinfo->nminnodes - 1]), depthinfo->minestimate) )
   {
      /* forget about node */
      --(depthinfo->nminnodes);
      --(eventhdlrdata->nrank1nodes);
   }

   /* increase the number of solved nodes at this depth */
   ++(depthinfo->nsolvednodes);

   /* decrease the counter for the number of open nodes */
   --eventhdlrdata->nnodesleft;
}

/** ensures sufficient size for depthInfo array */
static
SCIP_RETCODE ensureDepthInfoArraySize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   SCIP_NODE*            node                /**< node to be removed from the data structures of the event handler */
   )
{
   int nodedepth;
   int newsize;
   int oldsize;
   nodedepth = SCIPnodeGetDepth(node);
   oldsize = eventhdlrdata->maxdepth;
   newsize = oldsize;

   /* create depth info array with small initial size or enlarge the existing array if new node is deeper */
   if( oldsize == 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &eventhdlrdata->depthinfos, 10) );
      newsize = 10;
   }
   else if( nodedepth + 1 >= eventhdlrdata->maxdepth )
   {
      assert(nodedepth > 0);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &eventhdlrdata->depthinfos, oldsize, 2 * nodedepth) ); /*lint !e647*/
      newsize = 2 * nodedepth;
   }

   /* create the according depth information pointers */
   if( newsize > oldsize )
   {
      int c;

      for( c = oldsize; c < newsize; ++c )
      {
         SCIP_CALL( createDepthinfo(scip, &(eventhdlrdata->depthinfos[c])) );

      }

      eventhdlrdata->maxdepth = newsize;
   }
   assert(newsize > nodedepth);

   return SCIP_OKAY;
}

/** ensures the capacity of the event handler data structures and removes the current node */
static
SCIP_RETCODE releaseNodeInformation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   SCIP_NODE*            node                /**< node to be removed from the data structures of the event handler */
   )
{

   assert(scip != NULL);
   assert(node != NULL);
   assert(eventhdlrdata != NULL);

   /* ensure the depth info data structure can hold this node */
   SCIP_CALL( ensureDepthInfoArraySize(scip, eventhdlrdata, node) );


   /* in case that selected nodes were cut off in between two calls to this method, build data structures from scratch again */
   if( SCIPgetNDelayedCutoffs(scip) > eventhdlrdata->lastndelayedcutoffs || eventhdlrdata->newbestsol
         || eventhdlrdata->nnodesleft - 1 != SCIPgetNNodesLeft(scip) )
   {
      SCIP_CALL( recomputeNodeInformation(scip, eventhdlrdata) );

      eventhdlrdata->lastndelayedcutoffs = SCIPgetNDelayedCutoffs(scip);
   }
   else
   {
      /* remove the node from the data structures */
      releaseNodeFromDepthInfo(scip, eventhdlrdata, node);
   }

   assert(eventhdlrdata->nnodesleft == SCIPgetNNodesLeft(scip));

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** ensures correctness of counters by explicitly summing up all children, leaves, and siblings with small estimates */
static
int checkLeavesBelowIncumbent(
   SCIP* scip
   )
{
   SCIP_NODE** nodes;
   SCIP_RETCODE retcode;
   int nnodes;
   int n;
   SCIP_Real upperbound = SCIPgetUpperbound(scip);
   int nodesbelow = 0;

   /* compare children estimate and current upper bound */
   retcode = SCIPgetChildren(scip, &nodes, &nnodes);
   assert(retcode == SCIP_OKAY);

   for( n = 0; n < nnodes; ++n )
   {
      if( SCIPisLT(scip, SCIPnodeGetEstimate(nodes[n]), upperbound) )
         ++nodesbelow;
   }

   /* compare sibling estimate and current upper bound */
   retcode = SCIPgetSiblings(scip, &nodes, &nnodes);
   assert(retcode == SCIP_OKAY);

   for( n = 0; n < nnodes; ++n )
   {
      if( SCIPisLT(scip, SCIPnodeGetEstimate(nodes[n]), upperbound) )
         ++nodesbelow;
   }

   /* compare leaf node and current upper bound */
   retcode = SCIPgetLeaves(scip, &nodes, &nnodes);
   assert(retcode == SCIP_OKAY);

   for( n = 0; n < nnodes; ++n )
   {
      if( SCIPisLT(scip, SCIPnodeGetEstimate(nodes[n]), upperbound) )
         ++nodesbelow;
   }

   assert(nodesbelow <= SCIPgetNNodesLeft(scip));
   return nodesbelow;
}
#endif

/** get the point of the X axis for the regression according to the user choice of X type (time/nodes/iterations)*/
static
SCIP_Real getX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_Real x;

   switch( eventhdlrdata->logregression_xtype )
   {
   case 'l':
      /* get number of LP iterations so far */
      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVED )
         x = (SCIP_Real)SCIPgetNLPIterations(scip);
      else
         x = 1.0;
      break;
   case 'n':
      /* get total number of solving nodes so far */
      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVED )
         x = (SCIP_Real)SCIPgetNTotalNodes(scip);
      else
         x = 1.0;
      break;
   case 't':
      /* get solving time */
      x = SCIPgetSolvingTime(scip);
      break;
   default:
      x = 1.0;
      break;
   }

   /* prevent the calculation of logarithm too close to zero */
   x = MAX(x, .1);
   x = log(x);

   return x;
}





/** get axis intercept of current tangent to logarithmic regression curve */
static
SCIP_Real getCurrentRegressionTangentAxisIntercept(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data structure */
   )
{
   SCIP_REGRESSION* regression;
   SCIP_Real currentx;
   SCIP_Real regressionslope;

   assert(scip != NULL);
   assert(eventhdlrdata != NULL);

   regression = eventhdlrdata->regression;
   assert(regression != NULL);

   /* don't rely on too few (<= 2) observations */
   if( SCIPregressionGetNObservations(regression) <= 2 )
      return SCIPinfinity(scip);

   currentx = getX(scip, eventhdlrdata);
   regressionslope = SCIPregressionGetSlope(regression);

   return regressionslope * currentx + SCIPregressionGetIntercept(regression) - regressionslope;
}

/*
 * Local methods
 */

/** checks if rank-1 transition has been reached, that is, when all open nodes have a best-estimate higher than the best
 *  previously checked node at this depth
 */
static
SCIP_Bool checkRankOneTransition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   /* at least one solution is required for the transition */
   if( SCIPgetNSols(scip) > 0 )
      return (SCIPgetNNodes(scip) > eventhdlrdata->nodeoffset && getNRank1Nodes(scip) == 0);
   else
      return FALSE;
}

/** check if Best-Estimate criterion was reached, that is, when the active estimate is not better than the current incumbent solution */
static
SCIP_Bool checkEstimateCriterion(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   if( SCIPgetNSols(scip) > 0 )
      return ((SCIPgetNNodes(scip) > eventhdlrdata->nodeoffset) && (eventhdlrdata->nnodesbelowincumbent == 0));
   else
      return FALSE;
}

/** check if logarithmic phase transition has been reached.
 *
 *  the logarithmic phase transition is reached when the slope of the logarithmic primal progress (as a function of the number of
 *  LP iterations or solving nodes) becomes gentle. More concretely, we measure the slope by calculating the axis intercept of the tangent of
 *  the logarithmic primal progress. We then compare this axis intercept to the first and current primal bound and say that
 *  the logarithmic phase transition is reached as soon as the axis intercept passes the current primal bound so that the
 *  scalar becomes negative.
 *
 *  While it would be enough to directly compare the primal bound and the axis intercept of the
 *  tangent to check the criterion, the scalar allows for a continuous indicator how far the phase transition is still ahead
 */
static
SCIP_Bool checkLogCriterion(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real axisintercept = getCurrentRegressionTangentAxisIntercept(scip, eventhdlrdata);
      if( !SCIPisInfinity(scip, axisintercept) )
      {
         SCIP_Real primalbound;
         SCIP_Real lambda;
         SCIP_Real firstprimalbound = SCIPgetFirstPrimalBound(scip);

         primalbound = SCIPgetPrimalbound(scip);

         /* lambda is the scalar to describe the axis intercept as a linear combination of the current and the first primal bound
          * as intercept = pb_0 + lambda * (pb - pb_0) */
         lambda = (axisintercept - primalbound) / (firstprimalbound - primalbound);

         if( SCIPisNegative(scip, lambda) )
            return TRUE;
      }
   }
   return FALSE;
}

/** check if incumbent solution is nearly optimal; we allow a relative deviation of 10^-9 */
static
SCIP_Bool checkOptimalSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_Real referencevalue;
   SCIP_Real primalbound;

   referencevalue = eventhdlrdata->optimalvalue;
   primalbound = SCIPgetPrimalbound(scip);

   if(!SCIPisInfinity(scip, REALABS(primalbound)) && !SCIPisInfinity(scip, referencevalue) )
   {
      SCIP_Real max = MAX3(1.0, REALABS(primalbound), REALABS(referencevalue)); /*lint !e666*/

      if( EPSZ((primalbound - referencevalue)/max, 1e-9) )
         return TRUE;
   }
   return FALSE;
}

/** check if we are in the proof phase */
static
SCIP_Bool transitionPhase3(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   if( eventhdlrdata->solvingphase == SOLVINGPHASE_PROOF && !eventhdlrdata->fallback )
      return TRUE;

   /* check criterion based on selected transition method */
   switch( eventhdlrdata->transitionmethod )
   {
      case 'r':

         /* check rank-1 transition */
         if( checkRankOneTransition(scip, eventhdlrdata) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "reached rank-1 transition: nodes: %lld, rank-1: %d bound: %9.5g time: %.2f\n",
                  SCIPgetNNodes(scip), getNRank1Nodes(scip), SCIPgetPrimalbound(scip), SCIPgetSolvingTime(scip));
            return TRUE;
         }
         break;
      case 'o':

         /* cheat and use knowledge about optimal solution */
         if( checkOptimalSolution(scip, eventhdlrdata) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "optimal solution found: %lld, bound: %9.5g time: %.2f\n",
                  SCIPgetNNodes(scip), SCIPgetPrimalbound(scip), SCIPgetSolvingTime(scip));
            return TRUE;
         }
         break;
      case 'e':

         /* check best-estimate transition */
         if( checkEstimateCriterion(scip, eventhdlrdata) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "reached best-estimate transition: nodes: %lld, estimate: %d bound: %9.5g time: %.2f\n",
                  SCIPgetNNodes(scip), eventhdlrdata->nnodesbelowincumbent, SCIPgetPrimalbound(scip), SCIPgetSolvingTime(scip));
            return TRUE;
         }
         return FALSE;
      case 'l':

         /* check logarithmic transition */
         if( checkLogCriterion(scip, eventhdlrdata) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "reached a logarithmic phase transition: %.2f\n", SCIPgetSolvingTime(scip));
            return TRUE;
         }
         break;
      default:
         return FALSE;
   }

   return FALSE;
}

/* determine the solving phase: feasibility phase if no solution was found yet, otherwise improvement phase or proof phase
 * depending on whether selected transition criterion was already reached and fallback is active or not
 */
static
void determineSolvingPhase(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   /* without solution, we are in the feasibility phase */
   if( SCIPgetNSols(scip) == 0 )
      eventhdlrdata->solvingphase = SOLVINGPHASE_FEASIBILITY;
   else if( eventhdlrdata->solvingphase != SOLVINGPHASE_PROOF || eventhdlrdata->fallback )
      eventhdlrdata->solvingphase = SOLVINGPHASE_IMPROVEMENT;

   if( eventhdlrdata->solvingphase == SOLVINGPHASE_IMPROVEMENT && transitionPhase3(scip, eventhdlrdata) )
      eventhdlrdata->solvingphase = SOLVINGPHASE_PROOF;
}

/** changes parameters by using emphasis settings */
static
SCIP_RETCODE changeEmphasisParameters(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_PARAMEMPHASIS paramemphasis;

   /* choose the appropriate emphasis settings for the new solving phase */
   switch(eventhdlrdata->solvingphase)
   {
      case SOLVINGPHASE_FEASIBILITY:
         paramemphasis = SCIP_PARAMEMPHASIS_PHASEFEAS;
         break;
      case SOLVINGPHASE_IMPROVEMENT:
         paramemphasis = SCIP_PARAMEMPHASIS_PHASEIMPROVE;
         break;
      case SOLVINGPHASE_PROOF:
         paramemphasis = SCIP_PARAMEMPHASIS_PHASEPROOF;
         break;
      case SOLVINGPHASE_UNINITIALIZED:
      default:
         SCIPdebugMsg(scip, "Unknown solving phase: %d -> ABORT!\n ", eventhdlrdata->solvingphase);
         SCIPABORT();
         paramemphasis = SCIP_PARAMEMPHASIS_DEFAULT;
         break;
   }

   SCIP_CALL( SCIPsetEmphasis(scip, paramemphasis, FALSE) );

   return SCIP_OKAY;
}

/** change general solving strategy of SCIP depending on the phase by reading from settings file */
static
SCIP_RETCODE changeParametersUsingSettingsFiles(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   FILE* file;
   char* paramfilename = NULL;

   /* choose the settings file for the new solving phase */
   switch(eventhdlrdata->solvingphase)
   {
      case SOLVINGPHASE_FEASIBILITY:
         paramfilename = eventhdlrdata->feassetname;
         break;
      case SOLVINGPHASE_IMPROVEMENT:
         paramfilename = eventhdlrdata->improvesetname;
         break;
      case SOLVINGPHASE_PROOF:
         paramfilename = eventhdlrdata->proofsetname;
         break;
      case SOLVINGPHASE_UNINITIALIZED:
      default:
         SCIPdebugMsg(scip, "Unknown solving phase: %d -> ABORT!\n ", eventhdlrdata->solvingphase);
         return SCIP_INVALIDCALL;
   }

   assert(paramfilename != NULL);

   /* return if no there is no user-specified settings file for the current phase */
   if( strcmp(paramfilename, DEFAULT_SETNAME) == 0 )
      return SCIP_OKAY;

   file = fopen(paramfilename, "r");

   /* test if file could be found and print a warning if not */
   if( file == NULL )
   {
      SCIPwarningMessage(scip, "Parameter file <%s> not found--keeping settings as before.\n", paramfilename);
   }
   else
   {
      /* we can close the file */
      fclose(file);

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Reading parameters from file <%s>\n", paramfilename);

      SCIP_CALL( SCIPreadParams(scip, paramfilename) );
   }

   return SCIP_OKAY;
}

/** fix/unfix relevant solving parameters that should not accidentally be set to default values */
static
SCIP_RETCODE fixOrUnfixRelevantParameters(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   SCIP_Bool             fix                 /**< should the parameters be fixed (true) or unfixed? */
   )
{
   int p;
   const char* relevantparams[] = {
      "limits/time",
      "limits/nodes",
      "limits/totalnodes",
      "limits/stallnodes",
      "limits/memory",
      "limits/gap",
      "limits/absgap",
      "limits/solutions",
      "limits/bestsol",
      "limits/maxsol",
      "limits/maxorigsol",
      "limits/restarts",
      "limits/autorestartnodes",
      "limits/softtime",
      "solvingphases/enabled",
      "solvingphases/fallback",
      "solvingphases/interruptoptimal",
      "solvingphases/nodeoffset",
      "solvingphases/feassetname",
      "solvingphases/proofsetname",
      "solvingphases/optimalvalue",
      "solvingphases/improvesetname",
      "solvingphases/testmode",
      "solvingphases/transitionmethod",
      "solvingphases/useemphsettings",
      "solvingphases/userestart1to2",
      "solvingphases/userestart2to3",
      "solvingphases/xtype"
   };
   int nrelevantparams = 28;

   /* fix or unfix all specified limit parameters */
   for( p = 0; p < nrelevantparams; ++p )
   {
      if( fix )
      {
         SCIP_CALL( SCIPfixParam(scip, relevantparams[p]) );
      }
      else
      {
         SCIP_CALL( SCIPunfixParam(scip, relevantparams[p]) );
      }
   }

   /* fix or unfix all collected, non-default parameters after problem transformation */
   for( p = 0; p < eventhdlrdata->nnondefaultparams; ++p )
   {
      if( fix && ! SCIPparamIsFixed(eventhdlrdata->nondefaultparams[p]) )
      {
         SCIP_CALL( SCIPfixParam(scip, SCIPparamGetName(eventhdlrdata->nondefaultparams[p])) );
      }
      else if( ! fix && SCIPparamIsFixed(eventhdlrdata->nondefaultparams[p]) )
      {
         SCIP_CALL( SCIPunfixParam(scip, SCIPparamGetName(eventhdlrdata->nondefaultparams[p])) );
      }
   }

   return SCIP_OKAY;
}

/** change settings depending whether emphasis settings should be used, or settings files */
static
SCIP_RETCODE adaptSolverBehavior(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   /* fix relevant parameters such that they are not overwritten */
   SCIP_CALL( fixOrUnfixRelevantParameters(scip, eventhdlrdata, TRUE) );

   /* change settings using emphasis */
   if( eventhdlrdata->useemphsettings )
   {
      SCIP_CALL( changeEmphasisParameters(scip, eventhdlrdata) );
   }
   else
   {
      /* reset to default settings; this happens automatically when using emphasis settings */
      SCIP_CALL( SCIPsetEmphasis(scip, SCIP_PARAMEMPHASIS_DEFAULT, FALSE) );
   }

   /* read optional, phase-specific settings */
   SCIP_CALL( changeParametersUsingSettingsFiles(scip, eventhdlrdata) );

   /* unfix relevant parameters that have been fixed for changing emphasis */
   SCIP_CALL( fixOrUnfixRelevantParameters(scip, eventhdlrdata, FALSE) );

   return SCIP_OKAY;
}

/* apply the user-specified phase-based settings: A phase transition invokes the read of phase-specific settings from a file */
static
SCIP_RETCODE applySolvingPhase(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SOLVINGPHASE oldsolvingphase;
   SCIP_Bool restart;

   /* return immediately if we are in the proof phase */
   if( eventhdlrdata->solvingphase == SOLVINGPHASE_PROOF && !eventhdlrdata->fallback )
      return SCIP_OKAY;

   /* save current solving phase */
   oldsolvingphase = eventhdlrdata->solvingphase;

   /* determine current solving phase */
   determineSolvingPhase(scip, eventhdlrdata);


   /* nothing has changed */
   if( oldsolvingphase == eventhdlrdata->solvingphase )
      return SCIP_OKAY;


   /* check if the solving process should be interrupted when the current solution is optimal */
   if( eventhdlrdata->solvingphase == SOLVINGPHASE_PROOF && eventhdlrdata->transitionmethod == 'o' &&
         eventhdlrdata->interruptoptimal )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Solution is optimal. Calling user interruption.\n");

      /* we call interrupt solve but do not return yet because user-specified settings for the proof phase are applied first */
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   /* check if a restart should be performed after phase transition */
   if( eventhdlrdata->solvingphase == SOLVINGPHASE_IMPROVEMENT && eventhdlrdata->userestart1to2 )
      restart = TRUE;
   else if( eventhdlrdata->solvingphase == SOLVINGPHASE_PROOF && eventhdlrdata->userestart2to3 )
      restart = TRUE;
   else
      restart = FALSE;

   /* inform SCIP that a restart should be performed */
   if( restart )
   {
      SCIP_CALL( SCIPrestartSolve(scip) );
   }

   /* change general solving settings depending on solving strategy */
   SCIP_CALL( adaptSolverBehavior(scip, eventhdlrdata) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,"Changed solving phase to phase %d.\n", eventhdlrdata->solvingphase);

   return SCIP_OKAY;

}

/** update the logarithmic regression */
static
SCIP_RETCODE updateLogRegression(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< data of event handler */
   )
{
   SCIP_Real regressionx;
   SCIP_Real regressiony;

   regressionx = getX(scip, eventhdlrdata);
   regressiony = SCIPgetPrimalbound(scip);

   /* remove the last observation if it has been observed at the same x */
   if( SCIPisEQ(scip, eventhdlrdata->lastx, regressionx) )
   {
      SCIPregressionRemoveObservation(eventhdlrdata->regression, eventhdlrdata->lastx, eventhdlrdata->lasty);
   }

   /* add the new observation to the regression and save it if another update is necessary */
   SCIPregressionAddObservation(eventhdlrdata->regression, regressionx, regressiony);
   eventhdlrdata->lastx = regressionx;
   eventhdlrdata->lasty = regressiony;

   return SCIP_OKAY;
}

/** update data structures based on the event type caught */
static
SCIP_RETCODE updateDataStructures(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< data of event handler */
   SCIP_EVENTTYPE        eventtype           /**< type of the caught event */
   )
{
   SCIP_NODE** children;
   int nchildren;

   switch( eventtype )
   {
      /* store that a new best solution was found, but delay the update of node information until a node was solved */
      case SCIP_EVENTTYPE_BESTSOLFOUND:
         eventhdlrdata->newbestsol = TRUE;

         /* update logarithmic regression of solution process */
         SCIP_CALL( updateLogRegression(scip, eventhdlrdata) );

         break;

      /* release the focus node from the open node data structures */
      case SCIP_EVENTTYPE_NODEFOCUSED:
         assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

         SCIP_CALL( releaseNodeInformation(scip, eventhdlrdata, SCIPgetCurrentNode(scip)));
         assert(eventhdlrdata->nnodesbelowincumbent <= SCIPgetNNodesLeft(scip));

         break;

      /* store node information for child nodes */
      case SCIP_EVENTTYPE_NODEBRANCHED:
         assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

         /* if we lost track of exact number of open search nodes, we recompute node information from scratch */
         if( eventhdlrdata->newbestsol || eventhdlrdata->nnodesleft + SCIPgetNChildren(scip) != SCIPgetNNodesLeft(scip) )
         {
            SCIP_CALL( recomputeNodeInformation(scip, eventhdlrdata) );
            eventhdlrdata->newbestsol = FALSE;

            return SCIP_OKAY;
         }
         else
         {
            SCIP_CALL( SCIPgetChildren(scip, &children, &nchildren) );
            SCIP_CALL( addNodesInformation(scip, eventhdlrdata, children, nchildren) );
         }

         assert(eventhdlrdata->nnodesleft == SCIPgetNNodesLeft(scip));
         break;

      default:
         break;
   }

   /* ensure that required tree information was correctly computed; only available in solving stage and at the beginning
    * or end of a node solution process because we delay the recomputation of the node information)
    */
   assert(SCIPgetStage(scip) != SCIP_STAGE_SOLVING ||
          (eventtype == SCIP_EVENTTYPE_BESTSOLFOUND) ||
          (eventhdlrdata->nnodesleft == SCIPgetNNodesLeft(scip) && eventhdlrdata->nnodesbelowincumbent == checkLeavesBelowIncumbent(scip)));

   return SCIP_OKAY;
}

/** test all criteria whether they have been reached */
static
void testCriteria(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< data of event handler */
   )
{
   assert(scip != NULL);
   assert(eventhdlrdata != NULL);

   if( ! eventhdlrdata->logreached && checkLogCriterion(scip, eventhdlrdata) )
   {
      eventhdlrdata->logreached = TRUE;
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Log criterion reached after %lld nodes, %.2f sec.\n",
         SCIPgetNNodes(scip), SCIPgetSolvingTime(scip));
   }
   if( ! eventhdlrdata->rank1reached && checkRankOneTransition(scip, eventhdlrdata) )
   {
      eventhdlrdata->rank1reached = TRUE;
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Rank 1 criterion reached after %lld nodes, %.2f sec.\n",
         SCIPgetNNodes(scip), SCIPgetSolvingTime(scip));
   }

   if( ! eventhdlrdata->estimatereached && checkEstimateCriterion(scip, eventhdlrdata) )
   {
      eventhdlrdata->estimatereached = TRUE;
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Estimate criterion reached after %lld nodes, %.2f sec.\n",
         SCIPgetNNodes(scip), SCIPgetSolvingTime(scip));
   }

   if( ! eventhdlrdata->optimalreached && checkOptimalSolution(scip, eventhdlrdata) )
   {
      eventhdlrdata->optimalreached = TRUE;
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Optimum reached after %lld nodes, %.2f sec.\n",
         SCIPgetNNodes(scip), SCIPgetSolvingTime(scip));
   }
}

/*
 * Callback methods of event handler
 */

/** copy method for event handler (called when SCIP copies plugins) */
/* todo this code needs to stay disabled as long as the soft limit event handler is not copied, because we save
 * the soft time limit parameter but this will crash as soon as we are in a SCIP copy */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_EVENTCOPY(eventCopySolvingphase)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* call inclusion method of event handler */
   SCIP_CALL( SCIPincludeEventHdlrSolvingphase(scip) );

   return SCIP_OKAY;
}
#else
#define eventCopySolvingphase NULL
#endif

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeSolvingphase)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPregressionFree(&eventhdlrdata->regression);

   SCIPfreeBlockMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolSolvingphase)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   eventhdlrdata->depthinfos = NULL;
   eventhdlrdata->maxdepth = 0;
   eventhdlrdata->nnodesbelowincumbent = 0;
   eventhdlrdata->nnodesleft = 0;
   eventhdlrdata->nrank1nodes = 0;
   eventhdlrdata->lastndelayedcutoffs = SCIPgetNDelayedCutoffs(scip);
   eventhdlrdata->newbestsol = FALSE;

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolSolvingphase)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   /* free all data storage acquired during this branch-and-bound run */
   if( eventhdlrdata->maxdepth > 0 )
   {
      int c;

      /* free depth information */
      for( c = 0; c < eventhdlrdata->maxdepth; ++c )
      {
         SCIP_CALL( freeDepthinfo(scip, &(eventhdlrdata->depthinfos[c])) );
      }

      /* free depth information array */
      SCIPfreeBlockMemoryArray(scip, &eventhdlrdata->depthinfos, eventhdlrdata->maxdepth);
      eventhdlrdata->maxdepth = 0;
   }

   return SCIP_OKAY;
}

/** collects all parameters that are set to non-default values and stores them in eventhdlrdata */
static
SCIP_RETCODE collectNondefaultParams(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< data of event handler */
   )
{
   SCIP_PARAM** params;
   int nparams;
   int p;

   params = SCIPgetParams(scip);
   nparams = SCIPgetNParams(scip);

   eventhdlrdata->nnondefaultparams = 0;
   eventhdlrdata->nondefaultparams = NULL;
   eventhdlrdata->nondefaultparamssize = 0;

   /* loop over parameters and store the non-default ones */
   for( p = 0; p < nparams; ++p )
   {
      SCIP_PARAM* param = params[p];

      /* collect parameter if it is nondefault */
      if( ! SCIPparamIsDefault(param) )
      {
         if( eventhdlrdata->nnondefaultparams == 0 )
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &eventhdlrdata->nondefaultparams, 8) );
            eventhdlrdata->nondefaultparamssize = 8;
         }
         else if( eventhdlrdata->nnondefaultparams == eventhdlrdata->nondefaultparamssize )
         {
            eventhdlrdata->nondefaultparamssize *= 2;
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &eventhdlrdata->nondefaultparams, \
                  eventhdlrdata->nnondefaultparams, eventhdlrdata->nondefaultparamssize) );

         }

         eventhdlrdata->nondefaultparams[eventhdlrdata->nnondefaultparams++] = param;
      }
   }


   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitSolvingphase)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* initialize the solving phase */
   eventhdlrdata->solvingphase = SOLVINGPHASE_UNINITIALIZED;

   /* none of the transitions is reached yet */
   eventhdlrdata->optimalreached = FALSE;
   eventhdlrdata->logreached = FALSE;
   eventhdlrdata->rank1reached = FALSE;
   eventhdlrdata->estimatereached = FALSE;
   eventhdlrdata->nnondefaultparams = 0;
   eventhdlrdata->nondefaultparams = NULL;
   eventhdlrdata->nondefaultparamssize = 0;

   /* apply solving phase for the first time after problem was transformed to apply settings for the feasibility phase */
   if( eventhdlrdata->enabled )
   {
      /* collect non-default parameters */
      SCIP_CALL( collectNondefaultParams(scip, eventhdlrdata) );

      SCIP_CALL( applySolvingPhase(scip, eventhdlrdata) );
   }

   /* only start catching events if event handler is enabled or in test mode */
   if( eventhdlrdata->enabled || eventhdlrdata->testmode )
   {
      SCIP_CALL( SCIPcatchEvent(scip, EVENTHDLR_EVENT, eventhdlr, NULL, &eventhdlrdata->eventfilterpos) );
   }

   /* reset solving regression */
   SCIPregressionReset(eventhdlrdata->regression);
   eventhdlrdata->lastx = SCIP_INVALID;
   eventhdlrdata->lasty = SCIP_INVALID;

   return SCIP_OKAY;
}
/** deinitialization method of event handler (called before problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitSolvingphase)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* free collected, non-default parameters */
   SCIPfreeBlockMemoryArrayNull(scip, &eventhdlrdata->nondefaultparams, eventhdlrdata->nondefaultparamssize);

   return SCIP_OKAY;
}


/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecSolvingphase)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   eventtype = SCIPeventGetType(event);
   assert(eventtype & (EVENTHDLR_EVENT));
   assert(eventtype != SCIP_EVENTTYPE_NODEFOCUSED || SCIPeventGetNode(event) == SCIPgetCurrentNode(scip));


   /* update data structures depending on the event */
   SCIP_CALL( updateDataStructures(scip, eventhdlrdata, eventtype) );

   /* if the phase-based solver is enabled, we check if a phase transition occurred and alter the settings accordingly */
   if( eventhdlrdata->enabled )
   {
      SCIP_CALL( applySolvingPhase(scip, eventhdlrdata) );
   }


   /* in test mode, we check every transition criterion */
   if( eventhdlrdata->testmode )
   {
      testCriteria(scip, eventhdlrdata);
   }

   return SCIP_OKAY;
}

/*
 * displays that come with this event handler
 */

/* defines for the rank 1 node display */
#define DISP_NAME_NRANK1NODES         "nrank1nodes"
#define DISP_DESC_NRANK1NODES         "current number of rank1 nodes left"
#define DISP_HEAD_NRANK1NODES         "rank1"
#define DISP_WIDT_NRANK1NODES         7
#define DISP_PRIO_NRANK1NODES         40000
#define DISP_POSI_NRANK1NODES         500
#define DISP_STRI_NRANK1NODES         TRUE

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputNRank1Nodes)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_NRANK1NODES) == 0);
   assert(scip != NULL);

   /* ouput number of rank 1 nodes */
   SCIPdispInt(SCIPgetMessagehdlr(scip), file, getNRank1Nodes(scip), DISP_WIDT_NRANK1NODES);

   return SCIP_OKAY;
}

/* display for the number of nodes below the current incumbent */
#define DISP_NAME_NNODESBELOWINC         "nnodesbelowinc"
#define DISP_DESC_NNODESBELOWINC         "current number of nodes with an estimate better than the current incumbent"
#define DISP_HEAD_NNODESBELOWINC         "nbInc"
#define DISP_WIDT_NNODESBELOWINC         6
#define DISP_PRIO_NNODESBELOWINC         40000
#define DISP_POSI_NNODESBELOWINC         550
#define DISP_STRI_NNODESBELOWINC         TRUE

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputNnodesbelowinc)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_NNODESBELOWINC) == 0);
   assert(scip != NULL);

   /* display the number of nodes with an estimate below the the current incumbent */
   SCIPdispInt(SCIPgetMessagehdlr(scip), file, getNNodesBelowIncumbent(scip), DISP_WIDT_NNODESBELOWINC);

   return SCIP_OKAY;
}

/** creates event handler for Solvingphase event */
SCIP_RETCODE SCIPincludeEventHdlrSolvingphase(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create solving phase event handler data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocBlockMemory(scip, &eventhdlrdata) );
   assert(eventhdlrdata != NULL);

   eventhdlrdata->feassetname = NULL;
   eventhdlrdata->improvesetname = NULL;
   eventhdlrdata->proofsetname = NULL;

   eventhdlrdata->depthinfos = NULL;
   eventhdlrdata->maxdepth = 0;
   eventhdlrdata->eventfilterpos = -1;

   /* create a regression */
   eventhdlrdata->regression = NULL;
   SCIP_CALL( SCIPregressionCreate(&eventhdlrdata->regression) );

   eventhdlr = NULL;

   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecSolvingphase, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* include the new displays into scip */
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_NRANK1NODES, DISP_DESC_NRANK1NODES, DISP_HEAD_NRANK1NODES, SCIP_DISPSTATUS_OFF,
         NULL, NULL, NULL, NULL, NULL, NULL, dispOutputNRank1Nodes, NULL, DISP_WIDT_NRANK1NODES, DISP_PRIO_NRANK1NODES, DISP_POSI_NRANK1NODES,
         DISP_STRI_NRANK1NODES) );
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_NNODESBELOWINC, DISP_DESC_NNODESBELOWINC, DISP_HEAD_NNODESBELOWINC, SCIP_DISPSTATUS_OFF,
         NULL, NULL, NULL, NULL, NULL, NULL, dispOutputNnodesbelowinc, NULL, DISP_WIDT_NNODESBELOWINC, DISP_PRIO_NNODESBELOWINC, DISP_POSI_NNODESBELOWINC,
         DISP_STRI_NNODESBELOWINC) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopySolvingphase) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeSolvingphase) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitSolvingphase) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitSolvingphase) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolSolvingphase) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolSolvingphase) );

   /* add Solvingphase event handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, EVENTHDLR_NAME "s/enabled", "should the event handler adapt the solver behavior?",
         &eventhdlrdata->enabled, FALSE, DEFAULT_ENABLED, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, EVENTHDLR_NAME "s/testmode", "should the event handler test all phase transitions?",
         &eventhdlrdata->testmode, FALSE, DEFAULT_TESTMODE, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, EVENTHDLR_NAME "s/feassetname", "settings file for feasibility phase -- precedence over emphasis settings",
         &eventhdlrdata->feassetname, FALSE, DEFAULT_SETNAME, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, EVENTHDLR_NAME "s/improvesetname", "settings file for improvement phase -- precedence over emphasis settings",
         &eventhdlrdata->improvesetname, FALSE, DEFAULT_SETNAME, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, EVENTHDLR_NAME "s/proofsetname", "settings file for proof phase -- precedence over emphasis settings",
         &eventhdlrdata->proofsetname, FALSE, DEFAULT_SETNAME, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, EVENTHDLR_NAME "s/nodeoffset", "node offset for rank-1 and estimate transitions", &eventhdlrdata->nodeoffset,
         FALSE, DEFAULT_NODEOFFSET, 1L, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, EVENTHDLR_NAME "s/fallback", "should the event handler fall back from optimal phase?",
         &eventhdlrdata->fallback, FALSE, DEFAULT_FALLBACK, NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip ,EVENTHDLR_NAME "s/transitionmethod",
         "transition method: Possible options are 'e'stimate,'l'ogarithmic regression,'o'ptimal-value based,'r'ank-1",
         &eventhdlrdata->transitionmethod, FALSE, DEFAULT_TRANSITIONMETHOD, TRANSITIONMETHODS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, EVENTHDLR_NAME "s/interruptoptimal",
         "should the event handler interrupt the solving process after optimal solution was found?",
         &eventhdlrdata->interruptoptimal, FALSE, DEFAULT_INTERRUPTOPTIMAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, EVENTHDLR_NAME "s/userestart1to2",
         "should a restart be applied between the feasibility and improvement phase?",
         &eventhdlrdata->userestart1to2, FALSE, DEFAULT_USERESTART1TO2, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, EVENTHDLR_NAME "s/userestart2to3",
         "should a restart be applied between the improvement and the proof phase?",
         &eventhdlrdata->userestart2to3, FALSE, DEFAULT_USERESTART2TO3, NULL, NULL) );

   SCIP_CALL(SCIPaddRealParam(scip, EVENTHDLR_NAME "s/optimalvalue", "optimal solution value for problem",
         &eventhdlrdata->optimalvalue, FALSE, SCIP_INVALID, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   /* add parameter for logarithmic regression */
   SCIP_CALL( SCIPaddCharParam(scip, EVENTHDLR_NAME "s/xtype", "x-type for logarithmic regression - (t)ime, (n)odes, (l)p iterations",
        &eventhdlrdata->logregression_xtype, FALSE, DEFAULT_LOGREGRESSION_XTYPE, LOGREGRESSION_XTYPES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, EVENTHDLR_NAME "s/useemphsettings",
      "should emphasis settings for the solving phases be used, or settings files?",
      &eventhdlrdata->useemphsettings, FALSE, DEFAULT_USEEMPHSETTINGS, NULL, NULL) );

   return SCIP_OKAY;
}
