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

/**@file   prop_vbounds.c
 * @brief  variable upper and lower bound propagator
 * @author Stefan Heinz
 * @author Jens Schulz
 * @author Gerald Gamrath
 *
 * This propagator uses global bound information provided by SCIP to deduce global and local bound changes.
 * It can take into account
 * - implications (bound change following from specific value of a binary variable)
 * - cliques (set of binary variables, each with a corresponding value, of which at most one variable can get the value)
 * - variable lower/upper bounds (bounds of arbitrary variables that depend linearly on the value of another variable)
 *
 * The propagator does not look at a variable in whole, but at one point in time only handles one specific bound (lower
 * or upper) of a variable and deduces changes for lower or upper bounds of other variables. The concept is as follows:
 *
 * 1) Extract variable bound data
 *
 *    Implications and cliques are stored in a way such that given a variable and its new value, we can access all bound
 *    changes that can be deduced from setting the variable to that value. However, for variable bounds, this currently
 *    does not hold, they are only stored in the other direction, i.e. for a bound of a given variable, we have a list
 *    of all other bounds of variables that directly influence the bound of the given variable and a linear function
 *    describing how they do this.
 *    For the propagation, we need the other direction, thus we store it in the propagator data when the branch-and-bound
 *    solving process is about to begin.
 *
 * 2) Topological sorting of bounds of variable
 *
 *    We compute a topological order of the bounds of variables. This is needed to define an order in which we will
 *    regard bounds of variables in the propagation process in order to avoid unneccessarily regarding the same variable
 *    bound multiple times because it was changed in the meantime when propagating another bound of a variable.
 *    Therefore, we implictly regard a directed graph, in which each node corresponds to a bound of a variable and there
 *    exists a directed edge from one node to another, if the bound corresponding to the former node influences the
 *    bound corresponding to the latter node. This is done by iteratively running a DFS until all nodes were visited.
 *    Note that there might be cycles in the graph, which are randomly broken, so the order is only almost topological.
 *
 * 3) Collecting bound changes
 *
 *    For each bound of a variable, which can trigger bound changes of other variables, the propagator catches all
 *    events informing about a global change of the bound or a local tightening of the bound. The event handler
 *    then adds the bound of the variable to a priority queue, with the key in the priority queue corresponding
 *    to the position of the bound in the topological sort.
 *
 * 4) Propagating Bounds
 *
 *    As long as there are bounds contained in the priority queue, the propagator pops one bound from the queue, which
 *    is the one most at the beginning of the topological sort, so it should not be influenced by propagating other
 *    bounds currently contained in the queue. Starting at this bound, all implication, clique, and variable bound
 *    information is used to deduce tigther bounds for other variables and change the bounds, if a tighter one is found.
 *    These bound changes trigger an event that will lead to adding the corresponding bound to the priority queue,
 *    if it is not contained, yet. The process is iterated until the priority queue contains no more bounds.
 *
 * Additionally, the propagator analyzes the conflict/clique graph during presolving. It uses Tarjan's algorithm to
 * search for strongly connected components, for each of which all variables can be aggregated to one. Additionally,
 * it may detect invalid assignments of binary variables and fix the variable to the only possible value left.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdint.h>

#include "scip/prop_vbounds.h"

/**@name Propagator properties
 *
 * @{
 */

#define PROP_NAME              "vbounds"
#define PROP_DESC              "propagates variable upper and lower bounds"
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP | SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY           3000000 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */

#define PROP_PRESOL_PRIORITY     -90000 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOLTIMING       SCIP_PRESOLTIMING_MEDIUM | SCIP_PRESOLTIMING_EXHAUSTIVE
#define PROP_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */
/**@} */

/**@name Event handler properties
 *
 * @{
 */

#define EVENTHDLR_NAME         "vbounds"
#define EVENTHDLR_DESC         "bound change event handler for for vbounds propagator"

/**@} */

/**@name Default parameter values
 *
 * @{
 */

#define DEFAULT_USEBDWIDENING      TRUE      /**< should bound widening be used to initialize conflict analysis? */
#define DEFAULT_USEIMPLICS         FALSE     /**< should implications be propagated? */
#define DEFAULT_USECLIQUES         FALSE     /**< should cliques be propagated? */
#define DEFAULT_USEVBOUNDS         TRUE      /**< should variable bounds be propagated? */
#define DEFAULT_DOTOPOSORT         TRUE      /**< should the bounds be topologically sorted in advance? */
#define DEFAULT_SORTCLIQUES        FALSE     /**< should cliques be regarded for the topological sort? */
#define DEFAULT_DETECTCYCLES       FALSE     /**< should cycles in the variable bound graph be identified? */
#define DEFAULT_MINNEWCLIQUES      0.1       /**< minimum number of new cliques to trigger another clique table analysis */
#define DEFAULT_MAXCLIQUESMEDIUM   50.0      /**< maximum number of cliques per variable to run clique table analysis in
                                              *   medium presolving */
#define DEFAULT_MAXCLIQUESEXHAUSTIVE 100.0   /**< maximum number of cliques per variable to run clique table analysis in
                                              *   exhaustive presolving */

/**@} */

/**@name Propagator defines
 *
 * @{
 *
 * The propagator works on indices representing a bound of a variable. This index will be called bound index in the
 * following. For a given active variable with problem index i (note that active variables have problem indices
 * between 0 and nactivevariable - 1), the bound index of its lower bound is 2*i, the bound index of its upper
 * bound is 2*i + 1. The other way around, a given bound index i corresponds to the variable with problem index
 * i/2 (rounded down), and to the lower bound, if i is even, to the upper bound if i is odd.
 * The following macros can be used to convert bound index into variable problem index and boundtype and vice versa.
 */
#define getLbIndex(idx) (2*(idx))
#define getUbIndex(idx) (2*(idx)+1)
#define getVarIndex(idx) ((idx)/2)
#define getBoundtype(idx) (((idx) % 2 == 0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER)
#define isIndexLowerbound(idx) ((idx) % 2 == 0)
#define getBoundString(lower) ((lower) ? "lb" : "ub")
#define getBoundtypeString(type) ((type) == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper")
#define indexGetBoundString(idx) (getBoundString(isIndexLowerbound(idx)))
#define getOtherBoundIndex(idx) ((idx) + 1 - 2 * ((idx) % 2))

/**@} */

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for catching bound changes */
   SCIP_VAR**            vars;               /**< array containing all variable which are considered within the propagator */
   SCIP_HASHMAP*         varhashmap;         /**< hashmap mapping from variable to index in the vars array */
   int*                  topoorder;          /**< array mapping on the bounds of variables in topological order;
                                              *   or -1, if the bound that should be at that position has no outgoing
                                              *   implications, cliques, or vbounds;
                                              *   i.e., for i < j and topoorder[i] != -1 != topoorder[j], the variable
                                              *   and boundtype represented by index topoorder[i] are earlier in the
                                              *   topological order than those represented by index topoorder[j]
                                              */
   int**                 vboundboundedidx;   /**< array storing for each bound index the bound indices of all bounds
                                              *   influenced by this bound through variable bounds */
   SCIP_Real**           vboundcoefs;        /**< array storing for each bound index the coefficients in the variable
                                              *   bounds influencing the corresponding bound index stored in
                                              *   vboundboundedidx */
   SCIP_Real**           vboundconstants;    /**< array storing for each bound index the constants in the variable
                                              *   bounds influencing the corresponding bound index stored in
                                              *   vboundboundedidx */
   int*                  nvbounds;           /**< array storing for each bound index the number of vbounds stored */
   int*                  vboundsize;         /**< array with sizes of vbound arrays for the nodes */
   int                   nbounds;            /**< number of bounds of variables regarded (two times number of active variables) */
   int                   lastpresolncliques; /**< number of cliques created until the last call to the presolver */
   SCIP_PQUEUE*          propqueue;          /**< priority queue to handle the bounds of variables that were changed and have to be propagated */
   SCIP_Bool*            inqueue;            /**< boolean array to store whether a bound of a variable is already contained in propqueue */
   SCIP_Bool             initialized;        /**< was the data for propagation already initialized? */
   SCIP_Real             minnewcliques;      /**< minimum percentage of new cliques to trigger another clique table analysis */
   SCIP_Real             maxcliquesmedium;   /**< maximum number of cliques per variable to run clique table analysis in medium presolving */
   SCIP_Real             maxcliquesexhaustive;/**< maximum number of cliques per variable to run clique table analysis in exhaustive presolving */
   SCIP_Bool             usebdwidening;      /**< should bound widening be used to initialize conflict analysis? */
   SCIP_Bool             useimplics;         /**< should implications be propagated? */
   SCIP_Bool             usecliques;         /**< should cliques be propagated? */
   SCIP_Bool             usevbounds;         /**< should variable bounds be propagated? */
   SCIP_Bool             dotoposort;         /**< should the bounds be topologically sorted in advance? */
   SCIP_Bool             sortcliques;        /**< should cliques be regarded for the topological sort? */
   SCIP_Bool             detectcycles;       /**< should cycles in the variable bound graph be identified? */
};

/** inference information */
struct InferInfo
{
   union
   {
      struct
      {
         unsigned int    pos:31;             /**< position of the variable which forced that propagation */
         unsigned int    boundtype:1;        /**< bound type which was the reason (0: lower, 1: upper) */
      } asbits;
      int                asint;              /**< inference information as a single int value */
   } val;
};
typedef struct InferInfo INFERINFO;

/** converts an integer into an inference information */
static
INFERINFO intToInferInfo(
   int                   i                   /**< integer to convert */
   )
{
   INFERINFO inferinfo;

   inferinfo.val.asint = i;

   return inferinfo;
}

/** converts an inference information into an int */
static
int inferInfoToInt(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asint;
}

/** returns the propagation rule stored in the inference information */
static
SCIP_BOUNDTYPE inferInfoGetBoundtype(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   assert((SCIP_BOUNDTYPE)inferinfo.val.asbits.boundtype == SCIP_BOUNDTYPE_LOWER
      || (SCIP_BOUNDTYPE)inferinfo.val.asbits.boundtype == SCIP_BOUNDTYPE_UPPER);
   return (SCIP_BOUNDTYPE)inferinfo.val.asbits.boundtype;
}

/** returns the position stored in the inference information */
static
int inferInfoGetPos(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return (int) inferinfo.val.asbits.pos;
}

/** constructs an inference information out of a position of a variable and a boundtype */
static
INFERINFO getInferInfo(
   int                   pos,                /**< position of the variable which forced that propagation */
   SCIP_BOUNDTYPE        boundtype           /**< propagation rule that deduced the value */
   )
{
   INFERINFO inferinfo;

   assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);
   assert((int)boundtype >= 0 && (int)boundtype <= 1); /*lint !e685 !e568q*/
   assert(pos >= 0);

   inferinfo.val.asbits.pos = (unsigned int) pos; /*lint !e732*/
   inferinfo.val.asbits.boundtype = (unsigned int) boundtype; /*lint !e641*/

   return inferinfo;
}

/*
 * Local methods
 */

/* returns the lower bound index of a variable */
static
int varGetLbIndex(
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var                 /**< variable to get the index for */
   )
{
   assert(SCIPhashmapExists(propdata->varhashmap, var) == ((size_t)SCIPhashmapGetImage(propdata->varhashmap, var) > 0));

   return getLbIndex((int)(size_t)SCIPhashmapGetImage(propdata->varhashmap, var) - 1);
}

/* returns the upper bound index of a variable */
static
int varGetUbIndex(
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var                 /**< variable to get the index for */
   )
{
   assert(SCIPhashmapExists(propdata->varhashmap, var) == ((size_t)SCIPhashmapGetImage(propdata->varhashmap, var) > 0));

   return getUbIndex((int)(size_t)SCIPhashmapGetImage(propdata->varhashmap, var) - 1);
}

/** reset propagation data */
static
void resetPropdata(
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   propdata->vars = NULL;
   propdata->varhashmap = NULL;
   propdata->topoorder = NULL;
   propdata->vboundboundedidx = NULL;
   propdata->vboundcoefs = NULL;
   propdata->vboundconstants = NULL;
   propdata->nvbounds = NULL;
   propdata->vboundsize = NULL;
   propdata->nbounds = 0;
   propdata->initialized = FALSE;
}

/** catches events for variables */
static
SCIP_RETCODE catchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Bool lower;
   int nbounds;
   int v;
   int idx;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(propdata->vars != NULL);
   assert(propdata->topoorder != NULL);

   /* catch variable events according to computed eventtypes */
   eventhdlr = propdata->eventhdlr;
   assert(eventhdlr != NULL);

   vars = propdata->vars;
   nbounds = propdata->nbounds;

   /* setup events */
   for( v = 0; v < nbounds; ++v )
   {
      idx = propdata->topoorder[v];
      assert(idx >= 0 && idx < nbounds);

      var = vars[getVarIndex(idx)];
      lower = isIndexLowerbound(idx);

      /* if the bound does not influence another bound by implications, cliques, or vbounds,
       * we do not create an event and do not catch changes of the bound;
       * we mark this by setting the value in topoorder to -1
       */
      if( propdata->nvbounds[idx] == 0 && SCIPvarGetNImpls(var, lower) == 0 && SCIPvarGetNCliques(var, lower) == 0 )
      {
         propdata->topoorder[v] = -1;
         continue;
      }

      /* determine eventtype that we want to catch depending on boundtype of variable */
      if( lower )
         eventtype = SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_GLBCHANGED;
      else
         eventtype = SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_GUBCHANGED;

      SCIP_CALL( SCIPcatchVarEvent(scip, var, eventtype, eventhdlr, (SCIP_EVENTDATA*) (uintptr_t) v, NULL) ); /*lint !e571*/
   }

   return SCIP_OKAY;
}

/** drops events for variables */
static
SCIP_RETCODE dropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Bool lower;
   int nbounds;
   int v;
   int idx;

   assert(propdata != NULL);

   eventhdlr = propdata->eventhdlr;
   assert(eventhdlr != NULL);

   vars = propdata->vars;
   nbounds = propdata->nbounds;

   for( v = 0; v < nbounds; ++v )
   {
      idx = propdata->topoorder[v];

      if( idx == -1 )
         continue;

      assert(idx >= 0 && idx < nbounds);

      var = vars[getVarIndex(idx)];
      lower = isIndexLowerbound(idx);

      /* determine eventtype that we catch and now want to drop depending on boundtype of variable */
      if( lower )
         eventtype = SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_GLBCHANGED;
      else
         eventtype = SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_GUBCHANGED;

      SCIP_CALL( SCIPdropVarEvent(scip, var, eventtype, eventhdlr, (SCIP_EVENTDATA*) (uintptr_t) v, -1) ); /*lint !e571*/
   }

   return SCIP_OKAY;
}

#define INITMEMSIZE 5

/* adds a vbound to the propagator data to store it internally and allow forward propagation */
static
SCIP_RETCODE addVbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int                   startidx,           /**< index of bound of variable influencing the other variable */
   int                   endidx,             /**< index of bound of variable which is influenced */
   SCIP_Real             coef,               /**< coefficient in the variable bound */
   SCIP_Real             constant            /**< constant in the variable bound */
   )
{
   int nvbounds;

   assert(scip != NULL);
   assert(propdata != NULL);

   if( propdata->vboundsize[startidx] == 0 )
   {
      /* allocate memory for storing vbounds */
      propdata->vboundsize[startidx] = INITMEMSIZE;

      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundboundedidx[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundcoefs[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundconstants[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
   }
   else if( propdata->nvbounds[startidx] >= propdata->vboundsize[startidx] )
   {
      /* reallocate memory for storing vbounds */
      propdata->vboundsize[startidx] = SCIPcalcMemGrowSize(scip, propdata->nvbounds[startidx] + 1);
      assert(propdata->nvbounds[startidx] < propdata->vboundsize[startidx]);

      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundboundedidx[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundcoefs[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundconstants[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
   }

   nvbounds = propdata->nvbounds[startidx];
   propdata->vboundboundedidx[startidx][nvbounds] = endidx;
   propdata->vboundcoefs[startidx][nvbounds] = coef;
   propdata->vboundconstants[startidx][nvbounds] = constant;
   (propdata->nvbounds[startidx])++;

   return SCIP_OKAY;
}

/** comparison method for two indices in the topoorder array, preferring higher indices because the order is reverse
 *  topological
 */
static
SCIP_DECL_SORTPTRCOMP(compVarboundIndices)
{
   int idx1 = (int)(size_t)elem1;
   int idx2 = (int)(size_t)elem2;

   return idx2 - idx1;
}

/* extract bound changes or infeasibility information from a cycle in the variable bound graph detected during
 * depth-first search
 */
static
SCIP_RETCODE extractCycle(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int*                  dfsstack,           /**< array of size number of nodes to store the stack;
                                              *   only needed for performance reasons */
   int*                  stacknextedge,      /**< array storing the next edge to be visited in dfs for all nodes on the
                                              *   stack/in the current path; negative numbers represent a clique,
                                              *   positive numbers an implication (smaller numbers) or a variable bound */
   int                   stacksize,          /**< current stack size */
   SCIP_Bool             samebound,          /**< does the cycle contain the same bound twice or both bounds of the same variable? */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */

   )
{
   SCIP_VAR** vars;
   SCIP_VAR* currvar;
   SCIP_Bool currlower;

   SCIP_Real coef = 1.0;
   SCIP_Real constant = 0.0;
   SCIP_Bool islower;
   SCIP_Real newbound;
   int cycleidx;
   int startidx;
   int ntmpimpls;
   int j;
   int k;

   assert(scip != NULL);
   assert(propdata != NULL);

   vars = propdata->vars;

   /* the last element on the stack is the end-point of the cycle */
   cycleidx = dfsstack[stacksize - 1];

   /* the same bound of the variable was visited already earlier on the current path, so the start-point of the cycle
    * has the same index
    */
   if( samebound )
      startidx = cycleidx;
   /* the other bound of the variable was visited earlier on the current path, so the start-point of the cycle
    * has the index of the other bound
    */
   else
      startidx = getOtherBoundIndex(cycleidx);

   /* search for the start-point of the cycle; since the endpoint is at position stacksize - 1 we start at stacksize - 2
    * and go backwards
    */
   for( j = stacksize - 2; dfsstack[j] != startidx && j >= 0; --j ){};
   assert(j >= 0);

   for( ; j < stacksize - 1; ++j )
   {
      currvar = vars[getVarIndex(dfsstack[j])];
      currlower = isIndexLowerbound(dfsstack[j]);

      /* if we do not use implications, we assume the number of implications to be 0 (as we did before) */
      ntmpimpls = (propdata->useimplics ? SCIPvarGetNImpls(currvar, currlower) : 0);

      /* stacknextedge[j] <= 0 --> the last outgoing edge traversed during dfs starting from node dfsstack[j] was given
       * by a clique
       */
      if( stacknextedge[j] <= 0 )
      {
         SCIP_Bool nextlower = isIndexLowerbound(dfsstack[j+1]);
#if defined(SCIP_DEBUG) || defined(SCIP_MORE_DEBUG)
         SCIP_CLIQUE** tmpcliques = SCIPvarGetCliques(currvar, currlower);
         SCIP_VAR** cliquevars;
         SCIP_Bool* cliquevals;
         int ntmpcliques = SCIPvarGetNCliques(currvar, currlower);
         int ncliquevars;
         int v;
#endif
         /* there are four cases:
          * a) lb(x) -> ub(y)   ==>   clique(x,y,...)    ==>   y <= 1 - x
          * b) lb(x) -> lb(y)   ==>   clique(x,~y,...)   ==>   y >= x
          * c) ub(x) -> ub(y)   ==>   clique(~x,y,...)   ==>   y <= x
          * d) ub(x) -> lb(y)   ==>   clique(~x,~y,...)  ==>   y >= 1 - x
          *
          * in cases b) and c), coef=1.0 and constant=0.0; these are the cases where both nodes represent
          * the same type of bound
          * in cases a) and d), coef=-1.0 and constant=1.0; both nodes represent different types of bounds
          *
          * we do not need to change the overall coef and constant in cases b) and c), but for the others
          */
         if( currlower != nextlower )
         {
            coef = -coef;
            constant = -constant + 1.0;
         }

         /* since the coefficient and constant only depend on the type of bounds of the two nodes (see below), we do not
          * need to search for the variable in the clique; however, if debug output is enabled, we want to print the
          * clique, if more debugging is enabled, we explicitly check that the variable and bound we expect are in the
          * clique
          */
#if defined(SCIP_DEBUG) || defined(SCIP_MORE_DEBUG)
         if( stacknextedge[j] == 0 )
         {
            k = ntmpcliques - 1;
         }
         else
         {
            /* we processed at least one edge, so the next one should be -2 or smaller (we have a -1 offset) */
            assert(stacknextedge[j] <= -2);

            k = -stacknextedge[j] - 2;

            assert(k < ntmpcliques);
         }

         cliquevars = SCIPcliqueGetVars(tmpcliques[k]);
         cliquevals = SCIPcliqueGetValues(tmpcliques[k]);
         ncliquevars = SCIPcliqueGetNVars(tmpcliques[k]);
#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "clique: ");
         for( v = 0; v < ncliquevars; ++v )
         {
            SCIPdebugMsg(scip, "%s%s ", cliquevals[v] ? "" : "~", SCIPvarGetName(cliquevars[v]));
         }
         SCIPdebugMsg(scip, "(equation:%d)\n", SCIPcliqueIsEquation(tmpcliques[k]));
#endif
#ifdef SCIP_MORE_DEBUG
         for( v = 0; v < ncliquevars; ++v )
         {
            if( cliquevars[v] == vars[getVarIndex(dfsstack[j+1])] && cliquevals[v] == !nextlower )
               break;
         }
         assert(v < ncliquevars);
#endif

         SCIPdebugMsg(scip, "%s(%s) -- (*%g + %g)[clique(<%s%s>,<%s%s>,...)] --> %s(%s)\n",
            indexGetBoundString(dfsstack[j]), SCIPvarGetName(currvar),
            (currlower != nextlower ? -1.0 : 1.0),
            (currlower != nextlower ? 1.0 : 0.0),
            (currlower ? "" : "~"), SCIPvarGetName(currvar),
            (nextlower ? "~" : ""), SCIPvarGetName(vars[getVarIndex(dfsstack[j+1])]),
            indexGetBoundString(dfsstack[j+1]), SCIPvarGetName(currvar));
#endif
      }
      /* stacknextedge[j] > 0 --> the last outgoing edge traversed during dfs starting from node dfsstack[j] was given
       * by an implication or vbound. Implications are looked at first, so if stacknextedge[j] <= ntmpimpls, it comes
       * from an implication
       */
      else if( stacknextedge[j] <= ntmpimpls )
      {
#ifndef NDEBUG
         SCIP_VAR** implvars = SCIPvarGetImplVars(currvar, currlower);
#endif
         SCIP_BOUNDTYPE* impltypes = SCIPvarGetImplTypes(currvar, currlower);
         SCIP_Real* implbounds = SCIPvarGetImplBounds(currvar, currlower);
         SCIP_VAR* nextvar = vars[getVarIndex(dfsstack[j+1])];
         SCIP_Real newconstant;
         SCIP_Real newcoef;

         k = stacknextedge[j] - 1;

         assert(dfsstack[j+1] == (impltypes[k] == SCIP_BOUNDTYPE_LOWER ?
               varGetLbIndex(propdata, implvars[k]) : varGetUbIndex(propdata, implvars[k])));

         if( impltypes[k] == SCIP_BOUNDTYPE_LOWER )
         {
            newcoef = implbounds[k] - SCIPvarGetLbLocal(nextvar);

            if( currlower )
            {
               newconstant = SCIPvarGetLbLocal(nextvar);
            }
            else
            {
               newconstant = implbounds[k];
               newcoef *= -1.0;
            }
         }
         else
         {
            assert(impltypes[k] == SCIP_BOUNDTYPE_UPPER);

            newcoef = SCIPvarGetUbLocal(nextvar) - implbounds[k];

            if( currlower )
            {
               newconstant = SCIPvarGetUbLocal(nextvar);
               newcoef *= -1.0;
            }
            else
            {
               newconstant = implbounds[k];
            }
         }

         coef = coef * newcoef;
         constant = constant * newcoef + newconstant;

         SCIPdebugMsg(scip, "%s(%s) -- (*%g + %g) --> %s(%s)\n",
            indexGetBoundString(dfsstack[j]), SCIPvarGetName(vars[getVarIndex(dfsstack[j])]),
            newcoef, newconstant,
            indexGetBoundString(dfsstack[j+1]), SCIPvarGetName(vars[getVarIndex(dfsstack[j+1])]));
      }
      /* the edge was given by a variable bound relation */
      else
      {
         assert(stacknextedge[j] > ntmpimpls);

         k = stacknextedge[j] - ntmpimpls - 1;
         assert(k < propdata->nvbounds[dfsstack[j]]);
         assert(propdata->vboundboundedidx[dfsstack[j]][k] == dfsstack[j+1]);

         SCIPdebugMsg(scip, "%s(%s) -- (*%g + %g) --> %s(%s)\n",
            indexGetBoundString(dfsstack[j]), SCIPvarGetName(vars[getVarIndex(dfsstack[j])]),
            propdata->vboundcoefs[dfsstack[j]][k], propdata->vboundconstants[dfsstack[j]][k],
            indexGetBoundString(dfsstack[j+1]), SCIPvarGetName(vars[getVarIndex(dfsstack[j+1])]));

         coef = coef * propdata->vboundcoefs[dfsstack[j]][k];
         constant = constant * propdata->vboundcoefs[dfsstack[j]][k] + propdata->vboundconstants[dfsstack[j]][k];
      }
   }

   SCIPdebugMsg(scip, "==> %s(%s) -- (*%g + %g) --> %s(%s)\n",
      indexGetBoundString(startidx), SCIPvarGetName(vars[getVarIndex(startidx)]),
      coef, constant,
      indexGetBoundString(cycleidx), SCIPvarGetName(vars[getVarIndex(cycleidx)]));

   islower = isIndexLowerbound(cycleidx);

   /* we have a relation x <=/>= coef * x + constant now
    * (the relation depends on islower, i.e., whether the last node in the cycle is a lower or upper bound)
    * case 1) coef is 1.0 --> x cancels out and we have a statement 0 <=/>= constant.
    *         if we have a >= relation and constant is positive, we have a contradiction 0 >= constant
    *         if we have a <= relation and constant is negative, we have a contradiction 0 <= constant
    * case 2) coef != 1.0 --> we have a relation x - coef * x <=/>= constant
    *                                      <=> (1 - coef) * x <=/>= constant
    *         if coef < 1.0 this gives us x >= constant / (1 - coef) (if islower=TRUE)
    *                                  or x <= constant / (1 - coef) (if islower=FALSE)
    *         if coef > 1.0, the relation signs need to be switched.
    */
   if( SCIPisEQ(scip, coef, 1.0) )
   {
      if( islower && SCIPisFeasPositive(scip, constant) )
      {
         SCIPdebugMsg(scip, "-> infeasible aggregated variable bound relation 0 >= %g\n", constant);
         *infeasible = TRUE;
      }
      else if( !islower && SCIPisFeasNegative(scip, constant) )
      {
         SCIPdebugMsg(scip, "-> infeasible aggregated variable bound relation 0 <= %g\n", constant);
         *infeasible = TRUE;
      }
   }
   else
   {
      SCIP_Bool tightened;

      newbound = constant / (1.0 - coef);

      if( SCIPisGT(scip, coef, 1.0) )
         islower = !islower;

      if( islower )
      {
         SCIPdebugMsg(scip, "-> found new lower bound: <%s>[%g,%g] >= %g\n", SCIPvarGetName(vars[getVarIndex(cycleidx)]),
            SCIPvarGetLbLocal(vars[getVarIndex(cycleidx)]), SCIPvarGetUbLocal(vars[getVarIndex(cycleidx)]), newbound);
         SCIP_CALL( SCIPtightenVarLb(scip, vars[getVarIndex(cycleidx)], newbound, FALSE, infeasible, &tightened) );
      }
      else
      {
         SCIPdebugMsg(scip, "-> found new upper bound: <%s>[%g,%g] <= %g\n", SCIPvarGetName(vars[getVarIndex(cycleidx)]),
            SCIPvarGetLbLocal(vars[getVarIndex(cycleidx)]), SCIPvarGetUbLocal(vars[getVarIndex(cycleidx)]), newbound);
         SCIP_CALL( SCIPtightenVarUb(scip, vars[getVarIndex(cycleidx)], newbound, FALSE, infeasible, &tightened) );
      }

      if( tightened )
         SCIPdebugMsg(scip, "---> applied new bound\n");
   }

   return SCIP_OKAY;
}

#define VISITED 1
#define ACTIVE  2
/** performs depth-first-search in the implicitly given directed graph from the given start index */
static
SCIP_RETCODE dfs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int                   startnode,          /**< node to start the depth-first-search */
   int*                  visited,            /**< array to store for each node, whether it was already visited */
   int*                  dfsstack,           /**< array of size number of nodes to store the stack;
                                              *   only needed for performance reasons */
   int*                  stacknextedge,      /**< array of size number of nodes to store the next edge to be visited in
                                              *   dfs for all nodes on the stack/in the current path; only needed for
                                              *   performance reasons */
   int*                  dfsnodes,           /**< array of nodes that can be reached starting at startnode, in reverse
                                              *   dfs order */
   int*                  ndfsnodes,          /**< pointer to store number of nodes that can be reached starting at
                                              *   startnode */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* startvar;
   SCIP_Bool lower;
   int stacksize;
   int curridx;
   int nimpls;
   int idx;
   /* for cycle detection, we need to mark currently active nodes, otherwise we just mark them as visited */
   int visitedflag = (propdata->detectcycles ? ACTIVE : VISITED);

   assert(startnode >= 0);
   assert(startnode < propdata->nbounds);
   assert(visited != NULL);
   assert(visited[startnode] == 0);
   assert(dfsstack != NULL);
   assert(dfsnodes != NULL);
   assert(ndfsnodes != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   vars = propdata->vars;

   /* put start node on the stack */
   dfsstack[0] = startnode;
   stacknextedge[0] = 0;
   stacksize = 1;
   idx = -1;

   /* we run until no more bounds indices are on the stack, i.e. all changed bounds were propagated */
   while( stacksize > 0 )
   {
      /* get next node from stack */
      curridx = dfsstack[stacksize - 1];

      /* mark current node as visited */
      assert((visited[curridx] != 0) == (stacknextedge[stacksize - 1] != 0));
      visited[curridx] = visitedflag;

      startvar = vars[getVarIndex(curridx)];
      lower = isIndexLowerbound(curridx);

      /* if the variable was fixed in the meantime, it is a loose end in the variable bound graph */
      if( SCIPisFeasGE(scip, SCIPvarGetLbGlobal(startvar), SCIPvarGetUbGlobal(startvar)) )
         goto REMOVE;

      nimpls = 0;

      if( propdata->sortcliques && propdata->usecliques && stacknextedge[stacksize - 1] == 0 )
         stacknextedge[stacksize - 1] = -1;

      /* stacknextedge is negative, if the last visited edge from the current node belongs to a clique;
       * the index of the clique in the variable's clique list equals abs(stacknextedge) - 1
       */
      if( propdata->sortcliques && propdata->usecliques && stacknextedge[stacksize - 1] < 0 )
      {
         SCIP_CLIQUE** cliques;
         int ncliques;
         int j;
         int i;
         SCIP_Bool found;

         ncliques = SCIPvarGetNCliques(startvar, lower);
         cliques = SCIPvarGetCliques(startvar, lower);
         found = FALSE;

         assert(stacknextedge[stacksize - 1] == -1 || -stacknextedge[stacksize - 1] - 1 <= ncliques);

         /* iterate over all not yet handled cliques and search for an unvisited node */
         for( j = -stacknextedge[stacksize - 1] - 1; j < ncliques; ++j )
         {
            SCIP_VAR** cliquevars;
            SCIP_Bool* cliquevals;
            int ncliquevars;

            cliquevars = SCIPcliqueGetVars(cliques[j]);
            cliquevals = SCIPcliqueGetValues(cliques[j]);
            ncliquevars = SCIPcliqueGetNVars(cliques[j]);

            for( i = 0; i < ncliquevars; ++i )
            {
               if( cliquevars[i] == startvar )
                  continue;

               if( cliquevals[i] )
                  idx = varGetUbIndex(propdata, cliquevars[i]);
               else
                  idx = varGetLbIndex(propdata, cliquevars[i]);

               /* we reached a variable that was already visited on the active path, so we have a cycle in the variable
                * bound graph and try to extract valid bound changes from it or detect infeasibility
                */
               if( idx >= 0 && (visited[idx] == ACTIVE || visited[getOtherBoundIndex(idx)] == ACTIVE)
                  && !SCIPisFeasGE(scip, SCIPvarGetLbGlobal(cliquevars[i]), SCIPvarGetUbGlobal(cliquevars[i])) )
               {
                  SCIPdebugMsg(scip, "found cycle\n");

                  dfsstack[stacksize] = idx;
                  stacknextedge[stacksize - 1] = -j - 2;

                  SCIP_CALL( extractCycle(scip, propdata, dfsstack, stacknextedge, stacksize + 1,
                        visited[idx] == ACTIVE, infeasible) );

                  if( *infeasible )
                     break;
               }

               /* break when the first unvisited node is reached */
               if( idx >= 0 && !visited[idx] )
               {
                  found = TRUE;
                  break;
               }
            }
            if( found )
               break;
         }

         /* we stopped because we found an unhandled node and not because we reached the end of the list */
         if( found )
         {
            assert(idx >= 0);
            assert(!visited[idx]);
            assert(j < ncliques);

            SCIPdebugMsg(scip, "clique: %s(%s) -> %s(%s)\n", getBoundString(lower), SCIPvarGetName(startvar),
               indexGetBoundString(idx), SCIPvarGetName(vars[getVarIndex(idx)]));

            /* put the adjacent node onto the stack */
            dfsstack[stacksize] = idx;
            stacknextedge[stacksize] = 0;
            stacknextedge[stacksize - 1] = -j - 2;
            stacksize++;
            assert(stacksize <= propdata->nbounds);

            /* restart while loop, get next index from stack */
            continue;
         }
         else
         {
            /* we did not find an edge to an unhandled node given by a clique */
            stacknextedge[stacksize - 1] = 0;
         }
      }
      assert(stacknextedge[stacksize - 1] >= 0);

      /* go over edges given by implications */
      if( propdata->useimplics )
      {
         nimpls = SCIPvarGetNImpls(startvar, lower);

         if( stacknextedge[stacksize - 1] < nimpls )
         {
            SCIP_VAR** implvars;
            SCIP_BOUNDTYPE* impltypes;
            int* implids;
            int i;

            implvars = SCIPvarGetImplVars(startvar, lower);
            impltypes = SCIPvarGetImplTypes(startvar, lower);
            implids = SCIPvarGetImplIds(startvar, lower);

            for( i = stacknextedge[stacksize - 1]; i < nimpls; ++i )
            {
               /* it might happen that implications point to inactive variables (normally, those are removed when a
                * variable becomes inactive, but in some cases, it cannot be done), we have to ignore these variables
                */
               if( !SCIPvarIsActive(implvars[i]) )
                  continue;

               /* implication is just a shortcut, so we dont regard it now, because will later go the long way, anyway;
                * however, if we do regard cliques for the topological order, we use them to get a better order
                */
               if( propdata->usecliques && !propdata->sortcliques && implids[i] < 0 )
                  continue;

               idx = (impltypes[i] == SCIP_BOUNDTYPE_LOWER ?
                  varGetLbIndex(propdata, implvars[i]) : varGetUbIndex(propdata, implvars[i]));

               /* we reached a variable that was already visited on the active path, so we have a cycle in the variable
                * bound graph and try to extract valid bound changes from it or detect infeasibility
                */
               if( idx >= 0 && (visited[idx] == ACTIVE || visited[getOtherBoundIndex(idx)] == ACTIVE)
                  && !SCIPisFeasGE(scip, SCIPvarGetLbGlobal(implvars[i]), SCIPvarGetUbGlobal(implvars[i])) )
               {
                  SCIPdebugMsg(scip, "found cycle\n");

                  dfsstack[stacksize] = idx;
                  stacknextedge[stacksize - 1] = i + 1;

                  SCIP_CALL( extractCycle(scip, propdata, dfsstack, stacknextedge, stacksize + 1,
                        visited[idx] == ACTIVE, infeasible) );

                  if( *infeasible )
                     break;
               }

               /* break when the first unvisited node is reached */
               if( idx >= 0 && !visited[idx] )
                  break;
            }

            /* we stopped because we found an unhandled node and not because we reached the end of the list */
            if( i < nimpls )
            {
               assert(!visited[idx]);

               SCIPdebugMsg(scip, "impl: %s(%s) -> %s(%s)\n", getBoundString(lower), SCIPvarGetName(startvar),
                  indexGetBoundString(idx), SCIPvarGetName(vars[getVarIndex(idx)]));

               /* put the adjacent node onto the stack */
               dfsstack[stacksize] = idx;
               stacknextedge[stacksize] = 0;
               stacknextedge[stacksize - 1] = i + 1;
               stacksize++;
               assert(stacksize <= propdata->nbounds);

               /* restart while loop, get next index from stack */
               continue;
            }
            else
            {
               stacknextedge[stacksize - 1] = nimpls;
            }
         }
      }
      assert(stacknextedge[stacksize - 1] >= nimpls);

      /* go over edges corresponding to varbounds */
      if( propdata->usevbounds )
      {
         int nvbounds;
         int* vboundidx;
         int i;

         nvbounds = propdata->nvbounds[curridx];
         vboundidx = propdata->vboundboundedidx[curridx];

         /* iterate over all vbounds for the given bound */
         for( i = stacknextedge[stacksize - 1] - nimpls; i < nvbounds; ++i )
         {
            idx = vboundidx[i];
            assert(idx >= 0);

            if( (visited[idx] == ACTIVE || visited[getOtherBoundIndex(idx)] == ACTIVE)
               && !SCIPisFeasGE(scip, SCIPvarGetLbGlobal(vars[getVarIndex(idx)]), SCIPvarGetUbGlobal(vars[getVarIndex(idx)])) )
            {
               SCIPdebugMsg(scip, "found cycle\n");

               dfsstack[stacksize] = idx;
               stacknextedge[stacksize - 1] = nimpls + i + 1;

               /* we reached a variable that was already visited on the active path, so we have a cycle in the variable
                * bound graph and try to extract valid bound changes from it or detect infeasibility
                */
               SCIP_CALL( extractCycle(scip, propdata, dfsstack, stacknextedge, stacksize + 1,
                     visited[idx] == ACTIVE, infeasible) );

               if( *infeasible )
                  break;
            }

            /* break when the first unvisited node is reached */
            if( !visited[idx] )
               break;
         }

         if( *infeasible )
            break;

         /* we stopped because we found an unhandled node and not because we reached the end of the list */
         if( i < nvbounds )
         {
            assert(!visited[idx]);

            SCIPdebugMsg(scip, "vbound: %s(%s) -> %s(%s)\n", getBoundString(lower), SCIPvarGetName(startvar),
               indexGetBoundString(idx), SCIPvarGetName(vars[getVarIndex(idx)]));

            /* put the adjacent node onto the stack */
            dfsstack[stacksize] = idx;
            stacknextedge[stacksize] = 0;
            stacknextedge[stacksize - 1] = nimpls + i + 1;
            stacksize++;
            assert(stacksize <= propdata->nbounds);

            /* restart while loop, get next index from stack */
            continue;
         }

      }
   REMOVE:
      /* the current node was completely handled, remove it from stack */
      stacksize--;

      SCIPdebugMsg(scip, "topoorder[%d] = %s(%s)\n", *ndfsnodes, getBoundString(lower), SCIPvarGetName(startvar));

      /* store node in the sorted nodes array */
      dfsnodes[(*ndfsnodes)] = curridx;
      assert(visited[curridx] == visitedflag);
      visited[curridx] = VISITED;
      (*ndfsnodes)++;
   }

   return SCIP_OKAY;
}

/** sort the bounds of variables topologically */
static
SCIP_RETCODE topologicalSort(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   )
{
   int* dfsstack;
   int* stacknextedge;
   int* visited;
   int nsortednodes;
   int nbounds;
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(infeasible != NULL);

   nbounds = propdata->nbounds;

   SCIP_CALL( SCIPallocBufferArray(scip, &dfsstack, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stacknextedge, nbounds) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &visited, nbounds) );

   nsortednodes = 0;

   /* while there are unvisited nodes, run dfs starting from one of these nodes; the dfs orders are stored in the
    * topoorder array, later dfs calls are just appended after the stacks of previous dfs calls, which gives us a
    * reverse topological order
    */
   for( i = 0; i < nbounds && !(*infeasible); ++i )
   {
      if( !visited[i] )
      {
         SCIP_CALL( dfs(scip, propdata, i, visited, dfsstack, stacknextedge, propdata->topoorder, &nsortednodes, infeasible) );
      }
   }
   assert((nsortednodes == nbounds) || (*infeasible));

   SCIPfreeBufferArray(scip, &visited);
   SCIPfreeBufferArray(scip, &stacknextedge);
   SCIPfreeBufferArray(scip, &dfsstack);

   return SCIP_OKAY;
}

/** initializes the internal data for the variable bounds propagator */
static
SCIP_RETCODE initData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< vbounds propagator */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   int nvars;
   int nbounds;
   int startidx;
   int v;
   int n;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(infeasible != NULL);

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(!propdata->initialized);

   SCIPdebugMsg(scip, "initialize vbounds propagator for problem <%s>\n", SCIPgetProbName(scip));

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nbounds = 2 * nvars;

   *infeasible = FALSE;

   /* store size of the bounds of variables array */
   propdata->nbounds = nbounds;

   if( nbounds == 0 )
      return SCIP_OKAY;

   propdata->initialized = TRUE;

   /* prepare priority queue structure */
   SCIP_CALL( SCIPpqueueCreate(&propdata->propqueue, nvars, 2.0, compVarboundIndices) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->inqueue, nbounds) );
   BMSclearMemoryArray(propdata->inqueue, nbounds);

   /* we need to copy the variable since this array is the basis of the propagator and the corresponding variable array
    * within SCIP might change during the search
    */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &propdata->vars, vars, nvars) );
   SCIP_CALL( SCIPhashmapCreate(&propdata->varhashmap, SCIPblkmem(scip), nvars) );

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPhashmapInsert(propdata->varhashmap, propdata->vars[v], (void*)(size_t)(v + 1)) );
   }

   /* allocate memory for the arrays of the propdata */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->topoorder, nbounds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->vboundboundedidx, nbounds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->vboundcoefs, nbounds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->vboundconstants, nbounds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->nvbounds, nbounds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->vboundsize, nbounds) );
   BMSclearMemoryArray(propdata->vboundboundedidx, nbounds);
   BMSclearMemoryArray(propdata->vboundcoefs, nbounds);
   BMSclearMemoryArray(propdata->vboundconstants, nbounds);
   BMSclearMemoryArray(propdata->nvbounds, nbounds);
   BMSclearMemoryArray(propdata->vboundsize, nbounds);

   for( v = 0; v < nbounds; ++v )
   {
      propdata->topoorder[v] = v;
      propdata->vboundboundedidx[v] = NULL;
      propdata->vboundcoefs[v] = NULL;
      propdata->vboundconstants[v] = NULL;
      propdata->nvbounds[v] = 0;
      propdata->vboundsize[v] = 0;
   }

   /* collect information about varbounds */
   for( v = 0; v < nbounds; ++v )
   {
      SCIP_VAR** vbvars;
      SCIP_VAR* var;
      SCIP_Real* coefs;
      SCIP_Real* constants;
      SCIP_Bool lower;
      int nvbvars;

      var = vars[getVarIndex(v)];
      lower = isIndexLowerbound(v);

      /* get the variable bound informations for the current variable */
      if( lower )
      {
         vbvars = SCIPvarGetVlbVars(var);
         coefs = SCIPvarGetVlbCoefs(var);
         constants = SCIPvarGetVlbConstants(var);
         nvbvars = SCIPvarGetNVlbs(var);
      }
      else
      {
         vbvars = SCIPvarGetVubVars(var);
         coefs = SCIPvarGetVubCoefs(var);
         constants = SCIPvarGetVubConstants(var);
         nvbvars = SCIPvarGetNVubs(var);
      }

      /* loop over all variable bounds; a variable lower bound has the form: x >= b*y + d,
       * a variable upper bound the form x <= b*y + d */
      for( n = 0; n < nvbvars; ++n )
      {
         SCIP_VAR* vbvar;
         SCIP_Real coef;
         SCIP_Real constant;

         vbvar = vbvars[n];
         coef = coefs[n];
         constant = constants[n];
         assert(vbvar != NULL);

         /* transform variable bound variable to an active variable, if possible */
         SCIP_CALL( SCIPgetProbvarSum(scip, &vbvar, &coef, &constant) );
         assert(vbvar != NULL);

         if( !SCIPvarIsActive(vbvar) )
            continue;

         /* if the coefficient is positive, the type of bound is the same for the bounded and the bounding variable */
         if( SCIPisPositive(scip, coef) )
            startidx = (lower ? varGetLbIndex(propdata, vbvar) : varGetUbIndex(propdata, vbvar));
         else
            startidx = (lower ? varGetUbIndex(propdata, vbvar) : varGetLbIndex(propdata, vbvar));
         assert(startidx >= 0);

         /* If the vbvar is binary, the vbound should be stored as an implication already.
          * However, it might happen that vbvar was integer when the variable bound was added, but was converted
          * to a binary variable later during presolving when its upper bound was changed to 1. In this case,
          * the implication might not have been created.
          */
         if( SCIPvarGetType(vbvar) == SCIP_VARTYPE_BINARY
            && SCIPvarHasImplic(vbvar, isIndexLowerbound(startidx), var, getBoundtype(v)) )
         {
            SCIPdebugMsg(scip, "varbound <%s> %s %g * <%s> + %g not added to propagator data due to reverse implication\n",
               SCIPvarGetName(var), (lower ? ">=" : "<="), coef,
               SCIPvarGetName(vbvar), constant);
         }
         else
         {
            SCIP_CALL( addVbound(scip, propdata, startidx, v, coef, constant) );

            SCIPdebugMsg(scip, "varbound <%s> %s %g * <%s> + %g added to propagator data\n",
               SCIPvarGetName(var), (lower ? ">=" : "<="), coef,
               SCIPvarGetName(vbvar), constant);

         }
      }
   }

   /* sort the bounds topologically */
   if( propdata->dotoposort )
   {
      SCIP_CALL( topologicalSort(scip, propdata, infeasible) );
   }

   /* catch variable events */
   SCIP_CALL( catchEvents(scip, propdata) );

   return SCIP_OKAY;
}

/** resolves a propagation by adding the variable which implied that bound change */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable to be reported */
   SCIP_BOUNDTYPE        boundtype,          /**< bound to be reported */
   SCIP_BDCHGIDX*        bdchgidx            /**< the index of the bound change, representing the point of time where
                                              *   the change took place, or NULL for the current local bounds */
   )
{
   assert(propdata != NULL);
   assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

   SCIPdebugMsg(scip, " -> add %s bound of variable <%s> as reason\n",
      getBoundtypeString(boundtype), SCIPvarGetName(var));

   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
      break;
   case SCIP_BOUNDTYPE_UPPER:
      SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
      break;
   default:
      SCIPerrorMessage("invalid bound type <%d>\n", boundtype);
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   return SCIP_OKAY;
}

/** relaxes bound of give variable as long as the given inference bound still leads to a cutoff and add that bound
 *  change to the conflict set
 */
static
SCIP_RETCODE relaxVbdvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which the upper bound should be relaxed */
   SCIP_BOUNDTYPE        boundtype,          /**< boundtype used for the variable bound variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where
                                              *   the change took place, or NULL for the current local bounds */
   SCIP_Real             relaxedbd           /**< relaxed bound */
   )
{
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_CALL( SCIPaddConflictRelaxedLb(scip, var, bdchgidx, relaxedbd) );
   }
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      SCIP_CALL( SCIPaddConflictRelaxedUb(scip, var, bdchgidx, relaxedbd) );
   }

   return SCIP_OKAY;
}

/** compute the relaxed bound which is sufficient to propagate the inference lower bound of given variable */
static
SCIP_Real computeRelaxedLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which was propagated */
   SCIP_Real             inferlb,            /**< inference lower bound */
   SCIP_Real             coef,               /**< inference variable bound coefficient used */
   SCIP_Real             constant            /**< inference variable bound constant used */
   )
{
   SCIP_Real relaxedbd;

   if( SCIPvarIsIntegral(var) && inferlb < SCIPgetHugeValue(scip) * SCIPfeastol(scip) )
      relaxedbd = (inferlb - 1.0 + 2*SCIPfeastol(scip) - constant) / coef;
   else
      relaxedbd = (inferlb - constant) / coef;

   /* check the computed relaxed lower/upper bound is a proper reason for the inference bound which has to be explained */
   assert(SCIPisEQ(scip, inferlb, SCIPadjustedVarLb(scip, var, relaxedbd * coef + constant)));

   if( coef > 0.0 )
      relaxedbd += SCIPfeastol(scip);
   else
      relaxedbd -= SCIPfeastol(scip);

   return relaxedbd;
}

/** analyzes an infeasibility which was reached by updating the lower bound of the inference variable above its upper
 *  bound
 */
static
SCIP_RETCODE analyzeConflictLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             infervar,           /**< variable which led to a cutoff */
   SCIP_Real             inferlb,            /**< lower bound which led to infeasibility */
   SCIP_VAR*             vbdvar,             /**< variable which is the reason for the lower bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< bound which is the reason for the lower bound change */
   SCIP_Real             coef,               /**< inference variable bound coefficient used */
   SCIP_Real             constant,           /**< inference variable bound constant used */
   SCIP_Bool             canwide             /**< can bound widening be used (for vbounds) or not
                                              *   (for implications or cliques) */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);
   assert(infervar != NULL);
   assert(SCIPisEQ(scip, SCIPvarGetUbLocal(infervar), SCIPgetVarUbAtIndex(scip, infervar, NULL, FALSE)));
   assert(SCIPisEQ(scip, SCIPgetVarUbAtIndex(scip, infervar, NULL, TRUE), SCIPgetVarUbAtIndex(scip, infervar, NULL, FALSE)));
   assert(SCIPisGT(scip, inferlb, SCIPvarGetUbLocal(infervar)));
   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   /* check if conflict analysis is applicable */
   if( !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   if( canwide && propdata->usebdwidening )
   {
      SCIP_Real relaxedbd;
      SCIP_Real relaxedub;

      SCIPdebugMsg(scip, "try to create conflict using bound widening order: inference variable, variable bound variable\n");

      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

      /* adjust lower bound */
      inferlb = SCIPadjustedVarLb(scip, infervar, inferlb);

      /* compute a relaxed upper bound which would be sufficient to be still infeasible */
      if( SCIPvarIsIntegral(infervar) )
         relaxedub = inferlb - 1.0;
      else
         relaxedub = inferlb - 2*SCIPfeastol(scip);

      /* try to relax inference variable upper bound such that the infeasibility is still given */
      SCIP_CALL( SCIPaddConflictRelaxedUb(scip, infervar, NULL, relaxedub) );

      /* collect the upper bound which is reported to the conflict analysis */
      relaxedub = SCIPgetConflictVarUb(scip, infervar);

      /* adjust inference bound with respect to the upper bound reported to the conflict analysis */
      if( SCIPvarIsIntegral(infervar) )
         relaxedub = relaxedub + 1.0;
      else
         relaxedub = relaxedub + 2*SCIPfeastol(scip);

      /* compute the relaxed bound which is sufficient to propagate the inference lower bound of given variable */
      relaxedbd = computeRelaxedLowerbound(scip, infervar, relaxedub, coef, constant);

      /* try to relax variable bound variable */
      SCIP_CALL( relaxVbdvar(scip, vbdvar, boundtype, NULL, relaxedbd) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
   }
   else
   {
      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

      /* add upper bound of the variable for which we tried to change the lower bound */
      SCIP_CALL( SCIPaddConflictUb(scip, infervar, NULL) );

      /* add (correct) bound of the variable which let to the new lower bound */
      SCIP_CALL( resolvePropagation(scip, propdata, vbdvar, boundtype, NULL) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
   }

   return SCIP_OKAY;
}

/** compute the relaxed bound which is sufficient to propagate the inference upper bound of given variable */
static
SCIP_Real computeRelaxedUpperbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which was propagated */
   SCIP_Real             inferub,            /**< inference upper bound */
   SCIP_Real             coef,               /**< inference variable bound coefficient used */
   SCIP_Real             constant            /**< inference variable bound constant used */
   )
{
   SCIP_Real relaxedbd;

   if( SCIPvarIsIntegral(var) && inferub < SCIPgetHugeValue(scip) * SCIPfeastol(scip) )
      relaxedbd = (inferub + 1.0 - 2*SCIPfeastol(scip) - constant) / coef;
   else
      relaxedbd = (inferub - constant) / coef;

   /* check the computed relaxed lower/upper bound is a proper reason for the inference bound which has to be explained */
   assert(SCIPisEQ(scip, inferub, SCIPadjustedVarUb(scip, var, relaxedbd * coef + constant)));

   if( coef > 0.0 )
      relaxedbd -= SCIPfeastol(scip);
   else
      relaxedbd += SCIPfeastol(scip);

   return relaxedbd;
}

/** analyzes an infeasibility which was reached by updating the upper bound of the inference variable below its lower
 *  bound
 */
static
SCIP_RETCODE analyzeConflictUpperbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             infervar,           /**< variable which led to a cutoff */
   SCIP_Real             inferub,            /**< upper bound which led to infeasibility */
   SCIP_VAR*             vbdvar,             /**< variable which is the reason for the upper bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< bound which is the reason for the upper bound change */
   SCIP_Real             coef,               /**< inference variable bound coefficient used */
   SCIP_Real             constant,           /**< inference variable bound constant used */
   SCIP_Bool             canwide             /**< can bound widening be used (for vbounds) or not (for inplications or cliques) */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);
   assert(infervar != NULL);
   assert(SCIPisEQ(scip, SCIPvarGetLbLocal(infervar), SCIPgetVarLbAtIndex(scip, infervar, NULL, FALSE)));
   assert(SCIPisEQ(scip, SCIPgetVarLbAtIndex(scip, infervar, NULL, TRUE), SCIPgetVarLbAtIndex(scip, infervar, NULL, FALSE)));
   assert(SCIPisLT(scip, inferub, SCIPvarGetLbLocal(infervar)));
   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   /* check if conflict analysis is applicable */
   if( !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   if( canwide && propdata->usebdwidening )
   {
      SCIP_Real relaxedbd;
      SCIP_Real relaxedlb;

      SCIPdebugMsg(scip, "try to create conflict using bound widening order: inference variable, variable bound variable\n");

      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

      /* adjust upper bound */
      inferub = SCIPadjustedVarUb(scip, infervar, inferub);

      /* compute a relaxed lower bound which would be sufficient to be still infeasible */
      if( SCIPvarIsIntegral(infervar) )
         relaxedlb = inferub + 1.0;
      else
         relaxedlb = inferub + 2*SCIPfeastol(scip);

      /* try to relax inference variable lower bound such that the infeasibility is still given */
      SCIP_CALL( SCIPaddConflictRelaxedLb(scip, infervar, NULL, relaxedlb) );

      /* collect the lower bound which is reported to the conflict analysis */
      relaxedlb = SCIPgetConflictVarLb(scip, infervar);

      /* adjust inference bound with respect to the upper bound reported to the conflict analysis */
      if( SCIPvarIsIntegral(infervar) )
         relaxedlb = relaxedlb - 1.0;
      else
         relaxedlb = relaxedlb - 2*SCIPfeastol(scip);

      /* compute the relaxed bound which is sufficient to propagate the inference upper bound of given variable */
      relaxedbd = computeRelaxedUpperbound(scip, infervar, relaxedlb, coef, constant);

      /* try to relax variable bound variable */
      SCIP_CALL( relaxVbdvar(scip, vbdvar, boundtype, NULL, relaxedbd) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
   }
   else
   {
      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

      /* add lower bound of the variable for which we tried to change the upper bound */
      SCIP_CALL( SCIPaddConflictLb(scip, infervar, NULL) );

      /* add (correct) bound of the variable which let to the new upper bound */
      SCIP_CALL( resolvePropagation(scip, propdata, vbdvar, boundtype, NULL) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
   }

   return SCIP_OKAY;
}

/* tries to tighten the (global) lower bound of the given variable to the given new bound */
static
SCIP_RETCODE tightenVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< vbounds propagator */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable whose lower bound should be tightened */
   SCIP_Real             newlb,              /**< new lower bound for the variable */
   SCIP_Bool             global,             /**< is the bound globally valid? */
   SCIP_VAR*             vbdvar,             /**< variable which is the reason for the lower bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< bound which is the reason for the lower bound change */
   SCIP_Bool             force,              /**< should domain changes for continuous variables be forced */
   SCIP_Real             coef,               /**< coefficient in vbound constraint causing the propagation;
                                              *   or 0.0 if propagation is caused by clique or implication */
   SCIP_Real             constant,           /**< constant in vbound constraint causing the propagation;
                                              *   or 0.0 if propagation is caused by clique or implication */
   SCIP_Bool             canwide,            /**< can bound widening be used (for vbounds) or not (for inplications or cliques) */
   int*                  nchgbds,            /**< pointer to increase, if a bound was changed */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation */
   )
{
   INFERINFO inferinfo;
   SCIP_Real lb;
   SCIP_Bool tightened;
   SCIP_Bool infeasible;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(propdata != NULL);
   assert(var != NULL);
   assert(nchgbds != NULL);
   assert(result != NULL);

   lb = SCIPvarGetLbLocal(var);

   /* check that the new upper bound is better */
   if( (SCIPvarIsIntegral(var) && newlb - lb > 0.5) || (force && SCIPisGT(scip, newlb, lb)) )
      force = TRUE;
   else
      force = FALSE;

   /* try to tighten the lower bound */
   if( global )
   {
      SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newlb, force, &infeasible, &tightened) );
   }
   else
   {
      inferinfo = getInferInfo(boundtype == SCIP_BOUNDTYPE_LOWER ? varGetLbIndex(propdata, vbdvar) : varGetUbIndex(propdata, vbdvar), boundtype);

      SCIP_CALL( SCIPinferVarLbProp(scip, var, newlb, prop, inferInfoToInt(inferinfo), force, &infeasible, &tightened) );
   }

   if( infeasible )
   {
      /* the infeasible results comes from the fact that the new lower bound lies above the current upper bound */
      assert(SCIPisGT(scip, newlb, SCIPvarGetUbLocal(var)));
      assert(!global || SCIPisGT(scip, newlb, SCIPvarGetUbGlobal(var)));

      SCIPdebugMsg(scip, "tightening%s lower bound of variable <%s> to %g due the %s bound of variable <%s> led to infeasibility\n",
         (global ? " global" : ""), SCIPvarGetName(var), newlb, getBoundtypeString(boundtype), SCIPvarGetName(vbdvar));

      if( global )
      {
         /* cutoff the root node */
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
      }
      else
      {
         /* analyzes a infeasibility via conflict analysis */
         SCIP_CALL( analyzeConflictLowerbound(scip, propdata, var, newlb, vbdvar, boundtype, coef, constant, canwide) );
      }
      *result = SCIP_CUTOFF;
   }
   else if( tightened )
   {
      SCIPdebugMsg(scip, "tightened%s lower bound of variable <%s> to %g due the %s bound of variable <%s>\n",
         (global ? " global" : ""), SCIPvarGetName(var), newlb, getBoundtypeString(boundtype), SCIPvarGetName(vbdvar));
      (*nchgbds)++;
   }

   return SCIP_OKAY;
}

/* tries to tighten the (global) upper bound of the given variable to the given new bound */
static
SCIP_RETCODE tightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< vbounds propagator */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable whose upper bound should be tightened */
   SCIP_Real             newub,              /**< new upper bound of the variable */
   SCIP_Bool             global,             /**< is the bound globally valid? */
   SCIP_VAR*             vbdvar,             /**< variable which is the reason for the upper bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< bound which is the reason for the upper bound change */
   SCIP_Bool             force,              /**< should domain changes for continuous variables be forced */
   SCIP_Real             coef,               /**< coefficient in vbound constraint causing the propagation;
                                              *   or 0.0 if propagation is caused by clique or implication */
   SCIP_Real             constant,           /**< constant in vbound constraint causing the propagation;
                                              *   or 0.0 if propagation is caused by clique or implication */
   SCIP_Bool             canwide,            /**< can bound widening be used (for vbounds) or not (for inplications or cliques) */
   int*                  nchgbds,            /**< pointer to increase, if a bound was changed */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation */
   )
{
   INFERINFO inferinfo;
   SCIP_Real ub;
   SCIP_Bool tightened;
   SCIP_Bool infeasible;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(propdata != NULL);
   assert(var != NULL);
   assert(nchgbds != NULL);
   assert(result != NULL);

   ub = SCIPvarGetUbLocal(var);

   /* check that the new upper bound is better */
   if( (SCIPvarIsIntegral(var) && ub - newub > 0.5) || (force && SCIPisLT(scip, newub, ub)) )
      force = TRUE;
   else
      force = FALSE;

   /* try to tighten the upper bound */
   if( global )
   {
      SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newub, force, &infeasible, &tightened) );
   }
   else
   {
      inferinfo = getInferInfo(boundtype == SCIP_BOUNDTYPE_LOWER ? varGetLbIndex(propdata, vbdvar) : varGetUbIndex(propdata, vbdvar), boundtype);

      SCIP_CALL( SCIPinferVarUbProp(scip, var, newub, prop, inferInfoToInt(inferinfo), force, &infeasible, &tightened) );
   }

   if( infeasible )
   {
      /* the infeasible results comes from the fact that the new upper bound lies below the current lower bound */
      assert(SCIPisLT(scip, newub, SCIPvarGetLbLocal(var)));
      assert(!global || SCIPisLT(scip, newub, SCIPvarGetLbGlobal(var)));

      SCIPdebugMsg(scip, "tightening%s upper bound of variable <%s> to %g due the %s bound of variable <%s> led to infeasibility\n",
         (global ? " global" : ""), SCIPvarGetName(var), newub, getBoundtypeString(boundtype), SCIPvarGetName(vbdvar));

      if( global )
      {
         /* cutoff the root node */
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
      }
      else
      {
         /* analyzes a infeasibility via conflict analysis */
         SCIP_CALL( analyzeConflictUpperbound(scip, propdata, var, newub, vbdvar, boundtype, coef, constant, canwide) );
      }
      *result = SCIP_CUTOFF;
   }
   else if( tightened )
   {
      SCIPdebugMsg(scip, "tightened%s upper bound of variable <%s> to %g due the %s bound of variable <%s>\n",
         (global ? " global" : ""), SCIPvarGetName(var), newub, getBoundtypeString(boundtype), SCIPvarGetName(vbdvar));
      (*nchgbds)++;
   }

   return SCIP_OKAY;
}

/** performs propagation of variables lower and upper bounds, implications, and cliques */
static
SCIP_RETCODE propagateVbounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< vbounds propagator */
   SCIP_Bool             force,              /**< should domain changes for continuous variables be forced */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   SCIP_VAR* startvar;
   SCIP_BOUNDTYPE starttype;
   SCIP_Real startbound;
   SCIP_Real globalbound;
   int startpos;
   int topopos;
   int v;
   int n;
   int nchgbds;
   int nbounds;
   SCIP_Bool lower;
   SCIP_Bool global;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(result != NULL);

   (*result) = SCIP_DIDNOTRUN;

   /* we do not run the propagator in presolving, because we want to avoid doing the expensive creation of the graph twice */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
      return SCIP_OKAY;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* initialize propagator data needed for propagation, if not done yet */
   if( !propdata->initialized )
   {
      SCIP_Bool infeasible;

      SCIP_CALL( initData(scip, prop, &infeasible) );

      if( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }
   assert(propdata->nbounds == 0 || propdata->propqueue != NULL);

   vars = propdata->vars;
   nbounds = propdata->nbounds;

   if( nbounds == 0 )
      return SCIP_OKAY;

   /* propagate all variables if we are in repropagation */
   if( SCIPinRepropagation(scip) )
   {
      SCIP_VAR* var;
      int idx;

      for( v = nbounds - 1; v >= 0; --v )
      {
         idx = propdata->topoorder[v];
         if( idx != -1 && !propdata->inqueue[v] )
         {
            var = vars[getVarIndex(idx)];
            lower = isIndexLowerbound(idx);
            if( !SCIPvarIsBinary(var) || (lower && SCIPvarGetLbLocal(var) > 0.5)
                  || (!lower && SCIPvarGetUbLocal(var) < 0.5) )
            {
               SCIP_CALL( SCIPpqueueInsert(propdata->propqueue, (void*)(size_t)(v + 1)) ); /*lint !e571 !e776*/
               propdata->inqueue[v] = TRUE;
            }
         }
      }
   }

   /* return if no bound changes are in the priority queue (no changed bounds to handle since last propagation) */
   if( SCIPpqueueNElems(propdata->propqueue) == 0 )
   {
      (*result) = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   nchgbds = 0;

   SCIPdebugMsg(scip, "varbound propagator: %d elements in the propagation queue\n", SCIPpqueueNElems(propdata->propqueue));

   /* get variable bound of highest priority from priority queue and try to deduce bound changes for other variables;
    * the priority queue is ordered w.r.t the topological sort of the varbound graph
    */
   while( SCIPpqueueNElems(propdata->propqueue) > 0 )
   {
      topopos = ((int)(size_t)SCIPpqueueRemove(propdata->propqueue)) - 1;
      assert(propdata->inqueue[topopos]);
      startpos = propdata->topoorder[topopos];
      assert(startpos >= 0);
      propdata->inqueue[topopos] = FALSE;

      startvar = vars[getVarIndex(startpos)];
      starttype = getBoundtype(startpos);
      lower = (starttype == SCIP_BOUNDTYPE_LOWER);
      startbound = ( lower ? SCIPvarGetLbLocal(startvar) : SCIPvarGetUbLocal(startvar) );
      globalbound = ( lower ? SCIPvarGetLbGlobal(startvar) : SCIPvarGetUbGlobal(startvar));
      global = SCIPisEQ(scip, startbound, globalbound);

      SCIPdebugMsg(scip, "propagate new %s bound of %g of variable <%s>:\n",
         getBoundtypeString(starttype), startbound, SCIPvarGetName(startvar));

      /* there should be neither implications nor cliques for non-binary variables */
      assert(SCIPvarIsBinary(startvar) || SCIPvarGetNImpls(startvar, lower) == 0);
      assert(SCIPvarIsBinary(startvar) || SCIPvarGetNCliques(startvar, lower) == 0);

      if( SCIPvarIsBinary(startvar) )
      {
         /* we only propagate binary variables if the lower bound changed to 1.0 or the upper bound changed to 0.0 */
         if( lower != (startbound > 0.5) )
            continue;

         /* propagate implications */
         if( propdata->useimplics )
         {
            int nimplvars;

            /* if the lower bound of the startvar was changed, it was fixed to 1.0, otherwise it was fixed to 0.0;
             * get all implications for this varfixing
             */
            nimplvars = SCIPvarGetNImpls(startvar, lower);

            /* if there are implications for the varfixing, propagate them */
            if( nimplvars > 0 )
            {
               SCIP_VAR** implvars;
               SCIP_BOUNDTYPE* impltypes;
               SCIP_Real* implbounds;
               int* implids;

               implvars = SCIPvarGetImplVars(startvar, lower);
               impltypes = SCIPvarGetImplTypes(startvar, lower);
               implbounds = SCIPvarGetImplBounds(startvar, lower);
               implids = SCIPvarGetImplIds(startvar, lower);

               for( n = 0; n < nimplvars; ++n )
               {
                  /* implication is just a shortcut, so we do not propagate it now,
                   * because we will propagate the longer way, anyway
                   */
                  if( implids[n] < 0 )
                     continue;

                  /* it might happen that implications point to inactive variables (normally, those are removed when a
                   * variable becomes inactive, but in some cases, it cannot be done), we have to ignore these variables
                   */
                  if( !SCIPvarIsActive(implvars[n]) )
                     continue;

                  if( impltypes[n] == SCIP_BOUNDTYPE_LOWER )
                  {
                     SCIP_CALL( tightenVarLb(scip, prop, propdata, implvars[n], implbounds[n], global, startvar,
                           starttype, force, 0.0, 0.0, FALSE, &nchgbds, result) );
                  }
                  else
                  {
                     SCIP_CALL( tightenVarUb(scip, prop, propdata, implvars[n], implbounds[n], global, startvar,
                           starttype, force, 0.0, 0.0, FALSE, &nchgbds, result) );
                  }

                  if( *result == SCIP_CUTOFF )
                     return SCIP_OKAY;
               }
            }
         }

         /* propagate cliques */
         if( propdata->usecliques )
         {
            int ncliques;

            /* if the lower bound of the startvar was changed, it was fixed to 1.0, otherwise it was fixed to 0.0;
             * get all cliques for this varfixing
             */
            ncliques = SCIPvarGetNCliques(startvar, lower);

            /* if there are cliques for the varfixing, propagate them */
            if( ncliques > 0 )
            {
               SCIP_CLIQUE** cliques;
               int j;

               cliques = SCIPvarGetCliques(startvar, lower);

               for( j = 0; j < ncliques; ++j )
               {
                  SCIP_VAR** cliquevars;
                  SCIP_Bool* cliquevals;
                  int ncliquevars;

                  cliquevars = SCIPcliqueGetVars(cliques[j]);
                  cliquevals = SCIPcliqueGetValues(cliques[j]);
                  ncliquevars = SCIPcliqueGetNVars(cliques[j]);

                  /* fix all variables except for the startvar to the value which is not in the clique */
                  for( n = 0; n < ncliquevars; ++n )
                  {
                     if( cliquevars[n] == startvar )
                        continue;

                     /* try to tighten the bound */
                     if( cliquevals[n] )
                     {
                        /* unnegated variable is in clique, so it has to be fixed to 0.0 */
                        SCIP_CALL( tightenVarUb(scip, prop, propdata, cliquevars[n], 0.0, global, startvar, starttype,
                              force, 0.0, 0.0, FALSE, &nchgbds, result) );
                     }
                     else
                     {
                        /* negated variable is in clique, so it has to be fixed to 1.0 */
                        SCIP_CALL( tightenVarLb(scip, prop, propdata, cliquevars[n], 1.0, global, startvar, starttype,
                              force, 0.0, 0.0, FALSE, &nchgbds, result) );
                     }
                     if( *result == SCIP_CUTOFF )
                        return SCIP_OKAY;
                  }
               }
            }
         }
      }

      /* propagate vbounds */
      if( propdata->usevbounds )
      {
         SCIP_VAR* boundedvar;
         SCIP_Real newbound;
         SCIP_Real coef;
         SCIP_Real constant;

         /* iterate over all vbounds for the given bound */
         for( n = 0; n < propdata->nvbounds[startpos]; ++n )
         {
            boundedvar = vars[getVarIndex(propdata->vboundboundedidx[startpos][n])];
            coef = propdata->vboundcoefs[startpos][n];
            constant = propdata->vboundconstants[startpos][n];

            /* compute new bound */
            newbound = startbound * coef + constant;

            /* try to tighten the bound */
            if( isIndexLowerbound(propdata->vboundboundedidx[startpos][n]) )
            {
               SCIP_CALL( tightenVarLb(scip, prop, propdata, boundedvar, newbound, global, startvar, starttype, force,
                     coef, constant, TRUE, &nchgbds, result) );
            }
            else
            {
               SCIP_CALL( tightenVarUb(scip, prop, propdata, boundedvar, newbound, global, startvar, starttype, force,
                     coef, constant, TRUE, &nchgbds, result) );
            }

            if( *result == SCIP_CUTOFF )
               return SCIP_OKAY;
         }
      }
   }

   SCIPdebugMsg(scip, "tightened %d variable bounds\n", nchgbds);

   /* set the result depending on whether bound changes were found or not */
   if( nchgbds > 0 )
      (*result) = SCIP_REDUCEDDOM;
   else
      (*result) = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/**@name Callback methods of propagator
 *
 * @{
 */

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyVbounds)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of propagator */
   SCIP_CALL( SCIPincludePropVbounds(scip) );

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);

   SCIPfreeBlockMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}

/** presolving initialization method of propagator (called when presolving is about to begin) */
static
SCIP_DECL_PROPINITPRE(propInitpreVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->lastpresolncliques = 0;

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int v;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* free data stored for propagation */
   if( propdata->initialized )
   {
      /* drop all variable events */
      SCIP_CALL( dropEvents(scip, propdata) );

      /* release all variables */
      for( v = 0; v < propdata->nbounds; ++v )
      {
         /* free vbound data */
         if( propdata->vboundsize[v] > 0 )
         {
            SCIPfreeMemoryArray(scip, &propdata->vboundboundedidx[v]);
            SCIPfreeMemoryArray(scip, &propdata->vboundcoefs[v]);
            SCIPfreeMemoryArray(scip, &propdata->vboundconstants[v]);
         }
      }

      /* free priority queue */
      SCIPpqueueFree(&propdata->propqueue);

      /* free arrays */
      SCIPfreeBlockMemoryArray(scip, &propdata->vboundsize, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->nvbounds, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->vboundconstants, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->vboundcoefs, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->vboundboundedidx, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->inqueue, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->topoorder, propdata->nbounds);

      /* free variable array and hashmap */
      SCIPhashmapFree(&propdata->varhashmap);
      SCIPfreeBlockMemoryArray(scip, &propdata->vars, propdata->nbounds / 2);
   }

   /* reset propagation data */
   resetPropdata(propdata);

   return SCIP_OKAY;
}

/** performs Tarjan's algorithm for strongly connected components in the implicitly given directed implication graph
 *  from the given start index; each variable x is represented by two nodes lb(x) = 2*idx(x) and ub(x) = 2*idx(x)+1
 *  where lb(x) means that the lower bound of x should be changed, i.e., that x is fixed to 1, and vice versa.
 *
 *  The algorithm is an iterative version of Tarjans algorithm
 *  (see https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm)
 *  with some additional tweaks.
 *  Each clique x_1 + ... + x_k <= 1 is represented by k(k-1) arcs (lb(x_i),ub(x_j)), j != i.
 *  This quadratic number can blow up the running time of Tarjan's algorithm, which is linear in the number of
 *  nodes and arcs of the graph. However, it suffices to consider only k of these arcs during the course of the algorithm.
 *  To this end, when we first come to a node lb(x_i) of the clique, traverse all arcs (lb(x_i),ub(x_j)) for this particular i,
 *  and store that we entered the clique via lb(x_i). Next time we come to any node lb(x_i') of the clique, we know
 *  that the only arc pointing to an unvisited node is (lb(x_i'),ub(x_i)), all other edges can be disregarded.
 *  After that, we can disregard the clique for the further search.
 *  Additionally, we try to identify infeasible fixings for binary variables. Those can be given by a path
 *  from x=1 to x=0 (or vice versa) or if x=0 (or 1) implies both y=0 and y=1.
 */
static
SCIP_RETCODE tarjan(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   startnode,          /**< node to start the depth-first-search */
   int*                  startindex,         /**< next index to assign to a processed node */
   SCIP_Shortbool*       nodeonstack,        /**< array to store the whether a each node is on the stack */
   int*                  nodeindex,          /**< array to store the dfs index for each node */
   int*                  nodelowlink,        /**< array to store the lowlink for each node */
   SCIP_Shortbool*       nodeinfeasible,     /**< array to store whether the fixing of a node was detected to be infeasible */
   int*                  dfsstack,           /**< array of size number of nodes to store the stack */
   int*                  predstackidx,       /**< for each node on the stack: stack position of its predecessor in the Tarjan search */
   int*                  stacknextclique,    /**< array of size number of nodes to store the next clique to be regarded in
                                              *   the algorithm for all nodes on the stack */
   int*                  stacknextcliquevar, /**< array of size number of nodes to store the next variable in the next clique to be
                                              * regarded in the algorithm for all nodes on the stack */
   int*                  topoorder,          /**< array with reverse (almost) topological ordering of the nodes */
   int*                  nordered,           /**< number of ordered nodes (disconnected nodes are disregarded) */
   int*                  cliquefirstentry,   /**< node from which a clique was entered for the first time; needed because when
                                              *   entering the clique a second time, only the other bound corresponding to this node
                                              *   remains to be processed */
   int*                  cliquecurrentexit,  /**< for cliques which define an arc on the current path: target node of this arc */
   int*                  sccvars,            /**< array with all nontrivial strongly connected components in the graph */
   int*                  sccstarts,          /**< start indices of SCCs in sccvars array; one additional entry at the end
                                              *   to give length of used part of sccvars array */
   int*                  nsccs,              /**< pointer to store number of strongly connected components */
   int*                  infeasnodes,        /**< sparse array with node indices of infeasible nodes */
   int*                  ninfeasnodes,       /**< pointer to store the number of infeasible nodes */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* startvar;
   SCIP_Bool lower;
   int label = *startindex;
   int stacksize;
   int currstackidx;
   int curridx;
   int idx;

   assert(startnode >= 0);
   assert(startnode < 2 * (SCIPgetNVars(scip) - SCIPgetNContVars(scip)));
   assert(nodeindex != NULL);
   assert(nodeindex[startnode] == 0);
   assert(nodelowlink != NULL);
   assert(nodelowlink[startnode] == 0);
   assert(dfsstack != NULL);
   assert(stacknextclique != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   vars = SCIPgetVars(scip);

   /* put start node on the stack */
   dfsstack[0] = startnode;
   stacknextclique[0] = 0;
   stacknextcliquevar[0] = 0;
   predstackidx[0] = -1;
   stacksize = 1;
   idx = -1;
   currstackidx = 0;
#ifdef DEBUG_TARJAN
   SCIPdebugMsg(scip, "put %s(%s) on stack[%d]\n", indexGetBoundString(dfsstack[stacksize-1]),
      SCIPvarGetName(vars[getVarIndex(dfsstack[stacksize-1])]), stacksize-1);
#endif

   /* we run until no more bounds indices are on the stack, i.e., no further nodes are connected to the startnode */
   while( stacksize > 0 )
   {
      SCIP_CLIQUE** cliques;
      int ncliques;
      SCIP_Bool found;
      int clqidx = -1;
      int j;
      int i;

      /* get next node from stack */
      curridx = dfsstack[currstackidx];
      assert(nodelowlink[curridx] <= nodeindex[curridx]);

      startvar = vars[getVarIndex(curridx)];
      lower = isIndexLowerbound(curridx);

      /* mark current node as on stack, assign index and lowlink */
      if( nodeindex[curridx] == 0 )
      {
         assert(!nodeonstack[curridx]);
         assert(stacknextclique[currstackidx] == 0);
         assert(stacknextcliquevar[currstackidx] == 0);
         nodeonstack[curridx] = 1;
         nodeindex[curridx] = label;
         nodelowlink[curridx] = label;
         ++label;

#ifdef DEBUG_TARJAN
         {
            ncliques = SCIPvarGetNCliques(startvar, lower);
            cliques = SCIPvarGetCliques(startvar, lower);

            SCIPdebugMsg(scip, "variable %s(%s) has %d cliques\n", indexGetBoundString(curridx), SCIPvarGetName(startvar),
               ncliques);
            for( j = 0; j < ncliques; ++j )
            {
               SCIP_VAR** cliquevars;
               SCIP_Bool* cliquevals;
               int ncliquevars;

               clqidx = SCIPcliqueGetIndex(cliques[j]);
               cliquevars = SCIPcliqueGetVars(cliques[j]);
               cliquevals = SCIPcliqueGetValues(cliques[j]);
               ncliquevars = SCIPcliqueGetNVars(cliques[j]);

               SCIPdebugMsg(scip, "clique %d [%d vars, stacksize: %d]...\n", clqidx, ncliquevars, stacksize);
               for( int v = 0; v < ncliquevars; ++v )
                  SCIPdebugMsgPrint(scip, " %s<%s>", cliquevals[v] ? "" : "~", SCIPvarGetName(cliquevars[v]));
               SCIPdebugMsgPrint(scip,  "\n");
            }
         }
#endif
      }
      /* we just did a backtrack and still need to investigate some outgoing edges of the node;
       * however, we should have investigated some of the outgoing edges before
       */
      else
      {
         assert(stacknextclique[currstackidx] > 0 || stacknextcliquevar[currstackidx] > 0);
         assert(nodeindex[curridx] < label);
      }
      assert(stacknextclique[currstackidx] >= 0);

      ncliques = SCIPvarGetNCliques(startvar, lower);
      cliques = SCIPvarGetCliques(startvar, lower);
      found = FALSE;

      /* iterate over all not yet handled cliques and search for an unvisited node */
      for( j = stacknextclique[currstackidx]; j < ncliques; ++j )
      {
         SCIP_VAR** cliquevars;
         SCIP_Bool* cliquevals;
         int ncliquevars;

         clqidx = SCIPcliqueGetIndex(cliques[j]);
         cliquevars = SCIPcliqueGetVars(cliques[j]);
         cliquevals = SCIPcliqueGetValues(cliques[j]);
         ncliquevars = SCIPcliqueGetNVars(cliques[j]);

         /* we did not look at this clique before from the current node, i.e., we did not backtrack now from another
          * node which was reached via this clique
          */
         if( stacknextcliquevar[currstackidx] == 0 )
         {
#ifdef DEBUG_TARJAN
            SCIPdebugMsg(scip, "clique %d [%d vars, stacksize: %d]...\n", clqidx, ncliquevars, stacksize);
            for( int v = 0; v < ncliquevars; ++v )
               SCIPdebugMsgPrint(scip, " %s<%s>", cliquevals[v] ? "" : "~", SCIPvarGetName(cliquevars[v]));
            SCIPdebugMsgPrint(scip, "\n");
#endif
            /* the clique was not entered before, remember that we first entered it from curridx
             * (add 1 to distinguish it from 0 initialization)
             */
            if( cliquefirstentry[clqidx] == 0 )
            {
               cliquefirstentry[clqidx] = curridx + 1;
            }
            else
            {
               int cliquefirstentryidx = (cliquefirstentry[clqidx] > 0 ? cliquefirstentry[clqidx] : -cliquefirstentry[clqidx]) - 1;
               int infeasnode = -1;
               assert(cliquefirstentryidx != curridx);

               /* The node by which we entered the clique the first time is still on the stack, so there is a
                * way from that node to the node by which we are entering the clique right now.
                * Since these two assignments together violate the clique and the second assignment is implied by the first,
                * the first one is infeasible
                */
               if( nodeonstack[cliquefirstentryidx] && !nodeinfeasible[cliquefirstentryidx] )
               {
                  SCIPdebugMsg(scip, "infeasible assignment (1): %s(%s)\n", indexGetBoundString(cliquefirstentryidx),
                     SCIPvarGetName(vars[getVarIndex(cliquefirstentryidx)]));
                  infeasnode = cliquefirstentryidx;
               }
               /* the first entry point of the clique was also implied by the current startnode, so this node implies
                * two variables in the clique and is therefore infeasible
                */
               else if( nodeindex[cliquefirstentryidx] >= *startindex && !nodeinfeasible[startnode] )
               {
                  SCIPdebugMsg(scip, "infeasible assignment (2): %s(%s)\n", indexGetBoundString(startnode),
                     SCIPvarGetName(vars[getVarIndex(startnode)]));
                  infeasnode = startnode;
               }

               /* we identified an infeasibility */
               if( infeasnode >= 0 )
               {
                  /* both values are invalid for the variable, the whole problem is infeasible */
                  if( nodeinfeasible[getOtherBoundIndex(infeasnode)] )
                  {
                     *infeasible = TRUE;
                     return SCIP_OKAY;
                  }
                  infeasnodes[*ninfeasnodes] = infeasnode;
                  nodeinfeasible[infeasnode] = TRUE;
                  ++(*ninfeasnodes);

                  /* the last node by which the clique was exited is not the negation of the current node and still on
                   * the stack: update the lowlink of the current node
                   */
                  if( cliquecurrentexit[clqidx] > 0
                     && curridx != getOtherBoundIndex(cliquecurrentexit[clqidx] - 1)
                     && nodeonstack[cliquecurrentexit[clqidx] - 1]
                     && nodeindex[cliquecurrentexit[clqidx] - 1] < nodelowlink[curridx] )
                  {
                     nodelowlink[curridx] = nodeindex[cliquecurrentexit[clqidx] - 1];
                  }
               }
               /* clique is entered for the second time; there is only one edge left to investigate, namely the edge to
                * the negation of the first entry point
                */
               else if( cliquefirstentry[clqidx] > 0 )
               {
#ifdef DEBUG_TARJAN
                  SCIPdebugMsg(scip, "entering clique %d a second time\n", clqidx);
#endif
                  idx = getOtherBoundIndex(cliquefirstentry[clqidx] - 1);

                  /* node was not investigated yet, we found the next node to process */
                  if( nodeindex[idx] == 0 )
                     found = TRUE;
                  /* update lowlink if the node is on the stack */
                  else if( nodeonstack[idx] && nodeindex[idx] < nodelowlink[curridx] )
                     nodelowlink[curridx] = nodeindex[idx];

                  /* cliquefirstentry[clqidx] < 0 means that we entered the clique at least two times already */
                  cliquefirstentry[clqidx] = -cliquefirstentry[clqidx];
               }
               else
               {
#ifdef DEBUG_TARJAN
                  SCIPdebugMsg(scip, "skip clique %d: visited more than twice already!\n", clqidx);
#endif
               }
               stacknextcliquevar[currstackidx] = ncliquevars;
            }
         }

         /* iterate over variables in the clique; start where we stopped last time */
         for( i = stacknextcliquevar[currstackidx]; i < ncliquevars; ++i )
         {
            if( cliquevars[i] == startvar )
               continue;

            if( !SCIPvarIsActive(cliquevars[i]) )
               continue;

            if( cliquevals[i] )
               idx = getUbIndex(SCIPvarGetProbindex(cliquevars[i]));
            else
               idx = getLbIndex(SCIPvarGetProbindex(cliquevars[i]));
            assert(idx >= 0);

            /* break when the first unvisited node is reached */
            if( nodeindex[idx] == 0 )
            {
               assert(!nodeonstack[idx]);
               stacknextcliquevar[currstackidx] = i + 1;
               found = TRUE;
               break;
            }
            else if( nodeonstack[idx] && nodeindex[idx] < nodelowlink[curridx] )
            {
               nodelowlink[curridx] = nodeindex[idx];
            }
         }
         if( found )
         {
            if( stacknextcliquevar[currstackidx] < ncliquevars )
               stacknextclique[currstackidx] = j;
            else
            {
               stacknextclique[currstackidx] = j + 1;
               stacknextcliquevar[currstackidx] = 0;
            }
            break;
         }
         else
         {
            assert(i == ncliquevars);
            stacknextclique[currstackidx] = j + 1;
            stacknextcliquevar[currstackidx] = 0;
         }
      }
      assert(found || j == ncliques);
      assert(found || stacknextclique[currstackidx] == ncliques);

      /* we stopped because we found an unhandled node and not because we reached the end of the list */
      if( found )
      {
         int otheridx = getOtherBoundIndex(idx);
         int infeasnode = -1;

         assert(idx >= 0);
         assert(!nodeonstack[idx]);
         assert(j < ncliques);
         assert(clqidx >= 0);

         /* the negated node corresponding to the next node is already on the stack -> the negated assignment is
          * infeasible
          */
         if( nodeonstack[otheridx] && !nodeinfeasible[otheridx] )
         {
            SCIPdebugMsg(scip, "infeasible assignment (3): %s(%s)\n", indexGetBoundString(otheridx),
               SCIPvarGetName(vars[getVarIndex(otheridx)]));
            infeasnode = otheridx;
         }
         /* the negated node corresponding to the next node was reached from the same startnode -> the startnode is
          * infeasible
          */
         else if( nodeindex[otheridx] >= *startindex && !nodeinfeasible[startnode] )
         {
            SCIPdebugMsg(scip, "infeasible assignment (4): %s(%s)\n", indexGetBoundString(startnode),
               SCIPvarGetName(vars[getVarIndex(startnode)]));
            infeasnode = startnode;
         }
         /* treat infeasible case */
         if( infeasnode >= 0 )
         {
            if( nodeinfeasible[getOtherBoundIndex(infeasnode)] )
            {
               *infeasible = TRUE;
               return SCIP_OKAY;
            }
            infeasnodes[*ninfeasnodes] = infeasnode;
            nodeinfeasible[infeasnode] = TRUE;
            ++(*ninfeasnodes);
         }

         SCIPdebugMsg(scip, "clique: %s(%s) -> %s(%s)\n", getBoundString(lower), SCIPvarGetName(startvar),
            indexGetBoundString(idx), SCIPvarGetName(vars[getVarIndex(idx)]));

         /* put the adjacent node onto the stack */
         dfsstack[stacksize] = idx;
         stacknextclique[stacksize] = 0;
         stacknextcliquevar[stacksize] = 0;
         cliquecurrentexit[clqidx] = idx + 1;
         predstackidx[stacksize] = currstackidx;
         currstackidx = stacksize;
         stacksize++;
         assert(stacksize <= 2 * (SCIPgetNVars(scip) - SCIPgetNContVars(scip)));

#ifdef DEBUG_TARJAN
         SCIPdebugMsg(scip, "put %s(%s) on stack[%d]\n", indexGetBoundString(dfsstack[stacksize-1]), SCIPvarGetName(vars[getVarIndex(dfsstack[stacksize-1])]), stacksize-1);
#endif
         /* restart while loop, get next index from stack */
         continue;
      }
      assert(stacknextclique[currstackidx] == ncliques);

      /* no node with a smaller index can be reached from this node -> it is the root of a SCC,
       * consisting of all nodes above it on the stack, including the node itself
       */
      if( nodelowlink[curridx] == nodeindex[curridx] )
      {
         /* we are only interested in SCCs with more than one node */
         if( dfsstack[stacksize-1] != curridx )
         {
            int sccvarspos = sccstarts[*nsccs];

            SCIPdebugMsg(scip, "SCC:");

            /* store the SCC in sccvars */
            do{
               stacksize--;
               idx = dfsstack[stacksize];
               nodeonstack[idx] = 0;
               sccvars[sccvarspos] = idx;
               ++sccvarspos;
               SCIPdebugMsgPrint(scip, " %s(%s)", indexGetBoundString(idx), SCIPvarGetName(vars[getVarIndex(idx)]));
#ifdef DEBUG_TARJAN
               SCIPdebugMsg(scip, "remove %s(%s) from stack[%d]\n", indexGetBoundString(dfsstack[stacksize]), SCIPvarGetName(vars[getVarIndex(dfsstack[stacksize])]), stacksize);
#endif

            }
            while( idx != curridx );
            SCIPdebugMsgPrint(scip, "\n");
            ++(*nsccs);
            sccstarts[*nsccs] = sccvarspos;
         }
         /* trivial SCC: remove the single node from the stack, but don't store it as a SCC */
         else
         {
            stacksize--;
#ifdef DEBUG_TARJAN
            SCIPdebugMsg(scip, "remove %s(%s) from stack[%d]\n", indexGetBoundString(dfsstack[stacksize]), SCIPvarGetName(vars[getVarIndex(dfsstack[stacksize])]), stacksize);
#endif
            idx = dfsstack[stacksize];
            nodeonstack[idx] = 0;
            assert(nodeindex[idx] > 0);
         }
      }
      /* in a pure dfs, the node would now leave the stack, add it to the array of nodes in reverse topological order */
      if( topoorder != NULL && (stacksize > 0 || label > *startindex + 1) )
      {
         topoorder[*nordered] = curridx;
         ++(*nordered);
      }

      /* the current node was handled, backtrack */
      if( stacksize > 0 )
      {
         idx = dfsstack[predstackidx[currstackidx]];
         nodelowlink[idx] = MIN(nodelowlink[idx], nodelowlink[curridx]);
         currstackidx = predstackidx[currstackidx];
      }
   }

   *startindex = label;

   return SCIP_OKAY;
}


/** apply fixings and aggregations found by the clique graph analysis */
static
SCIP_RETCODE applyFixingsAndAggregations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of active variables */
   int*                  infeasnodes,        /**< sparse array with node indices of infeasible nodes */
   int                   ninfeasnodes,       /**< pointer to store the number of infeasible nodes */
   SCIP_Shortbool*       nodeinfeasible,     /**< array to store whether the fixing of a node was detected to be infeasible */
   int*                  sccvars,            /**< array with all nontrivial strongly connected components in the graph */
   int*                  sccstarts,          /**< start indices of SCCs in sccvars array; one additional entry at the end
                                              *   to give length of used part of sccvars array */
   int                   nsccs,              /**< pointer to store number of strongly connected components */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nfixedvars,         /**< pointer to number of fixed variables, increment when fixing another one */
   int*                  naggrvars,          /**< pointer to number of aggregated variables, increment when aggregating another one */
   SCIP_RESULT*          result              /**< pointer to store result of the call */
   )
{
   int i = 0;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(infeasible != NULL);

   /* for all infeasible node: fix variable to the other bound */
   if( !(*infeasible) && ninfeasnodes > 0 )
   {
      for( i = 0; i < ninfeasnodes; ++i )
      {
         SCIP_VAR* var = vars[getVarIndex(infeasnodes[i])];
         SCIP_Bool lower = isIndexLowerbound(infeasnodes[i]);
         SCIP_Bool fixed;

         assert(nodeinfeasible[infeasnodes[i]]);
         nodeinfeasible[infeasnodes[i]] = FALSE;

         SCIP_CALL( SCIPfixVar(scip, var, lower ? 0.0 : 1.0, infeasible, &fixed) );

         SCIPdebugMsg(scip, "fix <%s>[%d] to %g: inf=%d, fixed=%d\n",
            SCIPvarGetName(var), infeasnodes[i], lower ? 0.0 : 1.0, *infeasible, fixed);

         /* fixing was infeasible */
         if( *infeasible )
            break;

         /* increase fixing counter and update result pointer */
         if( fixed )
         {
            *result = SCIP_SUCCESS;
            ++(*nfixedvars);
         }
      }
   }
   assert((*infeasible) || i == ninfeasnodes);

   /* clear clean buffer array (if we did not enter the block above or stopped early due to an infeasibility) */
   for( ; i < ninfeasnodes; ++i )
   {
      assert(nodeinfeasible[infeasnodes[i]]);
      nodeinfeasible[infeasnodes[i]] = FALSE;
   }

   if( !(*infeasible) && nsccs > 0 )
   {
      /* for each strongly connected component: aggregate all variables to the first one */
      for( i = 0; i < nsccs; ++i )
      {
         SCIP_VAR* startvar;
         SCIP_Bool lower;
         SCIP_Bool aggregated;
         SCIP_Bool redundant;
         int v;

         assert(sccstarts[i] < sccstarts[i+1] - 1);

         /* get variable and boundtype for first node of the SCC */
         startvar = vars[getVarIndex(sccvars[sccstarts[i]])];
         lower = isIndexLowerbound(sccvars[sccstarts[i]]);

         for( v = sccstarts[i] + 1; v < sccstarts[i+1]; ++v )
         {
            /* aggregate variables: if both nodes represent the same bound, we have x=1 <=> y=1,
             * and thus aggregate x - y = 0; if both represent different bounds we have
             * x=1 <=> y=0, so we aggregate x + y = 1
             */
            SCIP_CALL( SCIPaggregateVars(scip, startvar, vars[getVarIndex(sccvars[v])], 1.0,
                  lower == isIndexLowerbound(sccvars[v]) ? -1.0 : 1.0,
                  lower == isIndexLowerbound(sccvars[v]) ? 0.0 : 1.0,
                  infeasible, &redundant, &aggregated) );

            SCIPdebugMsg(scip, "aggregate <%s> + %g <%s> = %g: inf=%d, red=%d, aggr=%d\n",
               SCIPvarGetName(startvar), lower == isIndexLowerbound(sccvars[v]) ? -1.0 : 1.0,
               SCIPvarGetName(vars[getVarIndex(sccvars[v])]), lower == isIndexLowerbound(sccvars[v]) ? 0.0 : 1.0,
               *infeasible, redundant, aggregated);

            /* aggregation was infeasible */
            if( *infeasible )
               break;

            /* increase aggregation counter and update result pointer */
            if( aggregated )
            {
               ++(*naggrvars);
               *result = SCIP_SUCCESS;
            }
         }
      }
   }

   return SCIP_OKAY;
}



/** presolving method of propagator: search for strongly connected components in the implication graph and
 *  aggregate all variables within a component; additionally, identifies infeasible variable assignments
 *  as a side product if a path from x=1 to x=0 (or vice versa) is found or x=1 implies both y=0 and y=1
 *  The identification of such assignments depends on the order in which variable bounds are processed;
 *  therefore, we are doing a second run with the bounds processed in (almost) topological order.
 */
static
SCIP_DECL_PROPPRESOL(propPresolVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** tmpvars;
   SCIP_VAR** vars;
   int* dfsstack;
   int* stacknextclique;
   int* stacknextcliquevar;
   int* nodeindex;
   int* nodelowlink;
   int* predstackidx;
   int* cliquefirstentry;
   int* cliquecurrentexit;
   int* topoorder;
   int* sccvars;
   int* sccstarts;
   int* infeasnodes;
   SCIP_Shortbool* nodeonstack;
   SCIP_Shortbool* nodeinfeasible;
   int ninfeasnodes;
   int nsccs;
   int nbounds;
   int nbinvars;
   int ncliques;
   int startindex = 1;
   int nordered = 0;
   int i;
   SCIP_Bool infeasible = FALSE;

   assert(scip != NULL);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   ncliques = SCIPgetNCliques(scip);

   *result = SCIP_DIDNOTRUN;

   if( ncliques < 2 )
      return SCIP_OKAY;

   /* too many cliques for medium presolving */
   if( presoltiming == SCIP_PRESOLTIMING_MEDIUM && ncliques > propdata->maxcliquesmedium * SCIPgetNBinVars(scip) )
      return SCIP_OKAY;

   /* too many cliques for medium presolving */
   if( ncliques > propdata->maxcliquesexhaustive * SCIPgetNBinVars(scip) )
      return SCIP_OKAY;

   /* only run if enough new cliques were created since the last successful call */
   if( SCIPgetNCliquesCreated(scip) < (1.0 + propdata->minnewcliques) * propdata->lastpresolncliques )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   nbinvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   nbounds = 2 * nbinvars;

   /* cleanup cliques, stop if this proved infeasibility already */
   SCIP_CALL( SCIPcleanupCliques(scip, &infeasible) );

   if( infeasible )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   tmpvars = SCIPgetVars(scip);

   /* duplicate variable array; needed to get the fixings right later */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, tmpvars, nbinvars) );

   SCIP_CALL( SCIPallocBufferArray(scip, &dfsstack, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stacknextclique, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stacknextcliquevar, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &predstackidx, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &topoorder, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sccvars, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sccstarts, nbinvars + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &infeasnodes, nbounds) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nodeindex, nbounds) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nodelowlink, nbounds) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &cliquefirstentry, ncliques) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &cliquecurrentexit, ncliques) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nodeonstack, nbounds) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &nodeinfeasible, nbounds) );
   sccstarts[0] = 0;
   nsccs = 0;
   ninfeasnodes = 0;

   /* while there are unvisited nodes, run Tarjan's algorithm starting from one of these nodes */
   for( i = 0; i < nbounds && !infeasible; ++i )
   {
      if( nodeindex[i] == 0 )
      {
         SCIP_CALL( tarjan(scip, i, &startindex, nodeonstack, nodeindex, nodelowlink, nodeinfeasible,
               dfsstack, predstackidx, stacknextclique, stacknextcliquevar, topoorder, &nordered,
               cliquefirstentry, cliquecurrentexit, sccvars, sccstarts, &nsccs,
               infeasnodes, &ninfeasnodes, &infeasible) );
      }
   }
   assert(nordered <= nbounds);

   /* aggregate all variables within a SCC and fix all variables for which one bounds was proven infeasible */
   if( ninfeasnodes > 0 || nsccs > 0 )
   {
      SCIP_CALL( applyFixingsAndAggregations(scip, vars, infeasnodes, ninfeasnodes, nodeinfeasible,
            sccvars, sccstarts, nsccs, &infeasible, nfixedvars, naggrvars, result) );
   }

   /* second round, now with topological order! */
   if( !infeasible && nordered > 0 )
   {
      SCIP_VAR** vars2;
      int nbounds2;

      assert(nordered > 1);

      /* we already fixed or aggregated some variables in the first run, so we better clean up the cliques */
      if( *result == SCIP_SUCCESS )
      {
         SCIP_CALL( SCIPcleanupCliques(scip, &infeasible) );

         if( infeasible )
            goto TERMINATE;
      }

      nbounds2 = 2 * (SCIPgetNVars(scip) - SCIPgetNContVars(scip));
      ncliques = SCIPgetNCliques(scip);

      SCIP_CALL( SCIPduplicateBufferArray(scip, &vars2, tmpvars, nbounds2/2) );

      /* clear arrays that should be initialized to 0 */
      BMSclearMemoryArray(nodeonstack, nbounds2);
      BMSclearMemoryArray(nodeindex, nbounds2);
      BMSclearMemoryArray(nodelowlink, nbounds2);
      BMSclearMemoryArray(cliquefirstentry, ncliques);
      BMSclearMemoryArray(cliquecurrentexit, ncliques);
      sccstarts[0] = 0;
      nsccs = 0;
      ninfeasnodes = 0;
      startindex = 1;

      /* while there are unvisited nodes, run Tarjan's algorithm starting from one of these nodes */
      for( i = nordered - 1; i >= 0  && !infeasible; --i )
      {
         int varindex;
         int startpos;
         assert(topoorder[i] < nbounds);

         /* variable of next node in topological order */
         varindex = SCIPvarGetProbindex(vars[getVarIndex(topoorder[i])]);

         /* variable was not fixed after the first run */
         if( varindex >= 0 )
         {
            startpos = isIndexLowerbound(topoorder[i]) ? getLbIndex(varindex) : getUbIndex(varindex);
            if( nodeindex[startpos] == 0 )
            {
               SCIP_CALL( tarjan(scip, startpos, &startindex, nodeonstack, nodeindex, nodelowlink, nodeinfeasible,
                     dfsstack, predstackidx, stacknextclique, stacknextcliquevar, NULL, NULL,
                     cliquefirstentry, cliquecurrentexit, sccvars, sccstarts, &nsccs,
                     infeasnodes, &ninfeasnodes, &infeasible) );
            }
         }
      }

      /* aggregate all variables within a SCC and fix all variables for which one bounds was proven infeasible */
      if( ninfeasnodes > 0 || nsccs > 0 )
      {
         SCIP_CALL( applyFixingsAndAggregations(scip, vars2, infeasnodes, ninfeasnodes, nodeinfeasible,
               sccvars, sccstarts, nsccs, &infeasible, nfixedvars, naggrvars, result) );
      }

      SCIPfreeBufferArray(scip, &vars2);
   }

 TERMINATE:
   if( infeasible )
      *result = SCIP_CUTOFF;
#ifndef NDEBUG
   for( i = 0; i < nbounds; ++i )
   {
      assert(nodeinfeasible[i] == FALSE);
   }
#endif
   SCIPfreeCleanBufferArray(scip, &nodeinfeasible);

   SCIPfreeBufferArray(scip, &cliquecurrentexit);
   SCIPfreeBufferArray(scip, &cliquefirstentry);

   SCIPfreeBufferArray(scip, &nodelowlink);
   SCIPfreeBufferArray(scip, &nodeindex);
   SCIPfreeBufferArray(scip, &nodeonstack);
   SCIPfreeBufferArray(scip, &infeasnodes);
   SCIPfreeBufferArray(scip, &sccstarts);
   SCIPfreeBufferArray(scip, &sccvars);
   SCIPfreeBufferArray(scip, &topoorder);
   SCIPfreeBufferArray(scip, &predstackidx);
   SCIPfreeBufferArray(scip, &stacknextcliquevar);
   SCIPfreeBufferArray(scip, &stacknextclique);
   SCIPfreeBufferArray(scip, &dfsstack);
   SCIPfreeBufferArray(scip, &vars);

   propdata->lastpresolncliques = SCIPgetNCliquesCreated(scip);

   return SCIP_OKAY;
}



/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecVbounds)
{  /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;

   /* perform variable lower and upper bound propagation */
   SCIP_CALL( propagateVbounds(scip, prop, FALSE, result) );

   assert((*result) == SCIP_CUTOFF || (*result) == SCIP_DIDNOTRUN
      || (*result) == SCIP_DIDNOTFIND || (*result) == SCIP_REDUCEDDOM);

   return SCIP_OKAY;
}

/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   SCIP_VAR* startvar;
   SCIP_BOUNDTYPE starttype;
   int pos;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   starttype = inferInfoGetBoundtype(intToInferInfo(inferinfo));
   pos = inferInfoGetPos(intToInferInfo(inferinfo));
   assert(pos >= 0);
   assert(pos < propdata->nbounds);

   vars = propdata->vars;
   assert(vars != NULL);
   startvar = vars[getVarIndex(pos)];
   assert(startvar != NULL);
   assert(startvar != infervar);

   SCIPdebugMsg(scip, "explain %s bound change of variable <%s>\n",
      getBoundtypeString(boundtype), SCIPvarGetName(infervar));

   if( !SCIPvarIsBinary(startvar) && propdata->usebdwidening )
   {
      int* vboundboundedidx;
      SCIP_Real constant;
      SCIP_Real coef;
      int inferidx;
      int nvbounds;
      int b;

      nvbounds = propdata->nvbounds[pos];
      vboundboundedidx = propdata->vboundboundedidx[pos];

      inferidx = boundtype == SCIP_BOUNDTYPE_LOWER ? varGetLbIndex(propdata, infervar) : varGetUbIndex(propdata, infervar);
      assert(inferidx >= 0);

      for( b = 0; b < nvbounds; ++b )
      {
         if( vboundboundedidx[b] == inferidx )
            break;
      }
      assert(b < nvbounds);

      coef = propdata->vboundcoefs[pos][b];
      constant = propdata->vboundconstants[pos][b];
      assert(!SCIPisZero(scip, coef));

      /* compute the relaxed bound which is sufficient to propagate the inference bound of given variable */
      if( boundtype == SCIP_BOUNDTYPE_LOWER )
         relaxedbd = computeRelaxedLowerbound(scip, infervar, relaxedbd, coef, constant);
      else
         relaxedbd = computeRelaxedUpperbound(scip, infervar, relaxedbd, coef, constant);

      /* try to relax variable bound variable */
      SCIP_CALL( relaxVbdvar(scip, startvar, starttype, bdchgidx, relaxedbd) );
   }
   else
   {
      SCIP_CALL( resolvePropagation(scip, propdata, startvar, starttype, bdchgidx) );
   }

   (*result) = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods of event handler
 *
 * @{
 */

/** execution method of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecVbound)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int idx;

   assert(eventhdlr != NULL);

   propdata = (SCIP_PROPDATA*)SCIPeventhdlrGetData(eventhdlr);
   assert(propdata != NULL);

   idx = (int) (size_t) eventdata;
   assert(idx >= 0);

   SCIPdebugMsg(scip, "eventexec (type=%llu): try to add sort index %d: %s(%s) to priority queue\n", SCIPeventGetType(event),
      idx, indexGetBoundString(propdata->topoorder[idx]),
      SCIPvarGetName(propdata->vars[getVarIndex(propdata->topoorder[idx])]));

   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_GUBCHANGED && SCIPvarIsBinary(SCIPeventGetVar(event))
      && SCIPeventGetNewbound(event) > 0.5 )
      return SCIP_OKAY;

   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_GLBCHANGED && SCIPvarIsBinary(SCIPeventGetVar(event))
      && SCIPeventGetNewbound(event) < 0.5 )
      return SCIP_OKAY;

   assert(getVarIndex(propdata->topoorder[idx]) < SCIPgetNVars(scip));
   assert(SCIPvarGetType(propdata->vars[getVarIndex(propdata->topoorder[idx])]) != SCIP_VARTYPE_BINARY
      || (isIndexLowerbound(propdata->topoorder[idx]) == (SCIPeventGetNewbound(event) > 0.5)));

   /* add the bound change to the propagation queue, if it is not already contained */
   if( !propdata->inqueue[idx] )
   {
      SCIP_CALL( SCIPpqueueInsert(propdata->propqueue, (void*)(size_t)(idx + 1)) ); /*lint !e571 !e776*/
      propdata->inqueue[idx] = TRUE;
   }
   assert(SCIPpqueueNElems(propdata->propqueue) > 0);

   return SCIP_OKAY;
}

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the vbounds propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropVbounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create vbounds propagator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );

   /* reset propagation data */
   resetPropdata(propdata);

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecVbounds, propdata) );
   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyVbounds) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeVbounds) );
   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitpreVbounds) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolVbounds) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropVbounds) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolVbounds, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS,
         PROP_PRESOLTIMING) );

   /* include event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &propdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecVbound, (SCIP_EVENTHDLRDATA*)propdata) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/usebdwidening", "should bound widening be used to initialize conflict analysis?",
         &propdata->usebdwidening, FALSE, DEFAULT_USEBDWIDENING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/useimplics", "should implications be propagated?",
         &propdata->useimplics, TRUE, DEFAULT_USEIMPLICS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/usecliques", "should cliques be propagated?",
         &propdata->usecliques, TRUE, DEFAULT_USECLIQUES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/usevbounds", "should vbounds be propagated?",
         &propdata->usevbounds, TRUE, DEFAULT_USEVBOUNDS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/dotoposort", "should the bounds be topologically sorted in advance?",
         &propdata->dotoposort, FALSE, DEFAULT_DOTOPOSORT, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/sortcliques", "should cliques be regarded for the topological sort?",
         &propdata->sortcliques, TRUE, DEFAULT_SORTCLIQUES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/detectcycles", "should cycles in the variable bound graph be identified?",
         &propdata->detectcycles, FALSE, DEFAULT_DETECTCYCLES, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "propagating/" PROP_NAME "/minnewcliques", "minimum percentage of new cliques to trigger another clique table analysis",
         &propdata->minnewcliques, FALSE, DEFAULT_MINNEWCLIQUES, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/maxcliquesmedium",
         "maximum number of cliques per variable to run clique table analysis in medium presolving",
         &propdata->maxcliquesmedium, FALSE, DEFAULT_MAXCLIQUESMEDIUM, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/maxcliquesexhaustive",
         "maximum number of cliques per variable to run clique table analysis in exhaustive presolving",
         &propdata->maxcliquesexhaustive, FALSE, DEFAULT_MAXCLIQUESEXHAUSTIVE, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** returns TRUE if the propagator has the status that all variable lower and upper bounds are propgated */
SCIP_Bool SCIPisPropagatedVbounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   prop = SCIPfindProp(scip, PROP_NAME);
   assert(prop != NULL);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   return (SCIPpqueueNElems(propdata->propqueue) == 0);
}

/** performs propagation of variables lower and upper bounds */
SCIP_RETCODE SCIPexecPropVbounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             force,              /**< should domain changes for continuous variables be forced */
   SCIP_RESULT*          result              /**< pointer to store result */
   )
{
   SCIP_PROP* prop;

   prop = SCIPfindProp(scip, PROP_NAME);
   assert(prop != NULL);

   /* perform variable lower and upper bound propagation */
   SCIP_CALL( propagateVbounds(scip, prop, force, result) );

   assert((*result) == SCIP_CUTOFF || (*result) == SCIP_DIDNOTRUN
      || (*result) == SCIP_DIDNOTFIND || (*result) == SCIP_REDUCEDDOM);

   return SCIP_OKAY;
}

/**@} */
