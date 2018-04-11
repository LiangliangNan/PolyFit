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

/**@file   reopt.c
 * @brief  data structures and methods for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/mem.h"
#include "scip/event.h"
#include "scip/scip.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/var.h"
#include "scip/lp.h"
#include "scip/misc.h"
#include "scip/reopt.h"
#include "scip/tree.h"
#include "scip/primal.h"
#include "scip/sepastore.h"
#include "scip/cutpool.h"
#include "scip/prob.h"
#include "scip/cons.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_linear.h"
#include "scip/clock.h"
#include "scip/heur_reoptsols.h"
#include "scip/history.h"
#include "blockmemshell/memory.h"

#define DEFAULT_MEM_VARAFTERDUAL    10
#define DEFAULT_MEM_VAR             10
#define DEFAULT_MEM_NODES         1000
#define DEFAULT_MEM_RUN            200
#define DEFAULT_MEM_DUALCONS        10

#define DEFAULT_RANDSEED            67

/* event handler properties */
#define EVENTHDLR_NAME         "Reopt"
#define EVENTHDLR_DESC         "node event handler for reoptimization"

/* ---------------- Callback methods of event handler ---------------- */

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecReopt)
{  /*lint --e{715}*/
   SCIP_NODE*          eventnode;
   SCIP_Real           oldbound;
   SCIP_Real           newbound;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(SCIPvarGetType(SCIPeventGetVar(event)) != SCIP_VARTYPE_CONTINUOUS);

   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   eventnode = SCIPgetCurrentNode(scip);
   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);

   assert( eventnode != NULL );

   /* skip if the node is not the focus nodes */
   if( SCIPnodeGetType(eventnode) != SCIP_NODETYPE_FOCUSNODE || SCIPnodeGetDepth(eventnode) != SCIPgetEffectiveRootDepth(scip) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "catch event for node %lld: <%s>: %g -> %g\n", SCIPnodeGetNumber(eventnode),
         SCIPvarGetName(SCIPeventGetVar(event)), SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));

   assert(SCIPisFeasLT(scip, newbound, oldbound) || SCIPisFeasGT(scip, newbound, oldbound));

   SCIP_CALL( SCIPaddReoptDualBndchg(scip, eventnode, SCIPeventGetVar(event), newbound, oldbound) );

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolReopt)
{
   SCIP_VAR** vars;
   int varnr;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   if( !SCIPisReoptEnabled(scip) )
      return SCIP_OKAY;

   vars = SCIPgetVars(scip);
   for(varnr = 0; varnr < SCIPgetNVars(scip); ++varnr)
   {
      if( SCIPvarGetType(vars[varnr]) != SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_CALL(SCIPcatchVarEvent(scip, vars[varnr], SCIP_EVENTTYPE_GBDCHANGED, eventhdlr, NULL, NULL));
      }
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolReopt)
{
   SCIP_VAR** vars;
   int varnr;
   assert(scip != NULL);

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   if( !SCIPisReoptEnabled(scip) )
      return SCIP_OKAY;

   vars = SCIPgetVars(scip);

   for(varnr = 0; varnr < SCIPgetNVars(scip); ++varnr)
   {
      if( SCIPvarGetType(vars[varnr]) == SCIP_VARTYPE_BINARY )
      {
         SCIP_CALL(SCIPdropVarEvent(scip, vars[varnr], SCIP_EVENTTYPE_GBDCHANGED , eventhdlr, NULL, -1));
      }
   }
   return SCIP_OKAY;
}

/* ---------------- Callback methods of reoptimization methods ---------------- */

/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that sols[pos] array can store at least num entries */
static
SCIP_RETCODE ensureSolsSize(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   num,                /**< minimum number of entries to store */
   int                   runidx              /**< run index for which the memory should checked */
   )
{
   assert(runidx >= 0);
   assert(runidx <= reopt->runsize);

   if( num > reopt->soltree->solssize[runidx] )
   {
      int newsize = SCIPsetCalcMemGrowSize(set, num + 1);

      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->soltree->sols[runidx],
            reopt->soltree->solssize[runidx], newsize) ); /*lint !e866 */

      reopt->soltree->solssize[runidx] = newsize;
   }
   assert(num <= reopt->soltree->solssize[runidx]);

   return SCIP_OKAY;
}

/** ensures, that sols array can store at least num entries */
static
SCIP_RETCODE ensureRunSize(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< gloabl SCIP settings */
   int                   num,                /**< minimum number of entries to store */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   if( num >= reopt->runsize )
   {
      int s;
      int newsize = SCIPsetCalcMemGrowSize(set, num+1);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->soltree->sols, reopt->runsize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->soltree->nsols, reopt->runsize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->soltree->solssize, reopt->runsize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->prevbestsols, reopt->runsize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->varhistory, reopt->runsize, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&reopt->objs, newsize) );

      for(s = reopt->runsize; s < newsize; s++)
      {
         reopt->varhistory[s] = NULL;
         reopt->prevbestsols[s] = NULL;
         reopt->objs[s] = NULL;
         reopt->soltree->solssize[s] = 0;
         reopt->soltree->nsols[s] = 0;
         reopt->soltree->sols[s] = NULL;
      }

      reopt->runsize = newsize;
   }
   assert(num < reopt->runsize);

   return SCIP_OKAY;
}

/** check the memory of the reoptimization tree and if necessary reallocate */
static
SCIP_RETCODE reopttreeCheckMemory(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopttree != NULL);
   assert(blkmem != NULL);

   if( SCIPqueueIsEmpty(reopttree->openids) )
   {
      int newsize;
      unsigned int id;

      assert(reopttree->nreoptnodes == (int)(reopttree->reoptnodessize));

      newsize = SCIPsetCalcMemGrowSize(set, (int)reopttree->reoptnodessize+1);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopttree->reoptnodes, reopttree->reoptnodessize, newsize) ); /*lint !e647*/

      for( id = reopttree->reoptnodessize; id < (unsigned int)newsize; id++ )
      {
         SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void*) (size_t) id) ); /*lint !e571*/
         reopttree->reoptnodes[id] = NULL;
      }

      reopttree->reoptnodessize = (unsigned int)newsize;
   }

   return SCIP_OKAY;
}

/** check allocated memory of a node within the reoptimization tree and if necessary reallocate */
static
SCIP_RETCODE reoptnodeCheckMemory(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   var_mem,            /**< memory for variables */
   int                   child_mem,          /**< memory for child nodes */
   int                   conss_mem           /**< memory for constraints */
   )
{
   int newsize;

   assert(reoptnode != NULL);
   assert(blkmem != NULL);
   assert(var_mem >= 0);
   assert(child_mem >= 0);
   assert(conss_mem >= 0);

   /* check allocated memory for variable and bound information */
   if( var_mem > 0 )
   {
      if( reoptnode->varssize == 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->vars, var_mem) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->varbounds, var_mem) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->varboundtypes, var_mem) );
         reoptnode->varssize = var_mem;
      }
      else if( reoptnode->varssize < var_mem )
      {
         newsize = SCIPsetCalcMemGrowSize(set, var_mem+1);
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->vars, reoptnode->varssize, newsize) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->varbounds, reoptnode->varssize, newsize) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->varboundtypes, reoptnode->varssize, newsize) );
         reoptnode->varssize = newsize;
      }
   }

   /* check allocated memory for child node information */
   if( child_mem > 0 )
   {
      if( reoptnode->allocchildmem == 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->childids, child_mem) );
         reoptnode->nchilds = 0;
         reoptnode->allocchildmem = child_mem;
      }
      else if( reoptnode->allocchildmem < child_mem )
      {
         newsize = SCIPsetCalcMemGrowSize(set, child_mem+1);
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->childids, reoptnode->allocchildmem, newsize) );
         reoptnode->allocchildmem = newsize;
      }
   }

   /* check allocated memory for add constraints */
   if( conss_mem > 0 )
   {
      if( reoptnode->consssize == 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->conss, conss_mem) );
         reoptnode->nconss = 0;
         reoptnode->consssize = conss_mem;
      }
      else if( reoptnode->consssize < conss_mem )
      {
         newsize = SCIPsetCalcMemGrowSize(set, conss_mem);
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->conss, reoptnode->consssize, newsize) );
         reoptnode->consssize = newsize;
      }
   }

   return SCIP_OKAY;
}

/*
 * local methods
 */

/** returns the number of stored solutions in the subtree induced by @p solnode */
static
int soltreeNInducedSols(
   SCIP_SOLNODE*         solnode             /**< node within the solution tree */
   )
{
   SCIP_SOLNODE* sibling;
   int nsols;

   assert(solnode != NULL);

   if( solnode->child == NULL && solnode->sol == NULL )
      return 0;
   if( solnode->child == NULL && solnode->sol != NULL )
      return 1;

   nsols = 0;
   sibling = solnode->child;

   /* traverse through the list */
   while( sibling != NULL )
   {
      nsols += soltreeNInducedSols(sibling);
      sibling = sibling->sibling;
   }

   return nsols;
}

/** returns the similarity of the objective functions of two given iterations */
static
SCIP_Real reoptSimilarity(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   obj1_id,            /**< id of one objective function */
   int                   obj2_id,            /**< id of the other objective function */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars               /**< number of problem variables */
   )
{
   SCIP_Real similarity;
   SCIP_Real norm_obj1;
   SCIP_Real norm_obj2;
   int v;

   assert(reopt != NULL);
   assert(vars != NULL);
   assert(nvars >= 0);

   similarity = 0.0;
   norm_obj1 = 0.0;
   norm_obj2 = 0.0;

   /* calculate similarity */
   for( v = 0; v < nvars; v++ )
   {
      SCIP_VAR* origvar;
      SCIP_VAR* transvar;
      SCIP_Real c1;
      SCIP_Real c2;
      SCIP_Real lb;
      SCIP_Real ub;

      origvar = vars[v];

      /* get the original variable */
      if( !SCIPvarIsOriginal(origvar) )
      {
         SCIP_RETCODE retcode;
         SCIP_Real constant = 0.0;
         SCIP_Real scalar = 1.0;

         retcode = SCIPvarGetOrigvarSum(&origvar, &scalar, &constant);

         if( retcode != SCIP_OKAY )
            return SCIP_INVALID;
      }
      assert(origvar != NULL && SCIPvarIsOriginal(origvar));

      /* get the transformed variable, we skip globally fixed variables */
      transvar = SCIPvarGetTransVar(origvar);
      assert(transvar != NULL);

      lb = SCIPvarGetLbLocal(transvar);
      ub = SCIPvarGetUbLocal(transvar);

      if( SCIPsetIsFeasLT(set, lb, ub) )
      {
         int probidx;

         probidx = SCIPvarGetIndex(origvar);
         assert(0 <= probidx && probidx < reopt->nobjvars);

         c1 = reopt->objs[obj1_id][probidx];
         c2 = reopt->objs[obj2_id][probidx];

         /* vector product */
         similarity += c1*c2;
         norm_obj1 += SQR(c1);
         norm_obj2 += SQR(c2);
      }
   }

   /* divide similarity by norms of the objective vectors */
   norm_obj1 = SQRT(norm_obj1);
   norm_obj2 = SQRT(norm_obj2);

   if( !SCIPsetIsZero(set, norm_obj1) && !SCIPsetIsZero(set, norm_obj2) )
      similarity /= (norm_obj1 * norm_obj2);

   /* make sure that we are between -1.0 und +1.0 */
   similarity = MAX(similarity, -1.0);
   similarity = MIN(similarity, 1.0);

   return similarity;
}

/** delete the given reoptimization node */
static
SCIP_RETCODE reoptnodeDelete(
   SCIP_REOPTNODE**      reoptnode,          /**< node of the reoptimization tree */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert((*reoptnode) != NULL );
   assert(blkmem != NULL );

   /* delete data for constraints */
   if( (*reoptnode)->consssize > 0 )
   {
      int c;

      assert((*reoptnode)->conss != NULL);

      for( c = 0; c < (*reoptnode)->nconss; c++ )
      {
         assert((*reoptnode)->conss[c] != NULL);
         assert((*reoptnode)->conss[c]->vals != NULL);
         assert((*reoptnode)->conss[c]->vars != NULL);

         BMSfreeBlockMemoryArrayNull(blkmem, &(*reoptnode)->conss[c]->boundtypes, (*reoptnode)->conss[c]->varssize);
         BMSfreeBlockMemoryArrayNull(blkmem, &(*reoptnode)->conss[c]->vals, (*reoptnode)->conss[c]->varssize);
         BMSfreeBlockMemoryArrayNull(blkmem, &(*reoptnode)->conss[c]->vars, (*reoptnode)->conss[c]->varssize);
         BMSfreeBlockMemory(blkmem, &(*reoptnode)->conss[c]); /*lint !e866*/
      }
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->conss, (*reoptnode)->consssize);
      (*reoptnode)->nconss = 0;
      (*reoptnode)->consssize = 0;
      (*reoptnode)->conss = NULL;
   }

   /* free list of children */
   if( (*reoptnode)->childids != NULL )
   {
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->childids, (*reoptnode)->allocchildmem);
      (*reoptnode)->nchilds = 0;
      (*reoptnode)->allocchildmem = 0;
      (*reoptnode)->childids = NULL;
   }

   /* delete dual constraint */
   if( (*reoptnode)->dualredscur != NULL )
   {
      assert((*reoptnode)->dualredscur->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualredscur->boundtypes, (*reoptnode)->dualredscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualredscur->vals, (*reoptnode)->dualredscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualredscur->vars, (*reoptnode)->dualredscur->varssize);
      BMSfreeBlockMemory(blkmem, &(*reoptnode)->dualredscur);
      (*reoptnode)->dualredscur = NULL;
   }

   /* delete dual constraint */
   if( (*reoptnode)->dualredsnex != NULL )
   {
      assert((*reoptnode)->dualredsnex->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualredsnex->boundtypes, (*reoptnode)->dualredsnex->varssize);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualredsnex->vals, (*reoptnode)->dualredsnex->varssize);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualredsnex->vars, (*reoptnode)->dualredsnex->varssize);
      BMSfreeBlockMemory(blkmem, &(*reoptnode)->dualredsnex);
      (*reoptnode)->dualredsnex = NULL;
   }

   /* free boundtypes */
   if ((*reoptnode)->varboundtypes != NULL )
   {
      assert((*reoptnode)->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->varboundtypes, (*reoptnode)->varssize);
      (*reoptnode)->varboundtypes = NULL;
   }

   /* free bounds */
   if ((*reoptnode)->varbounds != NULL )
   {
      assert((*reoptnode)->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->varbounds, (*reoptnode)->varssize);
      (*reoptnode)->varbounds = NULL;
   }

   /* free variables */
   if ((*reoptnode)->vars != NULL )
   {
      assert((*reoptnode)->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->vars, (*reoptnode)->varssize);
      (*reoptnode)->vars = NULL;
   }

   (*reoptnode)->varssize = 0;

   /* free afterdual-boundtypes */
   if ((*reoptnode)->afterdualvarboundtypes != NULL )
   {
      assert((*reoptnode)->afterdualvarssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->afterdualvarboundtypes, (*reoptnode)->afterdualvarssize);
      (*reoptnode)->afterdualvarboundtypes = NULL;
   }

   /* free afterdual-bounds */
   if ((*reoptnode)->afterdualvarbounds != NULL )
   {
      assert((*reoptnode)->afterdualvarssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->afterdualvarbounds, (*reoptnode)->afterdualvarssize);
      (*reoptnode)->afterdualvarbounds = NULL;
   }

   /* free afterdual-variables */
   if ((*reoptnode)->afterdualvars != NULL )
   {
      assert((*reoptnode)->afterdualvarssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->afterdualvars, (*reoptnode)->afterdualvarssize);
      (*reoptnode)->afterdualvars = NULL;
   }

   (*reoptnode)->afterdualvarssize = 0;

   BMSfreeBlockMemory(blkmem, reoptnode);
   (*reoptnode) = NULL;

   return SCIP_OKAY;
}

/** reset the given reoptimization node */
static
SCIP_RETCODE reoptnodeReset(
   SCIP_REOPTNODE*       reoptnode,          /**< reoptimization node */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reoptnode != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* remove and delete all constraints */
   if( reoptnode->nconss > 0 )
   {
      int c;

      assert(reoptnode->conss != NULL);
      assert(reoptnode->consssize > 0);

      for( c = 0; c < reoptnode->nconss; c++ )
      {
         if( !reoptnode->conss[c]->linear )
         {
            assert(reoptnode->conss[c]->boundtypes != NULL);
            BMSfreeBlockMemoryArray(blkmem, &reoptnode->conss[c]->boundtypes, reoptnode->conss[c]->varssize);
         }
         BMSfreeBlockMemoryArray(blkmem, &reoptnode->conss[c]->vals, reoptnode->conss[c]->varssize);
         BMSfreeBlockMemoryArray(blkmem, &reoptnode->conss[c]->vars, reoptnode->conss[c]->varssize);
         BMSfreeBlockMemory(blkmem, &reoptnode->conss[c]); /*lint !e866 */
      }
      reoptnode->nconss = 0;
   }

   /* remove all children */
   if( reoptnode->childids != NULL )
      reoptnode->nchilds = 0;

   /* delete dual constraint */
   if( reoptnode->dualredscur != NULL )
   {
      assert(reoptnode->dualredscur->varssize > 0);
      if( !reoptnode->dualredscur->linear )
      {
         assert(reoptnode->dualredscur->boundtypes != NULL);
         BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualredscur->boundtypes, reoptnode->dualredscur->varssize);
      }
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualredscur->vals, reoptnode->dualredscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualredscur->vars, reoptnode->dualredscur->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualredscur);
      reoptnode->dualredscur = NULL;
   }

   /* delete dual constraint */
   if( reoptnode->dualredsnex != NULL )
   {
      assert(reoptnode->dualredsnex->varssize > 0);
      if( !reoptnode->dualredsnex->linear )
      {
         assert(reoptnode->dualredsnex->boundtypes != NULL);
         BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualredsnex->boundtypes, reoptnode->dualredsnex->varssize);
      }
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualredsnex->vals, reoptnode->dualredsnex->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualredsnex->vars, reoptnode->dualredsnex->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualredsnex);
      reoptnode->dualredsnex = NULL;
   }

   reoptnode->parentID = 0;
   reoptnode->nvars = 0;
   reoptnode->nafterdualvars = 0;
   reoptnode->dualreds = FALSE;
   reoptnode->reopttype = (unsigned int)SCIP_REOPTTYPE_NONE;
   reoptnode->lowerbound = -SCIPsetInfinity(set);

   return SCIP_OKAY;
}

/** delete the node stored at position @p nodeID of the reoptimization tree */
static
SCIP_RETCODE reopttreeDeleteNode(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          id,                 /**< id of a node */
   SCIP_Bool             softreset           /**< delete at the end of the solving process */
   )
{
   assert(reopttree != NULL );
   assert(id < reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[id] != NULL );

   if( softreset )
   {
      SCIP_CALL( reoptnodeReset(reopttree->reoptnodes[id], set, blkmem) );
   }
   else
   {
      SCIP_CALL( reoptnodeDelete(&reopttree->reoptnodes[id], blkmem) );
   }

   assert(softreset || reopttree->reoptnodes[id] == NULL);
   assert(reopttree->reoptnodes[id] == NULL || reopttree->reoptnodes[id]->conss == NULL || reopttree->reoptnodes[id]->nconss == 0);
   assert(reopttree->reoptnodes[id] == NULL || reopttree->reoptnodes[id]->childids == NULL || reopttree->reoptnodes[id]->nchilds == 0);

   --reopttree->nreoptnodes;

   return SCIP_OKAY;
}

/** constructor of the solution tree */
static
SCIP_RETCODE createSolTree(
   SCIP_SOLTREE*         soltree,            /**< solution tree */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int s;

   assert(soltree != NULL);

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &soltree->sols, DEFAULT_MEM_RUN) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &soltree->nsols, DEFAULT_MEM_RUN) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &soltree->solssize, DEFAULT_MEM_RUN) );

   for( s = 0; s < DEFAULT_MEM_RUN; s++ )
   {
      soltree->nsols[s] = 0;
      soltree->solssize[s] = 0;
      soltree->sols[s] = NULL;
   }

   /* allocate the root node */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &soltree->root) );
   soltree->root->sol = NULL;
   soltree->root->value = SCIP_INVALID;
   soltree->root->updated = FALSE;
   soltree->root->father = NULL;
   soltree->root->child = NULL;
   soltree->root->sibling = NULL;

   return SCIP_OKAY;
}

/** free the given solution node */
static
SCIP_RETCODE soltreefreeNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< the primal */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SOLNODE**        solnode             /**< node within the solution tree */
   )
{
   SCIP_SOLNODE* child;
   SCIP_SOLNODE* sibling;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(primal != NULL || set->stage == SCIP_STAGE_INIT);
   assert(solnode != NULL);
   assert(blkmem != NULL);

   child = (*solnode)->child;

   /* traverse through the list and free recursive all subtree */
   while( child != NULL )
   {
      SCIP_CALL( soltreefreeNode(reopt, set, primal, blkmem, &child) );
      assert(child != NULL);

      sibling = child->sibling;
      BMSfreeBlockMemoryNull(blkmem, &child);
      child = sibling;
   }

   if( (*solnode)->sol != NULL )
   {
      assert(set->stage == SCIP_STAGE_PROBLEM);

      SCIP_CALL( SCIPsolFree(&(*solnode)->sol, blkmem, primal) );
   }

   return SCIP_OKAY;
}

/** free the solution tree */
static
SCIP_RETCODE freeSolTree(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          origprimal,         /**< the origprimal */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(reopt->soltree != NULL);
   assert(reopt->soltree->root != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* free all nodes recursive */
   SCIP_CALL( soltreefreeNode(reopt, set, origprimal, blkmem, &reopt->soltree->root) );
   BMSfreeBlockMemoryNull(blkmem, &reopt->soltree->root);

   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->sols, reopt->runsize);
   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->nsols, reopt->runsize);
   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->solssize, reopt->runsize);

   BMSfreeMemory(&reopt->soltree);

   return SCIP_OKAY;
}

/** creates and adds a solution node to the solution tree */
static
SCIP_RETCODE solnodeAddChild(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SOLNODE*         curnode,            /**< current node in the solution tree */
   SCIP_SOLNODE**        child,              /**< pointer to store the node representing the solution value */
   SCIP_VAR*             var,                /**< variable represented by this node */
   SCIP_Real             val,                /**< value the child shell represent */
   SCIP_Bool*            added               /**< TRUE iff we created a new node, i.e, we have not seen this solution so far */
   )
{
   SCIP_SOLNODE* solnode;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(curnode != NULL);
   assert(child != NULL && *child == NULL);
   assert(!SCIPsetIsInfinity(set, -val) && !SCIPsetIsInfinity(set, val));

   /* get the first node of the child node list */
   *child = curnode->child;

   /* this is the first solution in the subtree induced by the current node */
   if( *child == NULL )
   {
      assert(soltreeNInducedSols(curnode) == 0);

      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &solnode) );
      solnode->sol = NULL;
      solnode->updated = FALSE;
      solnode->father = curnode;
      solnode->child = NULL;
      solnode->sibling = NULL;
      solnode->value = val;
#ifndef NDEBUG
      assert(var != NULL);
      solnode->var = var;
#endif

      *added = TRUE;
      *child = solnode;

      curnode->child = *child;

#ifdef SCIP_MORE_DEBUG
      SCIPsetDebugMsg(set, "-> create new node %p: value=%g, sibling=%p\n", (void*) solnode, solnode->value,
            (void*) solnode->sibling);
#endif
   }
   else
   {
      /* we traverse through all children */
      while( *child != NULL )
      {
#ifdef SCIP_MORE_DEBUG
         SCIPsetDebugMsg(set, "-> check %p: father=%p, value=%g, sibling=%p\n", (void*) *child, (void*) (*child)->father,
               (*child)->value, (void*) (*child)->sibling);
#endif
         /* we found a node repesenting this solution value */
         if( SCIPsetIsEQ(set, val, (*child)->value) )
            break;

         /* we are at the end of the list */
         if( (*child)->sibling == NULL )
         {
            /* create a new solnode */
            SCIP_ALLOC( BMSallocBlockMemory(blkmem, &solnode) );
            solnode->sol = NULL;
            solnode->updated = FALSE;
            solnode->father = curnode;
            solnode->child = NULL;
            solnode->value = val;
#ifndef NDEBUG
            assert(var != NULL);
            solnode->var = var;
#endif
            *added = TRUE;

            /* we have to append the new node at the end of the list. but we have to check whether the insertion before
             * the current node would be correct. in that case, we switch the values, the child pointer, and the
             * solution
             */
            solnode->sibling = NULL;
            (*child)->sibling = solnode;

#ifdef SCIP_MORE_DEBUG
            SCIPsetDebugMsg(set, "-> create new node %p: value=%g, sibling=%p\n", (void*) solnode, solnode->value,
               (void*) solnode->sibling);
#endif
            /* the given value is lower than the current, insertion before the current node would be correct
             * in this case we do not have to change the child pointer
             */
            if( SCIPsetIsLT(set, val, (*child)->value) )
            {
#ifdef SCIP_MORE_DEBUG
               SCIPsetDebugMsg(set, "-> need to switch:\n");
               SCIPsetDebugMsg(set, "   before switching: node %p witch child=%p, sibling=%p, sol=%p, value=%g\n",
                  (void*) (*child), (void*) (*child)->child, (void*) (*child)->sibling, (void*) (*child)->sol,
                  (*child)->value);
               SCIPsetDebugMsg(set, "                     node %p witch child=%p, sibling=%p, sol=%p, value=%g\n",
                  (void*) solnode, (void*) solnode->child, (void*) solnode->sibling, (void*) solnode->sol,
                  solnode->value);
#endif
               /* switch child pointer */
               solnode->child = (*child)->child;
               (*child)->child = NULL;

               /* switch solution values */
               solnode->value = (*child)->value;
               (*child)->value = val;
               assert(SCIPsetIsLT(set, (*child)->value, solnode->value));

               /* switch solution pointer */
               solnode->sol = (*child)->sol;
               (*child)->sol = NULL;
#ifdef SCIP_MORE_DEBUG
               SCIPsetDebugMsg(set, "    after switching: node %p witch child=%p, sibling=%p, sol=%p, value=%g\n",
                  (void*) (*child), (void*) (*child)->child, (void*) (*child)->sibling, (void*) (*child)->sol,
                  (*child)->value);
               SCIPsetDebugMsg(set, "                     node %p witch child=%p, sibling=%p, sol=%p, value=%g\n",
                  (void*) solnode, (void*) solnode->child, (void*) solnode->sibling, (void*) solnode->sol,
                  solnode->value);
#endif
            }
            /* set the child pointer to the new created solnode */
            else
               (*child) = solnode;

            break;
         }

         /* the next sibling represents a solution value of larger size.
          * we insert a new node between the current child and the next sibling.
          */
         if( SCIPsetIsLT(set, val, (*child)->sibling->value) )
         {
            /* create a new solnode that points to the sibling of the current child */
            SCIP_ALLOC( BMSallocBlockMemory(blkmem, &solnode) );
            solnode->sol = NULL;
            solnode->updated = FALSE;
            solnode->father = curnode;
            solnode->child = NULL;
            solnode->sibling = (*child)->sibling;
            solnode->value = val;
#ifndef NDEBUG
            assert(var != NULL);
            solnode->var = var;
#endif
            *added = TRUE;

            /* change the poiter of the next sibling to the new node */
            (*child)->sibling = solnode;

            *child = solnode;
#ifdef SCIP_MORE_DEBUG
            SCIPsetDebugMsg(set, "-> create new node %p: value=%g, sibling=%p\n", (void*) solnode, solnode->value,
                  (void*) solnode->sibling);
#endif
            break;
         }

         /* go to the next sibling */
         *child = (*child)->sibling;
      }

#ifdef SCIP_DEBUG
      /* check whether the insert was correct and the list is increasing */
      solnode = curnode->child;
      assert(solnode != NULL);

      while( solnode->sibling != NULL )
      {
         assert(SCIPsetIsLT(set, solnode->value, solnode->sibling->value));
         solnode = solnode->sibling;
      }
#endif
   }
   return SCIP_OKAY;
}

/** add a solution to the solution tree */
static
SCIP_RETCODE soltreeAddSol(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal,         /**< orig primal */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            vars,               /**< array of original variables */
   SCIP_SOL*             sol,                /**< solution to add */
   SCIP_SOLNODE**        solnode,            /**< current solution node */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             bestsol,            /**< is the solution an optimal (best found) solution */
   SCIP_Bool*            added               /**< pointer to store the result */
   )
{
   SCIP_SOLNODE* cursolnode;
   SCIP_Bool purelp;
   int varid;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(origprimal != NULL);
   assert(blkmem != NULL);
   assert(vars != NULL);
   assert(sol != NULL);
   assert(solnode != NULL);

   cursolnode = reopt->soltree->root;
   *added = FALSE;
   purelp = TRUE;

   if( set->reopt_savesols > 0 )
   {
#ifdef MORE_DEBUG
      SCIPsetDebugMsg(set, "try to add solution found by <%s>\n", (SCIPsolGetHeur(sol) == NULL ?
            "relaxation" : SCIPheurGetName(SCIPsolGetHeur(sol))));
#endif

      for( varid = 0; varid < nvars; varid++ )
      {
         if( SCIPvarGetType(vars[varid]) != SCIP_VARTYPE_CONTINUOUS )
         {
            SCIP_SOLNODE* child;

            purelp = FALSE;
            child = NULL;
            SCIP_CALL( solnodeAddChild(set, blkmem, cursolnode, &child, vars[varid],
                  SCIPsolGetVal(sol, set, stat, vars[varid]), added) );
            assert(child != NULL);
            cursolnode = child;
         }
      }

      /* the solution was added or is an optimal solution */
      if( (*added || bestsol) && !purelp )
      {
         SCIP_SOL* copysol;

         assert(cursolnode->child == NULL);

         if( *added )
         {
            SCIP_CALL( SCIPsolCopy(&copysol, blkmem, set, stat, origprimal, sol) );
            cursolnode->sol = copysol;
         }
         else
            /* this is a pseudo add; we do not want to save this solution more than once, but we will link this solution
             * to the solution storage of this round
             */
            (*added) = TRUE;

         if( bestsol )
         {
            assert(reopt->prevbestsols != NULL);
            assert(cursolnode->sol != NULL);

            reopt->prevbestsols[reopt->run-1] = cursolnode->sol;
         }

         (*solnode) = cursolnode;
      }
   }

   return SCIP_OKAY;
}

/** reset all marks 'updated' to FALSE */
static
void soltreeResetMarks(
   SCIP_SOLNODE*         node                /**< node within the solution tree */
   )
{
   assert(node != NULL);

   if( node->child != NULL )
   {
      SCIP_SOLNODE* child;

      /* the node is no leaf */
      assert(node->sol == NULL);
      assert(!node->updated);

      child = node->child;

      /* traverse through the list of siblings */
      while( child != NULL )
      {
         soltreeResetMarks(child);
         child = child->sibling;
      }
   }
   else
   {
      /* the node is a leaf */
      assert(node->father != NULL);
      assert(node->sol != NULL);
      node->updated = FALSE;
   }
}

/** allocate memory for a node within the reoptimization tree */
static
SCIP_RETCODE createReoptnode(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          id                  /**< id of the node to create */
   )
{
   assert(reopttree != NULL );
   assert(id < reopttree->reoptnodessize);

   SCIPsetDebugMsg(set, "create a reoptnode at ID %u\n", id);

   if( reopttree->reoptnodes[id] == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopttree->reoptnodes[id]) ); /*lint !e866*/

      reopttree->reoptnodes[id]->conss = NULL;
      reopttree->reoptnodes[id]->nconss = 0;
      reopttree->reoptnodes[id]->consssize = 0;
      reopttree->reoptnodes[id]->childids = NULL;
      reopttree->reoptnodes[id]->allocchildmem = 0;
      reopttree->reoptnodes[id]->nchilds = 0;
      reopttree->reoptnodes[id]->nvars = 0;
      reopttree->reoptnodes[id]->nafterdualvars = 0;
      reopttree->reoptnodes[id]->parentID = 0;
      reopttree->reoptnodes[id]->dualreds = FALSE;
      reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_NONE;
      reopttree->reoptnodes[id]->varssize = 0;
      reopttree->reoptnodes[id]->afterdualvarssize = 0;
      reopttree->reoptnodes[id]->vars = NULL;
      reopttree->reoptnodes[id]->varbounds = NULL;
      reopttree->reoptnodes[id]->varboundtypes = NULL;
      reopttree->reoptnodes[id]->afterdualvars = NULL;
      reopttree->reoptnodes[id]->afterdualvarbounds = NULL;
      reopttree->reoptnodes[id]->afterdualvarboundtypes = NULL;
      reopttree->reoptnodes[id]->dualredscur = NULL;
      reopttree->reoptnodes[id]->dualredsnex = NULL;
      reopttree->reoptnodes[id]->lowerbound = -SCIPsetInfinity(set);
   }
   else
   {
      assert(reopttree->reoptnodes[id]->nvars == 0);
      assert(reopttree->reoptnodes[id]->nafterdualvars == 0);
      reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_NONE;
      reopttree->reoptnodes[id]->lowerbound = -SCIPsetInfinity(set);
   }

   /* increase the counter */
   ++reopttree->nreoptnodes;

   assert(reopttree->nreoptnodes + SCIPqueueNElems(reopttree->openids) == (int)reopttree->reoptnodessize);

   return SCIP_OKAY;
}

/** constructor of the reoptimization tree */
static
SCIP_RETCODE createReopttree(
   SCIP_REOPTTREE*       reopttree,          /**< pointer to the reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   unsigned int id;

   assert(reopttree != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* allocate memory */
   reopttree->reoptnodessize = DEFAULT_MEM_NODES;
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopttree->reoptnodes, reopttree->reoptnodessize) );

   /* initialize the queue of open IDs */
   SCIP_CALL( SCIPqueueCreate(&reopttree->openids, (int)reopttree->reoptnodessize, 2.0) );

   /* fill the queue, but reserve the 0 for the root */
   for( id = 1; id < reopttree->reoptnodessize; id++ )
   {
      reopttree->reoptnodes[id] = NULL;
      SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void*) (size_t) id) ); /*lint !e571*/
   }
   assert(SCIPqueueNElems(reopttree->openids) == (int)(reopttree->reoptnodessize)-1);

   reopttree->nreoptnodes = 0;
   reopttree->ntotalfeasnodes = 0;
   reopttree->nfeasnodes = 0;
   reopttree->ninfnodes = 0;
   reopttree->ntotalinfnodes= 0;
   reopttree->nprunednodes = 0;
   reopttree->ntotalprunednodes= 0;
   reopttree->ncutoffreoptnodes = 0;
   reopttree->ntotalcutoffreoptnodes = 0;

   /* initialize the root node */
   reopttree->reoptnodes[0] = NULL;
   SCIP_CALL( createReoptnode(reopttree, set, blkmem, 0) );

   return SCIP_OKAY;
}

/** clears the reopttree, e.g., to restart and solve the next problem from scratch */
static
SCIP_RETCODE clearReoptnodes(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             softreset           /**< delete nodes before exit the solving process */
   )
{
   unsigned int id;

   assert(reopttree != NULL );

   /* clear queue with open IDs */
   SCIPqueueClear(reopttree->openids);
   assert(SCIPqueueNElems(reopttree->openids) == 0);

   /* delete all data about nodes */
   for( id = 0; id < reopttree->reoptnodessize; id++ )
   {
      if( reopttree->reoptnodes[id] != NULL )
      {
         SCIP_CALL( reopttreeDeleteNode(reopttree, set, blkmem, id, softreset) );
         assert(reopttree->reoptnodes[id] == NULL || reopttree->reoptnodes[id]->nvars == 0);
      }

      if( id > 0 )
      {
         SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void* ) (size_t ) id) ); /*lint !e571*/
      }
   }
   assert(SCIPqueueNElems(reopttree->openids) == (int)(reopttree->reoptnodessize)-1);

   reopttree->nreoptnodes = 0;

   return SCIP_OKAY;
}

/** free the reoptimization tree */
static
SCIP_RETCODE freeReoptTree(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree data */
   SCIP_SET*             set,                /**< global SCIP settings  */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopttree != NULL);
   assert(blkmem != NULL);

   /* free nodes */
   SCIP_CALL( clearReoptnodes(reopttree, set, blkmem, FALSE) );

   /* free the data */
   BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes, reopttree->reoptnodessize);
   SCIPqueueFree(&reopttree->openids);

   /* free the tree itself */
   BMSfreeMemory(&reopttree);

   return SCIP_OKAY;
}

/** check memory for the constraint to handle bound changes based on dual information */
static
SCIP_RETCODE checkMemDualCons(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   size                /**< size which need to be allocated */
   )
{
   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(size > 0);

   if( reopt->dualreds == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopt->dualreds) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualreds->vars, size) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualreds->vals, size) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualreds->boundtypes, size) );
      reopt->dualreds->varssize = size;
      reopt->dualreds->nvars = 0;
   }
   else if( reopt->dualreds->varssize < size )
   {
      int newsize = SCIPsetCalcMemGrowSize(set, size+1);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->dualreds->vars, reopt->dualreds->varssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->dualreds->vals, reopt->dualreds->varssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->dualreds->boundtypes, reopt->dualreds->varssize, newsize) );
      reopt->dualreds->varssize = newsize;
   }

   return SCIP_OKAY;
}

/** check the memory to store global constraints */
static
SCIP_RETCODE checkMemGlbCons(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   mem                 /**< memory which has to be allocated */
   )
{
   int c;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(mem > 0);

   if( mem > 0 )
   {
      if( reopt->glbconss == NULL )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->glbconss, mem) );
         reopt->nglbconss = 0;
         reopt->allocmemglbconss = mem;

         for( c = 0; c < reopt->allocmemglbconss; c++ )
            reopt->glbconss[c] = NULL;

      }
      else if( reopt->allocmemglbconss < mem )
      {
         int newsize = SCIPsetCalcMemGrowSize(set, mem+1);

         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->glbconss, reopt->allocmemglbconss, newsize) );

         for( c = reopt->allocmemglbconss; c < newsize; c++ )
            reopt->glbconss[c] = NULL;

         reopt->allocmemglbconss = newsize;
      }
   }

   return SCIP_OKAY;
}

/** reactivate globally valid constraints that were deactivated and necessary to ensure correctness */
static
SCIP_RETCODE cleanActiveConss(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int nentries;
   int i;

   assert(reopt != NULL);
   assert(reopt->activeconss != NULL);

   nentries = SCIPhashmapGetNEntries(reopt->activeconss);

   /* loop over all entries of the hashmap and reactivate deactivated constraints */
   for( i = 0; i < nentries; i++ )
   {
      SCIP_CONS* cons;
      SCIP_HASHMAPENTRY* entry = SCIPhashmapGetEntry(reopt->activeconss, i);

      if( entry == NULL )
         continue;

      cons = (SCIP_CONS*)SCIPhashmapEntryGetImage(entry);
      assert(cons != NULL);

      SCIP_CALL( SCIPreleaseCons(set->scip, &cons) );
   }

   return SCIP_OKAY;
}

/** update the bound changes made by constraint propagations during current iteration; stop saving the bound changes if
 *  we reach a branching decision based on a dual information.
 */
static
SCIP_RETCODE updateConstraintPropagation(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int          id,                 /**< id of the node */
   SCIP_Bool*            transintoorig       /**< transform variables into originals */
   )
{
   int nvars;
   int nconsprops;
   int naddedbndchgs;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(0 < id && id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL );

   /* get the number of all stored constraint propagations */
   SCIPnodeGetNDomchg(node, NULL, &nconsprops, NULL);
   nvars = reopt->reopttree->reoptnodes[id]->nvars;

   if( nconsprops > 0 )
   {
      /* check the memory */
      SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[id], set, blkmem, nvars + nconsprops, 0, 0) );

      SCIPnodeGetConsProps(node,
            &reopt->reopttree->reoptnodes[id]->vars[nvars],
            &reopt->reopttree->reoptnodes[id]->varbounds[nvars],
            &reopt->reopttree->reoptnodes[id]->varboundtypes[nvars],
            &naddedbndchgs,
            reopt->reopttree->reoptnodes[id]->varssize-nvars);

      assert(nvars + naddedbndchgs <= reopt->reopttree->reoptnodes[id]->varssize);

      reopt->reopttree->reoptnodes[id]->nvars += naddedbndchgs;

      *transintoorig = TRUE;
   }

   return SCIP_OKAY;
}

/** save bound changes made after the first bound change based on dual information, e.g., mode by strong branching
 *
 *  This method can be used during reoptimization. If we want to reconstruct a node containing dual bound changes we
 *  have to split the node into the original one and at least one node representing the pruned part. All bound changes,
 *  i.e., (constraint) propagation, made after the first bound change based on dual information are still valid for
 *  the original node after changing the objective function. thus, we can store them for the following iterations.
 *
 *  It should be noted, that these bound changes will be found by (constraint) propagation methods anyway after changing
 *  the objective function. do not saving these information and find them again might be useful for conflict analysis.
 */
static
SCIP_RETCODE saveAfterDualBranchings(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int          id,                 /**< id of the node */
   SCIP_Bool*            transintoorig       /**< transform variables into originals */
   )
{
   int nbranchvars;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(0 < id && id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL );

   nbranchvars = 0;

   /* allocate memory */
   if (reopt->reopttree->reoptnodes[id]->afterdualvarssize == 0)
   {
      assert(reopt->reopttree->reoptnodes[id]->afterdualvars == NULL );
      assert(reopt->reopttree->reoptnodes[id]->afterdualvarbounds == NULL );
      assert(reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes == NULL );

      /* allocate block memory for node information */
      reopt->reopttree->reoptnodes[id]->afterdualvarssize = DEFAULT_MEM_VARAFTERDUAL;
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvars), \
            reopt->reopttree->reoptnodes[id]->afterdualvarssize) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvarbounds), \
            reopt->reopttree->reoptnodes[id]->afterdualvarssize) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes), \
            reopt->reopttree->reoptnodes[id]->afterdualvarssize) );
   }

   assert(reopt->reopttree->reoptnodes[id]->afterdualvarssize > 0);
   assert(reopt->reopttree->reoptnodes[id]->nafterdualvars >= 0);

   SCIPnodeGetBdChgsAfterDual(node,
         reopt->reopttree->reoptnodes[id]->afterdualvars,
         reopt->reopttree->reoptnodes[id]->afterdualvarbounds,
         reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes,
         reopt->reopttree->reoptnodes[id]->nafterdualvars,
         &nbranchvars,
         reopt->reopttree->reoptnodes[id]->afterdualvarssize);

   if( nbranchvars > reopt->reopttree->reoptnodes[id]->afterdualvarssize )
   {
      int newsize = SCIPsetCalcMemGrowSize(set, nbranchvars+1);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvars), \
            reopt->reopttree->reoptnodes[id]->afterdualvarssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvarbounds), \
            reopt->reopttree->reoptnodes[id]->afterdualvarssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes), \
            reopt->reopttree->reoptnodes[id]->afterdualvarssize, newsize) );
      reopt->reopttree->reoptnodes[id]->afterdualvarssize = newsize;

      SCIPnodeGetBdChgsAfterDual(node,
            reopt->reopttree->reoptnodes[id]->afterdualvars,
            reopt->reopttree->reoptnodes[id]->afterdualvarbounds,
            reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes,
            reopt->reopttree->reoptnodes[id]->nafterdualvars,
            &nbranchvars,
            reopt->reopttree->reoptnodes[id]->afterdualvarssize);
   }

   /* the stored variables of this node need to be transformed into the original space */
   if( nbranchvars > 0 )
      *transintoorig = TRUE;

   SCIPsetDebugMsg(set, " -> save %d bound changes after dual reductions\n", nbranchvars);

   assert(nbranchvars <= reopt->reopttree->reoptnodes[id]->afterdualvarssize); /* this should be the case */

   reopt->reopttree->reoptnodes[id]->nafterdualvars = nbranchvars;

   return SCIP_OKAY;
}

/** store cuts that are active in the current LP */
static
SCIP_RETCODE storeCuts(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LP*              lp,                 /**< current LP */
   unsigned int          id                  /**< id in the reopttree */
   )
{
   SCIP_ROW** lprows;
   int nlprows;
   int r;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(blkmem != NULL);

   lprows = SCIPlpGetRows(lp);
   nlprows = SCIPlpGetNRows(lp);

   for( r = 0; r < nlprows; r++ )
   {
      /* we can break if we reach the first row that is not part of the current LP */
      if( SCIProwGetLPPos(lprows[r]) == -1 )
         break;

      /* currently we only want to store cuts generated by a seperator */
      if( SCIProwGetOrigintype(lprows[r]) == SCIP_ROWORIGINTYPE_SEPA && SCIProwGetAge(lprows[r]) <= set->reopt_maxcutage )
      {
         SCIP_VAR** cutvars;
         SCIP_COL** cols;
         SCIP_Real* cutvals;
         SCIP_Real lhs;
         SCIP_Real rhs;
         int ncutvars;
         int c;

         ncutvars = SCIProwGetNLPNonz(lprows[r]);
         lhs = SCIProwGetLhs(lprows[r]);
         rhs = SCIProwGetRhs(lprows[r]);

         /* subtract row constant */
         if( !SCIPsetIsInfinity(set, -lhs) )
            lhs -= SCIProwGetConstant(lprows[r]);
         if( !SCIPsetIsInfinity(set, rhs) )
            rhs -= SCIProwGetConstant(lprows[r]);

         cutvals = SCIProwGetVals(lprows[r]);
         cols = SCIProwGetCols(lprows[r]);

         SCIP_CALL( SCIPsetAllocBufferArray(set, &cutvars, ncutvars) );

         for( c = 0; c < ncutvars; c++ )
         {
            SCIP_Real constant;
            SCIP_Real scalar;

            cutvars[c] = SCIPcolGetVar(cols[c]);
            assert(cutvars[c] != NULL);

            constant = 0.0;
            scalar = 1.0;

            SCIP_CALL( SCIPvarGetOrigvarSum(&cutvars[c], &scalar, &constant) );
            assert(cutvars[c] != NULL);
            assert(!SCIPsetIsZero(set, scalar));

            /* subtract constant from sides */
            if( !SCIPsetIsZero(set, constant) && !SCIPsetIsInfinity(set, -lhs) )
               lhs -= constant;
            if( !SCIPsetIsZero(set, constant) && !SCIPsetIsInfinity(set, rhs) )
               rhs -= constant;

            cutvals[c] = cutvals[c]/scalar;
         }

         /* add cut as a linear constraint */
         SCIP_CALL( SCIPreoptnodeAddCons(reopt->reopttree->reoptnodes[id], set, blkmem, cutvars, cutvals, NULL,
               lhs, rhs, ncutvars, REOPT_CONSTYPE_CUT, TRUE) );

         SCIPsetFreeBufferArray(set, &cutvars);
      }
   }

   return SCIP_OKAY;
}

/** transform variable and bounds back to the original space */
static
SCIP_RETCODE transformIntoOrig(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   unsigned int          id                  /**< id of the node */
   )
{
   int varnr;

   assert(reopt != NULL );
   assert(0 < id && id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL );

   /* transform branching variables and bound changes applied before the first dual reduction */
   for( varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nvars; varnr++ )
   {
      SCIP_Real constant = 0.0;
      SCIP_Real scalar = 1.0;

      if( !SCIPvarIsOriginal(reopt->reopttree->reoptnodes[id]->vars[varnr]) )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&reopt->reopttree->reoptnodes[id]->vars[varnr], &scalar, &constant)) ;
         reopt->reopttree->reoptnodes[id]->varbounds[varnr] = (reopt->reopttree->reoptnodes[id]->varbounds[varnr] - constant) / scalar;
      }
      assert(SCIPvarIsOriginal(reopt->reopttree->reoptnodes[id]->vars[varnr]));
   }

   /* transform bound changes affected by dual reduction */
   for( varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nafterdualvars; varnr++ )
   {
      SCIP_Real constant = 0.0;
      SCIP_Real scalar = 1.0;

      if( !SCIPvarIsOriginal(reopt->reopttree->reoptnodes[id]->afterdualvars[varnr]) )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&reopt->reopttree->reoptnodes[id]->afterdualvars[varnr], &scalar, &constant)) ;
         reopt->reopttree->reoptnodes[id]->afterdualvarbounds[varnr]
            = (reopt->reopttree->reoptnodes[id]->afterdualvarbounds[varnr] - constant) / scalar;
      }
      assert(SCIPvarIsOriginal(reopt->reopttree->reoptnodes[id]->afterdualvars[varnr]));
   }

   return SCIP_OKAY;
}

/** search the next node along the root path that was saved by reoptimization */
static
SCIP_RETCODE getLastSavedNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_NODE**           parent,             /**< parent node within the search tree */
   unsigned int*         parentid,           /**< id of the parent node */
   int*                  nbndchgs            /**< number of bound changes */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(reopt->reopttree->reoptnodes != NULL);

   (*nbndchgs) = 0;
   (*parent) = node;

   /* look for a saved parent along the root-path */
   while( SCIPnodeGetDepth(*parent) != 0 )
   {
      int nbranchings = 0;
      int nconsprop = 0;

      if( set->reopt_saveconsprop )
         SCIPnodeGetNDomchg((*parent), &nbranchings, &nconsprop, NULL);
      else
         SCIPnodeGetNDomchg((*parent), &nbranchings, NULL, NULL);

      (*nbndchgs) = (*nbndchgs) + nbranchings + nconsprop;
      (*parent) = SCIPnodeGetParent(*parent);
      (*parentid) = SCIPnodeGetReoptID(*parent);

      if( SCIPnodeGetDepth(*parent) == 0)
      {
         (*parentid) = 0;
         break;
      }
      else if( SCIPnodeGetReopttype((*parent)) >= SCIP_REOPTTYPE_TRANSIT )
      {
         /* this is a special case: due to re-propagation the node could be already deleted. We need to reset reoptid
          * and reopttype and continue upto we have found the last stored node
          */
         if( reopt->reopttree->reoptnodes[*parentid] == NULL )
         {
            SCIPnodeSetReoptID(*parent, 0);
            SCIPnodeSetReopttype(*parent, SCIP_REOPTTYPE_NONE);
         }
         else
         {
            assert(reopt->reopttree->reoptnodes[*parentid] != NULL);
            assert(SCIPnodeGetReoptID((*parent)) < reopt->reopttree->reoptnodessize);
            assert((*parentid) && (*parentid) < reopt->reopttree->reoptnodessize);
            break;
         }
      }
   }

   return SCIP_OKAY;
}

/** adds the id @p childid to the array of child nodes of @p parentid */
static
SCIP_RETCODE reoptAddChild(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          parentid,           /**< id of the parent node */
   unsigned int          childid             /**< id of the child node */
   )
{
   int nchilds;

   assert(reopttree != NULL);
   assert(blkmem != NULL);
   assert(parentid < (unsigned int)reopttree->reoptnodessize);
   assert(childid < (unsigned int)reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[parentid] != NULL);

   nchilds = reopttree->reoptnodes[parentid]->nchilds;

   /* ensure that the array is large enough */
   SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[parentid], set, blkmem, 0, nchilds+1, 0) );
   assert(reopttree->reoptnodes[parentid]->allocchildmem > nchilds);

   /* add the child */
   reopttree->reoptnodes[parentid]->childids[nchilds] = childid;
   ++reopttree->reoptnodes[parentid]->nchilds;

   SCIPsetDebugMsg(set, "add ID %u as a child of ID %u.\n", childid, parentid);

   return SCIP_OKAY;
}

/** move all children to the next node (along the root path) stored in the reoptimization tree */
static
SCIP_RETCODE moveChildrenUp(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          nodeid,             /**< id of the node */
   unsigned int          parentid            /**< id of the parent node */
   )
{
   unsigned int childid;
   int varnr;
   int nvars;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(0 < nodeid && nodeid < reopt->reopttree->reoptnodessize);
   assert(parentid < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[nodeid]->childids != NULL);

   /* ensure that enough memory at the parentID is available */
   SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[parentid], set, blkmem, 0,
         reopt->reopttree->reoptnodes[parentid]->nchilds + reopt->reopttree->reoptnodes[nodeid]->nchilds, 0) );

   while( reopt->reopttree->reoptnodes[nodeid]->nchilds > 0 )
   {
      int nchilds;

      nchilds = reopt->reopttree->reoptnodes[nodeid]->nchilds;
      childid = reopt->reopttree->reoptnodes[nodeid]->childids[nchilds-1];
      assert(0 < childid && childid < reopt->reopttree->reoptnodessize);

      /* check the memory */
      SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[childid], set, blkmem,
            reopt->reopttree->reoptnodes[childid]->nvars + reopt->reopttree->reoptnodes[nodeid]->nvars, 0, 0) );
      assert(reopt->reopttree->reoptnodes[childid]->varssize >= reopt->reopttree->reoptnodes[childid]->nvars
         + reopt->reopttree->reoptnodes[nodeid]->nvars);

      /* save branching information */
      for( varnr = 0; varnr < reopt->reopttree->reoptnodes[nodeid]->nvars; varnr++ )
      {
         nvars = reopt->reopttree->reoptnodes[childid]->nvars;
         reopt->reopttree->reoptnodes[childid]->vars[nvars] = reopt->reopttree->reoptnodes[nodeid]->vars[varnr];
         reopt->reopttree->reoptnodes[childid]->varbounds[nvars] = reopt->reopttree->reoptnodes[nodeid]->varbounds[varnr];
         reopt->reopttree->reoptnodes[childid]->varboundtypes[nvars] = reopt->reopttree->reoptnodes[nodeid]->varboundtypes[varnr];
         ++reopt->reopttree->reoptnodes[childid]->nvars;
      }

      /* update the ID of the parent node */
      reopt->reopttree->reoptnodes[childid]->parentID = parentid;

      /* insert the node as a child */
      SCIP_CALL( reoptAddChild(reopt->reopttree, set, blkmem, parentid, childid) );

      /* reduce the number of child nodes by 1 */
      --reopt->reopttree->reoptnodes[nodeid]->nchilds;
   }

   return SCIP_OKAY;
}

/** delete all nodes in the subtree induced by nodeID */
static
SCIP_RETCODE deleteChildrenBelow(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          id,                 /**< id of the node */
   SCIP_Bool             delnodeitself,      /**< should the node be deleted after deleting the induced subtree? */
   SCIP_Bool             exitsolve           /**< will the solving process end after deletion */
   )
{
   assert(reopttree != NULL );
   assert(blkmem != NULL);
   assert(id < reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[id] != NULL);

   /* delete all children below */
   if( reopttree->reoptnodes[id]->childids != NULL && reopttree->reoptnodes[id]->nchilds > 0 )
   {
      SCIPsetDebugMsg(set, "-> delete subtree induced by ID %u (hard remove = %u)\n", id, exitsolve);

      while( reopttree->reoptnodes[id]->nchilds > 0 )
      {
         int nchilds;
         unsigned int childid;

         nchilds = reopttree->reoptnodes[id]->nchilds;
         childid = reopttree->reoptnodes[id]->childids[nchilds-1];
         assert(0 < childid && childid < reopttree->reoptnodessize);

         SCIP_CALL( deleteChildrenBelow(reopttree, set, blkmem, childid, TRUE, exitsolve) );

         --reopttree->reoptnodes[id]->nchilds;
      }
   }

   /* delete node data*/
   if( delnodeitself )
   {
      SCIP_CALL( reopttreeDeleteNode(reopttree, set, blkmem, id, exitsolve) );
      SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void*) (size_t) id) );
   }

   return SCIP_OKAY;
}

/** replaces a reoptimization nodes by its stored child nodes */
static
SCIP_RETCODE shrinkNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int          id,                 /**< id of the node */
   SCIP_Bool*            shrank,             /**< pointer to store if the node was shrank */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_REOPTNODE** reoptnodes;

   assert(reopt != NULL);
   assert(node != NULL);
   assert(id < reopt->reopttree->reoptnodessize);

   reoptnodes = reopt->reopttree->reoptnodes;
   assert(reoptnodes != NULL);
   assert(reoptnodes[id] != NULL);

   if( reoptnodes[id]->childids != NULL && reoptnodes[id]->nchilds > 0 )
   {
      int ndomchgs = 0;
      unsigned int parentid = 0;
      SCIP_NODE* parent = NULL;

      SCIP_CALL( getLastSavedNode(reopt, set, node, &parent, &parentid, &ndomchgs) );

      assert(parentid != id);
      assert(reoptnodes[parentid] != NULL );
      assert(reoptnodes[parentid]->childids != NULL && reoptnodes[parentid]->nchilds);

      /* check if we want move all children to the next saved node above
       * we want to shrink the path if either
       * - the maximal number of bound changes fix and the number of bound changes is
       *   less than the given threshold set->reopt_maxdiffofnodes
       * or
       * - the number is calculated dynamically and the number of bound changes
       *   is less than log2(SCIPgetNBinVars - (#vars of parent))
       * */
      if( ndomchgs <= set->reopt_maxdiffofnodes )
      {
         int c;

         SCIPsetDebugMsg(set, " -> shrink node %lld at ID %u, replaced by %d child nodes.\n", SCIPnodeGetNumber(node),
            id, reoptnodes[id]->nchilds);

         /* copy the references of child nodes to the parent*/
         SCIP_CALL( moveChildrenUp(reopt, set, blkmem, id, parentid) );

         /* delete the current node */
         c = 0;
         while( reoptnodes[parentid]->childids[c] != id )
         {
            ++c;
            assert(c < reoptnodes[parentid]->nchilds);
         }

         assert(reoptnodes[parentid]->childids[c] == id);

         /* replace the childid at position c by the last one */
         reoptnodes[parentid]->childids[c] = reoptnodes[parentid]->childids[reoptnodes[parentid]->nchilds-1];
         --reoptnodes[parentid]->nchilds;

         SCIP_CALL( reopttreeDeleteNode(reopt->reopttree, set, blkmem, id, TRUE) );
         SCIP_CALL( SCIPqueueInsert(reopt->reopttree->openids, (void*) (size_t) id) );

         *shrank = TRUE;

         /* set the reopttype to none */
         SCIPnodeSetReopttype(node, SCIP_REOPTTYPE_NONE);
      }
   }

   return SCIP_OKAY;
}

/** change all reopttypes in the subtree induced by @p nodeID */
static
SCIP_RETCODE changeReopttypeOfSubtree(
   SCIP_REOPTTREE*       reopttree,          /**< reopttree */
   unsigned int          id,                 /**< id of the node */
   SCIP_REOPTTYPE        reopttype           /**< reopttype */
   )
{
   assert(reopttree != NULL);
   assert(id < reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[id] != NULL);

   if( reopttree->reoptnodes[id]->childids != NULL && reopttree->reoptnodes[id]->nchilds > 0 )
   {
      unsigned int childid;
      int nchildids;
      int seenids = 0;

      nchildids = reopttree->reoptnodes[id]->nchilds;

      while( seenids < nchildids )
      {
         /* get childID */
         childid = reopttree->reoptnodes[id]->childids[seenids];
         assert(childid < reopttree->reoptnodessize);
         assert(reopttree->reoptnodes[childid] != NULL);

         /* change the reopttype of the node iff the node is neither infeasible nor induces an
          * infeasible subtree and if the node contains no bound changes based on dual decisions
          */
         if( reopttree->reoptnodes[childid]->reopttype != SCIP_REOPTTYPE_STRBRANCHED
            && reopttree->reoptnodes[childid]->reopttype != SCIP_REOPTTYPE_INFSUBTREE ) /*lint !e641*/
            reopttree->reoptnodes[childid]->reopttype = reopttype; /*lint !e641*/

         /* change reopttype of subtree */
         SCIP_CALL( changeReopttypeOfSubtree(reopttree, childid, reopttype) );

         ++seenids;
      }
   }

   return SCIP_OKAY;
}

/** delete the constraint handling dual information for the current iteration and replace it with the dual constraint
 *  for the next iteration
 */
static
SCIP_RETCODE reoptnodeUpdateDualConss(
   SCIP_REOPTNODE*       reoptnode,          /**< reoptimization node */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reoptnode != NULL);
   assert(blkmem != NULL);

   if( reoptnode->dualredscur != NULL )
   {
      SCIPdebugMessage("reset dual information (current run)\n");

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualredscur->boundtypes, reoptnode->dualredscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualredscur->vals, reoptnode->dualredscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualredscur->vars, reoptnode->dualredscur->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualredscur);
      reoptnode->dualredscur = NULL;
   }

   if( reoptnode->dualredsnex != NULL )
   {
      SCIPdebugMessage("set dual information of next run to current run\n");
      reoptnode->dualredscur = reoptnode->dualredsnex;
      reoptnode->dualredsnex = NULL;
   }

   reoptnode->dualreds = (reoptnode->dualredscur != NULL ? TRUE : FALSE);

   return SCIP_OKAY;
}

/** calculates a (local) similarity of a given node and returns if the subproblem should be solved from scratch */
static
SCIP_RETCODE reoptCheckLocalRestart(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_VAR**            transvars,          /**< transformed variables */
   int                   ntransvars,         /**< number of transformed variables */
   SCIP_Bool*            localrestart        /**< pointer to store if we want to restart solving the (sub)problem */
   )
{
   unsigned int id;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(transvars != NULL);

   /* node == NULL is equivalent to node == root, this case should be handled by SCIPreoptCheckReopt */
   assert(node != NULL);

   *localrestart = FALSE;

   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* set the id to -1 if the node is not part of the reoptimization tree */
   if( SCIPnodeGetDepth(node) > 0 && id == 0 )
      return SCIP_OKAY;

   if( set->reopt_objsimdelay > -1 )
   {
      SCIP_Real sim = 0.0;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real oldcoef;
      SCIP_Real newcoef;
      int v;
      int idx;

      if( id == 0 )
         reopt->nlocrestarts = 0;

      /* since the stored objective functions are already normalize the dot-product is equivalent to the similarity */
      for( v = 0; v < ntransvars; v++ )
      {
         lb = SCIPvarGetLbLocal(transvars[v]);
         ub = SCIPvarGetUbLocal(transvars[v]);

         /* skip already fixed variables */
         if( SCIPsetIsFeasLT(set, lb, ub) )
         {
            idx = SCIPvarGetProbindex(transvars[v]);
            assert(0 <= idx && idx < ntransvars);

            oldcoef = SCIPreoptGetOldObjCoef(reopt, reopt->run-1, idx);
            newcoef = SCIPreoptGetOldObjCoef(reopt, reopt->run, idx);

            sim += (oldcoef * newcoef);
         }
      }

      /* delete the stored subtree and information about bound changes
       * based on dual information */
      if( SCIPsetIsLT(set, sim, set->reopt_objsimdelay) )
      {
         /* set the flag */
         *localrestart = TRUE;

         ++reopt->nlocrestarts;
         ++reopt->ntotallocrestarts;

         /* delete the stored subtree */
         SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, id, FALSE, FALSE) );

         /* delete the stored constraints; we do this twice in a row because we want to delete both constraints */
         SCIP_CALL( reoptnodeUpdateDualConss(reopt->reopttree->reoptnodes[id], blkmem) );
         SCIP_CALL( reoptnodeUpdateDualConss(reopt->reopttree->reoptnodes[id], blkmem) );
      }

      SCIPsetDebugMsg(set, " -> local similarity: %.4f%s\n", sim, *localrestart ? " (solve subproblem from scratch)" : "");
   }

   return SCIP_OKAY;
}

/** save ancestor branching information up to the next stored node */
static
SCIP_RETCODE saveAncestorBranchings(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the branch and bound tree */
   SCIP_NODE*            parent,             /**< parent node */
   unsigned int          id,                 /**< id of the node */
   unsigned int          parentid            /**< id of the parent node */
   )
{
   int nbranchvars;

   assert(reopttree != NULL );
   assert(node != NULL );
   assert(parent != NULL );
   assert(1 <= id && id < reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[id] != NULL );
   assert(parentid < reopttree->reoptnodessize);
   assert(parentid == 0 || reopttree->reoptnodes[parentid] != NULL ); /* if the root is the next saved node, the nodedata can be NULL */

   SCIPsetDebugMsg(set, " -> save ancestor branchings\n");

   /* allocate memory */
   if (reopttree->reoptnodes[id]->varssize == 0)
   {
      assert(reopttree->reoptnodes[id]->vars == NULL );
      assert(reopttree->reoptnodes[id]->varbounds == NULL );
      assert(reopttree->reoptnodes[id]->varboundtypes == NULL );

      /* allocate memory for node information */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], set, blkmem, DEFAULT_MEM_VAR, 0, 0) );
   }

   assert(reopttree->reoptnodes[id]->varssize > 0);
   assert(reopttree->reoptnodes[id]->nvars == 0);

   SCIPnodeGetAncestorBranchingsPart(node, parent,
         reopttree->reoptnodes[id]->vars,
         reopttree->reoptnodes[id]->varbounds,
         reopttree->reoptnodes[id]->varboundtypes,
         &nbranchvars,
         reopttree->reoptnodes[id]->varssize);

   if( nbranchvars >  reopttree->reoptnodes[id]->varssize )
   {
      /* reallocate memory */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], set, blkmem, nbranchvars, 0, 0) );

      SCIPnodeGetAncestorBranchingsPart(node, parent,
            reopttree->reoptnodes[id]->vars,
            reopttree->reoptnodes[id]->varbounds,
            reopttree->reoptnodes[id]->varboundtypes,
            &nbranchvars,
            reopttree->reoptnodes[id]->varssize);
   }

   assert(nbranchvars <= reopttree->reoptnodes[id]->varssize); /* this should be the case */

   reopttree->reoptnodes[id]->nvars = nbranchvars;

   assert(nbranchvars <= reopttree->reoptnodes[id]->varssize);
   assert(reopttree->reoptnodes[id]->vars != NULL );

   return SCIP_OKAY;
}

/** transform a constraint with linear representation into reoptimization constraint data */
static
SCIP_RETCODE saveConsLinear(
   SCIP_REOPTCONSDATA*   reoptconsdata,      /**< reoptimization constraint data */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS*            cons,               /**< linear constraint that should be stored */
   SCIP_Bool*            success             /**< pointer to store the success */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool allocbuffervals;
   int v;

   assert(reoptconsdata != NULL);
   assert(cons != NULL);

   *success = FALSE;
   allocbuffervals = FALSE;
   reoptconsdata->linear = TRUE;

   vars = NULL;
   vals = NULL;
   SCIP_CALL( SCIPconsGetNVars(cons, set, &reoptconsdata->nvars, success) );
   assert(*success);

   /* allocate memory for variables and values; boundtypes are not needed */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptconsdata->vars, reoptconsdata->nvars) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptconsdata->vals, reoptconsdata->nvars) );
   reoptconsdata->varssize = reoptconsdata->nvars;

   /* only needed for bounddisjuction constraints, thus we set them to NULL to avoid compiler warnings */
   reoptconsdata->boundtypes = NULL;

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);

   /* get all variables, values, and sides */
   if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
   {
      vars = SCIPgetVarsLinear(NULL, cons);
      vals = SCIPgetValsLinear(NULL, cons);
      reoptconsdata->lhs = SCIPgetLhsLinear(NULL, cons);
      reoptconsdata->rhs = SCIPgetRhsLinear(NULL, cons);
   }
   else if( strcmp(SCIPconshdlrGetName(conshdlr), "logicor") == 0 )
   {
      vars = SCIPgetVarsLogicor(NULL, cons);

      /* initialize values to 1.0 */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, reoptconsdata->nvars) );
      allocbuffervals = TRUE;

      for( v = 0; v < reoptconsdata->nvars; v++ )
         vals[v] = 1.0;

      reoptconsdata->lhs = 1.0;
      reoptconsdata->rhs = SCIPsetInfinity(set);
   }
   else if( strcmp(SCIPconshdlrGetName(conshdlr), "setppc") == 0 )
   {
      vars = SCIPgetVarsSetppc(NULL, cons);

      /* initialize values to 1.0 */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, reoptconsdata->nvars) );
      allocbuffervals = TRUE;

      for( v = 0; v < reoptconsdata->nvars; v++ )
         vals[v] = 1.0;

      switch( SCIPgetTypeSetppc(NULL, cons) ) {
      case SCIP_SETPPCTYPE_PARTITIONING:
         reoptconsdata->lhs = 1.0;
         reoptconsdata->rhs = 1.0;
         break;
      case SCIP_SETPPCTYPE_PACKING:
         reoptconsdata->lhs = -SCIPsetInfinity(set);
         reoptconsdata->rhs = 1.0;
         break;
      case SCIP_SETPPCTYPE_COVERING:
         reoptconsdata->lhs = 1.0;
         reoptconsdata->rhs = SCIPsetInfinity(set);
         break;
      default:
         *success = FALSE;
         return SCIP_OKAY;
      }
   }
   else
   {
      assert(strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 || strcmp(SCIPconshdlrGetName(conshdlr), "logicor") == 0
         || strcmp(SCIPconshdlrGetName(conshdlr), "setppc") == 0);

      SCIPerrorMessage("Cannot handle constraints of type <%s> in saveConsLinear.\n", SCIPconshdlrGetName(conshdlr));
      return SCIP_INVALIDDATA;
   }
   assert(vars != NULL);
   assert(vals != NULL);

   /* transform all variables into the original space */
   for( v = 0; v < reoptconsdata->nvars; v++ )
   {
      SCIP_Real constant = 0.0;
      SCIP_Real scalar = 1.0;

      assert(vars[v] != NULL);

      reoptconsdata->vars[v] = vars[v];
      reoptconsdata->vals[v] = vals[v];

      SCIP_CALL( SCIPvarGetOrigvarSum(&reoptconsdata->vars[v], &scalar, &constant) );
      assert(!SCIPsetIsZero(set, scalar));

      assert(!SCIPsetIsInfinity(set, REALABS(reoptconsdata->vals[v])));
      reoptconsdata->vals[v] *= scalar;

      if( !SCIPsetIsZero(set, constant) && !SCIPsetIsInfinity(set, -reoptconsdata->lhs) )
         reoptconsdata->lhs -= constant;
      if( !SCIPsetIsZero(set, constant) && !SCIPsetIsInfinity(set, reoptconsdata->rhs) )
         reoptconsdata->rhs -= constant;
   }

   /* free buffer if needed */
   if( allocbuffervals )
   {
      SCIPsetFreeBufferArray(set, &vals);
   }

   return SCIP_OKAY;
}

/** transform a bounddisjunction constraint into reoptimization constraint data */
static
SCIP_RETCODE saveConsBounddisjuction(
   SCIP_REOPTCONSDATA*   reoptconsdata,      /**< reoptimization constraint data */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS*            cons,               /**< bounddisjuction constraint that should be stored */
   SCIP_Bool*            success             /**< pointer to store the success */
   )
{
   SCIP_VAR** vars;
   SCIP_CONSHDLR* conshdlr;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Real* bounds;
   int v;

   assert(reoptconsdata != NULL);
   assert(cons != NULL);

   *success = FALSE;
   reoptconsdata->linear = FALSE;

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "bounddisjunction") == 0);

   if( strcmp(SCIPconshdlrGetName(conshdlr), "bounddisjunction") != 0 )
   {
      SCIPerrorMessage("Cannot handle constraints of type <%s> in saveConsBounddisjuction.\n",
         SCIPconshdlrGetName(conshdlr));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPconsGetNVars(cons, set, &reoptconsdata->nvars, success) );
   assert(*success);

   /* allocate memory for variables and values; boundtypes are not needed */
   vars = SCIPgetVarsBounddisjunction(NULL, cons);
   bounds = SCIPgetBoundsBounddisjunction(NULL, cons);
   boundtypes = SCIPgetBoundtypesBounddisjunction(NULL, cons);
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptconsdata->vars, vars, reoptconsdata->nvars) );
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptconsdata->vals, bounds, reoptconsdata->nvars) );
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptconsdata->boundtypes, boundtypes, reoptconsdata->nvars) );
   reoptconsdata->varssize = reoptconsdata->nvars;
   reoptconsdata->lhs = SCIP_UNKNOWN;
   reoptconsdata->rhs = SCIP_UNKNOWN;

   /* transform all variables into the original space */
   for( v = 0; v < reoptconsdata->nvars; v++ )
   {
      SCIP_Real constant = 0.0;
      SCIP_Real scalar = 1.0;

      assert(reoptconsdata->vars[v] != NULL);

      SCIP_CALL( SCIPvarGetOrigvarSum(&reoptconsdata->vars[v], &scalar, &constant) );
      assert(!SCIPsetIsZero(set, scalar));

      assert(!SCIPsetIsInfinity(set, REALABS(reoptconsdata->vals[v])));
      reoptconsdata->vals[v] -= constant;
      reoptconsdata->vals[v] *= scalar;

      /* due to multipling with a negative scalar the relation need to be changed */
      if( SCIPsetIsNegative(set, scalar) )
         reoptconsdata->boundtypes[v] = (SCIP_BOUNDTYPE)(SCIP_BOUNDTYPE_UPPER - reoptconsdata->boundtypes[v]); /*lint !e656*/
   }

   return SCIP_OKAY;
}

/** save additional all constraints that were additionally added to @p node */
static
SCIP_RETCODE saveLocalConssData(
   SCIP_REOPTTREE*       reopttree,          /**< reopttree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the branch and bound tree */
   unsigned int          id                  /**< id of the node*/
   )
{
   SCIP_CONS** addedcons;
   int naddedconss;
   int addedconsssize;
   int nconss;
   int c;

   assert(node != NULL );
   assert(reopttree != NULL);
   assert(id < reopttree->reoptnodessize);

   /* save the added pseudo-constraint */
   if( SCIPnodeGetNAddedConss(node) > 0 )
   {
      addedconsssize = SCIPnodeGetNAddedConss(node);

      SCIPsetDebugMsg(set, " -> save %d locally added constraints\n", addedconsssize);

      /* get memory */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &addedcons, addedconsssize) );
      SCIPnodeGetAddedConss(node, addedcons, &naddedconss, addedconsssize);

      nconss = reopttree->reoptnodes[id]->nconss;

      /* check memory for added constraints */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], set, blkmem, 0, 0, naddedconss) );

      /* since the first nconss are already stored in the data structure, we skip them */
      for( c = nconss; c < naddedconss; c++ )
      {
         SCIP_CONSHDLR* conshdlr;
         SCIP_Bool islinear;
         SCIP_Bool success;

         conshdlr = SCIPconsGetHdlr(addedcons[c]);

         /* check whether the constraint has a linear representation */
         islinear = (strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0
               || strcmp(SCIPconshdlrGetName(conshdlr), "logicor") == 0
               || strcmp(SCIPconshdlrGetName(conshdlr), "setppc") == 0);

         SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopttree->reoptnodes[id]->conss[c]) ); /*lint !e866*/

         success = FALSE;

         /* the constraint has a linear representation */
         if( islinear )
         {
            SCIP_CALL( saveConsLinear(reopttree->reoptnodes[id]->conss[c], set, blkmem, addedcons[c], &success) );
            assert(success);

            /* increase the counter for added constraints */
            ++reopttree->reoptnodes[id]->nconss;
         }
         else
         {
            assert(strcmp(SCIPconshdlrGetName(conshdlr), "bounddisjunction") == 0);
            SCIP_CALL( saveConsBounddisjuction(reopttree->reoptnodes[id]->conss[c], set, blkmem, addedcons[c], &success) );
            assert(success);

            /* increase the counter for added constraints */
            ++reopttree->reoptnodes[id]->nconss;
         }
         assert(reopttree->reoptnodes[id]->conss[c]->nvars > 0);

         if( strcmp("reopt_inf", SCIPconsGetName(addedcons[c])) == 0 )
            reopttree->reoptnodes[id]->conss[c]->constype = REOPT_CONSTYPE_INFSUBTREE;
         else if( strcmp("reopt_dual", SCIPconsGetName(addedcons[c])) == 0 )
            reopttree->reoptnodes[id]->conss[c]->constype = REOPT_CONSTYPE_DUALREDS;
         else
            reopttree->reoptnodes[id]->conss[c]->constype = REOPT_CONSTYPE_UNKNOWN;
      }

      assert(reopttree->reoptnodes[id]->nconss == naddedconss);
      SCIPsetFreeBufferArray(set, &addedcons);
   }

   return SCIP_OKAY;
}

/** collect all bound changes based on dual information
 *
 *  If the bound changes are global, all information are already stored because they were caught by the event handler.
 *  otherwise, we have to use SCIPnodeGetDualBoundchgs.
 *
 *  Afterwards, we check if the constraint will be added in the next iteration or after splitting the node.
 */
static
SCIP_RETCODE collectDualInformation(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int          id,                 /**< id of the node */
   SCIP_REOPTTYPE        reopttype           /**< reopttype */
   )
{
   SCIP_Bool cons_is_next = TRUE;
   int nbndchgs;
   int v;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id]->dualreds);
   assert(node != NULL);
   assert(blkmem != NULL);

   /* first case, all bound changes were global */
   if( reopt->currentnode == SCIPnodeGetNumber(node) && reopt->dualreds != NULL && reopt->dualreds->nvars > 0 )
   {
      nbndchgs = reopt->dualreds->nvars;
   }
   else
   {
      assert(reopt->currentnode == SCIPnodeGetNumber(node));

      /* get the number of bound changes based on dual information */
      nbndchgs = SCIPnodeGetNDualBndchgs(node);

      /* ensure that enough memory is allocated */
      SCIP_CALL( checkMemDualCons(reopt, set, blkmem, nbndchgs) );

      /* collect the bound changes */
      SCIPnodeGetDualBoundchgs(node, reopt->dualreds->vars, reopt->dualreds->vals, reopt->dualreds->boundtypes,
            &nbndchgs, reopt->dualreds->varssize);
      assert(nbndchgs <= reopt->dualreds->varssize);

      reopt->dualreds->nvars = nbndchgs;
      reopt->dualreds->linear = FALSE;

      /* transform the variables into the original space */
      for( v = 0; v < nbndchgs; v++ )
      {
         SCIP_Real constant = 0.0;
         SCIP_Real scalar = 1.0;

         SCIP_CALL( SCIPvarGetOrigvarSum(&reopt->dualreds->vars[v], &scalar, &constant) );
         reopt->dualreds->vals[v] = (reopt->dualreds->vals[v] - constant) / scalar;

         assert(SCIPvarIsOriginal(reopt->dualreds->vars[v]));
      }
   }

   assert(nbndchgs > 0);

   /* due to the strong branching initialization it can be possible that two
    * constraints handling dual information are stored at the same time.
    * During reoptimizing a node we add the constraint stored at dualredscur only,
    * i.e, if dualredscur is not NULL, we need to store the constraint for the next
    * iteration at dualredsnex because the constraint stored at dualredscur is needed
    * to split the constraint in the current iteration.
    */
   if( reopt->reopttree->reoptnodes[id]->dualredscur != NULL )
   {
      assert(reopt->reopttree->reoptnodes[id]->dualredsnex == NULL);
      cons_is_next = FALSE;
   }
   assert((cons_is_next && reopt->reopttree->reoptnodes[id]->dualredscur == NULL)
       || (!cons_is_next && reopt->reopttree->reoptnodes[id]->dualredsnex == NULL));

   /* the constraint will be added next */
   if( cons_is_next )
   {
      assert(reopt->reopttree->reoptnodes[id]->dualredscur == NULL);
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopt->reopttree->reoptnodes[id]->dualredscur) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualredscur->vars, \
            reopt->dualreds->vars, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualredscur->vals, \
            reopt->dualreds->vals, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualredscur->boundtypes, \
            reopt->dualreds->boundtypes, nbndchgs) );

      reopt->reopttree->reoptnodes[id]->dualredscur->nvars = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualredscur->varssize = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualredscur->lhs = 1.0;
      reopt->reopttree->reoptnodes[id]->dualredscur->rhs = SCIPsetInfinity(set);
      reopt->reopttree->reoptnodes[id]->dualredscur->constype = (reopttype == SCIP_REOPTTYPE_STRBRANCHED ?
         REOPT_CONSTYPE_DUALREDS : REOPT_CONSTYPE_INFSUBTREE);
      reopt->reopttree->reoptnodes[id]->dualredscur->linear = FALSE;

      SCIPsetDebugMsg(set, " -> save dual information of type 1: node %lld, nvars %d, constype %d\n",
            SCIPnodeGetNumber(node), reopt->reopttree->reoptnodes[id]->dualredscur->nvars,
            reopt->reopttree->reoptnodes[id]->dualredscur->constype);
   }
   /* the constraint will be added after next */
   else
   {
      assert(reopt->reopttree->reoptnodes[id]->dualredsnex == NULL);
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopt->reopttree->reoptnodes[id]->dualredsnex) );
      reopt->reopttree->reoptnodes[id]->dualredsnex->nvars = -1;

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualredsnex->vars, \
            reopt->dualreds->vars, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualredsnex->vals, \
            reopt->dualreds->vals, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualredsnex->boundtypes, \
            reopt->dualreds->boundtypes, nbndchgs) );
      reopt->reopttree->reoptnodes[id]->dualredsnex->nvars = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualredsnex->varssize = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualredsnex->lhs = 1.0;
      reopt->reopttree->reoptnodes[id]->dualredsnex->rhs = SCIPsetInfinity(set);
      reopt->reopttree->reoptnodes[id]->dualredsnex->constype = (reopttype == SCIP_REOPTTYPE_STRBRANCHED ?
         REOPT_CONSTYPE_DUALREDS : REOPT_CONSTYPE_INFSUBTREE);

      SCIPsetDebugMsg(set, " -> save dual information of type 2: node %lld, nvars %d, constype %d\n",
         SCIPnodeGetNumber(node), reopt->reopttree->reoptnodes[id]->dualredsnex->nvars,
         reopt->reopttree->reoptnodes[id]->dualredsnex->constype);
   }

   return SCIP_OKAY;
}

/** adds a node of the branch and bound tree to the reoptimization tree */
static
SCIP_RETCODE addNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< current node */
   SCIP_REOPTTYPE        reopttype,          /**< reason for storing the node*/
   SCIP_Bool             saveafterdual,      /**< save branching decisions after the first dual */
   SCIP_Bool             isrootnode,         /**< node is the root node */
   SCIP_Real             lowerbound          /**< lower bound of the node */
   )
{
   SCIP_NODE* parent = NULL;
   SCIP_Bool shrank = FALSE;
   unsigned int id;
   unsigned int parentid = 0;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);

   if( set->reopt_maxsavednodes == 0 )
      return SCIP_OKAY;

   assert(reopttype == SCIP_REOPTTYPE_TRANSIT
       || reopttype == SCIP_REOPTTYPE_INFSUBTREE
       || reopttype == SCIP_REOPTTYPE_STRBRANCHED
       || reopttype == SCIP_REOPTTYPE_LOGICORNODE
       || reopttype == SCIP_REOPTTYPE_LEAF
       || reopttype == SCIP_REOPTTYPE_PRUNED
       || reopttype == SCIP_REOPTTYPE_FEASIBLE);

   /* start clock */
   SCIPclockStart(reopt->savingtime, set);

   /* the node was created by reoptimization, i.e., we need to update the
    * stored data */
   if( SCIPnodeGetReoptID(node) >= 1 )
   {
      SCIP_Bool transintoorig;

      assert(reopttype != SCIP_REOPTTYPE_LEAF);
      assert(!isrootnode);

      id = SCIPnodeGetReoptID(node);
      assert(id < reopt->reopttree->reoptnodessize);

      /* this is a special case:
       *   due to re-propagation of the an anchester node it can happen that we try to update a node that was created by
       *   reoptimization and already removed by deleteChildrenBelow. In this case we do not want to save the current
       *   node
       */
      if( reopt->reopttree->reoptnodes[id] == NULL )
      {
         parent = SCIPnodeGetParent(node);
         assert(parent != NULL);

         parentid = SCIPnodeGetReoptID(parent);

         /* traverse along the branching path until reaching a node that is part of the reoptimization tree or the root node */
         while( SCIPnodeGetDepth(parent) > 0 && reopt->reopttree->reoptnodes[parentid] == NULL )
         {
            /* the parent node is not part of the reoptimization, reset the reoptid and reopttype of the parent node */
            SCIPnodeSetReoptID(parent, 0);
            SCIPnodeSetReopttype(parent, SCIP_REOPTTYPE_NONE);

            parent = SCIPnodeGetParent(parent);
            assert(parent != NULL);

            parentid = SCIPnodeGetReoptID(parent);
         }

         /* the anchestor node has to be part of the reoptimization tree. either the parent is the root itself or
          * marked to be a leaf, pruned or feasible
          */
         assert(reopt->reopttree->reoptnodes[parentid] != NULL);
         assert(parentid == 0
             || reopt->reopttree->reoptnodes[parentid]->reopttype == SCIP_REOPTTYPE_FEASIBLE
             || reopt->reopttree->reoptnodes[parentid]->reopttype == SCIP_REOPTTYPE_INFSUBTREE
             || reopt->reopttree->reoptnodes[parentid]->reopttype == SCIP_REOPTTYPE_LEAF
             || reopt->reopttree->reoptnodes[parentid]->reopttype == SCIP_REOPTTYPE_PRUNED); /*lint !e641*/

         SCIPsetDebugMsg(set, " -> skip saving\n");
         SCIPnodeSetReoptID(node, 0);
         SCIPnodeSetReopttype(node, SCIP_REOPTTYPE_NONE);

         /* stop clock */
         SCIPclockStop(reopt->savingtime, set);

         return SCIP_OKAY;
      }

      SCIPsetDebugMsg(set, "update node %lld at ID %u:\n", SCIPnodeGetNumber(node), id);

      transintoorig = FALSE;

      /* store separated cuts */
      if( set->reopt_usecuts )
      {
         SCIP_CALL( storeCuts(reopt, set, blkmem, lp, id) );
      }

      /* save primal bound changes made after the first dual bound change */
      if( saveafterdual )
      {
         assert(reopttype == SCIP_REOPTTYPE_STRBRANCHED);
         SCIP_CALL( saveAfterDualBranchings(reopt, set, blkmem, node, id, &transintoorig) );
      }

      /* update constraint propagations */
      if( set->reopt_saveconsprop )
      {
         SCIP_CALL( updateConstraintPropagation(reopt, set, blkmem, node, id, &transintoorig) );
      }

      /* ensure that all variables describing the branching path are original */
      if( transintoorig )
      {
         SCIP_CALL( transformIntoOrig(reopt, id) );
      }

      /* update the lowerbound if the new lower bound is finite */
      if( !SCIPsetIsInfinity(set, REALABS(lowerbound)) )
         reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;
      SCIPsetDebugMsg(set, " -> reopttype: %u, lowerbound: %g\n", reopttype, reopt->reopttree->reoptnodes[id]->lowerbound);

#ifdef SCIP_MORE_DEBUG
      {
         int varnr;
         SCIPsetDebugMsg(set, " -> saved variables:\n");
         for (varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nvars; varnr++)
         {
            SCIPsetDebugMsg(set, "  <%s> %s %g\n", SCIPvarGetName(reopt->reopttree->reoptnodes[id]->vars[varnr]),
               reopt->reopttree->reoptnodes[id]->varboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
               "=>" : "<=", reopt->reopttree->reoptnodes[id]->varbounds[varnr]);
         }
         for (varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nafterdualvars; varnr++)
         {
            int varnr;
            SCIPsetDebugMsg(set, " -> saved variables:\n");
            for (varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nvars; varnr++)
            {
               SCIPsetDebugMsg(set, "  <%s> %s %g\n", SCIPvarGetName(reopt->reopttree->reoptnodes[id]->vars[varnr]),
                  reopt->reopttree->reoptnodes[id]->varboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                  "=>" : "<=", reopt->reopttree->reoptnodes[id]->varbounds[varnr]);
            }
            for (varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nafterdualvars; varnr++)
            {
               SCIPsetDebugMsg(set, "  <%s> %s %g (after dual red.)\n", SCIPvarGetName(reopt->reopttree->reoptnodes[id]->afterdualvars[varnr]),
                  reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                  "=>" : "<=", reopt->reopttree->reoptnodes[id]->afterdualvarbounds[varnr]);
            }
         }
      }
#endif

      /* update LPI state */
      switch( reopttype )
      {
      case SCIP_REOPTTYPE_TRANSIT:
         if( set->reopt_shrinkinner )
         {
            SCIP_CALL( shrinkNode(reopt, set, node, id, &shrank, blkmem) );
         }
         goto TRANSIT;

      case SCIP_REOPTTYPE_LOGICORNODE:
      case SCIP_REOPTTYPE_LEAF:
         goto TRANSIT;

      case SCIP_REOPTTYPE_INFSUBTREE:
         /* delete the whole subtree induced be the current node */
         SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, id, FALSE, FALSE) );
         goto PSEUDO;

      case SCIP_REOPTTYPE_STRBRANCHED:
         goto PSEUDO;

      case SCIP_REOPTTYPE_FEASIBLE:
         /* delete the subtree */
         if( set->reopt_reducetofrontier )
         {
            SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, id, FALSE, FALSE) );
            SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );
         }
         /* dive through all children and change the reopttype to PRUNED */
         else
         {
            SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, id, SCIP_REOPTTYPE_PRUNED) );
         }
         goto FEASIBLE;

      case SCIP_REOPTTYPE_PRUNED:
         /* delete the subtree */
         if( set->reopt_reducetofrontier )
         {
            SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, id, FALSE, FALSE) );
            SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );
         }
         /* dive through all children and change the reopttype to LEAF */
         else
         {
            SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, id, SCIP_REOPTTYPE_PRUNED) );
         }

         /* increase number of reoptimized nodes that could be pruned */
         ++reopt->reopttree->ncutoffreoptnodes;
         ++reopt->reopttree->ntotalcutoffreoptnodes;

         goto PRUNED;

      default:
         break;
      } /*lint !e788*/

      /* stop clock */
      SCIPclockStart(reopt->savingtime, set);

      return SCIP_OKAY;
   }

   /* get new IDs */
   SCIP_CALL( reopttreeCheckMemory(reopt->reopttree, set, blkmem) );

   /* the current node is the root node */
   if( isrootnode )
   {
      id = 0;

      /* save local constraints
       * note: currently, there will be no constraint to save because all global constraints are added by calling
       * SCIPprobAddCons.
       */
      if (SCIPnodeGetNAddedConss(node) >= 1)
      {
         assert(reopt->reopttree->reoptnodes[id]->nconss == 0);

         SCIP_CALL( saveLocalConssData(reopt->reopttree, set, blkmem, node, id) );
      }

      /* store separated cuts
       * note: we need to call this after saveLocalConssData to be sure that the local conss array is ordered, first all
       * local constraints, then cuts
       */
      if( set->reopt_usecuts )
      {
         SCIP_CALL( storeCuts(reopt, set, blkmem, lp, id) );
      }

      switch( reopttype )
      {
      case SCIP_REOPTTYPE_TRANSIT:
         /* ensure that no dual constraints are stored */
         SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );

         /* update the lowerbound */
         if( !SCIPsetIsInfinity(set, REALABS(lowerbound)) )
            reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;

         goto TRANSIT;

      case SCIP_REOPTTYPE_INFSUBTREE:
      case SCIP_REOPTTYPE_STRBRANCHED:
         reopt->reopttree->reoptnodes[0]->reopttype = (unsigned int)reopttype;
         reopt->reopttree->reoptnodes[0]->dualreds = TRUE;
         reopt->reopttree->reoptnodes[0]->nvars = 0;

         if( reopttype == SCIP_REOPTTYPE_INFSUBTREE )
         {
            /* delete the whole subtree induced be the current node */
            SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, 0, FALSE, FALSE) );
         }

         /* update the lowerbound */
         if( !SCIPsetIsInfinity(set, REALABS(lowerbound)) )
            reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;

         SCIPsetDebugMsg(set, "update node %d at ID %d:\n", 1, 0);
         SCIPsetDebugMsg(set, " -> nvars: 0, ncons: 0, parentID: -, reopttype: %u, lowerbound: %g\n", reopttype,
               reopt->reopttree->reoptnodes[id]->lowerbound);

         goto PSEUDO;

      case SCIP_REOPTTYPE_FEASIBLE:
         ++reopt->reopttree->ntotalfeasnodes;
         ++reopt->reopttree->nfeasnodes;
         reopt->reopttree->reoptnodes[0]->reopttype = (unsigned int)SCIP_REOPTTYPE_FEASIBLE;
         reopt->reopttree->reoptnodes[0]->dualreds = FALSE;

         if( reopt->reopttree->reoptnodes[0]->childids != NULL && reopt->reopttree->reoptnodes[0]->nchilds > 0 )
         {
            /* delete the subtree */
            if( set->reopt_reducetofrontier )
            {
               SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, 0, FALSE, FALSE) );
               SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );
            }
            /* dive through all children and change the reopttype to LEAF */
            else
            {
               SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, 0, SCIP_REOPTTYPE_PRUNED) );
            }
         }
         else
            SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );

         /* update the lowerbound */
         if( !SCIPsetIsInfinity(set, REALABS(lowerbound)) )
            reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;

         SCIPsetDebugMsg(set, "update node %d at ID %d:\n", 1, 0);
         SCIPsetDebugMsg(set, " -> nvars: 0, ncons: 0, parentID: -, reopttype: %u, lowerbound: %g\n", reopttype,
               reopt->reopttree->reoptnodes[id]->lowerbound);

         break;

      case SCIP_REOPTTYPE_PRUNED:
         ++reopt->reopttree->nprunednodes;
         ++reopt->reopttree->ntotalprunednodes;
         reopt->reopttree->reoptnodes[0]->reopttype = (unsigned int)SCIP_REOPTTYPE_PRUNED;
         reopt->reopttree->reoptnodes[0]->dualreds = FALSE;

         if( reopt->reopttree->reoptnodes[0]->childids != NULL && reopt->reopttree->reoptnodes[0]->nchilds > 0 )
         {
            /* delete the subtree */
            if( set->reopt_reducetofrontier )
            {
               SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, 0, FALSE, FALSE) );
               SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );
            }
            /* dive through all children and change the reopttype to LEAF */
            else
            {
               SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, 0, SCIP_REOPTTYPE_PRUNED) );
            }
         }
         else
            SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );

         /* update the lowerbound if it was not set */
         if( !SCIPsetIsInfinity(set, REALABS(lowerbound)) )
            reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;

         SCIPsetDebugMsg(set, "update node %d at ID %d:\n", 1, 0);
         SCIPsetDebugMsg(set, " -> nvars: 0, ncons: 0, parentID: -, reopttype: %u, lowerbound:%g \n", reopttype,
               reopt->reopttree->reoptnodes[id]->lowerbound);

         break;

      default:
         assert(reopttype == SCIP_REOPTTYPE_TRANSIT
            || reopttype == SCIP_REOPTTYPE_INFSUBTREE
            || reopttype == SCIP_REOPTTYPE_STRBRANCHED
            || reopttype == SCIP_REOPTTYPE_PRUNED
            || reopttype == SCIP_REOPTTYPE_FEASIBLE);
         break;
      }/*lint !e788*/

      /* reset the information of dual bound changes */
      reopt->currentnode = -1;
      if( reopt->dualreds != NULL )
         reopt->dualreds->nvars = 0;

      /* stop clock */
      SCIPclockStop(reopt->savingtime, set);

      return SCIP_OKAY;
   }
   else
   {
      int nbndchgdiff;
      SCIP_Bool transintoorig;

      SCIPsetDebugMsg(set, "try to add node #%lld to the reopttree\n", SCIPnodeGetNumber(node));
      SCIPsetDebugMsg(set, " -> reopttype = %u\n", reopttype);

      /* check if we really want to save this node:
       *   1. save the node if reopttype is at least SCIP_REOPTTYPE_INFSUBTREE
       *   2. save the node if the number of bound changes of this node
       *      and the last saved node is at least a given number n
       */

      /* get the ID of the last saved node or 0 for the root */
      SCIP_CALL( getLastSavedNode(reopt, set, node, &parent, &parentid, &nbndchgdiff) );

      if( (reopttype < SCIP_REOPTTYPE_INFSUBTREE && nbndchgdiff <= set->reopt_maxdiffofnodes)
        || reopt->reopttree->reoptnodes[parentid]->reopttype >= SCIP_REOPTTYPE_LEAF ) /*lint !e641*/
      {
         SCIPsetDebugMsg(set, " -> skip saving\n");

         /* stop clock */
         SCIPclockStop(reopt->savingtime, set);

         return SCIP_OKAY;
      }

      /* check if there are free slots to store the node */
      SCIP_CALL( reopttreeCheckMemory(reopt->reopttree, set, blkmem) );

      id = (unsigned int) (size_t) SCIPqueueRemove(reopt->reopttree->openids);

      SCIPsetDebugMsg(set, " -> save at ID %u\n", id);

      assert(reopt->reopttree->reoptnodes[id] == NULL
         || (reopt->reopttree->reoptnodes[id]->nvars == 0 && reopt->reopttree->reoptnodes[id]->nconss == 0));
      assert(id >= 1 && id < reopt->reopttree->reoptnodessize);
      assert(!isrootnode);

      /* get memory for nodedata */
      assert(reopt->reopttree->reoptnodes[id] == NULL || reopt->reopttree->reoptnodes[id]->nvars == 0);
      SCIP_CALL( createReoptnode(reopt->reopttree, set, blkmem, id) );
      reopt->reopttree->reoptnodes[id]->parentID = parentid;

      assert(parent != NULL );
      assert((SCIPnodeGetDepth(parent) == 0 && parentid == 0) || (SCIPnodeGetDepth(parent) >= 1 && parentid > 0));
      assert(id >= 1);

      /* create the array of "child nodes" if they not exist */
      if( reopt->reopttree->reoptnodes[parentid]->childids == NULL
         || reopt->reopttree->reoptnodes[parentid]->allocchildmem == 0 )
      {
         SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[parentid], set, blkmem, 0, 2, 0) );
      }

      /* add the new node as a "child node" of the last saved reoptminization node */
      SCIP_CALL( reoptAddChild(reopt->reopttree, set, blkmem, parentid, id) );

      /* save branching path */
      SCIP_CALL( saveAncestorBranchings(reopt->reopttree, set, blkmem, node, parent, id, parentid) );

      /* save bound changes after some dual reduction */
      if( saveafterdual )
      {
         assert(reopttype == SCIP_REOPTTYPE_STRBRANCHED);
         SCIP_CALL( saveAfterDualBranchings(reopt, set, blkmem, node, id, &transintoorig) );
      }
      else
      {
         SCIPsetDebugMsg(set, " -> skip saving bound changes after dual reductions.\n");
      }

      /* transform all bounds of branched variables and ensure that they are original. */
      SCIP_CALL( transformIntoOrig(reopt, id) );

      /* save pseudo-constraints (if one exists) */
      if (SCIPnodeGetNAddedConss(node) >= 1)
      {
         assert(reopt->reopttree->reoptnodes[id]->nconss == 0);

         SCIP_CALL( saveLocalConssData(reopt->reopttree, set, blkmem, node, id) );
      }

      /* store separated cuts
       * note: we need to call this after saveLocalConssData to be sure that the local conss array is ordered, first all
       * local constraints, then cuts
       */
      if( set->reopt_usecuts )
      {
         SCIP_CALL( storeCuts(reopt, set, blkmem, lp, id) );
      }

      /* update the lowerbound if it was not set */
      if( !SCIPsetIsInfinity(set, REALABS(lowerbound)) )
         reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;

      /* set ID */
      SCIPnodeSetReoptID(node, id);

      /* set the REOPTTYPE */
      SCIPnodeSetReopttype(node, reopttype);

      SCIPsetDebugMsg(set, "save node #%lld successful\n", SCIPnodeGetNumber(node));
      SCIPsetDebugMsg(set, " -> nvars: %d, ncons: %d, parentID: %u, reopttype: %d, lowerbound: %g\n",
         reopt->reopttree->reoptnodes[id]->nvars + reopt->reopttree->reoptnodes[id]->nafterdualvars,
         reopt->reopttree->reoptnodes[id]->nconss, reopt->reopttree->reoptnodes[id]->parentID,
         reopttype, reopt->reopttree->reoptnodes[id]->lowerbound);
#ifdef SCIP_MORE_DEBUG
      {
         int varnr;
         for (varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nvars; varnr++)
         {
            int varnr;
            for (varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nvars; varnr++)
            {
               SCIPsetDebugMsg(set, "  <%s> %s %g\n", SCIPvarGetName(reopt->reopttree->reoptnodes[id]->vars[varnr]),
                  reopt->reopttree->reoptnodes[id]->varboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                  "=>" : "<=", reopt->reopttree->reoptnodes[id]->varbounds[varnr]);
            }
            for (varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nafterdualvars; varnr++)
            {
               SCIPsetDebugMsg(set, "  <%s> %s %g (after dual red.)\n",
                  SCIPvarGetName(reopt->reopttree->reoptnodes[id]->afterdualvars[varnr]),
                  reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                  "=>" : "<=", reopt->reopttree->reoptnodes[id]->afterdualvarbounds[varnr]);
            }
         }
      }
#endif
   }

   switch( reopttype )
   {
   case SCIP_REOPTTYPE_TRANSIT:
   case SCIP_REOPTTYPE_LOGICORNODE:
   case SCIP_REOPTTYPE_LEAF:
  TRANSIT:
      if( !shrank )
         reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)reopttype;
      else
      {
         SCIPnodeSetReoptID(node, 0);
         SCIPnodeSetReopttype(node, SCIP_REOPTTYPE_NONE);
      }
      break;

   case SCIP_REOPTTYPE_INFSUBTREE:
   case SCIP_REOPTTYPE_STRBRANCHED:
  PSEUDO:
      assert(reopt->currentnode == SCIPnodeGetNumber(node));

      reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)reopttype;
      reopt->reopttree->reoptnodes[id]->dualreds = TRUE;

      /* get all the dual information and decide if the constraint need
       * to be added next or after next */
      SCIP_CALL( collectDualInformation(reopt, set, blkmem, node, id, reopttype) );

      break;

   case SCIP_REOPTTYPE_FEASIBLE:
  FEASIBLE:
      reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_FEASIBLE;
      reopt->reopttree->reoptnodes[id]->dualreds = FALSE;
      ++reopt->reopttree->nfeasnodes;
      ++reopt->reopttree->ntotalfeasnodes;

      break;

   case SCIP_REOPTTYPE_PRUNED:
  PRUNED:
      reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_PRUNED;
      reopt->reopttree->reoptnodes[id]->dualreds = FALSE;
      ++reopt->reopttree->nprunednodes;
      ++reopt->reopttree->ntotalprunednodes;

      break;

   default:
      assert(reopttype == SCIP_REOPTTYPE_TRANSIT
         || reopttype == SCIP_REOPTTYPE_LOGICORNODE
         || reopttype == SCIP_REOPTTYPE_LEAF
         || reopttype == SCIP_REOPTTYPE_INFSUBTREE
         || reopttype == SCIP_REOPTTYPE_STRBRANCHED
         || reopttype == SCIP_REOPTTYPE_FEASIBLE
         || reopttype == SCIP_REOPTTYPE_PRUNED);
      break;
   } /*lint !e788*/

   /* stop clock */
   SCIPclockStop(reopt->savingtime, set);

   /* reset the information of dual bound changes */
   reopt->currentnode = -1;
   if( reopt->dualreds != NULL )
      reopt->dualreds->nvars = 0;

   return SCIP_OKAY;
}

/** delete the stored information about dual bound changes of the last focused node */
static
void deleteLastDualBndchgs(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   if( reopt->dualreds != NULL && reopt->dualreds->nvars > 0 )
   {
      SCIPdebugMessage("delete %d dual variable information about node %lld\n", reopt->dualreds->nvars,
            reopt->currentnode);
      reopt->dualreds->nvars = 0;
      reopt->currentnode = -1;
   }
}

/** delete the stored constraints that dual information at the given reoptimization node */
static
SCIP_RETCODE reoptnodeResetDualConss(
   SCIP_REOPTNODE*       reoptnode,          /**< reoptimization node */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reoptnode != NULL);
   assert(blkmem != NULL);

   if( reoptnode->dualredscur != NULL )
   {
      SCIP_REOPTCONSDATA* reoptconsdata;

      SCIPdebugMessage("reset dual information (current run)\n");

      reoptconsdata = reoptnode->dualredscur;

      BMSfreeBlockMemoryArray(blkmem, &reoptconsdata->boundtypes, reoptconsdata->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptconsdata->vals, reoptconsdata->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptconsdata->vars, reoptconsdata->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualredscur);
      reoptnode->dualredscur = NULL;
   }

   if( reoptnode->dualredsnex != NULL )
   {
      SCIP_REOPTCONSDATA* reoptconsdata;

      SCIPdebugMessage("reset dual information (next run)\n");

      reoptconsdata = reoptnode->dualredsnex;

      BMSfreeBlockMemoryArray(blkmem, &reoptconsdata->boundtypes, reoptconsdata->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptconsdata->vals, reoptconsdata->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptconsdata->vars, reoptconsdata->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualredsnex);
      reoptnode->dualredsnex = NULL;
   }

   reoptnode->dualreds = FALSE;

   return SCIP_OKAY;
}


/** transform given set of variables, bounds and boundtypes into a global cut.
 *
 *  @note: boundtypes can be NULL if all variables are binary or a MIP solution should be separated.
 *  @note: continuous variables will be skiped if boundtypes is NULL
 */
static
SCIP_RETCODE addGlobalCut(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< variables of the cut */
   SCIP_Real*            vals,               /**< values of the cut */
   SCIP_BOUNDTYPE*       boundtypes,         /**< bounds of the cut */
   int                   nvars,              /**< number of variables in the cut */
   int                   nbinvars,           /**< number of binary variables */
   int                   nintvars            /**< number of integer variables */
   )
{
   SCIP_REOPTCONSDATA* reoptconsdata;
   int nglbconss;
   int nvarsadded;
   int v;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(nbinvars + nintvars == nvars);

   nvarsadded = 0;

   /* check whether we have enough memory allocated */
   SCIP_CALL( checkMemGlbCons(reopt, set, blkmem, 10) );
   nglbconss = reopt->nglbconss;
   reoptconsdata = NULL;

   if( reopt->glbconss[nglbconss] == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopt->glbconss[nglbconss]) ); /*lint !e866*/
      reoptconsdata = reopt->glbconss[nglbconss];

      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptconsdata->vars, (int)(nbinvars+2*nintvars)) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptconsdata->vals, (int)(nbinvars+2*nintvars)) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptconsdata->boundtypes, (int)(nbinvars+2*nintvars)) );
      reoptconsdata->varssize = (int)(nbinvars+2*nintvars);
      reoptconsdata->nvars = 0;
   }
   else
   {
      assert(reopt->glbconss[nglbconss]->nvars == 0);
      assert(reopt->glbconss[nglbconss]->varssize > 0);

      reoptconsdata = reopt->glbconss[nglbconss];

      if( reoptconsdata->varssize < nbinvars+2*nintvars )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptconsdata->vars, reoptconsdata->varssize, \
               (int)(nbinvars+2*nintvars)) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptconsdata->vals, reoptconsdata->varssize, \
               (int)(nbinvars+2*nintvars)) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptconsdata->boundtypes, reoptconsdata->varssize, \
               (int)(nbinvars+2*nintvars)) );
         reoptconsdata->varssize = (int)(nbinvars+2*nintvars);
      }
   }
   assert(reoptconsdata != NULL);

   reoptconsdata->lhs = 1.0;
   reoptconsdata->rhs = SCIPsetInfinity(set);
   reoptconsdata->linear = FALSE;
   reoptconsdata->constype = REOPT_CONSTYPE_CUT;

   for( v = 0; v < nvars; v++ )
   {
      assert(nvarsadded < reoptconsdata->varssize);
      assert(vars[v] != NULL);
      assert(SCIPvarIsOriginal(vars[v]));
      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsIntegral(set, vals[v]));

      /* if no boundtypes are given we skip continuous variables, otherwise we would add trivial clauses:
       * a)       x <= ub
       * b) lb <= x
       * c) (x <= val) or (x >= val)
       */
      if( boundtypes == NULL && SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY )
      {
         reoptconsdata->vars[nvarsadded] = vars[v];

         if( SCIPsetIsEQ(set, vals[v], 1.0) )
         {
            assert(boundtypes == NULL || boundtypes[v] == SCIP_BOUNDTYPE_LOWER);
            reoptconsdata->vals[nvarsadded] = 0.0;
            reoptconsdata->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_UPPER;
         }
         else
         {
            assert(SCIPsetIsEQ(set, vals[v], 0.0));
            assert(boundtypes == NULL || boundtypes[v] == SCIP_BOUNDTYPE_UPPER);
            reoptconsdata->vals[nvarsadded] = 1.0;
            reoptconsdata->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_LOWER;
         }
         ++nvarsadded;
      }
      else if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS)
      {
         assert(boundtypes != NULL);

         reoptconsdata->vals[nvarsadded] = vals[v];
         reoptconsdata->boundtypes[nvarsadded] = (boundtypes[v] == SCIP_BOUNDTYPE_LOWER ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
         ++nvarsadded;
      }
      else
      {
         SCIP_Real roundedval;
         SCIP_Real ubglb;
         SCIP_Real lbglb;

         assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_IMPLINT);

         reoptconsdata->vars[nvarsadded] = vars[v];

         ubglb = SCIPvarGetUbGlobal(vars[v]);
         lbglb = SCIPvarGetLbGlobal(vars[v]);

         /* case 1  :      x == val == ub -> x <= ub-1
          * case 2  :      x == val == lb -> x >= lb+1
          * case 3.1:      x <= val <  ub -> x >= y+1
          * case 3.2:      x >= val >  lb -> x <= y-1
          * case 4  : lb < x == val <  ub -> (x <= y-1) or (x >= y+1)
          */

         /* case 1 */
         if( SCIPsetIsEQ(set, vals[v], ubglb) )
         {
            assert(boundtypes == NULL || boundtypes[v] == SCIP_BOUNDTYPE_LOWER);
            reoptconsdata->vals[nvarsadded] = ubglb - 1.0;
            reoptconsdata->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_UPPER;
            ++nvarsadded;
         }
         /* case 2 */
         else if( SCIPsetIsEQ(set, vals[v], lbglb) )
         {
            assert(boundtypes == NULL || boundtypes[v] == SCIP_BOUNDTYPE_UPPER);
            reoptconsdata->vals[nvarsadded] = lbglb + 1.0;
            reoptconsdata->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_LOWER;
            ++nvarsadded;
         }
         else if( boundtypes != NULL )
         {
            /* we round the solution value to get a 'clean' bound */
            assert(SCIPsetIsIntegral(set, vals[v]));
            roundedval = SCIPsetRound(set, vals[v]);

            /* case 3.1 */
            if( boundtypes[v] == SCIP_BOUNDTYPE_UPPER )
            {
               reoptconsdata->vals[nvarsadded] = roundedval + 1.0;
               reoptconsdata->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_LOWER;
               ++nvarsadded;
            }
            /* case 3.2 */
            else
            {
               assert(boundtypes[v] == SCIP_BOUNDTYPE_LOWER);
               reoptconsdata->vals[nvarsadded] = roundedval - 1.0;
               reoptconsdata->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_UPPER;
               ++nvarsadded;
            }
         }
         /* case 4: in this case we have to add two clauses: (x <= val-1) and (x >= val+1) */
         else
         {
            /* we round the solution value to get a 'clean' bound */
            assert(SCIPsetIsIntegral(set, vals[v]));
            roundedval = SCIPsetRound(set, vals[v]);

            /* first clause: x <= val-1 */
            reoptconsdata->vals[nvarsadded] = roundedval - 1.0;
            reoptconsdata->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_UPPER;
            ++nvarsadded;

            /* second clause:  x >= val+1 */
            reoptconsdata->vars[nvarsadded] = vars[v];
            reoptconsdata->vals[nvarsadded] = roundedval + 1.0;
            reoptconsdata->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_LOWER;
            ++nvarsadded;
         }
      }
   }
   assert(nvars <= nvarsadded);
   assert(nvarsadded == nbinvars + 2 * nintvars);

   reoptconsdata->nvars = nvarsadded;
   ++reopt->nglbconss;

   return SCIP_OKAY;
}

/** generate a global constraint to separate an infeasible subtree */
static
SCIP_RETCODE saveGlobalCons(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   REOPT_CONSTYPE        consttype           /**< reopttype of the constraint */
   )
{
   assert(reopt != NULL);
   assert(node != NULL);

   if( consttype == REOPT_CONSTYPE_INFSUBTREE )
   {
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_BOUNDTYPE* boundtypes;
      int allocmem;
      int nbranchvars;
      int nbinvars;
      int nintvars;
      int v;

      /* allocate memory to store the infeasible path */
      allocmem = SCIPnodeGetDepth(node);
      SCIP_CALL( SCIPsetAllocBufferArray(set, &vars, allocmem) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, allocmem) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &boundtypes, allocmem) );

      /* get the branching path */
      SCIPnodeGetAncestorBranchings(node, vars, vals, boundtypes, &nbranchvars, allocmem);

      if( allocmem < nbranchvars )
      {
         SCIP_CALL( SCIPsetReallocBufferArray(set, &vars, nbranchvars) );
         SCIP_CALL( SCIPsetReallocBufferArray(set, &vals, nbranchvars) );
         SCIP_CALL( SCIPsetReallocBufferArray(set, &boundtypes, nbranchvars) );
         allocmem = nbranchvars;

         SCIPnodeGetAncestorBranchings(node, vars, vals, boundtypes, &nbranchvars, allocmem);
      }

      /* we count the number of binary and (impl) integer variables */
      nbinvars = 0;
      nintvars = 0;
      for( v = 0; v < nbranchvars; v++ )
      {
         if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY )
            ++nbinvars;
         if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_IMPLINT )
            ++nintvars;
      }
      assert(nbinvars + nintvars == nbranchvars);

      SCIP_CALL( addGlobalCut(reopt, blkmem, set, vars, vals, boundtypes, nbranchvars, nbinvars, nintvars) );
      assert(!reopt->glbconss[reopt->nglbconss - 1]->linear);

      /* free buffer */
      SCIPsetFreeBufferArray(set, &boundtypes);
      SCIPsetFreeBufferArray(set, &vals);
      SCIPsetFreeBufferArray(set, &vars);
   }

   return SCIP_OKAY;
}


/** move all id of child nodes from reoptimization node stored at @p id1 to the node stored at @p id2 */
static
SCIP_RETCODE reoptMoveIDs(
   SCIP_REOPTTREE*       reopttree,          /**< reopttree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          id1,                /**< source id */
   unsigned int          id2                 /**< target id */
   )
{
   int c;
   int nchilds_id1;
   int nchilds_id2;

   assert(reopttree != NULL);
   assert(blkmem != NULL);
   assert(id1 < reopttree->reoptnodessize);
   assert(id2 < reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[id1] != NULL);
   assert(reopttree->reoptnodes[id2] != NULL);

   nchilds_id1 = reopttree->reoptnodes[id1]->nchilds;
   nchilds_id2 = reopttree->reoptnodes[id2]->nchilds;

   /* ensure that the array storing the child id's is large enough */
   SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id2], set, blkmem, 0, nchilds_id1+nchilds_id2, 0) );
   assert(reopttree->reoptnodes[id2]->allocchildmem >= nchilds_id1+nchilds_id2);

   SCIPsetDebugMsg(set, "move %d IDs: %u -> %u\n", nchilds_id1, id1, id2);

   /* move the ids */
   for( c = 0; c < nchilds_id1; c++ )
   {
#ifdef SCIP_DEBUG
      {
         /* check that no id is added twice */
         int k;
         for( k = 0; k < nchilds_id2; k++ )
            assert(reopttree->reoptnodes[id2]->childids[k] != reopttree->reoptnodes[id1]->childids[c]);
      }
#endif

      reopttree->reoptnodes[id2]->childids[nchilds_id2+c] = reopttree->reoptnodes[id1]->childids[c];
   }

   /* update the number of childs */
   reopttree->reoptnodes[id1]->nchilds = 0;
   reopttree->reoptnodes[id2]->nchilds += nchilds_id1;

   return SCIP_OKAY;
}

/** change all bound changes along the root path */
static
SCIP_RETCODE changeAncestorBranchings(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< search tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the branch and bound tree */
   unsigned int          id,                 /**< id of stored node */
   SCIP_Bool             afterdualintobranching /**< convert all bound changes made directly after the first bound
                                                 *   changes based on dual information into normal branchings
                                                 */
   )
{
   SCIP_REOPTTREE* reopttree;
   SCIP_REOPTNODE* reoptnode;
   int v;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(cliquetable != NULL);
   assert(node != NULL);
   assert(blkmem != NULL);

   reopttree = reopt->reopttree;
   assert(reopttree != NULL);
   assert(id < reopttree->reoptnodessize);

   reoptnode = reopttree->reoptnodes[id];
   assert(reoptnode != NULL);

   /* copy memory to ensure that only original variables are saved */
   if( reoptnode->nvars == 0 && reoptnode->nafterdualvars == 0)
      return SCIP_OKAY;

   /* change the bounds along the branching path */
   for( v = 0; v < reoptnode->nvars; v++ )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_BOUNDTYPE boundtype;
      SCIP_Real oldlb;
      SCIP_Real oldub;
      SCIP_Real newbound;

      var = reoptnode->vars[v];
      val = reoptnode->varbounds[v];
      boundtype = reoptnode->varboundtypes[v];

      assert(SCIPvarIsOriginal(var));
      SCIP_CALL( SCIPvarGetProbvarBound(&var, &val, &boundtype) );
      assert(SCIPvarIsTransformed(var));
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

      oldlb = SCIPvarGetLbLocal(var);
      oldub = SCIPvarGetUbLocal(var);
      newbound = val;

      assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsGT(set, newbound, oldlb) && SCIPsetIsFeasLE(set, newbound, oldub) )
      {
         SCIPvarAdjustLb(var, set, &newbound);

         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      }
      else if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsLT(set, newbound, oldub) && SCIPsetIsFeasGE(set, newbound, oldlb) )
      {
         SCIPvarAdjustUb(var, set, &newbound);

         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      }
#ifdef SCIP_MORE_DEBUG
      SCIPsetDebugMsg(set, "  (path) <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=", newbound);
#endif
   }

   if( afterdualintobranching && reoptnode->nafterdualvars > 0 )
   {
      /* check the memory to convert this bound changes into 'normal' */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], set, blkmem,
            reoptnode->nvars + reoptnode->nafterdualvars, 0, 0) );

      /* change the bounds */
      for( v = 0; v < reoptnode->nafterdualvars; v++ )
      {
         SCIP_VAR* var;
         SCIP_Real val;
         SCIP_BOUNDTYPE boundtype;
         SCIP_Bool bndchgd;
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Real newbound;

         var = reoptnode->afterdualvars[v];
         val = reoptnode->afterdualvarbounds[v];
         boundtype = reoptnode->afterdualvarboundtypes[v];

         assert(SCIPvarIsOriginal(var));
         SCIP_CALL( SCIPvarGetProbvarBound(&var, &val, &boundtype) );
         assert(SCIPvarIsTransformed(var));
         assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

         bndchgd = FALSE;

         oldlb = SCIPvarGetLbLocal(var);
         oldub = SCIPvarGetUbLocal(var);
         newbound = val;

         if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsGT(set, newbound, oldlb) && SCIPsetIsFeasLE(set, newbound, oldub) )
         {
            SCIPvarAdjustLb(var, set, &newbound);
            SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );

            bndchgd = TRUE;
         }
         else if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsLT(set, newbound, oldub) && SCIPsetIsFeasGE(set, newbound, oldlb) )
         {
            SCIPvarAdjustUb(var, set, &newbound);
            SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );

            bndchgd = TRUE;
         }

         assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

#ifdef SCIP_MORE_DEBUG
         SCIPsetDebugMsg(set, "   (prop) <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=", newbound);
#endif
         if( bndchgd )
         {
            int nvars;

            nvars = reoptnode->nvars;
            reoptnode->vars[nvars] = reoptnode->afterdualvars[v];
            reoptnode->varbounds[nvars] = reoptnode->afterdualvarbounds[v];
            reoptnode->varboundtypes[nvars] = reoptnode->afterdualvarboundtypes[v];
            ++reoptnode->nvars;
         }
      }

      /* free the afterdualvars, -bounds, and -boundtypes */
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->afterdualvarboundtypes, reoptnode->afterdualvarssize);
      reoptnode->afterdualvarboundtypes = NULL;

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->afterdualvarbounds, reoptnode->afterdualvarssize);
      reoptnode->afterdualvarbounds = NULL;

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->afterdualvars, reoptnode->afterdualvarssize);
      reoptnode->afterdualvars = NULL;

      reoptnode->nafterdualvars = 0;
      reoptnode->afterdualvarssize = 0;
   }

   return SCIP_OKAY;
}

/** add a constraint to ensure that at least one variable bound gets different */
static
SCIP_RETCODE addSplitcons(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< search tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_NODE*            node,               /**< node corresponding to the pruned part */
   unsigned int          id                  /**< id of stored node */
   )
{
   SCIP_CONS* cons;
   char name[SCIP_MAXSTRLEN];
   int v;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);
   assert(reopt->reopttree->reoptnodes[id]->dualreds);
   assert(reopt->reopttree->reoptnodes[id]->dualredscur != NULL);
   assert(scip != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(node != NULL);

   assert(reopt->reopttree->reoptnodes[id]->dualredscur->constype == REOPT_CONSTYPE_DUALREDS
         || reopt->reopttree->reoptnodes[id]->dualredscur->constype == REOPT_CONSTYPE_INFSUBTREE);

#ifndef NDEBUG
   if( reopt->reopttree->reoptnodes[id]->dualredscur->constype == REOPT_CONSTYPE_DUALREDS )
      SCIPsetDebugMsg(set, " create a split-node #%lld\n", SCIPnodeGetNumber(node));
   else
      SCIPsetDebugMsg(set, " separate an infeasible subtree\n");
#endif

   /* if the constraint consists of exactly one variable it can be interpreted
    * as a normal branching step, i.e., we can fix the variable to the negated bound */
   if( reopt->reopttree->reoptnodes[id]->dualredscur->nvars == 1 )
   {
      SCIP_REOPTCONSDATA* reoptconsdata;
      SCIP_VAR* var;
      SCIP_BOUNDTYPE boundtype;
      SCIP_Real oldlb;
      SCIP_Real oldub;
      SCIP_Real newbound;

      reoptconsdata = reopt->reopttree->reoptnodes[id]->dualredscur;
      assert(!reoptconsdata->linear);
      assert(reoptconsdata->vars != NULL);
      assert(reoptconsdata->vals != NULL);
      assert(reoptconsdata->boundtypes != NULL);

      var = reoptconsdata->vars[0];
      newbound = reoptconsdata->vals[0];
      boundtype = reoptconsdata->boundtypes[0];

      assert(SCIPvarIsOriginal(var));
      SCIP_CALL( SCIPvarGetProbvarBound(&var, &newbound, &boundtype) );
      assert(SCIPvarIsTransformed(var));

      oldlb = SCIPvarGetLbLocal(var);
      oldub = SCIPvarGetUbLocal(var);

      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         newbound = reoptconsdata->vals[0] - 1.0;
         assert(SCIPisLE(scip, newbound, oldub));
      }
      else
      {
         newbound = reoptconsdata->vals[0] + 1.0;
         assert(SCIPisGE(scip, newbound, oldlb));
      }
      boundtype = (SCIP_BOUNDTYPE) (1 - (int)boundtype);
      assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsGT(set, newbound, oldlb) && SCIPsetIsFeasLE(set, newbound, oldub) )
      {
         SCIPvarAdjustLb(var, set, &newbound);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      }
      else if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsLT(set, newbound, oldub) && SCIPsetIsFeasGE(set, newbound, oldlb) )
      {
         SCIPvarAdjustUb(var, set, &newbound);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      }

      SCIPsetDebugMsg(set, "  -> constraint consists of only one variable: <%s> %s %g\n", SCIPvarGetName(var),
            boundtype == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=", newbound);
   }
   else
   {
      SCIP_REOPTCONSDATA* reoptconsdata;
      SCIP_VAR** consvars;
      SCIP_Real consval;
      SCIP_BOUNDTYPE consboundtype;
      int nbinvars = 0;
      int nintvars = 0;
      int ncontvars = 0;

      reoptconsdata = reopt->reopttree->reoptnodes[id]->dualredscur;
      assert(!reoptconsdata->linear);
      assert(reoptconsdata->vars != NULL);
      assert(reoptconsdata->vals != NULL);
      assert(reoptconsdata->boundtypes != NULL);

      /* allocate buffer */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, reoptconsdata->nvars) );

      /* count number of binary, integer, and continuous variables */
      for( v = 0; v < reoptconsdata->nvars; v++ )
      {
         switch ( SCIPvarGetType(reoptconsdata->vars[v]) ) {
         case SCIP_VARTYPE_BINARY:
            ++nbinvars;
            break;
         case SCIP_VARTYPE_IMPLINT:
         case SCIP_VARTYPE_INTEGER:
            if( SCIPisEQ(scip, SCIPvarGetLbLocal(reoptconsdata->vars[v]), 0.0)
               && SCIPisEQ(scip, SCIPvarGetUbLocal(reoptconsdata->vars[v]), 1.0) )
               ++nbinvars;
            else
               ++nintvars;
            break;
         case SCIP_VARTYPE_CONTINUOUS:
            ++ncontvars;
            break;
         default:
            SCIPerrorMessage("Variable <%s> has to be either binary, (implied) integer, or continuous.\n",
               SCIPvarGetName(reoptconsdata->vars[v]));
            return SCIP_INVALIDDATA;
         }
      }

      if( reoptconsdata->constype == REOPT_CONSTYPE_INFSUBTREE )
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "reopt_inf");
      else
      {
         assert(reoptconsdata->constype == REOPT_CONSTYPE_DUALREDS);
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "reopt_dual");
      }

      /* case 1: all variables are binary. we use a logic-or constraint. */
      if( reoptconsdata->nvars == nbinvars )
      {
         for( v = 0; v < reoptconsdata->nvars; v++ )
         {
            consvars[v] = reoptconsdata->vars[v];
            consval = reoptconsdata->vals[v];
            consboundtype = SCIPsetIsFeasEQ(set, consval, 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;

           assert(SCIPvarIsOriginal(consvars[v]));
           SCIP_CALL( SCIPvarGetProbvarBound(&consvars[v], &consval, &consboundtype) );
           assert(SCIPvarIsTransformed(consvars[v]));
           assert(SCIPvarGetStatus(consvars[v]) != SCIP_VARSTATUS_MULTAGGR);

           if ( SCIPsetIsFeasEQ(set, consval, 1.0) )
           {
              SCIP_CALL( SCIPvarNegate(consvars[v], blkmem, set, stat, &consvars[v]) );
              assert(SCIPvarIsNegated(consvars[v]));
           }
         }

         SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, reoptconsdata->nvars, consvars,
               FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
      }
      /* case 2: at least one variable is integer or continuous. we use a bounddisjunction constraint. */
      else
      {
         SCIP_Real* consvals;
         SCIP_BOUNDTYPE* consboundtypes;

         assert(nintvars > 0 || ncontvars > 0);

         /* alloc buffer memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, reoptconsdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consboundtypes, reoptconsdata->nvars) );

         /* iterate over all variable and transform them */
         for( v = 0; v < reoptconsdata->nvars; v++ )
         {
            consvars[v] = reoptconsdata->vars[v];
            consvals[v] = reoptconsdata->vals[v];
            consboundtypes[v] = reoptconsdata->boundtypes[v];

            /* we have to switch the bounds.
             * case 1: integer variable with bound x <= u is transformed to u+1 <= x
             *                                 and l <= x is transformed to   x <= l-1
             * case 2: continuous variable with bound x <= u is transformed to u <= x
             *                                    and l <= x is transformed to x <= l
             */
            if( SCIPvarGetType(consvars[v]) == SCIP_VARTYPE_BINARY
             || SCIPvarGetType(consvars[v]) == SCIP_VARTYPE_INTEGER
             || SCIPvarGetType(consvars[v]) == SCIP_VARTYPE_IMPLINT )
            {
               if( consboundtypes[v] == SCIP_BOUNDTYPE_UPPER )
               {
                  consvals[v] += 1.0;
                  assert(SCIPsetIsLE(set, consvals[v], SCIPvarGetUbGlobal(consvars[v])));
               }
               else
               {
                  consvals[v] -= 1.0;
                  assert(SCIPsetIsGE(set, consvals[v], SCIPvarGetLbGlobal(consvars[v])));
               }
            }

            consboundtypes[v] = (SCIP_BOUNDTYPE)(1 - consboundtypes[v]); /*lint !e641*/

            assert(SCIPvarIsOriginal(consvars[v]));
            SCIP_CALL( SCIPvarGetProbvarBound(&consvars[v], &consvals[v], &consboundtypes[v]) );
            assert(SCIPvarIsTransformed(consvars[v]));
            assert(SCIPvarGetStatus(consvars[v]) != SCIP_VARSTATUS_MULTAGGR);
         }

         /* create the constraints and add them to the corresponding nodes */
         SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, reoptconsdata->nvars, consvars, consboundtypes,
               consvals, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

         /* free buffer memory */
         SCIPfreeBufferArray(scip, &consboundtypes);
         SCIPfreeBufferArray(scip, &consvals);
      }

      SCIPsetDebugMsg(set, " -> add constraint in node #%lld:\n", SCIPnodeGetNumber(node));
#ifdef SCIP_DEBUG_CONSS
      SCIPdebugPrintCons(scip, cons, NULL);
#endif

      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* free buffer */
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** fix all bounds ad stored in dualredscur at the given node @p node_fix */
static
SCIP_RETCODE fixBounds(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< search tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node corresponding to the fixed part */
   unsigned int          id,                 /**< id of stored node */
   SCIP_Bool             updatedualconss     /**< update constraint representing dual bound changes */
   )
{
   SCIP_REOPTTREE* reopttree;
   SCIP_REOPTNODE* reoptnode;
   int v;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(cliquetable != NULL);
   assert(node != NULL);
   assert(blkmem != NULL);

   reopttree = reopt->reopttree;
   assert(reopttree != NULL);
   assert(0 < id && id < reopttree->reoptnodessize);

   reoptnode = reopttree->reoptnodes[id];
   assert(reoptnode != NULL);
   assert(reoptnode->dualreds);
   assert(reoptnode->dualredscur != NULL);

   /* ensure that the arrays to store the bound changes are large enough */
   SCIP_CALL( reoptnodeCheckMemory(reoptnode, set, blkmem, reoptnode->nvars + reoptnode->dualredscur->nvars, 0, 0) );

   for( v = 0; v < reoptnode->dualredscur->nvars; v++ )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_BOUNDTYPE boundtype;
      SCIP_Bool bndchgd;

      var = reoptnode->dualredscur->vars[v];
      val = reoptnode->dualredscur->vals[v];
      boundtype = reoptnode->dualredscur->boundtypes[v];

      SCIP_CALL(SCIPvarGetProbvarBound(&var, &val, &boundtype));
      assert(SCIPvarIsTransformedOrigvar(var));

      bndchgd = FALSE;

      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsGT(set, val, SCIPvarGetLbLocal(var))
         && SCIPsetIsFeasLE(set, val, SCIPvarGetUbLocal(var)) )
      {
         SCIPvarAdjustLb(var, set, &val);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, val, SCIP_BOUNDTYPE_LOWER, FALSE) );

         bndchgd = TRUE;
      }
      else if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsLT(set, val, SCIPvarGetUbLocal(var))
         && SCIPsetIsFeasGE(set, val, SCIPvarGetLbLocal(var)) )
      {
         SCIPvarAdjustUb(var, set, &val);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, val, SCIP_BOUNDTYPE_UPPER, FALSE) );

         bndchgd = TRUE;
      }
      else if( boundtype != SCIP_BOUNDTYPE_LOWER && boundtype != SCIP_BOUNDTYPE_UPPER )
      {
         SCIPerrorMessage("** Unknown boundtype: %d **\n", boundtype);
         return SCIP_INVALIDDATA;
      }
#ifdef SCIP_MORE_DEBUG
      SCIPsetDebugMsg(set, "  (dual) <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", val);
#endif
      /* add variable and bound to branching path information, because we don't want to delete this data */
      if( bndchgd )
      {
         int pos;
         SCIP_Real constant;
         SCIP_Real scalar;

         pos = reoptnode->nvars;

         reoptnode->vars[pos] = var;
         scalar = 1.0;
         constant = 0.0;
         SCIP_CALL( SCIPvarGetOrigvarSum(&reoptnode->vars[pos], &scalar, &constant) );
         assert(SCIPvarIsOriginal(reoptnode->vars[pos]));

         reoptnode->varbounds[pos] = reoptnode->dualredscur->vals[v];
         reoptnode->varboundtypes[pos] = (SCIPsetIsFeasEQ(set, reoptnode->varbounds[pos], 0.0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
         ++reoptnode->nvars;
      }
   }

   if( updatedualconss )
   {
      /* delete dualredscur and move dualredsnex -> dualredscur */
      SCIP_CALL( reoptnodeUpdateDualConss(reoptnode, blkmem) );
   }

   return SCIP_OKAY;
}

/** fix all bounds corresponding to dual bound changes in a previous iteration in the fashion of interdiction branching;
 *  keep the first negbndchg-1 bound changes as stored in dualredscur and negate the negbndchg-th bound.
 */
static
SCIP_RETCODE fixInterdiction(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< search tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< child node */
   unsigned int          id,                 /**< id of the node */
   int*                  perm,               /**< array of permuted indices */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_Real*            vals,               /**< bounds */
   SCIP_BOUNDTYPE*       boundtypes,         /**< boundtypes */
   int                   nvars,              /**< number of variables */
   int                   negbndchg           /**< index of the variable that should negated */
   )
{
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_BOUNDTYPE boundtype;
   int nbndchgs;
   int v;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(cliquetable != NULL);
   assert(node != NULL);
   assert(perm != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(boundtypes != NULL);
   assert(nvars >= 0);
   assert(blkmem != NULL);
   assert(0 < id && id < reopt->reopttree->reoptnodessize);

#ifndef NDEBUG
   {
      SCIP_REOPTTREE* reopttree;
      SCIP_REOPTNODE* reoptnode;

      reopttree = reopt->reopttree;
      assert(reopttree != NULL);

      reoptnode = reopttree->reoptnodes[id];
      assert(reoptnode != NULL);
      assert(reoptnode->dualreds);
   }
#endif

   nbndchgs = MIN(negbndchg, nvars);

   /* change the first nbndchg-1 bounds as stored in dualredscur and negate the negbndchg-th bound */
   for( v = 0; v < nbndchgs; v++ )
   {
      var = vars[perm[v]];
      val = vals[perm[v]];
      boundtype = boundtypes[perm[v]];

      SCIP_CALL(SCIPvarGetProbvarBound(&var, &val, &boundtype));
      assert(SCIPvarIsTransformedOrigvar(var));

      /* negate the last bound change */
      if( v == nbndchgs-1 )
      {
         boundtype = (SCIP_BOUNDTYPE)(SCIP_BOUNDTYPE_UPPER - boundtype); /*lint !e656*/
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && boundtype == SCIP_BOUNDTYPE_UPPER )
            val = val - 1.0;
         else if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && boundtype == SCIP_BOUNDTYPE_LOWER )
            val = val + 1.0;
      }

      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsGT(set, val, SCIPvarGetLbLocal(var))
         && SCIPsetIsFeasLE(set, val, SCIPvarGetUbLocal(var)) )
      {
         SCIPvarAdjustLb(var, set, &val);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, val, SCIP_BOUNDTYPE_LOWER, FALSE) );
      }
      else if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsLT(set, val, SCIPvarGetUbLocal(var))
         && SCIPsetIsFeasGE(set, val, SCIPvarGetLbLocal(var)) )
      {
         SCIPvarAdjustUb(var, set, &val);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, val, SCIP_BOUNDTYPE_UPPER, FALSE) );
      }
      else if( boundtype != SCIP_BOUNDTYPE_LOWER && boundtype != SCIP_BOUNDTYPE_UPPER )
      {
         SCIPerrorMessage("** Unknown boundtype: %d **\n", boundtype);
         return SCIP_INVALIDDATA;
      }
#ifdef SCIP_MORE_DEBUG
      SCIPsetDebugMsg(set, "  (dual) <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", val);
#endif
   }

   return SCIP_OKAY;
}

/** add all constraints stored at @p id to the given nodes @p node_fix and @p node_cons */
static
SCIP_RETCODE addLocalConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the branch and bound tree*/
   unsigned int          id                  /**< id of stored node */
   )
{
   int c;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(0 < id && id < reopt->reopttree->reoptnodessize);

   if( reopt->reopttree->reoptnodes[id]->nconss == 0 )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, " -> add %d constraint(s) to node #%lld:\n", reopt->reopttree->reoptnodes[id]->nconss,
      SCIPnodeGetNumber(node));

   for( c = 0; c < reopt->reopttree->reoptnodes[id]->nconss; c++ )
   {
      SCIP_CONS* cons;
      SCIP_REOPTCONSDATA* reoptconsdata;

      reoptconsdata = reopt->reopttree->reoptnodes[id]->conss[c];
      assert(reoptconsdata != NULL);
      assert(reoptconsdata->nvars > 0);
      assert(reoptconsdata->varssize >= reoptconsdata->nvars);

      if( reoptconsdata->constype == REOPT_CONSTYPE_CUT )
         continue;

      if( reoptconsdata->constype == REOPT_CONSTYPE_INFSUBTREE )
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "reopt_inf");
      else if( reoptconsdata->constype == REOPT_CONSTYPE_DUALREDS )
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "reopt_dual");
      else
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "reopt_unkn");

      if( reoptconsdata->linear )
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, reoptconsdata->nvars, reoptconsdata->vars, reoptconsdata->vals,
            reoptconsdata->lhs, reoptconsdata->rhs, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
      }
      else
      {
         assert(reoptconsdata->boundtypes != NULL);
         SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, reoptconsdata->nvars, reoptconsdata->vars, reoptconsdata->boundtypes,
            reoptconsdata->vals, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
      }
#ifdef SCIP_DEBUG_CONSS
      SCIPdebugPrintCons(scip, cons, NULL);
#endif
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}

/** reset the internal statistics at the beginning of a new iteration */
static
void resetStats(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   reopt->lastbranched = -1;
   reopt->currentnode = -1;
   reopt->lastseennode = -1;
   reopt->reopttree->nfeasnodes = 0;
   reopt->reopttree->ninfnodes = 0;
   reopt->reopttree->nprunednodes = 0;
   reopt->reopttree->ncutoffreoptnodes = 0;
}

/** check the stored bound changes of all child nodes for redundancy and infeasibility
 *
 *  Due to strongbranching initialization at node stored at @p id it can happen, that some bound changes stored in the
 *  child nodes of the reoptimization node stored at @p id become redundant or make the subproblem infeasible. in this
 *  method we remove all redundant bound changes and delete infeasible child nodes.
 */
static
SCIP_RETCODE dryBranch(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool*            runagain,           /**< pointer to store of this method should run again */
   unsigned int          id                  /**< id of stored node */
   )
{
   SCIP_REOPTNODE* reoptnode;
   unsigned int* cutoffchilds;
   int ncutoffchilds = 0;
   unsigned int* redchilds;
   int nredchilds = 0;
   int c;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes != NULL);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   reoptnode = reopt->reopttree->reoptnodes[id];

   *runagain = FALSE;

   SCIPsetDebugMsg(set, "start dry branching of node at ID %u\n", id);

   /* allocate buffer arrays */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cutoffchilds, reoptnode->nchilds) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &redchilds, reoptnode->nchilds) );

   /* iterate over all child nodes and check each bound changes
    * for redundancy and conflict */
   for( c = 0; c < reoptnode->nchilds; c++ )
   {
      SCIP_REOPTNODE* child;
      SCIP_Bool cutoff;
      SCIP_Bool redundant;
      int* redundantvars;
      int nredundantvars;
      int v;
      unsigned int childid;

      cutoff = FALSE;
      redundant = FALSE;
      nredundantvars = 0;

      childid = reoptnode->childids[c];
      assert(childid < reopt->reopttree->reoptnodessize);
      child = reopt->reopttree->reoptnodes[childid];
      assert(child != NULL);
#ifdef SCIP_MORE_DEBUG
      SCIPsetDebugMsg(set, "-> check child at ID %d (%d vars, %d conss):\n", childid, child->nvars, child->nconss);
#endif
      if( child->nvars > 0 )
      {
         /* allocate buffer memory to store the redundant variables */
         SCIP_CALL( SCIPsetAllocBufferArray(set, &redundantvars, child->nvars) );

         for( v = 0; v < child->nvars && !cutoff; v++ )
         {
            SCIP_VAR* transvar;
            SCIP_Real transval;
            SCIP_BOUNDTYPE transbndtype;
            SCIP_Real ub;
            SCIP_Real lb;

            transvar = child->vars[v];
            transval = child->varbounds[v];
            transbndtype = child->varboundtypes[v];

            /* transform into the transformed space */
            SCIP_CALL( SCIPvarGetProbvarBound(&transvar, &transval, &transbndtype) );

            lb = SCIPvarGetLbLocal(transvar);
            ub = SCIPvarGetUbLocal(transvar);

            /* check for infeasibility */
            if( SCIPsetIsFeasEQ(set, lb, ub) && !SCIPsetIsFeasEQ(set, lb, transval) )
            {
               SCIPsetDebugMsg(set, " -> <%s> is fixed to %g, can not change bound to %g -> cutoff\n",
                  SCIPvarGetName(transvar), lb, transval);

               cutoff = TRUE;
               break;
            }

            /* check for redundancy */
            if( SCIPsetIsFeasEQ(set, lb, ub) && SCIPsetIsFeasEQ(set, lb, transval) )
            {
               SCIPsetDebugMsg(set, " -> <%s> is already fixed to %g -> redundant bound change\n",
                  SCIPvarGetName(transvar), lb);

               redundantvars[nredundantvars] = v;
               ++nredundantvars;
            }
         }

         if( !cutoff && nredundantvars > 0 )
         {
            for( v = 0; v < nredundantvars; v++ )
            {
               /* replace the redundant variable by the last stored variable */
               child->vars[redundantvars[v]] = child->vars[child->nvars-1];
               child->varbounds[redundantvars[v]] = child->varbounds[child->nvars-1];
               child->varboundtypes[redundantvars[v]] = child->varboundtypes[child->nvars-1];
               --child->nvars;
            }
         }

         /* free buffer memory */
         SCIPsetFreeBufferArray(set, &redundantvars);
      }
      else if( child->nconss == 0 )
      {
         redundant = TRUE;
         SCIPsetDebugMsg(set, " -> redundant node found.\n");
      }

      if( cutoff )
      {
         cutoffchilds[ncutoffchilds] = childid;
         ++ncutoffchilds;
      }
      else if( redundant )
      {
         redchilds[nredchilds] = childid;
         ++nredchilds;
      }
   }

   SCIPsetDebugMsg(set, "-> found %d redundant and %d infeasible nodes\n", nredchilds, ncutoffchilds);

   /* delete all nodes that can be cut off */
   while( ncutoffchilds > 0 )
   {
      /* delete the node and the induced subtree */
      SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, cutoffchilds[ncutoffchilds-1], TRUE, TRUE) );

      /* find the position in the childid array */
      c = 0;
      while( reoptnode->childids[c] != cutoffchilds[ncutoffchilds-1] && c < reoptnode->nchilds )
         ++c;
      assert(reoptnode->childids[c] == cutoffchilds[ncutoffchilds-1]);

      /* replace the ID at position c by the last ID */
      reoptnode->childids[c] = reoptnode->childids[reoptnode->nchilds-1];
      --reoptnode->nchilds;

      /* decrease the number of nodes to cutoff */
      --ncutoffchilds;
   }

   /* replace all redundant nodes their child nodes or cutoff the node if it is a leaf */
   while( nredchilds > 0 )
   {
      /* find the position in the childid array */
      c = 0;
      while( reoptnode->childids[c] != redchilds[nredchilds-1] && c < reoptnode->nchilds )
         ++c;
      assert(reoptnode->childids[c] == redchilds[nredchilds-1]);

      /* the node is a leaf and we can cutoff them  */
      if( reopt->reopttree->reoptnodes[redchilds[nredchilds-1]]->nchilds == 0 )
      {
         /* delete the node and the induced subtree */
         SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, redchilds[nredchilds-1], TRUE, TRUE) );

         /* replace the ID at position c by the last ID */
         reoptnode->childids[c] = reoptnode->childids[reoptnode->nchilds-1];
         --reoptnode->nchilds;

         /* decrease the number of redundant nodes */
         --nredchilds;
      }
      else
      {
         int cc;
         int ncc;

         /* replace the ID at position c by the last ID */
         reoptnode->childids[c] = reoptnode->childids[reoptnode->nchilds-1];
         --reoptnode->nchilds;

         ncc = reopt->reopttree->reoptnodes[redchilds[nredchilds-1]]->nchilds;

         /* check the memory */
         SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[id], set, blkmem, 0, reoptnode->nchilds+ncc, 0) );

         /* add all IDs of child nodes to the current node */
         for( cc = 0; cc < ncc; cc++ )
         {
            reoptnode->childids[reoptnode->nchilds] = reopt->reopttree->reoptnodes[redchilds[nredchilds-1]]->childids[cc];
            ++reoptnode->nchilds;
         }

         /* delete the redundant node */
         SCIP_CALL( reopttreeDeleteNode(reopt->reopttree, set, blkmem, redchilds[nredchilds-1], TRUE) );
         SCIP_CALL( SCIPqueueInsert(reopt->reopttree->openids, (void*) (size_t) redchilds[nredchilds-1]) );

         /* decrease the number of redundant nodes */
         --nredchilds;

         /* update the flag to rerun this method */
         *runagain = TRUE;
      }
   }

   /* free buffer arrays */
   SCIPsetFreeBufferArray(set, &redchilds);
   SCIPsetFreeBufferArray(set, &cutoffchilds);

   return SCIP_OKAY;
}

/** return the number of all nodes in the subtree induced by the reoptimization node stored at @p id */
static
int reopttreeGetNNodes(
   SCIP_REOPTTREE*       reopttree,          /**< reopttree */
   unsigned int          id                  /**< id of stored node */
   )
{
   int nnodes = 0;
   int i;

   assert(reopttree != NULL);
   assert(id < reopttree->reoptnodessize);

   for( i = 0; i < reopttree->reoptnodes[id]->nchilds; i++ )
      nnodes += reopttreeGetNNodes(reopttree, reopttree->reoptnodes[id]->childids[i]);

   return nnodes + 1;
}

/** returns the number of leaf nodes of the induced subtree */
static
int reoptGetNLeaves(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   unsigned int          id                  /**< id of stored node */
   )
{
   int i;
   int nleaves = 0;

   assert(reopt != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   /* iterate over all child nods and check whether they are leaves or not */
   for( i = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++ )
   {
      unsigned int childid;

      childid = reopt->reopttree->reoptnodes[id]->childids[i];
      assert(childid < reopt->reopttree->reoptnodessize);

      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
         ++nleaves;
      else
         nleaves += reoptGetNLeaves(reopt, childid);
   }

   return nleaves;
}

/** returns all leaves of the subtree induced by the node stored at @p id*/
static
SCIP_RETCODE reoptGetLeaves(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure*/
   unsigned int          id,                 /**< id of stored node */
   unsigned int*         leaves,             /**< array of leave nodes */
   int                   leavessize,         /**< size of leaves array */
   int*                  nleaves             /**< pointer to store the number of leave nodes */
   )
{
   int i;
   int l;

   assert(reopt != NULL);
   assert(leavessize > 0 && leaves != NULL);
   assert((*nleaves) >= 0);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   for( i = 0, l = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++ )
   {
      unsigned int childid;

      assert(*nleaves <= leavessize);

      childid = reopt->reopttree->reoptnodes[id]->childids[i];
      assert(childid < reopt->reopttree->reoptnodessize);

      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
      {
         leaves[l] = reopt->reopttree->reoptnodes[id]->childids[i];
         ++l;
         ++(*nleaves);
      }
      else
      {
         int nleaves2 = 0;

         SCIP_CALL( reoptGetLeaves(reopt, childid, &leaves[l], leavessize - l, &nleaves2) );
         l += nleaves2;
         (*nleaves) += nleaves2;
      }
   }

   return SCIP_OKAY;
}

/** after restarting the reoptimization and an after compressing the search tree we have to delete all stored information */
static
SCIP_RETCODE reoptResetTree(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             softreset           /**< mark the nodes to overwriteable (TRUE) or delete them completely (FALSE) */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* clear the tree */
   SCIP_CALL( clearReoptnodes(reopt->reopttree, set, blkmem, softreset) );
   assert(reopt->reopttree->nreoptnodes == 0);

   /* reset the dual constraint */
   if( reopt->dualreds != NULL )
      reopt->dualreds->nvars = 0;

   reopt->currentnode = -1;

   return SCIP_OKAY;
}

/** restart the reoptimization by removing all stored information about nodes and increase the number of restarts */
static
SCIP_RETCODE reoptRestart(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* clear the tree */
   SCIP_CALL( reoptResetTree(reopt, set, blkmem, FALSE) );
   assert(reopt->reopttree->nreoptnodes == 0);

   /* allocate memory for the root node */
   SCIP_CALL( createReoptnode(reopt->reopttree, set, blkmem, 0) );

   reopt->nglbrestarts += 1;

   if( reopt->firstrestart == -1 )
      reopt->firstrestart = reopt->run;

   reopt->lastrestart = reopt->run;

   return SCIP_OKAY;
}

/** save the new objective function */
static
SCIP_RETCODE reoptSaveNewObj(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            origvars,           /**< original problem variables */
   int                   norigvars           /**< number of original problem variables */
   )
{
   int probidx;
   int v;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(origvars != NULL);
   assert(norigvars >= 0);

   /* check memory */
   SCIP_CALL( ensureRunSize(reopt, set, reopt->run, blkmem) );

   /* get memory and check whether we have to resize all previous objectives */
   if( reopt->nobjvars < norigvars )
   {
      int i;
      for( i = 0; i < reopt->run-1; i++ )
      {
         SCIP_ALLOC( BMSreallocMemoryArray(&reopt->objs[i], norigvars) ); /*lint !e866*/
         for( v = reopt->nobjvars-1; v < norigvars; v++ )
            reopt->objs[i][v] = 0.0;
      }
      reopt->nobjvars = norigvars;
   }
   SCIP_ALLOC( BMSallocClearMemoryArray(&reopt->objs[reopt->run-1], reopt->nobjvars) ); /*lint !e866*/

   /* save coefficients */
   for( v = 0; v < norigvars; v++ )
   {
      assert(SCIPvarIsOriginal(origvars[v]));

      probidx = SCIPvarGetIndex(origvars[v]);

      /* it can happen that the index is greater than the number of problem variables,
       * i.e., not all created variables were added
       */
      if( probidx >= reopt->nobjvars )
      {
         int i;
         int j;
         int newsize = SCIPsetCalcMemGrowSize(set, probidx+1);
         for( i = 0; i < reopt->run; i++ )
         {
            SCIP_ALLOC( BMSreallocMemoryArray(&reopt->objs[i], newsize) ); /*lint !e866*/
            for( j = reopt->nobjvars; j < newsize; j++ )
               reopt->objs[i][j] = 0.0;
         }
         reopt->nobjvars = newsize;
      }
      assert(0 <= probidx && probidx < reopt->nobjvars);

      reopt->objs[reopt->run-1][probidx] = SCIPvarGetObj(origvars[v]);

      /* update flag to remember if the objective function has changed */
      if( !reopt->objhaschanged && reopt->run >= 2
          && ! SCIPsetIsEQ(set, reopt->objs[reopt->run-2][probidx], reopt->objs[reopt->run-1][probidx]) )
         reopt->objhaschanged = TRUE;

      /* mark this objective as the first non empty */
      if( reopt->firstobj == -1 && reopt->objs[reopt->run-1][probidx] != 0 )
         reopt->firstobj = reopt->run-1;
   }

   /* calculate similarity to last objective */
   if( reopt->run-1 >= 1 )
   {
      /* calculate similarity to last objective */
      reopt->simtolastobj = reoptSimilarity(reopt, set, reopt->run-1, reopt->run-2, origvars, norigvars);

      if( reopt->simtolastobj == SCIP_INVALID )  /*lint !e777*/
         return SCIP_INVALIDRESULT;

      SCIPverbMessage(set->scip, SCIP_VERBLEVEL_HIGH, NULL, "new objective has similarity of %g compared to previous.\n",
         reopt->simtolastobj);
   }

   SCIPsetDebugMsg(set, "saved obj for run %d.\n", reopt->run);

   return SCIP_OKAY;
}

/** orders the variable by inference score */
static
SCIP_RETCODE getInferenceOrder(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int*                  perm,               /**< array of indices that need to be permuted */
   SCIP_VAR**            vars,               /**< variable array to permute */
   SCIP_Real*            bounds,             /**< bound array to permute in the same order */
   SCIP_BOUNDTYPE*       boundtypes,         /**< boundtype array to permute in the same order */
   int                   nvars               /**< number of variables */
   )
{
   SCIP_Real* infscore;
   int v;

   assert(set != NULL);
   assert(perm != NULL);
   assert(vars != NULL);
   assert(bounds != NULL);
   assert(boundtypes != NULL);
   assert(nvars >= 0);

   /* allocate buffer for the scores */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &infscore, nvars) );

   for( v = 0; v < nvars; v++ )
   {
      if( boundtypes[v] == SCIP_BOUNDTYPE_UPPER )
      {
         infscore[v] = 0.75 * SCIPvarGetAvgInferences(vars[v], stat, SCIP_BRANCHDIR_UPWARDS)
            + 0.25 * SCIPvarGetAvgInferences(vars[v], stat, SCIP_BRANCHDIR_DOWNWARDS);
      }
      else
      {
         infscore[v] = 0.25 * SCIPvarGetAvgInferences(vars[v], stat, SCIP_BRANCHDIR_UPWARDS)
               + 0.75 * SCIPvarGetAvgInferences(vars[v], stat, SCIP_BRANCHDIR_DOWNWARDS);
      }
   }

   /* permute indices by inference score */
   SCIPsortDownRealInt(infscore, perm, nvars);

   /* free buffer */
   SCIPsetFreeBufferArray(set, &infscore);

   return SCIP_OKAY;
}

/** create a global constraint to separate the given solution */
static
SCIP_RETCODE separateSolution(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_VAR**            vars,               /**< array of original problem variables */
   int                   nvars               /**< number of original problem variables */
   )
{
   SCIP_VAR** origvars;
   SCIP_Real* vals;
   int nintvars;
   int nbinvars;
   int v;
   int w;

   assert(reopt != NULL);
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(vars != NULL);
   assert(nvars != 0);
   assert(SCIPsolIsOriginal(sol));

   /* allocate buffer memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &origvars, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, nvars) );

   nbinvars = 0;
   nintvars = 0;

   /* get the solution values of the variables */
   for( v = 0, w = 0; v < nvars; v++ )
   {
      assert(SCIPvarIsOriginal(vars[v]));
      assert(nbinvars + nintvars == w);

      /* we do not want to create cuts for continous variables */
      if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY )
         ++nbinvars;
      if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_IMPLINT )
         ++nintvars;

      origvars[v] = vars[v];
      assert(origvars[v] != NULL);
      assert(SCIPvarIsOriginal(origvars[v]));

      vals[w] = SCIPsolGetVal(sol, set, stat, origvars[v]);
      ++w;
   }

   SCIP_CALL( addGlobalCut(reopt, blkmem, set, origvars, vals, NULL, w, nbinvars, nintvars) );

   /* free buffer memory */
   SCIPsetFreeBufferArray(set, &vals);
   SCIPsetFreeBufferArray(set, &origvars);

   return SCIP_OKAY;
}

/*
 * public methods
 */

/* ---------------- methods of general reoptimization ---------------- */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPreoptGetNRestartsGlobal
#undef SCIPreoptGetNRestartsLocal
#undef SCIPreoptGetNTotalRestartsLocal
#undef SCIPreoptGetFirstRestarts
#undef SCIPreoptGetLastRestarts
#undef SCIPreoptGetNFeasNodes
#undef SCIPreoptGetNTotalFeasNodes
#undef SCIPreoptGetNPrunedNodes
#undef SCIPreoptGetNTotalPrunedNodes
#undef SCIPreoptGetNCutoffReoptnodes
#undef SCIPreoptGetNTotalCutoffReoptnodes
#undef SCIPreoptGetNInfNodes
#undef SCIPreoptGetNTotalInfNodes
#undef SCIPreoptGetNInfSubtrees


/** returns the number of global restarts */
int SCIPreoptGetNRestartsGlobal(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->nglbrestarts;
}

/** returns the number of local restarts in the current run */
int SCIPreoptGetNRestartsLocal(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->nlocrestarts;
}

/** returns the number of local restarts over all runs */
int SCIPreoptGetNTotalRestartsLocal(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->ntotallocrestarts;
}

/** returns the number of iteration with the first global restarts */
int SCIPreoptGetFirstRestarts(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->firstrestart;
}

/** returns the number of iteration with the last global restarts */
int SCIPreoptGetLastRestarts(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->lastrestart;
}

/** returns the number of stored nodes providing an improving feasible LP solution in the current run */
int SCIPreoptGetNFeasNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->nfeasnodes;
}

/** returns the number of stored nodes providing an improving feasible LP solution over all runs */
int SCIPreoptGetNTotalFeasNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ntotalfeasnodes;
}

/** returns the number of stored nodes that exceeded the cutoff bound in the current run */
int SCIPreoptGetNPrunedNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->nprunednodes;
}

/** returns the number of stored nodes that exceeded the cutoff bound over all runs */
int SCIPreoptGetNTotalPrunedNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ntotalprunednodes;
}

/** rerturns the number of reoptimized nodes that were cutoff in the same iteration in the current run */
int SCIPreoptGetNCutoffReoptnodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ncutoffreoptnodes;
}

/** rerturns the number of reoptimized nodes that were cutoff in the same iteration over all runs */
int SCIPreoptGetNTotalCutoffReoptnodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ntotalcutoffreoptnodes;
}

/** returns the number of stored nodes with an infeasible LP in the current run */
int SCIPreoptGetNInfNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ninfnodes;
}

/** returns the number of stored nodes with an infeasible LP over all runs */
int SCIPreoptGetNTotalInfNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ntotalinfnodes;
}

/** constructor for the reoptimization data */
SCIP_RETCODE SCIPreoptCreate(
   SCIP_REOPT**          reopt,              /**< pointer to reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   int i;

   assert(reopt != NULL);

   SCIP_ALLOC( BMSallocMemory(reopt) );
   (*reopt)->runsize = DEFAULT_MEM_RUN;
   (*reopt)->run = 0;
   (*reopt)->simtolastobj = -2.0;
   (*reopt)->simtofirstobj = -2.0;
   (*reopt)->firstobj = -1;
   (*reopt)->currentnode = -1;
   (*reopt)->lastbranched = -1;
   (*reopt)->dualreds = NULL;
   (*reopt)->glbconss = NULL;
   (*reopt)->nglbconss = 0;
   (*reopt)->allocmemglbconss = 0;
   (*reopt)->ncheckedsols = 0;
   (*reopt)->nimprovingsols = 0;
   (*reopt)->noptsolsbyreoptsol = 0;
   (*reopt)->nglbrestarts = 0;
   (*reopt)->nlocrestarts = 0;
   (*reopt)->ntotallocrestarts = 0;
   (*reopt)->firstrestart = -1;
   (*reopt)->lastrestart = 0;
   (*reopt)->nobjvars = 0;
   (*reopt)->objhaschanged = FALSE;
   (*reopt)->consadded = FALSE;
   (*reopt)->addedconss = NULL;
   (*reopt)->naddedconss = 0;
   (*reopt)->addedconsssize = 0;
   (*reopt)->glblb = NULL;
   (*reopt)->glbub = NULL;
   (*reopt)->activeconss = NULL;

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*reopt)->varhistory, (*reopt)->runsize) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*reopt)->prevbestsols, (*reopt)->runsize) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*reopt)->objs, (*reopt)->runsize) );

   for( i = 0; i < (*reopt)->runsize; i++ )
   {
      (*reopt)->objs[i] = NULL;
      (*reopt)->prevbestsols[i] = NULL;
      (*reopt)->varhistory[i] = NULL;
   }

   /* clocks */
   SCIP_CALL( SCIPclockCreate(&(*reopt)->savingtime, SCIP_CLOCKTYPE_DEFAULT) );

   /* create and initialize SCIP_SOLTREE */
   SCIP_ALLOC( BMSallocMemory(&(*reopt)->soltree) );
   SCIP_CALL( createSolTree((*reopt)->soltree, blkmem) );

   /* create and initialize SCIP_REOPTTREE */
   SCIP_ALLOC( BMSallocMemory(&(*reopt)->reopttree) );
   SCIP_CALL( createReopttree((*reopt)->reopttree, set, blkmem) );

   /* create a random number generator */
   SCIP_CALL( SCIPrandomCreate(&(*reopt)->randnumgen, blkmem, (unsigned int)SCIPsetInitializeRandomSeed(set, DEFAULT_RANDSEED)) );

   /* create event handler for node events */
   eventhdlr = NULL;

   /* include event handler into SCIP */
   SCIP_CALL( SCIPeventhdlrCreate(&eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, NULL, NULL, NULL, NULL, eventInitsolReopt,
         eventExitsolReopt, NULL, eventExecReopt, NULL) );
   SCIP_CALL( SCIPsetIncludeEventhdlr(set, eventhdlr) );
   assert(eventhdlr != NULL);

   return SCIP_OKAY;
}

/* release all variables and constraints captured during reoptimization */
SCIP_RETCODE SCIPreoptReleaseData(
   SCIP_REOPT*           reopt,              /**< pointer to reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   /* release all added constraints and free the data */
   if( reopt->addedconss != NULL )
   {
      int c;
      for( c = 0; c < reopt->naddedconss; c++)
      {
         assert(reopt->addedconss[c] != NULL);

         SCIP_CALL( SCIPconsRelease(&reopt->addedconss[c], blkmem, set) );
      }

      BMSfreeBlockMemoryArray(blkmem, &reopt->addedconss, reopt->addedconsssize);
   }

   SCIP_CALL( cleanActiveConss(reopt, set) );

   return SCIP_OKAY;
}

/** frees reoptimization data */
SCIP_RETCODE SCIPreoptFree(
   SCIP_REOPT**          reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          origprimal,         /**< original primal */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(*reopt != NULL);
   assert(set != NULL);
   assert(origprimal != NULL || set->stage == SCIP_STAGE_INIT);
   assert(blkmem != NULL);

   /* free random number generator */
   SCIPrandomFree(&(*reopt)->randnumgen, blkmem);

   /* free reopttree */
   SCIP_CALL( freeReoptTree((*reopt)->reopttree, set, blkmem) );

   /* free solutions */
   if( set->stage >= SCIP_STAGE_PROBLEM )
   {
      int p;
      for( p = (*reopt)->run-1; p >= 0; p-- )
      {
         if( (*reopt)->soltree->sols[p] != NULL )
         {
            BMSfreeBlockMemoryArray(blkmem, &(*reopt)->soltree->sols[p], (*reopt)->soltree->solssize[p]); /*lint !e866*/
            (*reopt)->soltree->sols[p] = NULL;
         }

         /* we have to free all optimal solution separatly, because those solutions are not stored in the
          * solution reopt_sepabestsol = TRUE
          */
         if( set->reopt_sepabestsol && (*reopt)->prevbestsols[p] != NULL )
         {
            SCIP_CALL( SCIPsolFree(&(*reopt)->prevbestsols[p], blkmem, origprimal) );
         }

         if( (*reopt)->objs[p] != NULL )
         {
            BMSfreeMemoryArray(&(*reopt)->objs[p]);
         }
      }
   }

   /* free solution tree */
   SCIP_CALL( freeSolTree((*reopt), set, origprimal, blkmem) );

   if( (*reopt)->dualreds != NULL )
   {
      if( (*reopt)->dualreds->varssize > 0 )
      {
         assert(!(*reopt)->dualreds->linear);

         BMSfreeBlockMemoryArray(blkmem, &(*reopt)->dualreds->boundtypes, (*reopt)->dualreds->varssize);
         BMSfreeBlockMemoryArray(blkmem, &(*reopt)->dualreds->vals, (*reopt)->dualreds->varssize);
         BMSfreeBlockMemoryArray(blkmem, &(*reopt)->dualreds->vars, (*reopt)->dualreds->varssize);
         BMSfreeBlockMemory(blkmem, &(*reopt)->dualreds);
         (*reopt)->dualreds = NULL;
      }
   }

   if( (*reopt)->glbconss != NULL && (*reopt)->allocmemglbconss > 0 )
   {
      int c;

      /* free all constraint */
      for( c = 0; c < (*reopt)->allocmemglbconss; c++ )
      {
         if( (*reopt)->glbconss[c] != NULL )
         {
            if( (*reopt)->glbconss[c]->varssize > 0 )
            {
               BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss[c]->boundtypes, (*reopt)->glbconss[c]->varssize);
               BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss[c]->vals, (*reopt)->glbconss[c]->varssize);
               BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss[c]->vars, (*reopt)->glbconss[c]->varssize);
               (*reopt)->glbconss[c]->varssize = 0;
            }
            BMSfreeBlockMemory(blkmem, &(*reopt)->glbconss[c]); /*lint !e866*/
            --(*reopt)->nglbconss;
         }

      }
      assert((*reopt)->nglbconss == 0);

      BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss, (*reopt)->allocmemglbconss);
      (*reopt)->allocmemglbconss = 0;
   }

   /* clocks */
   SCIPclockFree(&(*reopt)->savingtime);

   SCIPhashmapFree(&(*reopt)->activeconss);
   (*reopt)->activeconss = NULL;

   if( (*reopt)->glblb != NULL )
   {
      SCIPhashmapFree(&(*reopt)->glblb);
      SCIPhashmapFree(&(*reopt)->glbub);
      (*reopt)->glblb = NULL;
      (*reopt)->glbub = NULL;
   }
   else
      assert((*reopt)->glbub == NULL);

   BMSfreeBlockMemoryArray(blkmem, &(*reopt)->varhistory, (*reopt)->runsize);
   BMSfreeBlockMemoryArray(blkmem, &(*reopt)->prevbestsols, (*reopt)->runsize);
   BMSfreeMemoryArray(&(*reopt)->objs);
   BMSfreeMemory(reopt);

   return SCIP_OKAY;
}

/** returns the number of constraints added by the reoptimization plug-in */
int SCIPreoptGetNAddedConss(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   )
{
   unsigned int id;

   assert(reopt != NULL);
   assert(node != NULL);

   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* set the id to -1 if the node is not part of the reoptimization tree */
   if( SCIPnodeGetDepth(node) > 0 && id == 0 )
      return SCIPnodeGetNAddedConss(node);

   if( id >= 1 && reopt->reopttree->reoptnodes[id]->nconss > 0 )
      return MAX(SCIPnodeGetNAddedConss(node), reopt->reopttree->reoptnodes[id]->nconss); /*lint !e666*/
   else
      return SCIPnodeGetNAddedConss(node);
}

/** add a solution to the solution tree */
SCIP_RETCODE SCIPreoptAddSol(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal,         /**< original primal */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SOL*             sol,                /**< solution to add */
   SCIP_Bool             bestsol,            /**< is the current solution an optimal solution? */
   SCIP_Bool*            added,              /**< pointer to store the information if the soltion was added */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of variables */
   int                   run                 /**< number of the current run (1,2,...) */
   )
{
   SCIP_SOLNODE* solnode = NULL;
   SCIP_HEUR* heur;
   int insertpos;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(sol != NULL);
   assert(run > 0);

   assert(reopt->soltree->sols[run-1] != NULL);

   /* if the solution was found by reoptsols the solutions is already stored */
   heur = SCIPsolGetHeur(sol);
   if( heur != NULL && strcmp(SCIPheurGetName(heur), "reoptsols") == 0 && bestsol )
      ++reopt->noptsolsbyreoptsol;
   else if( bestsol )
      reopt->noptsolsbyreoptsol = 0;

   /* check memory */
   SCIP_CALL( ensureSolsSize(reopt, set, blkmem, reopt->soltree->nsols[run-1]+1, run-1) );

   /* add solution to solution tree */
   SCIP_CALL( soltreeAddSol(reopt, set, stat, origprimal, blkmem, vars, sol, &solnode, nvars, bestsol, added) );

   if( (*added) )
   {
      assert(solnode != NULL);

      /* add solution */
      insertpos = reopt->soltree->nsols[run-1];
      reopt->soltree->sols[run-1][insertpos] = solnode;
      ++reopt->soltree->nsols[run-1];
      assert(reopt->soltree->nsols[run-1] <= set->reopt_savesols);
   }

   return SCIP_OKAY;
}

/** we want to store the optimal solution of each run in a separate array */
SCIP_RETCODE SCIPreoptAddOptSol(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SOL*             sol,                /**< solution to add */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal,         /**< original primal */
   SCIP_VAR**            vars,               /**< original problem variables */
   int                   nvars               /**< number of original problem variables */
   )
{
   /* cppcheck-suppress unassignedVariable */
   SCIP_SOL* solcopy;

   assert(reopt != NULL);
   assert(reopt->run-1 >= 0);
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(origprimal != NULL);

   SCIP_CALL( SCIPsolCopy(&solcopy, blkmem, set, stat, origprimal, sol) );
   reopt->prevbestsols[reopt->run-1] = solcopy;

   /* store a global constraint that cutsoff the solution */
   if( set->reopt_sepabestsol )
   {
      SCIP_CALL( separateSolution(reopt, blkmem, set, stat, sol, vars, nvars) );
   }

   return SCIP_OKAY;
}

/** add a new iteration after changing the objective function */
SCIP_RETCODE SCIPreoptAddRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data sturcture */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            origvars,           /**< original problem variables */
   int                   norigvars,          /**< number of original variables */
   int                   size                /**< number of expected solutions */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem !=  NULL);
   assert(origvars != NULL);

   /* increase number of runs */
   ++reopt->run;

   /* check memory */
   SCIP_CALL( ensureRunSize(reopt, set, reopt->run, blkmem) );

   /* allocate memory */
   reopt->soltree->solssize[reopt->run-1] = size;
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->soltree->sols[reopt->run-1], size) ); /*lint !e866*/

   /* reset flag */
   reopt->objhaschanged = FALSE;

   /* save the objective function */
   SCIP_CALL( reoptSaveNewObj(reopt, set, blkmem, origvars, norigvars) );

   resetStats(reopt);

   return SCIP_OKAY;
}

/** get the number of checked solutions during the reoptimization process */
int SCIPreoptGetNCheckedSols(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->ncheckedsols;
}

/** update the number of checked solutions during the reoptimization process */
void SCIPreoptAddNCheckedSols(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   ncheckedsols        /**< number of updated solutions */
   )
{
   assert(reopt != NULL);

   reopt->ncheckedsols += ncheckedsols;
}

/** get the number of checked solutions during the reoptimization process */
int SCIPreoptGetNImprovingSols(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->nimprovingsols;
}

/** update the number of checked solutions during the reoptimization process */
void SCIPreoptAddNImprovingSols(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   nimprovingsols      /**< number of improving solutions */
   )
{
   assert(reopt != NULL);

   reopt->nimprovingsols += nimprovingsols;
}

/** returns number of solutions stored in the solution tree of a given run */
int SCIPreoptGetNSolsRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   run                 /**< number of the run (1,2,..) */
   )
{
   assert(reopt != NULL);
   assert(0 < run && run <= reopt->runsize);

   if( reopt->soltree->sols[run-1] == NULL )
      return 0;
   else
      return reopt->soltree->nsols[run-1];
}

/** returns number of all solutions of all runs */
int SCIPreoptGetNSols(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   int nsols = 0;
   int r;

   assert(reopt != NULL);

   for( r = 0; r < reopt->run; r++)
      nsols += reopt->soltree->nsols[r];

   return nsols;
}

/** return the stored solutions of a given run */
SCIP_RETCODE SCIPreoptGetSolsRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   run,                /**< number of the run (1,2,...) */
   SCIP_SOL**            sols,               /**< array of solutions to fill */
   int                   solssize,           /**< length of the array */
   int*                  nsols               /**< pointer to store the number of added solutions */
   )
{
   int s;

   assert(reopt != NULL);
   assert(run > 0 && run <= reopt->run);
   assert(sols != NULL);

   assert(solssize > 0);
   assert(nsols != NULL);
   *nsols = 0;

   for( s = 0; s < reopt->soltree->nsols[run-1]; s++ )
   {
      if( !reopt->soltree->sols[run-1][s]->updated )
         ++(*nsols);
   }

   if( solssize < (*nsols) )
      return SCIP_OKAY;

   (*nsols) = 0;
   for( s = 0; s < reopt->soltree->nsols[run-1]; s++ )
   {
      if( !reopt->soltree->sols[run-1][s]->updated )
      {
         sols[*nsols] = reopt->soltree->sols[run-1][s]->sol;
         reopt->soltree->sols[run-1][s]->updated = TRUE;
         ++(*nsols);
      }
   }

   return SCIP_OKAY;
}

/** returns the number of saved solutions overall runs */
int SCIPreoptGetNSavedSols(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   int nsavedsols = 0;

   assert(reopt != NULL);
   assert(reopt->soltree->root != NULL);

   if( reopt->soltree->root->child != NULL )
      nsavedsols = soltreeNInducedSols(reopt->soltree->root);

   return nsavedsols;
}

/** check if the reoptimization process should be (locally) restarted.
 *
 *  First, we check whether the current node is the root node, e.g., node == NULL. in this case, we do not need to calculate
 *  the similarity again. we trigger a restart if
 *    1. the objective function has changed too much
 *    2. the number of stored nodes is exceeded
 *    3. the last n optimal solutions were found by heur_reoptsols (in this case, the stored tree was only needed to
 *       prove the optimality and this can be probably faster by solving from scratch)
 *
 *  If the current node is different to the root node we calculate the local similarity, i.e., exclude all variable
 *  that are already fixed by bounding.
 */
SCIP_RETCODE SCIPreoptCheckRestart(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< current node of the branch and bound tree (or NULL) */
   SCIP_VAR**            transvars,          /**< transformed problem variables */
   int                   ntransvars,         /**< number of transformed problem variables */
   SCIP_Bool*            restart             /**< pointer to store if the reoptimization process should be restarted */
   )
{
   SCIP_Real sim = 1.0;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(transvars != NULL);
   assert(ntransvars >= 0);
   assert(restart != NULL);

   *restart = FALSE;

   /* check if the whole reoptimization process should start from scratch */
   if( node == NULL )
   {
      /* compute the similarity to the objective function of the first run after restarting */
      if( reopt->run > 1 && set->reopt_objsimdelay > -1.0 )
      {
         sim = reoptSimilarity(reopt, set, reopt->run-1, MAX(0, reopt->lastrestart-1), transvars, ntransvars);

         if( sim == SCIP_INVALID )  /*lint !e777*/
            return SCIP_INVALIDRESULT;
      }

      /* check similarity */
      if( SCIPsetIsFeasLT(set, sim, set->reopt_objsimdelay) )
      {
         SCIPsetDebugMsg(set, "-> restart reoptimization (objective functions are not similar enough)\n");
         *restart = TRUE;
      }
      /* check size of the reoptimization tree */
      else if( reopt->reopttree->nreoptnodes > set->reopt_maxsavednodes )
      {
         SCIPsetDebugMsg(set, "-> restart reoptimization (node limit reached)\n");
         *restart = TRUE;
      }
      /* check if the tree was only needed to prove optimality */
      else if( reopt->noptsolsbyreoptsol >= set->reopt_forceheurrestart )
      {
         SCIPsetDebugMsg(set, "-> restart reoptimization (found last %d optimal solutions by <reoptsols>)\n",
               reopt->noptsolsbyreoptsol);
         reopt->noptsolsbyreoptsol = 0;
         *restart = TRUE;
      }

      if( *restart )
      {
         /* trigger a restart */
         SCIP_CALL( reoptRestart(reopt, set, blkmem) );
      }
   }
   /* check for a local restart, ie, start the solving process of an inner node from scatch */
   else
   {
      SCIP_CALL( reoptCheckLocalRestart(reopt, set, blkmem, node, transvars, ntransvars, restart) );
   }
   return SCIP_OKAY;
}

/** returns the similarity to the previous objective function, if no exist return -2.0 */
SCIP_Real SCIPreoptGetSimToPrevious(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);
   return reopt->simtolastobj;
}

/** returns the similarity to the first objective different to the zero-function function, if no exist return -2.0 */
SCIP_Real SCIPreoptGetSimToFirst(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);
   return reopt->simtofirstobj;
}

/** return the similarity between two of objective functions of two given runs */
SCIP_Real SCIPreoptGetSimilarity(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   run1,               /**< number of the first run */
   int                   run2,               /**< number of the second run */
   SCIP_VAR**            origvars,           /**< original problem variables */
   int                   norigvars           /**< number of original problem variables */
   )
{
   assert(reopt != NULL);
   assert(run1 > 0 && run1 <= reopt->run);
   assert(run2 > 0 && run2 <= reopt->run);
   assert(origvars != NULL);
   assert(norigvars >= 0);

   return reoptSimilarity(reopt, set, run1-1, run2-1, origvars, norigvars);
}

/** returns the best solution of the last run */
SCIP_SOL* SCIPreoptGetLastBestSol(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);
   assert(reopt->prevbestsols != NULL);

   if( reopt->run-2 < 0 )
      return NULL;
   else
      return reopt->prevbestsols[reopt->run-2];
}

/** returns the node of the reoptimization tree corresponding to the unique @p id */
SCIP_REOPTNODE* SCIPreoptGetReoptnode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   unsigned int          id                  /**< unique id */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   return reopt->reopttree->reoptnodes[id];
}

/** returns the coefficient of variable with index @p idx in run @p run */
SCIP_Real SCIPreoptGetOldObjCoef(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   run,                /**< number of the run (1,2,...) */
   int                   idx                 /**< index of original variable */
   )
{
   assert(reopt != NULL);
   assert(0 < run && run <= reopt->runsize);

   return reopt->objs[run-1][idx];
}

/** return the best solution of a given run.
 *
 *  @note the returned solution is part of the original space.
 */
SCIP_SOL* SCIPreoptGetBestSolRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   run                 /**< number of the run (1,2,...) */
   )
{
   assert(reopt != NULL);
   assert(0 < run && run <= reopt->run);

   return reopt->prevbestsols[run-1];
}

/** reset solving specific parameters */
SCIP_RETCODE SCIPreoptReset(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int c;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* clean addedconss array */
   for( c = 0; c < reopt->naddedconss; c++)
   {
      SCIP_CONS* cons;

      cons = reopt->addedconss[c];
      assert(cons != NULL);

      SCIP_CALL( SCIPconsRelease(&cons, blkmem, set) );
      reopt->addedconss[c] = NULL;
   }

   reopt->naddedconss = 0;
   reopt->consadded = FALSE;
   reopt->objhaschanged = FALSE;

   return SCIP_OKAY;
}

/** reset marks of stored solutions to not updated */
void SCIPreoptResetSolMarks(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   SCIP_SOLNODE* child;

   assert(reopt != NULL);
   assert(reopt->soltree != NULL);
   assert(reopt->soltree->root != NULL);

   child = reopt->soltree->root->child;

   /* traverse through the list */
   while( child != NULL )
   {
      soltreeResetMarks(child);
      child = child->sibling;
   }
}

/** returns the number of stored nodes in the subtree induced by @p node */
int SCIPreoptGetNNodes(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   )
{
   unsigned int id;

   assert(reopt != NULL);

   if( node == NULL || SCIPnodeGetDepth(node) == 0 )
      return reopt->reopttree->nreoptnodes;

   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* set the id to -1 if the node is not part of the reoptimization tree */
   if( SCIPnodeGetDepth(node) > 0 && id == 0 )
      return 0;

   assert(0 < id && id < reopt->reopttree->reoptnodessize);

   return reopttreeGetNNodes(reopt->reopttree, id);
}

/* ---------------- methods of general reoptimization nodes ---------------- */

/** In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPreoptnodeGetNVars
#undef SCIPreoptnodeGetNConss
#undef SCIPreoptnodeGetNDualBoundChgs
#undef SCIPreoptnodeGetNChildren
#undef SCIPreoptnodeGetLowerbound
#undef SCIPreoptnodeGetType

/** returns the number of bound changes stored in the reopttree at ID id */
int SCIPreoptnodeGetNVars(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reopttree */
   )
{
   assert(reoptnode != NULL);

   return reoptnode->nvars + reoptnode->nafterdualvars;
}

/** returns the number of bound changes at the node stored at ID id */
int SCIPreoptnodeGetNConss(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   )
{
   assert(reoptnode != NULL);

   return reoptnode->nconss;
}

/** returns the number of stored bound changes based on dual information in the reopttree at ID id */
int SCIPreoptnodeGetNDualBoundChgs(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   )
{
   assert(reoptnode != NULL);

   if( reoptnode->dualredscur == NULL )
      return 0;
   else
      return reoptnode->dualredscur->nvars;
}

/** returns the number of child nodes of @p reoptnode */
int SCIPreoptnodeGetNChildren(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   )
{
   assert(reoptnode != NULL);

   return reoptnode->nchilds;
}

/** return the lower bound stored at @p ID id */
SCIP_Real SCIPreoptnodeGetLowerbound(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   )
{
   assert(reoptnode != NULL);

   return reoptnode->lowerbound;
}

/** returns the type of the @p reoptnode */
SCIP_REOPTTYPE SCIPreoptnodeGetType(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   )
{
   assert(reoptnode != NULL);

   return (SCIP_REOPTTYPE)reoptnode->reopttype;
}

/** returns all added constraints at ID id */
void SCIPreoptnodeGetConss(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR***           vars,               /**< 2-dim array of variables */
   SCIP_Real**           bounds,             /**< 2-dim array of bounds */
   SCIP_BOUNDTYPE**      boundtypes,         /**< 2-dim array of boundtypes */
   int                   mem,                /**< allocated memory for constraints */
   int*                  nconss,             /**< pointer to store the number of constraints */
   int*                  nvars               /**< pointer to store the number of variables */
   )
{
   int c;

   assert(reoptnode != NULL);
   assert(vars != NULL);
   assert(bounds != NULL);
   assert(boundtypes != NULL);
   assert(nvars != NULL);
   assert(nconss != NULL);

   (*nconss) = reoptnode->nconss;

   if( mem < *nconss )
      return;

   for( c = 0; c < *nconss; c++ )
   {
      assert(vars[c] != NULL);
      assert(bounds[c] != NULL);

      vars[c] = reoptnode->conss[c]->vars;
      bounds[c] = reoptnode->conss[c]->vals;
      boundtypes[c] = reoptnode->conss[c]->boundtypes;
      nvars[c] = reoptnode->conss[c]->nvars;
   }
}

/** set the parent id */
void SCIPreoptnodeSetParentID(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   unsigned int          parentid            /**< id of the parent node */
   )
{
   assert(reoptnode != NULL);
   assert(parentid <= 536870911); /* id can be at most 2^29 - 1 */

   reoptnode->parentID = parentid;
}

/** returns the number of leaf nodes of the subtree induced by @p node (of the whole tree if node == NULL) */
int SCIPreoptGetNLeaves(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node                /**< node of the search tree (or NULL) */
   )
{
   int nleaves = 0;
   unsigned int id;
   int i;

   assert(reopt != NULL);

   id = (node == NULL) ? 0 : SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* return if the node is not part of the reoptimization tree */
   if( node != NULL && SCIPnodeGetDepth(node) > 0 && id == 0 )
      return nleaves;

   for( i = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++ )
   {
      unsigned int childid;

      childid = reopt->reopttree->reoptnodes[id]->childids[i]; /*lint !e713*/
      assert(childid < reopt->reopttree->reoptnodessize);

      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
         ++nleaves;
      else
         nleaves += reoptGetNLeaves(reopt, childid);
   }

   return nleaves;
}

/** save information that given node is infeasible */
SCIP_RETCODE SCIPreoptAddInfNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node                /**< node of the search tree */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);

   if( set->reopt_sepaglbinfsubtrees )
   {
      SCIP_CALL( saveGlobalCons(reopt, set, blkmem, node, REOPT_CONSTYPE_CUT) );
   }

   ++reopt->reopttree->ninfnodes;
   ++reopt->reopttree->ntotalinfnodes;

   return SCIP_OKAY;
}

/** check the reason for cut off a node and if necessary store the node */
SCIP_RETCODE SCIPreoptCheckCutoff(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_EVENTTYPE        eventtype,          /**< eventtype */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPSOLSTAT        lpsolstat,          /**< solution status of the LP */
   SCIP_Bool             isrootnode,         /**< the node is the root */
   SCIP_Bool             isfocusnode,        /**< the node is the current focus node */
   SCIP_Real             lowerbound,         /**< lower bound of the node */
   int                   effectiverootdepth  /**< effective root depth */
   )
{
   SCIP_Bool strongbranched;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(lp != NULL);
   assert(node != NULL);
   assert(eventtype == SCIP_EVENTTYPE_NODEBRANCHED || eventtype == SCIP_EVENTTYPE_NODEFEASIBLE || eventtype == SCIP_EVENTTYPE_NODEINFEASIBLE);

   if( reopt->lastseennode == SCIPnodeGetNumber(node) )
      return SCIP_OKAY;

   /* we do not want to store probing node */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   reopt->lastseennode = SCIPnodeGetNumber(node);

   SCIPsetDebugMsg(set, "catch event %" SCIP_EVENTTYPE_FORMAT " for node %lld (type:%d)\n", eventtype, SCIPnodeGetNumber(node), SCIPnodeGetType(node));

   /* case 1: the current node is the root node
    * we can skip if the root is (in)feasible or branched w/o bound
    * changes based on dual information.
    *
    * case 2: we need to store the current node if it contains
    * bound changes based on dual information or is a leave node
    */
   if( isrootnode )
   {
      if( SCIPreoptGetNDualBndchgs(reopt, node) > 0 )
      {
         goto CHECK;
      }
      else if( eventtype == SCIP_EVENTTYPE_NODEBRANCHED )
      {
         /* store or update the information */
         SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_TRANSIT, FALSE, isrootnode, lowerbound) );
      }
      else if( eventtype == SCIP_EVENTTYPE_NODEFEASIBLE )
      {
         /* delete saved dual information which would lead to split the node in a further iteration */
         SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );

         /* store or update the information */
         SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_FEASIBLE, FALSE, isrootnode, lowerbound) );
      }
      else if( eventtype == SCIP_EVENTTYPE_NODEINFEASIBLE )
      {
         /* delete saved dual information which would lead to split the node in a further iteration */
         SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );

         if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT )
         {
            SCIP_Real cutoffbound = SCIPlpGetCutoffbound(lp);
            lowerbound = MIN(lowerbound, cutoffbound);
         }

         /* store or update the information */
         SCIP_CALL( addNode(reopt, set, lp, blkmem, node, reopt->currentnode == 1 ? SCIP_REOPTTYPE_INFSUBTREE : SCIP_REOPTTYPE_PRUNED, FALSE,
               isrootnode, lowerbound) );
      }

      assert(reopt->currentnode == -1);
      assert(reopt->dualreds == NULL || reopt->dualreds->nvars == 0);

      return SCIP_OKAY;
   }

  CHECK:

   if( effectiverootdepth == SCIPnodeGetDepth(node) )
      strongbranched = SCIPreoptGetNDualBndchgs(reopt, node) > 0 ? TRUE : FALSE;
   else
      strongbranched = SCIPnodeGetNDualBndchgs(node) > 0 ? TRUE : FALSE;

   SCIPsetDebugMsg(set, "check the reason of cutoff for node %lld:\n", SCIPnodeGetNumber(node));
   SCIPsetDebugMsg(set, " -> focusnode       : %s\n", isfocusnode ? "yes" : "no");
   SCIPsetDebugMsg(set, " -> depth           : %d (eff. %d)\n", SCIPnodeGetDepth(node), effectiverootdepth);
   SCIPsetDebugMsg(set, " -> strong branched : %s\n", strongbranched ? "yes" : "no");
   SCIPsetDebugMsg(set, " -> LP lpsolstat    : %d\n", lpsolstat);

   switch( eventtype )
   {
   case SCIP_EVENTTYPE_NODEFEASIBLE:
      /* current node has to be the eventnode */
      assert(isfocusnode);

      SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_FEASIBLE);

      /* delete strong branching information of some exists */
      deleteLastDualBndchgs(reopt);

      SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_FEASIBLE, FALSE, isrootnode, lowerbound) );
      break;

   case SCIP_EVENTTYPE_NODEINFEASIBLE:
      /* We have to check if the current node is the event node.
       * if the current node is not the event node, we have to save this node, else we have to
       * look at LP lpsolstat and decide.
       */
      if( isfocusnode )
      {
         /* An after-branch heuristic says NODEINFEASIBLE, maybe the cutoff bound is reached.
          * because the node is already branched we have all children and can delete this node.
          */
         if( SCIPnodeGetNumber(node) == reopt->lastbranched )
         {
            deleteLastDualBndchgs(reopt);
            break;
         }

         /* If the node is strong branched, we possibly detect an infeasible subtree;
          * otherwise, the whole node is either infeasible or exceeds the cutoff bound.
          */
         if( strongbranched )
         {
            /* 1. the LP is infeasible: the (sub-)node is infeasible and can be discarded
             *    because the LP proves infeasibility. We have to store an infeasible subtree separated by a constraint.
             * 2. the LP exceeds the objective limit or was not solved, we have to store the node and can delete the
             *    strong branching information
             */
            if( lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE )
            {
               /* add a dummy variable, because the bound changes were not global in the sense of effective root depth */
               if( SCIPnodeGetDepth(node) > effectiverootdepth )
               {
                  SCIP_CALL( SCIPreoptAddDualBndchg(reopt, set, blkmem, node, NULL, 0.0, 1.0) );
               }

               SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_INFSUBTREE);
               SCIPsetDebugMsg(set, " -> new constype    : %d\n", REOPT_CONSTYPE_INFSUBTREE);

               /* save the node as a strong branched node */
               SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_INFSUBTREE, FALSE, isrootnode, lowerbound) );
            }
            else
            {
               assert(SCIP_LPSOLSTAT_OBJLIMIT || SCIP_LPSOLSTAT_OPTIMAL || SCIP_LPSOLSTAT_NOTSOLVED);

               /* delete strong branching information if some exists */
               deleteLastDualBndchgs(reopt);

               SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_PRUNED);
               SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_PRUNED, FALSE, isrootnode, lowerbound) );
            }
         }
         else
         {
            /* 1. the LP is infeasible: the whole node is infeasible and can be discarded
             * 2. the LP was not solved or exceeds the objective limit, we have to store the node
             */
            if( lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE )
            {
               SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_INFSUBTREE);
               SCIP_CALL( SCIPreoptAddInfNode(reopt, set, blkmem, node) );
            }
            else
            {
               assert(lpsolstat == SCIP_LPSOLSTAT_NOTSOLVED || lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT
                  || lpsolstat == SCIP_LPSOLSTAT_OPTIMAL);

               if( SCIPreoptGetNAddedConss(reopt, node) > 0 )
               {
                  SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_LOGICORNODE);
                  SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_LOGICORNODE, FALSE, isrootnode, lowerbound) );
               }
               else
               {
                  SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_PRUNED);
                  SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_PRUNED, FALSE, isrootnode, lowerbound) );
               }
            }
         }
      }
      else
      {
         SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_PRUNED);

         /* if the node was created by branch_nodereopt, nothing happens */
         SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_PRUNED, FALSE, isrootnode, lowerbound) );

      }
      break;

   case SCIP_EVENTTYPE_NODEBRANCHED:
      /* current node has to be the eventnode */
      assert(isfocusnode);

      reopt->lastbranched = SCIPnodeGetNumber(node);

      /* we have to check the depth of the current node. if the depth is equal to the effective
       * root depth, then all information about bound changes based on dual information already exists,
       * else we have to look at the domchg-data-structure.
       */
      if (SCIPnodeGetDepth(node) == effectiverootdepth)
      {
         /* Save the node if there are added constraints, because this means the node is a copy create by the
          * reoptimization plug-in and contains at least one logic-or-constraint */
         if( strongbranched )
         {
            SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_STRBRANCHED);
            SCIPsetDebugMsg(set, " -> new constype    : %d\n", REOPT_CONSTYPE_DUALREDS);
            SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_STRBRANCHED, FALSE, isrootnode, lowerbound) );
         }
         else if( SCIPreoptGetNAddedConss(reopt, node) > 0 )
         {
            SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_LOGICORNODE);
            SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_LOGICORNODE, FALSE, isrootnode, lowerbound) );
         }
         else
         {
            SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_TRANSIT);
            SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_TRANSIT, FALSE, isrootnode, lowerbound) );
         }
      }
      else
      {
         /* we only branch on binary variables and var == NULL indicates memory allocation w/o saving information.
          *
          * we have to do this in the following order:
          * 1) all bound-changes are local, thats way we have to mark the node to include bound changes based
          *    on dual information.
          * 2) save or update the node.
          */
         if( strongbranched )
         {
            SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_STRBRANCHED);
            SCIPsetDebugMsg(set, " -> new constype    : %d\n", REOPT_CONSTYPE_DUALREDS);
            SCIP_CALL( SCIPreoptAddDualBndchg(reopt, set, blkmem, node, NULL, 0.0, 1.0) );
            SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_STRBRANCHED, FALSE, isrootnode, lowerbound) );
         }
         else if( SCIPreoptGetNAddedConss(reopt, node) > 0 )
         {
            SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_LOGICORNODE);
            SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_LOGICORNODE, FALSE, isrootnode, lowerbound) );
         }
         else
         {
            SCIPsetDebugMsg(set, " -> new reopttype   : %d\n", SCIP_REOPTTYPE_TRANSIT);
            SCIP_CALL( addNode(reopt, set, lp, blkmem, node, SCIP_REOPTTYPE_TRANSIT, FALSE, isrootnode, lowerbound) );
         }
      }
      break;

   default:
      break;
   }

   assert(reopt->currentnode == -1);
   assert(reopt->dualreds == NULL || reopt->dualreds->nvars == 0);

   return SCIP_OKAY; /*lint !e438*/
}

/** store bound change based on dual information */
SCIP_RETCODE SCIPreoptAddDualBndchg(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             newval,             /**< new bound */
   SCIP_Real             oldval              /**< old bound */
   )
{
   SCIP_Real constant = 0.0;
   SCIP_Real scalar = 1.0;

   assert(reopt != NULL);
   assert(node != NULL);

   /* If var == NULL, we save all information by calling SCIPreoptNodeFinished().
    * In that case, all bound changes were not global and we can find them within the
    * domchg data structure.
    * Otherwise, we allocate memory and store the information.
    */
   if( var != NULL )
   {
      SCIP_BOUNDTYPE boundtype;
      int resizelength;
      int allocmem;

      if( SCIPsetFindBranchrule(set, "relpscost") != NULL )
      {
         SCIP_CALL( SCIPsetGetIntParam(set, "branching/relpscost/maxlookahead", &resizelength) );
      }
      else
         resizelength = 1;

      if( reopt->dualreds == NULL || reopt->dualreds->varssize == 0 )
         allocmem = DEFAULT_MEM_DUALCONS;
      else
         allocmem = reopt->dualreds->nvars + resizelength;

      /* allocate memory of necessary */
      SCIP_CALL( checkMemDualCons(reopt, set, blkmem, allocmem) );

      assert(reopt->dualreds->varssize > 0);
      assert(reopt->dualreds->nvars >= 0);
      assert(reopt->currentnode == -1 || reopt->dualreds->nvars > 0);
      assert((reopt->dualreds->nvars > 0 && reopt->currentnode == SCIPnodeGetNumber(node))
           || reopt->dualreds->nvars == 0);

      reopt->currentnode = SCIPnodeGetNumber(node);

      /* transform into the original space and then save the bound change */
      SCIP_CALL(SCIPvarGetOrigvarSum(&var, &scalar, &constant));
      newval = (newval - constant) / scalar;
      oldval = (oldval - constant) / scalar;

      assert(SCIPvarIsOriginal(var));

      if( SCIPsetIsEQ(set, oldval, newval) )
      {
         SCIPerrorMessage("cannot store equal bounds: old = %g, new = %g\n", oldval, newval);
         return SCIP_INVALIDDATA;
      }

      if( SCIPsetIsLT(set, newval, oldval) )
         boundtype = SCIP_BOUNDTYPE_UPPER;
      else
         boundtype = SCIP_BOUNDTYPE_LOWER;

      reopt->dualreds->vars[reopt->dualreds->nvars] = var;
      reopt->dualreds->vals[reopt->dualreds->nvars] = newval;
      reopt->dualreds->boundtypes[reopt->dualreds->nvars] = boundtype;
      ++reopt->dualreds->nvars;

      SCIPsetDebugMsg(set, ">> store %s bound change of <%s>: %g -> %g\n",
         (boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper"), SCIPvarGetName(var), oldval, newval);

      reopt->dualreds->linear = FALSE;
   }
   else
   {
      assert(reopt->currentnode == -1);
      assert(reopt->dualreds == NULL || reopt->dualreds->nvars == 0);

      reopt->currentnode = SCIPnodeGetNumber(node);
   }

   return SCIP_OKAY;
}

/** returns the number of bound changes based on dual information */
int SCIPreoptGetNDualBndchgs(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   )
{
   int ndualbndchgs = 0;

   assert(reopt != NULL);
   assert(node != NULL);

   if( SCIPnodeGetNumber(node) == reopt->currentnode )
   {
      assert(reopt->dualreds != NULL);
      ndualbndchgs = reopt->dualreds->nvars;
   }

   return ndualbndchgs;
}

/** returns the child nodes of @p node that need to be reoptimized next or NULL if @p node is a leaf */
SCIP_RETCODE SCIPreoptGetChildIDs(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int*         childs,             /**< array to store the child ids */
   int                   childssize,         /**< size of the childs array */
   int*                  nchilds             /**< pointer to store the number of child nodes */
   )
{
   SCIP_Bool runagain;
   unsigned int id;

   assert(reopt != NULL);
   assert(childssize > 0 && childs != NULL);
   assert(nchilds != NULL);

   (*nchilds) = 0;

   if( node == NULL )
      id = 0;
   else
      id = SCIPnodeGetReoptID(node);

   assert(id >= 1 || SCIPnodeGetDepth(node) == 0);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   /* check if there are redundant bound changes or infeasible nodes */
   runagain = TRUE;
   while( runagain && reopt->reopttree->reoptnodes[id]->nchilds > 0 )
   {
      SCIP_CALL( dryBranch(reopt, set, blkmem, &runagain, id) );
   }

   /* return the list of child nodes if some exists; otherwise return NULL */
   if( reopt->reopttree->reoptnodes[id]->childids != NULL && reopt->reopttree->reoptnodes[id]->nchilds > 0 )
   {
      int c;

      (*nchilds) = reopt->reopttree->reoptnodes[id]->nchilds;

      if( childssize < *nchilds )
         return SCIP_OKAY;

      for( c = 0; c < *nchilds; c++ )
         childs[c] = reopt->reopttree->reoptnodes[id]->childids[c];
   }

   return SCIP_OKAY;
}

/** returns all leaves of the subtree induced by @p node */
SCIP_RETCODE SCIPreoptGetLeaves(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int*         leaves,             /**< array to the the ids */
   int                   leavessize,         /**< size of leaves array */
   int*                  nleaves             /**< pointer to store the number of leave node */
   )
{
   unsigned int id;
   int i;

   assert(reopt != NULL);
   assert(leavessize > 0 && leaves != NULL);
   assert((*nleaves) >= 0);

   /* if the given node is we start from the root */
   if( node == NULL )
      id = 0;
   else
      id = SCIPnodeGetReoptID(node);

   /* return if the node is not part of the reoptimization tree */
   if( id == 0 && node != NULL )
   {
      (*nleaves) = 0;
      return SCIP_OKAY;
   }

   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   for( i = 0; i < leavessize; i++ )
      leaves[i] = 0;

   /* we traverse through all child nodes of the given node an collect all leave nodes of the subtrees induced by them */
   for( i = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++ )
   {
      unsigned int childid;

      assert(*nleaves + 1 <= leavessize);

      childid = reopt->reopttree->reoptnodes[id]->childids[i];
      assert(childid < reopt->reopttree->reoptnodessize);

      /* the node is already a leave */
      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
      {
         leaves[(*nleaves)] = reopt->reopttree->reoptnodes[id]->childids[i];
         ++(*nleaves);
      }
      /* go into the tree induced by the current child node */
      else
      {
         int nleaves2 = 0;

         SCIP_CALL( reoptGetLeaves(reopt, childid, &leaves[*nleaves], leavessize - (*nleaves), &nleaves2) );
         (*nleaves) += nleaves2;
      }
   }

   return SCIP_OKAY;
}

/** add all unprocessed nodes to the reoptimization tree */
SCIP_RETCODE SCIPreoptSaveOpenNodes(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE**           leaves,             /**< array of open leave nodes */
   int                   nleaves,            /**< number of open leave nodes */
   SCIP_NODE**           childs,             /**< array of open children nodes */
   int                   nchilds,            /**< number of open leave nodes */
   SCIP_NODE**           siblings,           /**< array of open sibling nodes */
   int                   nsiblings           /**< number of open leave nodes */
   )
{
   int n;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(nleaves >= 0);
   assert(nleaves == 0 || leaves != NULL);
   assert(nchilds >= 0);
   assert(nchilds == 0 || childs != NULL);
   assert(nsiblings >= 0);
   assert(nsiblings == 0 || siblings != NULL);

   SCIPsetDebugMsg(set, "save unprocessed nodes (%d leaves, %d children, %d siblings)\n", nleaves, nchilds, nsiblings);

   /* save open leaves */
   for( n = 0; n < nleaves; n++ )
   {
      SCIP_CALL( addNode(reopt, set, lp, blkmem, leaves[n], SCIP_REOPTTYPE_PRUNED, FALSE, FALSE,
            SCIPnodeGetLowerbound(leaves[n])) );
   }

   /* save open children */
   for( n = 0; n < nchilds; n++ )
   {
      SCIP_CALL( addNode(reopt, set, lp, blkmem, childs[n], SCIP_REOPTTYPE_PRUNED, FALSE, FALSE,
            SCIPnodeGetLowerbound(childs[n])) );
   }

   /* save open siblings */
   for( n = 0; n < nsiblings; n++ )
   {
      SCIP_CALL( addNode(reopt, set, lp, blkmem, siblings[n], SCIP_REOPTTYPE_PRUNED, FALSE, FALSE,
            SCIPnodeGetLowerbound(siblings[n])) );
   }

   return SCIP_OKAY;
}

/** merges the variable history of the current run with the stored history */
SCIP_RETCODE SCIPreoptMergeVarHistory(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR**            vars,               /**< original problem variables */
   int                   nvars               /**< number of original problem variables */
   )
{
   SCIP_VAR* transvar;
   SCIP_Real avginference[2];
   SCIP_Real avgcutoff[2];
   SCIP_Real bestsim;
   int bestrun;
   int idx;
   int d;
   int r;
   int v;

   assert(reopt != NULL);
   assert(stat != NULL);
   assert(nvars >= 0);

   if( !set->reopt_storevarhistory )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "start merging variable histories:\n");

   bestrun = reopt->run-2;
   bestsim = reopt->simtolastobj;

   /* find the run with the most similar objective */
   for( r = reopt->run-3; r >= 0 && reopt->objhaschanged && set->reopt_usepscost; r-- )
   {
      SCIP_Real sim;
      sim = reoptSimilarity(reopt, set, r, reopt->run-1, vars, nvars);

      if( sim == SCIP_INVALID )  /*lint !e777*/
         return SCIP_INVALIDRESULT;

      if( SCIPsetIsGT(set, sim, bestsim) )
      {
         bestsim = sim;
         bestrun = r;
      }
   }
   SCIPverbMessage(set->scip, SCIP_VERBLEVEL_NORMAL, NULL, "run %d has best similarity=%g\n", bestrun, bestsim);

   /* iterate through all variables and scale the histories */
   for( v = 0; v < nvars; v++ )
   {
      assert(SCIPvarIsOriginal(vars[v]));

      transvar = SCIPvarGetTransVar(vars[v]);
      assert(transvar != NULL);

      /* skip variable that are not active */
      if( !SCIPvarIsActive(transvar) )
         continue;

      idx = SCIPvarGetIndex(vars[v]);
      assert(0 <= idx && idx <= nvars);

      /* set the updated history for both directions */
      for( d = 0; d <= 1; d++ )
      {
         if( set->reopt_usepscost && !SCIPsetIsZero(set, reopt->varhistory[bestrun][idx]->pscostcount[d])
            && SCIPsetIsGT(set, bestsim, 0.985) ) /* 0.985 is a magic number determined in some experiments */
         {
            transvar->history->pscostcount[d] = 1.0;
            transvar->history->pscostweightedmean[d] = reopt->varhistory[bestrun][idx]->pscostweightedmean[d];
            transvar->history->pscostvariance[d] = 0.0;
            SCIPsetDebugMsg(set, "-> <%s> pscosts %4s: count=%g weightedmean=%g variance=%g\n", SCIPvarGetName(transvar),
               (d == 0 ? "down" : "up"), transvar->history->pscostcount[d], transvar->history->pscostweightedmean[d],
               transvar->history->pscostvariance[d]);
         }

         SCIPhistoryIncNBranchings(transvar->history, (SCIP_BRANCHDIR)d, 1);

         /* inference score */
         avginference[d] = SCIPhistoryGetAvgInferences(reopt->varhistory[reopt->run-2][idx], (SCIP_BRANCHDIR)d);
         SCIPhistoryIncInferenceSum(transvar->history, (SCIP_BRANCHDIR)d, avginference[d]);

         /* cutoff score */
         avgcutoff[d] = SCIPhistoryGetAvgCutoffs(reopt->varhistory[reopt->run-2][idx], (SCIP_BRANCHDIR)d);
         SCIPhistoryIncCutoffSum(transvar->history, (SCIP_BRANCHDIR)d, avgcutoff[d]);

         SCIPsetDebugMsg(set, "-> <%s> %4s scores: inf=%g cutoff=%g\n", SCIPvarGetName(transvar),
            (d == 0 ? "down" : "up"), avginference[d], avgcutoff[d]);
      }
   }

   return SCIP_OKAY;
}

/** updates the variable history */
SCIP_RETCODE SCIPreoptUpdateVarHistory(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            vars,               /**< original variable array */
   int                   nvars               /**< number of original variables */
   )
{
   int v;

   assert(reopt != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(nvars >= 0);

   if( !set->reopt_storevarhistory )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "updating variable history\n");

   if( reopt->varhistory[reopt->run-1] == NULL )
   {
      /* allocate memory */
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->varhistory[reopt->run-1], nvars) );

      for( v = 0; v < nvars; v++ )
      {
         SCIP_CALL( SCIPhistoryCreate(&(reopt->varhistory[reopt->run-1][v]), blkmem) );
      }
   }

   /* update the history and scale them */
   for( v = 0; v < nvars; v++ )
   {
      SCIP_VAR* transvar;
      int idx;

      assert(SCIPvarIsOriginal(vars[v]));
      idx = SCIPvarGetIndex(vars[v]);
      assert(idx >= 0 && idx < nvars);

      transvar = SCIPvarGetTransVar(vars[v]);
      assert(transvar != NULL);

      if( !SCIPvarIsActive(transvar) )
         continue;

      /* we store the complete history */
      SCIPhistoryReset(reopt->varhistory[reopt->run-1][idx]);
      SCIPhistoryUnite(reopt->varhistory[reopt->run-1][idx], transvar->history, FALSE);
   }

   return SCIP_OKAY;
}

/** reset the complete tree and set the given search frontier */
SCIP_RETCODE SCIPreoptApplyCompression(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPTNODE**      representatives,    /**< array of representatives */
   int                   nrepresentatives,   /**< number of representatives */
   SCIP_Bool*            success             /**< pointer to store if the method was successful */
   )
{
   SCIP_REOPTTREE* reopttree;
   unsigned int id;
   int r;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(representatives != NULL);
   assert(nrepresentatives > 0);

   reopttree = reopt->reopttree;

   /* reset the current search tree */
   SCIP_CALL( reoptResetTree(reopt, set, blkmem, FALSE) );
   assert(reopttree->nreoptnodes == 0);

   /* create a new root node */
   id = 0;
   SCIP_CALL( createReoptnode(reopttree, set, blkmem, id) );

   /* set the reopttype */
   reopttree->reoptnodes[0]->reopttype = (unsigned int)SCIP_REOPTTYPE_TRANSIT;

   /* add all representatives */
   for( r = 0; r < nrepresentatives; r++ )
   {
      /* get an empty slot*/
      id = (unsigned int) (size_t) SCIPqueueRemove(reopttree->openids);
      assert(1 <= id && id < reopttree->reoptnodessize);
      assert(reopttree->reoptnodes[id] == NULL);

      SCIP_CALL( createReoptnode(reopttree, set, blkmem, id) );
      assert(reopttree->reoptnodes[id] != NULL);

      /* set the new node
       * 1. copy all variables, bounds, and boundtypes
       * 2. copy all constraints
       * 3. set the parent relation
       */
      if( representatives[r]->nvars > 0 )
      {
         int v;

         assert(representatives[r]->nvars <= representatives[r]->varssize);

         for( v = 0; v < representatives[r]->nvars; v++ )
         {
            SCIP_CALL( SCIPreoptnodeAddBndchg(reopttree->reoptnodes[id], set, blkmem, representatives[r]->vars[v],
                  representatives[r]->varbounds[v], representatives[r]->varboundtypes[v]) );
         }
      }

      if( representatives[r]->nconss > 0 )
      {
         int c;

         assert(representatives[r]->nconss <= representatives[r]->consssize);

         for( c = 0; c < representatives[r]->nconss; c++ )
         {
            SCIP_CALL( SCIPreoptnodeAddCons(reopttree->reoptnodes[id], set, blkmem, representatives[r]->conss[c]->vars,
                  representatives[r]->conss[c]->vals, representatives[r]->conss[c]->boundtypes,
                  representatives[r]->conss[c]->lhs, representatives[r]->conss[c]->rhs,
                  representatives[r]->conss[c]->nvars, representatives[r]->conss[c]->constype,
                  representatives[r]->conss[c]->linear) );
         }
      }

      reopttree->reoptnodes[id]->parentID = representatives[r]->parentID; /*lint !e732*/

      assert(reopttree->reoptnodes[id]->parentID == 0);
      assert(reopttree->reoptnodes[id]->nvars >= 0);
      assert(reopttree->reoptnodes[id]->nvars <= reopttree->reoptnodes[id]->varssize);
      assert(reopttree->reoptnodes[id]->nconss >= 0);

      /* set the reopttype */
      if( reopttree->reoptnodes[id]->nconss == 0 )
         reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_LEAF;
      else
         reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_LOGICORNODE;

      /* add the representative as a child of the root */
      SCIP_CALL( reoptAddChild(reopttree, set, blkmem, 0, id) );
   }

   SCIPsetDebugMsg(set, "-> new tree consists of %d nodes, the root has %d child nodes.\n",
         reopttree->nreoptnodes, reopttree->reoptnodes[0]->nchilds);

   (*success) = TRUE;

   return SCIP_OKAY;
}

/** transforms a set of dual reductions into a linear constraint */
static
SCIP_RETCODE transformDualredsToLinear(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPTCONSDATA*   consdata,           /**< reoptimization constraint data that should represent to set of solutions
                                               *  pruned by the dual reductions
                                               */
   SCIP_REOPTCONSDATA*   dualreds            /**< set of dual reductions */
   )
{
   int v;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(consdata != NULL);
   assert(dualreds != NULL);

   /* we have to transform the set of bound changes into a linear constraint */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &consdata->vars, dualreds->vars, dualreds->nvars) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &consdata->vals, dualreds->nvars) );
   consdata->boundtypes = NULL;

   consdata->varssize = dualreds->nvars;
   consdata->nvars = dualreds->nvars;
   consdata->constype = REOPT_CONSTYPE_DUALREDS;
   consdata->linear = TRUE;

   /* set lhs and rhs */
   consdata->lhs = 1.0;
   consdata->rhs = SCIPsetInfinity(set);

   for( v = 0; v < consdata->nvars; v++ )
   {
      assert(consdata->vars[v] != NULL);

      /* the bound is 0.0, the variable has to appear with a coefficient +1.0 in the constraint, sides do not change */
      if( SCIPsetIsEQ(set, dualreds->vals[v], 0.0) )
      {
         assert(dualreds->boundtypes[v] == SCIP_BOUNDTYPE_UPPER);
         consdata->vals[v] = 1.0;
      }
      /* the bound is 1.0, the variable has to appear with a coefficient -1.0 in the constraint, we subtract -1.0 from lhs
       *   logicor:       sum x_i + ~y_i    >= 1
       *           <==>   sum x_i + (1-y_i) >= 1
       *           <==>   sum x_i - y_i     >= 0
       */
      else
      {
         assert(SCIPsetIsEQ(set, dualreds->vals[v], 1.0));
         assert(dualreds->boundtypes[v] == SCIP_BOUNDTYPE_LOWER);

         consdata->vals[v] = -1.0;
         consdata->lhs -= 1.0;
      }
   }

   return SCIP_OKAY;
}


/** transforms a set of dual reductions into a bounddisjuction constraint */
static
SCIP_RETCODE transformDualredsToBounddisjunction(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPTCONSDATA*   consdata,           /**< reoptimization constraint data that should represent to set of solutions
                                               *  pruned by the dual reductions
                                               */
   SCIP_REOPTCONSDATA*   dualreds            /**< set of dual reductions */
   )
{
   int v;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(consdata != NULL);
   assert(dualreds != NULL);

   /* we have to transform the set of bound changes into a linear constraint */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &consdata->vars, dualreds->vars, dualreds->nvars) );
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &consdata->vals, dualreds->vals, dualreds->nvars) );
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &consdata->boundtypes, dualreds->boundtypes, dualreds->nvars) );

   consdata->varssize = dualreds->nvars;
   consdata->nvars = dualreds->nvars;
   consdata->constype = REOPT_CONSTYPE_DUALREDS;
   consdata->linear = FALSE;

   /* set lhs and rhs */
   consdata->lhs = SCIP_UNKNOWN;
   consdata->rhs = SCIP_UNKNOWN;

   for( v = 0; v < consdata->nvars; v++ )
   {
      SCIP_Real glbbd;

      assert(consdata->vars[v] != NULL);

      /* we do the followung to transformations:
       * (a) x <= val   ==>   (x >= val+1)
       * (b) x >= val   ==>   (x <= val-1)
       */
      if( consdata->boundtypes[v] == SCIP_BOUNDTYPE_UPPER )
      {
         glbbd = SCIPvarGetUbGlobal(consdata->vars[v]);
         consdata->vals[v] = MIN(consdata->vals[v]+1.0, glbbd);
      }
      else
      {
         assert(dualreds->boundtypes[v] == SCIP_BOUNDTYPE_LOWER);
         glbbd = SCIPvarGetLbGlobal(consdata->vars[v]);
         consdata->vals[v] = MAX(glbbd, consdata->vals[v]-1.0);
      }
      consdata->boundtypes[v] = (SCIP_BOUNDTYPE)(SCIP_BOUNDTYPE_UPPER - consdata->boundtypes[v]); /*lint !e656*/
   }

   return SCIP_OKAY;
}

/** splits the root into several nodes and moves the child nodes of the root to one of the created nodes */
SCIP_RETCODE SCIPreoptSplitRoot(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int*                  ncreatedchilds,     /**< pointer to store the number of created nodes */
   int*                  naddedconss         /**< pointer to store the number added constraints */
   )
{
   SCIP_REOPTTREE* reopttree;
   SCIP_REOPTNODE** reoptnodes;
   SCIP_REOPTCONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* bounds;
   SCIP_BOUNDTYPE* boundtypes;
   int* perm = NULL;
   unsigned int id;
   int nbndchgs;
   int nchilds;
   int nvars = 0;
   int v;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);

   reopttree = reopt->reopttree;
   assert(reopttree != NULL);

   reoptnodes = reopttree->reoptnodes;
   assert(reoptnodes != NULL);
   assert(reoptnodes[0] != NULL);
   assert(reoptnodes[0]->dualreds);
   assert(reoptnodes[0]->reopttype == (unsigned int)SCIP_REOPTTYPE_STRBRANCHED);

   nchilds = reoptnodes[0]->nchilds;

   assert(reoptnodes[0]->dualredscur != NULL);
   nbndchgs = reoptnodes[0]->dualredscur->nvars;

   (*ncreatedchilds) = 0;
   (*naddedconss) = 0;

   /* create a node with all variables fixed, i.e., reconstruct the root of the last iteration */

   /* ensure that two free slots are available  */
   SCIP_CALL( reopttreeCheckMemory(reopttree, set, blkmem) );
   id = (unsigned int) (size_t) SCIPqueueRemove(reopttree->openids);

   assert(0 < id && id < reopt->reopttree->reoptnodessize);
   assert(reoptnodes[id] == NULL || reoptnodes[id]->nvars == 0);

   /*   1. create the node
    *   2. add all bound changes
    *   3. move all child nodes to id
    *   4. add id as a child of the root node
    */
   SCIP_CALL( createReoptnode(reopttree, set, blkmem, id) );
   reoptnodes[id]->parentID = 0;
   reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_TRANSIT;

   /* check memory */
   SCIP_CALL( reoptnodeCheckMemory(reoptnodes[id], set, blkmem, nbndchgs, nchilds, 0) );
   assert(reoptnodes[id]->varssize >= nbndchgs);
   assert(reoptnodes[id]->nvars == 0);
   assert(reoptnodes[id]->vars != NULL);
   assert(reoptnodes[id]->varbounds != NULL);
   assert(reoptnodes[id]->varboundtypes != NULL);

   /* create a permutation array */
   if( !set->reopt_usesplitcons )
   {
      assert(perm == NULL);
      SCIP_CALL( SCIPsetAllocBufferArray(set, &perm, nbndchgs) );
   }

   /* copy bounds */
   for( v = 0; v < nbndchgs; v++ )
   {
      reoptnodes[id]->vars[v] = reoptnodes[0]->dualredscur->vars[v];
      reoptnodes[id]->varbounds[v] = reoptnodes[0]->dualredscur->vals[v];
      reoptnodes[id]->varboundtypes[v] = reoptnodes[0]->dualredscur->boundtypes[v];
      ++reoptnodes[id]->nvars;

      /* fill a permutation array */
      if( !set->reopt_usesplitcons )
         perm[v] = v;   /*lint !e613*/
   }
   assert(reoptnodes[id]->nvars == reoptnodes[0]->dualredscur->nvars);

   /* move the children */
   SCIP_CALL( reoptMoveIDs(reopttree, set, blkmem, 0, id) );
   assert(reoptnodes[0]->nchilds == 0);

   /* add the new reoptimization node as a child of the root node */
   SCIP_CALL( reoptAddChild(reopttree, set, blkmem, 0, id) );

   ++(*ncreatedchilds);

   if( set->reopt_usesplitcons )
   {
      int nbinvars = 0;
      int nintvars = 0;
      int ncontvars = 0;

      assert(*ncreatedchilds == 1);

      /* ensure that there is a free slots */
      SCIP_CALL( reopttreeCheckMemory(reopttree, set, blkmem) );
      id = (unsigned int) (size_t) SCIPqueueRemove(reopttree->openids);
      assert(0 < id && id < reopt->reopttree->reoptnodessize);

      /* 1. create the node
       * 2. add the constraint to ensure that at least one
       *    variable gets different
       * 3. add id as a child of the root node
       */
      SCIP_CALL( createReoptnode(reopttree, set, blkmem, id) );
      reoptnodes[id]->parentID = 0;
      reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_LOGICORNODE;

      /* check memory for added constraints */
      SCIP_CALL( reoptnodeCheckMemory(reoptnodes[id], set, blkmem, 0, 0, 1) );

      /* create the constraint */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reoptnodes[id]->conss[0]) );
      consdata = reoptnodes[id]->conss[0];

      /* count number of binary, integer, and continuous varibales */
      for( v = 0; v < nbndchgs; v++ )
      {
         switch( SCIPvarGetType(reoptnodes[0]->dualredscur->vars[v]) ) {
         case SCIP_VARTYPE_BINARY:
            ++nbinvars;
            break;
         case SCIP_VARTYPE_INTEGER:
         case SCIP_VARTYPE_IMPLINT:
            ++nintvars;
            break;
         case SCIP_VARTYPE_CONTINUOUS:
            ++ncontvars;
            break;
         default:
            SCIPerrorMessage("Cannot handle vartype %d\n", SCIPvarGetType(reoptnodes[0]->dualredscur->vars[v]));
            return SCIP_INVALIDDATA;
         }
      }

      /* we create a linear constraint, since all variables are binary */
      if( nbinvars == nbndchgs )
      {
         SCIP_CALL( transformDualredsToLinear(reopt, set, blkmem, consdata, reoptnodes[0]->dualredscur) );
      }
      /* we create a bounddisjunction constraint, since at least one variable is (implicit) integer or continuous */
      else
      {
         assert(nintvars > 0 || ncontvars > 0);
         SCIP_CALL( transformDualredsToBounddisjunction(reopt, set, blkmem, consdata, reoptnodes[0]->dualredscur) );
      }
      ++reoptnodes[id]->nconss;

      /* add id as a child of the root node */
      SCIP_CALL( reoptAddChild(reopttree, set, blkmem, 0, id) );
      ++(*ncreatedchilds);

      ++(*naddedconss);
   }
   else
   {
      int c;

      assert(*ncreatedchilds == 1);
      assert(perm != NULL);

      vars = reoptnodes[0]->dualredscur->vars;
      bounds = reoptnodes[0]->dualredscur->vals;
      boundtypes = reoptnodes[0]->dualredscur->boundtypes;
      nvars = reoptnodes[0]->dualredscur->nvars;
      assert(perm[0] == 0 && perm[nvars-1] == nvars-1);

      /* calculate the order of the variables */
      switch (set->reopt_varorderinterdiction)
      {
         /* default order */
         case 'd':
            break;

         /* inference order */
         case 'i':
            SCIP_CALL( getInferenceOrder(set, stat, perm, vars, bounds, boundtypes, nvars) );
            break;

         /* random order */
         case 'r':
            SCIPrandomPermuteIntArray(reopt->randnumgen, perm, 0, nvars-1);
            break;

         default:
            return SCIP_INVALIDDATA;
      }

      /* create nvars nodes in the fashion of interdiction branching */
      for( c = 0; c < nvars; c++ )
      {
         /* ensure that two free slots are available  */
         SCIP_CALL( reopttreeCheckMemory(reopttree, set, blkmem) );
         id = (unsigned int) (size_t) SCIPqueueRemove(reopttree->openids);

         assert(0 < id && id < reopt->reopttree->reoptnodessize);
         assert(reoptnodes[id] == NULL || reoptnodes[id]->nvars == 0);

         /*   1. create the node
          *   2. fix the first v bound changes to vals[v] and v+1 to vals[v] +/- 1 (depending on the bound- and vartype)
          *   4. add the ID id as a child of the root node
          */
         SCIP_CALL( createReoptnode(reopttree, set, blkmem, id) );
         reoptnodes[id]->parentID = 0;
         reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_TRANSIT;

         /* check memory */
         SCIP_CALL( reoptnodeCheckMemory(reoptnodes[id], set, blkmem, c+1, 0, 0) );
         assert(reoptnodes[id]->varssize >= perm[c]+1);
         assert(reoptnodes[id]->nvars == 0);
         assert(reoptnodes[id]->vars != NULL);
         assert(reoptnodes[id]->varbounds != NULL);
         assert(reoptnodes[id]->varboundtypes != NULL);

         /* the permutation is the identity */
         if( set->reopt_varorderinterdiction == 'd' )
         {
            /* copy first c bound changes */
            for( v = 0; v < c; v++ )
            {
               reoptnodes[id]->vars[v] = vars[v];
               reoptnodes[id]->varbounds[v] = bounds[v];
               reoptnodes[id]->varboundtypes[v] = boundtypes[v];
            }
         }
         else
         {
            /* copy first c bound changes */
            for( v = 0; v < c; v++ )
            {
               reoptnodes[id]->vars[v] = vars[perm[v]];
               reoptnodes[id]->varbounds[v] = bounds[perm[v]];
               reoptnodes[id]->varboundtypes[v] = boundtypes[perm[v]];
            }
         }
         reoptnodes[id]->nvars += c;

         /* set bound change v+1 (= c) to vals[v] +/- 1 (depending on the bound- and vartype) */
         assert(v == c);
         reoptnodes[id]->vars[c] = vars[perm[c]];
         reoptnodes[id]->varbounds[c] = bounds[perm[c]];
         if( SCIPvarGetType(vars[perm[c]]) != SCIP_VARTYPE_CONTINUOUS )
         {
            if( boundtypes[perm[c]] == SCIP_BOUNDTYPE_LOWER )
               reoptnodes[id]->varbounds[c] -= 1.0;
            else
               reoptnodes[id]->varbounds[c] += 1.0;
         }
         reoptnodes[id]->varboundtypes[c] = (boundtypes[perm[c]] == SCIP_BOUNDTYPE_UPPER ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER);
         ++reoptnodes[id]->nvars;

         /* add dummy1 as a child of the root node */
         SCIP_CALL( reoptAddChild(reopttree, set, blkmem, 0, id) );

         ++(*ncreatedchilds);
      }

      assert(*ncreatedchilds == nvars+1);

      SCIPsetFreeBufferArray(set, &perm);
      perm = NULL;
   }
   assert(perm == NULL);

   /* free the current dualredscur and assign dualredsnex */
   assert(reoptnodes[0]->dualredscur->vars != NULL);
   assert(reoptnodes[0]->dualredscur->vals != NULL);
   assert(reoptnodes[0]->dualredscur->boundtypes != NULL);

   /* free the current dualredscur and assign dualredsnex */
   SCIP_CALL( reoptnodeUpdateDualConss(reoptnodes[0], blkmem) );

   /* change the reopttype of the root node */
   SCIPnodeSetReopttype(SCIPtreeGetRootNode(tree), SCIP_REOPTTYPE_TRANSIT);

   return SCIP_OKAY;
}

/** reset the stored information abound bound changes based on dual information */
SCIP_RETCODE SCIPreoptResetDualBndchgs(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node,               /**< node of the search tree */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   unsigned int id;

   assert(reopt != NULL);
   assert(node != NULL);

   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* return if the node is not part of the reoptimization tree */
   if( SCIPnodeGetDepth(node) > 0 && id == 0 )
      return SCIP_OKAY;

   /* reset the dual constraint */
   SCIP_CALL( reoptnodeResetDualConss(reopt->reopttree->reoptnodes[id], blkmem) );

   return SCIP_OKAY;
}

/** return the branching path stored of the given node in the reoptimization tree */
void SCIPreoptnodeGetPath(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR**            vars,               /**< array for variables */
   SCIP_Real*            vals,               /**< array for values */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array for bound types */
   int                   varssize,           /**< size of arrays vars, vals, and boundtypes */
   int*                  nbndchgs,           /**< pointer to store the number of bound changes */
   int*                  nbndchgsafterdual   /**< pointer to store the number of bound changes applied after
                                              *  the first dual reduction at the given node */
   )
{
   int v;
   int nvars2;
   int nafterdualvars2;

   assert(reopt != NULL);
   assert(reoptnode != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(boundtypes != NULL);

   (*nbndchgs) = reoptnode->nvars;
   (*nbndchgsafterdual) = reoptnode->nafterdualvars;

   /* return if the size of the given array is not large enough */
   if( varssize == 0 || varssize < *nbndchgs + *nbndchgsafterdual )
      return;

   /* add all bound changes made by branching (including dual reductions) */
   for( v = 0; v < *nbndchgs; v++ )
   {
      vars[v] = reoptnode->vars[v];
      vals[v] = reoptnode->varbounds[v];
      boundtypes[v] = reoptnode->varboundtypes[v];
   }

   /* add all bound changes made applied after a dual reduction */
   for( ; v < *nbndchgs + *nbndchgsafterdual; v++ )
   {
      vars[v] = reoptnode->afterdualvars[v-(*nbndchgs)];
      vals[v] = reoptnode->afterdualvarbounds[v-(*nbndchgs)];
      boundtypes[v] = reoptnode->afterdualvarboundtypes[v-(*nbndchgs)];
   }

   /* go along the root path within the reoptimization tree */
   if( reoptnode->parentID != 0 )
   {
      SCIP_REOPTNODE* parent;

      parent = reopt->reopttree->reoptnodes[reoptnode->parentID];
      SCIPreoptnodeGetPath(reopt, parent, &vars[v], &vals[v], &boundtypes[v], varssize, &nvars2, &nafterdualvars2);

      (*nbndchgs) += nvars2;
      (*nbndchgsafterdual) += nafterdualvars2;
   }
}

/** delete a node stored in the reoptimization tree */
SCIP_RETCODE SCIPreoptDeleteNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   unsigned int          id,                 /**< id of a stored node */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);
   assert(blkmem != NULL);

   SCIP_CALL( reopttreeDeleteNode(reopt->reopttree, set, blkmem, id, TRUE) );
   SCIP_CALL( SCIPqueueInsert(reopt->reopttree->openids, (void*) (size_t) id) );

   return SCIP_OKAY;
}

/** reactivate the given @p reoptnode and split them into several nodes if necessary */
SCIP_RETCODE SCIPreoptApply(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branching tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree to reactivate */
   unsigned int          id,                 /**< id of the node to reactivate */
   SCIP_Real             estimate,           /**< estimate of the child nodes that should be created */
   SCIP_NODE**           childnodes,         /**< array to store the created child nodes */
   int*                  ncreatedchilds,     /**< pointer to store number of created child nodes */
   int*                  naddedconss,        /**< pointer to store number of generated constraints */
   int                   childnodessize,     /**< available size of childnodes array */
   SCIP_Bool*            success             /**< pointer store the result */
   )
{
   assert(reopt != NULL);
   assert(scip != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(cliquetable != NULL);
   assert(blkmem != NULL);
   assert(reoptnode != NULL);
   assert(childnodes != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(success != NULL);

   SCIPsetDebugMsg(set, "reactivating node at id %u:\n", id);

   *success = FALSE;

   /* check if we need to split the node */
   if( reoptnode->reopttype == (unsigned int)SCIP_REOPTTYPE_STRBRANCHED
      || reoptnode->reopttype == (unsigned int)SCIP_REOPTTYPE_INFSUBTREE )
   {
      int c;

      assert(reoptnode->dualreds);

      /* we want use a constraint to split the node into two disjoint node */
      if( set->reopt_usesplitcons )
      {
         if( reoptnode->reopttype == (unsigned int)SCIP_REOPTTYPE_INFSUBTREE )
         {
            assert(reoptnode->dualredscur != NULL);
            assert(reoptnode->dualredscur->constype == REOPT_CONSTYPE_INFSUBTREE);
            (*ncreatedchilds) = 1;
         }
         else
         {
            assert(reoptnode->dualredscur != NULL);
            assert(reoptnode->dualredscur->constype == REOPT_CONSTYPE_DUALREDS);
            (*ncreatedchilds) = 2;
         }

         /* in both cases we add exactly one constraint */
         (*naddedconss) = 1;

         if( childnodessize < *ncreatedchilds )
            return SCIP_OKAY;

         /* generate the nodes */
         for( c = 0; c < *ncreatedchilds; c++ )
         {
            /* create the child node */
            SCIP_CALL( SCIPnodeCreateChild(&childnodes[c], blkmem, set, stat, tree, 1.0, estimate) );

            /* change all bounds; convert the bound changes after the first based on dual reductions into branching
             * for second node only. if we generate only one node, i.e., the pruned part, we do not need this
             * changes anyway.
             */
            SCIP_CALL( changeAncestorBranchings(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue,
                  cliquetable, blkmem, childnodes[c], id, c == 1) );

            /* add all local constraints */
            SCIP_CALL( addLocalConss(scip, reopt, set, stat, blkmem, childnodes[c], id) );

            /* we can use the old lowerbound if the objective function has not changed */
            if( !reopt->objhaschanged && SCIPsetIsGT(set, reopt->reopttree->reoptnodes[id]->lowerbound, estimate) )
               SCIPnodeSetEstimate(childnodes[c], set, reopt->reopttree->reoptnodes[id]->lowerbound);

            if( c == 0 )
            {
               /* in both cases the node generated first represents the pruned is currently not part of the reoptimization tree */
               SCIPnodeSetReopttype(childnodes[c], SCIP_REOPTTYPE_NONE);

               /* add the constraint to the node */
               assert(reopt->reopttree->reoptnodes[id]->dualredscur != NULL);
               SCIP_CALL( addSplitcons(reopt, scip, set, stat, blkmem, transprob, origprob, tree, lp, branchcand,
                     eventqueue, cliquetable, childnodes[c], id) );

               /* fixBounds() does the same, but in this case we go not into it */
               if( reoptnode->dualredscur->constype == REOPT_CONSTYPE_INFSUBTREE )
               {
                  assert(reoptnode->dualredscur->nvars > 0);
                  assert(reoptnode->dualredscur->varssize > 0);

                  /* delete dualredscur and move dualredsnex -> dualredscur */
                  SCIP_CALL( reoptnodeUpdateDualConss(reoptnode, blkmem) );
               }

               /* the added constraint could be deleted due to propagation, thus, we store the node in the reoptimization
                * tree. the node has to stored anyway, because of the constraint representing the dual reductions
                */
               SCIP_CALL( addNode(reopt, set, lp, blkmem, childnodes[c], SCIP_REOPTTYPE_LOGICORNODE, FALSE, FALSE,
                     -SCIPsetInfinity(set)) );
            }
            else
            {
               /* if we reach this lines of code, the current node represents the original node including all bound
                * changes based in dual information.
                */
               assert(reoptnode->dualredscur->constype == REOPT_CONSTYPE_DUALREDS);
               if( reoptnode->nconss == 0 )
                  SCIPnodeSetReopttype(childnodes[c], SCIP_REOPTTYPE_TRANSIT);
               else
                  SCIPnodeSetReopttype(childnodes[c], SCIP_REOPTTYPE_LOGICORNODE);

               /* fix all bound changes based on dual information and convert them into branchings */
               assert(reopt->reopttree->reoptnodes[id]->dualredscur != NULL);
               SCIP_CALL( fixBounds(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue, cliquetable,
                     blkmem, childnodes[c], id, TRUE) );

               /* set the unique id the id of the original node */
               SCIPnodeSetReoptID(childnodes[c], id);
            }
         }

         /* reset the stored dual constraints */
         SCIP_CALL( reoptnodeUpdateDualConss(reopt->reopttree->reoptnodes[id], blkmem) );

         /* set the reoptimization type */
         if( reopt->reopttree->reoptnodes[id]->dualreds )
            reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_STRBRANCHED;
         else
            reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_TRANSIT;

         *success = TRUE;
      }
      else
      {
         SCIP_VAR** vars;
         SCIP_Real* bounds;
         SCIP_BOUNDTYPE* boundtypes;
         int* perm = NULL;
         int nvars;

         vars = reoptnode->dualredscur->vars;
         bounds = reoptnode->dualredscur->vals;
         boundtypes = reoptnode->dualredscur->boundtypes;
         nvars = reoptnode->dualredscur->nvars;

         *ncreatedchilds = nvars+1;
         *naddedconss = 0;

         /* check if there is enough memory allocated */
         if( childnodessize < *ncreatedchilds )
            return SCIP_OKAY;

         /* create and fill permutation array */
         SCIP_CALL( SCIPsetAllocBufferArray(set, &perm, nvars) );
         for( c = 0; c < nvars; c++ )
            perm[c] = c;

         /* calculate the order of the variables */
         switch (set->reopt_varorderinterdiction)
         {
            /* default order */
            case 'd':
               break;

            /* inference order */
            case 'i':
               SCIP_CALL( getInferenceOrder(set, stat, perm, vars, bounds, boundtypes, nvars) );
               break;

            /* random order */
            case 'r':
               SCIPrandomPermuteIntArray(reopt->randnumgen, perm, 0, nvars-1);
               break;

            default:
               return SCIP_INVALIDDATA;
         }

         assert(reopt->reopttree->reoptnodes[id] != NULL);
         reoptnode = reopt->reopttree->reoptnodes[id];

         /* enough that the node need to split */
         assert(reoptnode->dualreds);

         /* iterate over all nodes and change the necessary bounds (nodes[0] corresponds to the original one)
          * we need to do this in the reverse order because we want to transform the bound changes based on dual information
          * into branching decisions at nodes[0].
          */
         for( c = nvars; c >= 0; c-- )
         {
            /* create the child node */
            SCIP_CALL( SCIPnodeCreateChild(&childnodes[c], blkmem, set, stat, tree, 1.0, estimate) );

#ifdef SCIP_MORE_DEBUG
            SCIPsetDebugMsg(set, " change bounds at node %lld\n", SCIPnodeGetNumber(childnodes[c]));
#endif

            /* change all bounds */
            SCIP_CALL( changeAncestorBranchings(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue,
                  cliquetable, blkmem, childnodes[c], id, FALSE) );

            /* reconstruct the original node and the pruned part, respectively */
            if( c == 0 )
            {
               /* fix bound changes based on dual information and convert all these bound changes to normal bound changes */
               SCIP_CALL( fixBounds(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue, cliquetable,
                     blkmem, childnodes[c], id, TRUE) );

               /* set the reopttype of the node */
               SCIPnodeSetReopttype(childnodes[c], SCIP_REOPTTYPE_TRANSIT);

               /* set the unique id */
               SCIPnodeSetReoptID(childnodes[c], id);
            }
            else
            {
               /* fix the first c bound changes and negate the (c+1)th */
               SCIP_CALL( fixInterdiction(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue, cliquetable,
                     blkmem, childnodes[c], id, perm, vars, bounds, boundtypes, nvars, c) );
            }

            /* add all local constraints */
            SCIP_CALL( addLocalConss(scip, reopt, set, stat, blkmem, childnodes[c], id) );

            /* we can use the old lowerbound if the objective function has not changed */
            if( !reopt->objhaschanged && SCIPsetIsGT(set, reopt->reopttree->reoptnodes[id]->lowerbound, estimate) )
               SCIPnodeSetEstimate(childnodes[c], set, reopt->reopttree->reoptnodes[id]->lowerbound);
         }

         /* free buffer array */
         SCIPsetFreeBufferArray(set, &perm);

         /* reset the stored dual constraints */
         SCIP_CALL( reoptnodeUpdateDualConss(reopt->reopttree->reoptnodes[id], blkmem) );

         /* set the reoptimization type to transit */
         if( reopt->reopttree->reoptnodes[id]->dualreds )
            reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_STRBRANCHED;
         else
            reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_TRANSIT;

         *success = TRUE;
      }
   }
   else
   {
      /* we need the create exactly one node to reconstruct the node itself and no additional constraint */
      (*ncreatedchilds) = 1;
      (*naddedconss) = 0;

      if( childnodessize < *ncreatedchilds )
         return SCIP_OKAY;

      /* create the child node */
      SCIP_CALL( SCIPnodeCreateChild(&childnodes[0], blkmem, set, stat, tree, 1.0, estimate) );

      /* change all bounds */
      assert(reoptnode->nafterdualvars == 0);
      SCIP_CALL( changeAncestorBranchings(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue,
            cliquetable, blkmem, childnodes[0], id, FALSE) );

      /* add all local constraints */
      SCIP_CALL( addLocalConss(scip, reopt, set, stat, blkmem, childnodes[0], id) );

      /* we can use the old lowerbound if the objective function has not changed */
      if( !reopt->objhaschanged && SCIPsetIsGT(set, reopt->reopttree->reoptnodes[id]->lowerbound, estimate) )
         SCIPnodeSetEstimate(childnodes[0], set, reopt->reopttree->reoptnodes[id]->lowerbound);

      /* set the reopttype */
      assert(reoptnode->reopttype != (unsigned int)SCIP_REOPTTYPE_INFSUBTREE
          && reoptnode->reopttype != (unsigned int)SCIP_REOPTTYPE_STRBRANCHED);
      SCIPnodeSetReopttype(childnodes[0], (SCIP_REOPTTYPE)reoptnode->reopttype);

      /* set the unique id */
      SCIPnodeSetReoptID(childnodes[0], id);

      *success = TRUE;
   }

   return SCIP_OKAY;
}

/** returns the time needed to store the nodes for reoptimization */
SCIP_Real SCIPreoptGetSavingtime(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return SCIPclockGetTime(reopt->savingtime);
}

/** add the stored constraints globally to the problem */
SCIP_RETCODE SCIPreoptApplyGlbConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   char name[SCIP_MAXSTRLEN];
   int c;

   assert(scip != NULL);
   assert(reopt != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);

   if( reopt->glbconss == NULL || reopt->nglbconss == 0 )
      return SCIP_OKAY;

   for( c = reopt->nglbconss-1; c >= 0; c-- )
   {
      SCIP_CONS* cons;
      SCIP_VAR** consvars;
      int nbinvars;
      int nintvars;
      int v;

      assert(reopt->glbconss[c] != NULL);
      assert(reopt->glbconss[c]->nvars > 0);

      cons = NULL;
      consvars = NULL;
      nbinvars = 0;
      nintvars = 0;

      /* check if we can use a logic-or or if we have to use a bounddisjuction constraint */
      for( v = 0; v < reopt->glbconss[c]->nvars; v++ )
      {
         if( SCIPvarGetType(reopt->glbconss[c]->vars[v]) == SCIP_VARTYPE_BINARY )
            ++nbinvars;
         else if( SCIPvarGetType(reopt->glbconss[c]->vars[v]) == SCIP_VARTYPE_INTEGER
               || SCIPvarGetType(reopt->glbconss[c]->vars[v]) == SCIP_VARTYPE_IMPLINT )
            ++nintvars;
         else
         {
            SCIPerrorMessage("Expected variable type binary or (impl.) integer for variable <%s> in global constraint at pos. %d.\n",
                  SCIPvarGetName(reopt->glbconss[c]->vars[v]), c);
            return SCIP_INVALIDDATA;
         }
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "glb_%s_%d_%d", reopt->glbconss[c]->constype == REOPT_CONSTYPE_CUT ? "cut" : "inf", reopt->run, c);

      /* @todo use active representatives */

      /* all variables are binary, we can create a logic-or constraint */
      if( nbinvars == reopt->glbconss[c]->nvars )
      {
         SCIPsetDebugMsg(set, "-> add logic-or constraints with %d binvars\n", nbinvars);

         /* allocate buffer */
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, reopt->glbconss[c]->nvars) );

         for( v = 0; v < reopt->glbconss[c]->nvars; v++ )
         {
            consvars[v] = reopt->glbconss[c]->vars[v];
            assert(SCIPvarIsOriginal(consvars[v]));

            /* negate the variable if it was fixed to 1 */
            if( SCIPsetIsFeasEQ(set, reopt->glbconss[c]->vals[v], 0.0) )
            {
               assert(reopt->glbconss[c]->boundtypes[v] == SCIP_BOUNDTYPE_UPPER);
               SCIP_CALL( SCIPvarNegate(consvars[v], blkmem, set, stat, &consvars[v]) );
            }
         }

         /* create the logic-or constraint */
         SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, reopt->glbconss[c]->nvars,
               consvars, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         /* free buffer */
         SCIPfreeBufferArray(scip, &consvars);
      }
      /* not all variables are binary, we need a bounddisjunction constraint */
      else
      {
         assert(reopt->glbconss[c]->nvars == nbinvars + 2*nintvars);

         SCIPsetDebugMsg(set, "-> add bounddisjuction constraints with %d binvars, %d intvars\n", nbinvars, (int) (2*nintvars));

         /* create the bounddisjuction constraint */
         SCIP_CALL( SCIPcreateConsBasicBounddisjunction(scip, &cons, name, reopt->glbconss[c]->nvars, reopt->glbconss[c]->vars,
               reopt->glbconss[c]->boundtypes, reopt->glbconss[c]->vals) );
      }

#ifdef SCIP_DEBUG_CONSS
      SCIPdebugPrintCons(scip, cons, NULL);
#endif

      SCIP_CALL( SCIPaddCons(scip, cons) );

      /* remember the constraint for re-activation */
      assert(!SCIPhashmapExists(reopt->activeconss, (void*)cons));
      SCIP_CALL( SCIPhashmapInsert(reopt->activeconss, (void*)cons, (void*)cons) );

      /* don't release the constraint because we would need to capture the constraint anyway */

      /* mark the constraint as empty */
      reopt->glbconss[c]->nvars = 0;
   }

   SCIPsetDebugMsg(set, "added %d gobal constraints\n", reopt->nglbconss);

   /* reset number of global constraints */
   reopt->nglbconss = 0;

   return SCIP_OKAY;
}

/** add the stored cuts to the separation storage */
SCIP_RETCODE SCIPreoptApplyCuts(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node,               /**< current focus node */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cutpool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_Bool             root                /**< bool whether the current node is the root */
   )
{
   SCIP_REOPTNODE* reoptnode;
   SCIP_Bool infeasible;
   unsigned int id;
   int ncuts;
   int c;

   assert(reopt != NULL);
   assert(node != NULL);
   assert(sepastore != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(eventqueue != NULL);
   assert(eventfilter != NULL);
   assert(lp != NULL);

   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* skip nodes that are node part of the reoptimization tree */
   if( id == 0 && SCIPnodeGetDepth(node) > 0 )
      return SCIP_OKAY;

   reoptnode = reopt->reopttree->reoptnodes[id];
   assert(reoptnode != NULL);

   ncuts = 0;
   for( c = reoptnode->nconss-1; c >= 0; c-- )
   {
      SCIP_REOPTCONSDATA* cons;

      cons = reoptnode->conss[c];
      assert(cons != NULL);

      if( cons->constype == REOPT_CONSTYPE_CUT )
      {
         SCIP_ROW* cut;
         SCIP_COL** cols;
         SCIP_Real* vals;
         char cutname[SCIP_MAXSTRLEN];
         int ncols;
         int v;

         SCIP_CALL( SCIPsetAllocBufferArray(set, &cols, cons->nvars) );
         SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, cons->nvars) );

         ncols = 0;
         for( v = 0; v < cons->nvars; v++ )
         {
            SCIP_VAR* transvar;

            assert(SCIPvarIsOriginal(cons->vars[v]));

            transvar = SCIPvarGetTransVar(cons->vars[v]);
            assert(transvar != NULL);
            assert(SCIPvarGetStatus(transvar) == SCIP_VARSTATUS_COLUMN);

            vals[ncols] = cons->vals[v];
            cols[ncols] = SCIPvarGetCol(transvar);
            assert(cols[ncols] != NULL);

            ++ncols;
         }
         assert(ncols == cons->nvars);

         (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "reoptcut_%d_%d", id, ncuts);
         infeasible = FALSE;

         if( id == 0 )
         {
            SCIP_CALL( SCIProwCreate(&cut, blkmem, set, stat, lp, cutname, ncols, cols, vals, cons->lhs, cons->rhs,
                  SCIP_ROWORIGINTYPE_REOPT, NULL, FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcutpoolAddRow(cutpool, blkmem, set, stat, lp, cut) );

            SCIPsetDebugMsg(set, "add cut <%s> of size %d to cutpool, [lhs, rhs] = [%g,%g] to node %lld\n", cutname,
               ncols, cons->lhs, cons->rhs, SCIPnodeGetNumber(node));
         }
         else
         {
            SCIP_CALL( SCIProwCreate(&cut, blkmem, set, stat, lp, cutname, ncols, cols, vals, cons->lhs, cons->rhs,
                  SCIP_ROWORIGINTYPE_REOPT, NULL, TRUE, TRUE, TRUE) );
            SCIP_CALL( SCIPsepastoreAddCut(sepastore, blkmem, set, stat, eventqueue, eventfilter, lp, cut, FALSE, root,
                  &infeasible) );

            SCIPsetDebugMsg(set, "add cut <%s> of size %d to sepastore, [lhs, rhs] = [%g,%g] to node %lld\n", cutname,
               ncols, cons->lhs, cons->rhs, SCIPnodeGetNumber(node));
         }

         SCIP_CALL( SCIProwRelease(&cut, blkmem, set, lp) );

         if( infeasible )
            SCIPsetDebugMsg(set, "cut %d stored at node %llu (id: %u) is infeasible.\n", c, SCIPnodeGetNumber(node), id);
         else
            ++ncuts;

         SCIPsetFreeBufferArray(set, &vals);
         SCIPsetFreeBufferArray(set, &cols);

         BMSfreeBlockMemoryArrayNull(blkmem, &reoptnode->conss[c]->boundtypes, reoptnode->conss[c]->varssize);
         BMSfreeBlockMemoryArray(blkmem, &reoptnode->conss[c]->vals, reoptnode->conss[c]->varssize);
         BMSfreeBlockMemoryArray(blkmem, &reoptnode->conss[c]->vars, reoptnode->conss[c]->varssize);
         BMSfreeBlockMemory(blkmem, &reoptnode->conss[c]); /*lint !e866*/
         --reoptnode->nconss;
      }
      else
      {
#ifndef NDEBUG
         int i;
         for( i = c-1; i >= 0; i-- )
            assert(reoptnode->conss[i]->constype != REOPT_CONSTYPE_CUT);
#endif
         break;
      }
   }

   return SCIP_OKAY;
}

/** check if the LP of the given node should be solved or not */
SCIP_Bool SCIPreoptGetSolveLP(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node of the current search tree */
   )
{
   unsigned int id;

   assert(reopt != NULL);
   assert(node != NULL);

   /* get the ID */
   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* return if the node is not part of the reoptimization tree */
   if( SCIPnodeGetDepth(node) > 0 && id == 0 )
      return TRUE;

   /* return always true if the parameter is set to 1.0 */
   if( SCIPsetIsGE(set, set->reopt_objsimrootlp, 1.0) )
      return TRUE;

   /* current node is the root */
   if( id == 0 )
   {
      if( reopt->reopttree->reoptnodes[0]->nchilds > 0 )
      {
         /* the objective function has changed only slightly */
         if( SCIPsetIsGE(set, reopt->simtolastobj, set->reopt_objsimrootlp) )
            return FALSE;
      }
   }
   else
   {
      /* solve node LP if the node type is greater or equal to solvelp or there were too many bound changes at the current node */
      if( reopt->reopttree->reoptnodes[id]->nvars < set->reopt_solvelpdiff && (int) SCIPnodeGetReopttype(node) < set->reopt_solvelp )
      {
         assert(reopt->reopttree->reoptnodes[id]->nchilds > 0);
         return FALSE;
      }
   }

   return TRUE;
}

/** initialize an empty node */
void SCIPreoptnodeInit(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(reoptnode != NULL);
   assert(set != NULL);

   reoptnode->conss = NULL;
   reoptnode->nconss = 0;
   reoptnode->consssize = 0;
   reoptnode->childids = NULL;
   reoptnode->allocchildmem = 0;
   reoptnode->nchilds = 0;
   reoptnode->nvars = 0;
   reoptnode->nafterdualvars = 0;
   reoptnode->parentID = 0;
   reoptnode->dualreds = FALSE;
   reoptnode->reopttype = (unsigned int)SCIP_REOPTTYPE_NONE;
   reoptnode->varssize = 0;
   reoptnode->afterdualvarssize = 0;
   reoptnode->vars = NULL;
   reoptnode->varbounds = NULL;
   reoptnode->varboundtypes = NULL;
   reoptnode->afterdualvars = NULL;
   reoptnode->afterdualvarbounds = NULL;
   reoptnode->afterdualvarboundtypes = NULL;
   reoptnode->dualredscur = NULL;
   reoptnode->dualredsnex = NULL;
   reoptnode->lowerbound = -SCIPsetInfinity(set);
}

/** reset the given reoptimization node */
SCIP_RETCODE SCIPreoptnodeReset(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPTNODE*       reoptnode           /**< reoptimization node */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(reoptnode != NULL);

   SCIP_CALL( reoptnodeReset(reoptnode, set, blkmem) );

   return SCIP_OKAY;
}

/** delete the given reoptimization node */
SCIP_RETCODE SCIPreoptnodeDelete(
   SCIP_REOPTNODE**      reoptnode,          /**< pointer of reoptnode */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reoptnode != NULL);
   assert(blkmem != NULL);

   SCIP_CALL( reoptnodeDelete(reoptnode, blkmem) );

   return SCIP_OKAY;
}

/** add a variable to a given reoptnode */
SCIP_RETCODE SCIPreoptnodeAddBndchg(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             val,                /**< value of the variable */
   SCIP_BOUNDTYPE        boundtype           /**< boundtype of the variable */
   )
{
   int nvars;

   assert(reoptnode != NULL);
   assert(var != NULL);
   assert(blkmem != NULL);

   nvars = reoptnode->nvars;

   SCIP_CALL( reoptnodeCheckMemory(reoptnode, set, blkmem, nvars + 1, 0, 0) );

   reoptnode->vars[nvars] = var;
   reoptnode->varbounds[nvars] = val;
   reoptnode->varboundtypes[nvars] = boundtype;
   ++reoptnode->nvars;

   return SCIP_OKAY;
}

/** add a constraint to a given reoptnode */
SCIP_RETCODE SCIPreoptnodeAddCons(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            vars,               /**< variables which are part of the constraint */
   SCIP_Real*            bounds,             /**< bounds of the variables */
   SCIP_BOUNDTYPE*       boundtypes,         /**< boundtypes of the variables (or NULL is the constraint is a cut) */
   SCIP_Real             lhs,                /**< lhs of the constraint */
   SCIP_Real             rhs,                /**< rhs of the constraint */
   int                   nvars,              /**< number of variables */
   REOPT_CONSTYPE        constype,           /**< type of the constraint */
   SCIP_Bool             linear              /**< the given constraint has a linear representation */
   )
{
   int nconss;

   assert(reoptnode != NULL);
   assert(set != NULL);
   assert(vars != NULL);
   assert(bounds != NULL);
   assert(REOPT_CONSTYPE_CUT || boundtypes != NULL);
   assert(nvars > 0);
   assert(blkmem != NULL);

   /* the constraint can be interpreted as a normal bound change */
   if( nvars == 1 )
   {
      assert(constype == REOPT_CONSTYPE_DUALREDS || constype == REOPT_CONSTYPE_INFSUBTREE);

      SCIPsetDebugMsg(set, "-> constraint has size 1 -> save as normal bound change.\n");

      if( SCIPvarGetType(vars[0]) == SCIP_VARTYPE_BINARY )
      {
         SCIP_CALL( SCIPreoptnodeAddBndchg(reoptnode, set, blkmem, vars[0], 1-bounds[0],
               1-bounds[0] == 1 ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER) );
      }
      else
      {
         SCIP_Real newbound;
         SCIP_BOUNDTYPE newboundtype;

         assert(SCIPvarGetType(vars[0]) == SCIP_VARTYPE_INTEGER);

         if( boundtypes[0] == SCIP_BOUNDTYPE_UPPER )
         {
            newbound = bounds[0] + 1.0;
            assert(SCIPsetIsLE(set, newbound, SCIPvarGetUbLocal(vars[0])));

            newboundtype = SCIP_BOUNDTYPE_LOWER;
         }
         else
         {
            newbound = bounds[0] - 1.0;
            assert(SCIPsetIsGE(set, newbound, SCIPvarGetLbLocal(vars[0])));

            newboundtype = SCIP_BOUNDTYPE_UPPER;
         }

         SCIP_CALL( SCIPreoptnodeAddBndchg(reoptnode, set, blkmem, vars[0], newbound, newboundtype) );
      }
   }
   else
   {
      nconss = reoptnode->nconss;

      SCIP_CALL( reoptnodeCheckMemory(reoptnode, set, blkmem, 0, 0, nconss+1) );

      /* create the constraint */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reoptnode->conss[nconss]) ); /*lint !e866*/
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptnode->conss[nconss]->vars, vars, nvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptnode->conss[nconss]->vals, bounds, nvars) );
      if( boundtypes != NULL )
      {
         assert(!linear);
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptnode->conss[nconss]->boundtypes, boundtypes, nvars) );
      }
      else
         reoptnode->conss[nconss]->boundtypes = NULL;

      reoptnode->conss[nconss]->varssize = nvars;
      reoptnode->conss[nconss]->nvars = nvars;
      reoptnode->conss[nconss]->lhs = lhs;
      reoptnode->conss[nconss]->rhs = rhs;
      reoptnode->conss[nconss]->constype = constype;
      reoptnode->conss[nconss]->linear = linear;
      ++reoptnode->nconss;
   }
   return SCIP_OKAY;
}

/** add a constraint to the reoptimization data structure */
SCIP_RETCODE SCIPreoptAddCons(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(cons != NULL);

   /* check memory */
   if( reopt->addedconsssize == 0 )
   {
      assert(reopt->addedconss == NULL);

      reopt->addedconsssize = 10;
      SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &reopt->addedconss, reopt->addedconsssize) );
   }
   else if( reopt->naddedconss == reopt->addedconsssize )
   {
      int newsize = SCIPsetCalcMemGrowSize(set, reopt->addedconsssize+1);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->addedconss, reopt->addedconsssize, newsize) );

      /* clear the array */
      BMSclearMemoryArray(&reopt->addedconss[reopt->addedconsssize], newsize - reopt->addedconsssize); /*lint !e866 */

      reopt->addedconsssize = newsize;
   }
   assert(reopt->naddedconss < reopt->addedconsssize);
   assert(reopt->addedconss[reopt->naddedconss] == NULL);

   reopt->addedconss[reopt->naddedconss] = cons;
   reopt->consadded = TRUE;
   ++reopt->naddedconss;

   /* capture the constraint */
   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

/** save global lower and upper bounds
 *
 *  @note this method should only be called once, i.e., after fishing presolving of the first problem
 */
SCIP_RETCODE SCIPreoptSaveGlobalBounds(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(reopt != NULL);
   assert(transprob != NULL);
   assert(reopt->glblb == NULL && reopt->glbub == NULL);

   nvars = SCIPprobGetNVars(transprob);
   vars = SCIPprobGetVars(transprob);

   /* create hashmaps */
   SCIP_CALL( SCIPhashmapCreate(&reopt->glbub, blkmem, nvars) );
   SCIP_CALL( SCIPhashmapCreate(&reopt->glblb, blkmem, nvars) );

   /* store the global bounds */
   for( i = 0; i < nvars; i++ )
   {
      assert(!SCIPhashmapExists(reopt->glblb, (void*)vars[i]));
      assert(!SCIPhashmapExists(reopt->glbub, (void*)vars[i]));

      SCIP_CALL( SCIPhashmapInsertReal(reopt->glblb, (void*)vars[i], SCIPvarGetLbGlobal(vars[i])) );
      SCIP_CALL( SCIPhashmapInsertReal(reopt->glbub, (void*)vars[i], SCIPvarGetUbGlobal(vars[i])) );
   }

   return SCIP_OKAY;
}

/** save active constraints
 *
 *  @note this method can only called once, i.e., after fishing presolving of the first problem
 */
SCIP_RETCODE SCIPreoptSaveActiveConss(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_CONS** conss;
   int nconss;
   int i;

   assert(reopt != NULL);
   assert(transprob != NULL);
   assert(reopt->activeconss == NULL);

   conss = transprob->conss;
   nconss = transprob->nconss;

   /* create hashmap */
   SCIP_CALL( SCIPhashmapCreate(&reopt->activeconss, blkmem, nconss) );

   for( i = 0; i < nconss; i++ )
   {
      assert(SCIPconsIsActive(conss[i]));
      assert(!SCIPhashmapExists(reopt->activeconss, (void*)conss[i]));

      SCIPconsCapture(conss[i]);
      SCIP_CALL( SCIPhashmapInsert(reopt->activeconss, (void*)conss[i], (void*)conss[i]) );
   }

   return SCIP_OKAY;
}

/** installs global lower and upper bounds */
SCIP_RETCODE SCIPreoptInstallBounds(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(reopt != NULL);
   assert(transprob != NULL);
   assert(reopt->glblb != NULL && reopt->glbub != NULL);
   assert(SCIPprobIsTransformed(transprob));

   nvars = SCIPprobGetNVars(transprob);
   vars = SCIPprobGetVars(transprob);

   /* install global lower and upper bounds */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      assert(SCIPhashmapExists(reopt->glblb, (void*)vars[i]));
      assert(SCIPhashmapExists(reopt->glbub, (void*)vars[i]));

      lb = SCIPhashmapGetImageReal(reopt->glblb, (void*)vars[i]);
      ub = SCIPhashmapGetImageReal(reopt->glbub, (void*)vars[i]);
      assert(lb < SCIP_INVALID && ub < SCIP_INVALID);

      /* reset the global bounds back */
      SCIP_CALL( SCIPvarChgLbGlobal(vars[i], blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, lb) );
      SCIP_CALL( SCIPvarChgUbGlobal(vars[i], blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, ub) );

      /* reset the local bounds back */
      SCIP_CALL( SCIPvarChgLbLocal(vars[i], blkmem, set, stat, lp, branchcand, eventqueue, lb) );
      SCIP_CALL( SCIPvarChgUbLocal(vars[i], blkmem, set, stat, lp, branchcand, eventqueue, ub) );
   }

   return SCIP_OKAY;
}

/** reactivate globally valid constraints that were deactivated and necessary to ensure correctness */
SCIP_RETCODE SCIPreoptResetActiveConss(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic SCIP statistics */
   )
{
   int nentries;
   int i;

   assert(reopt != NULL);
   assert(reopt->activeconss != NULL);

   nentries = SCIPhashmapGetNEntries(reopt->activeconss);

   /* loop over all entries of the hashmap and reactivate deactivated constraints */
   for( i = 0; i < nentries; i++ )
   {
      SCIP_CONS* cons;
      SCIP_HASHMAPENTRY* entry = SCIPhashmapGetEntry(reopt->activeconss, i);

      if( entry == NULL )
         continue;

      cons = (SCIP_CONS*)SCIPhashmapEntryGetImage(entry);
      assert(cons != NULL);

      /* it can happen that the constraint got globally deleted */
      if( SCIPconsIsDeleted(cons) )
         cons->deleted = FALSE;

      /* to ensure that the constraint will be added to all the data structures we need to deactivate the
       * constraint first.
       */
      if( SCIPconsIsActive(cons) )
      {
         SCIP_CALL( SCIPconsDeactivate(cons, set, stat) );
      }
      SCIP_CALL( SCIPconsActivate(cons, set, stat, -1, TRUE) );
   }

   return SCIP_OKAY;
}

/** returns whether a constraint is necessary to ensure correctness and cannot be deleted */
SCIP_Bool SCIPreoptConsCanBeDeleted(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_CONS*            cons                /**< problem constraint */
   )
{
   assert(reopt != NULL);
   assert(cons != NULL);

   /* the hashmap is not initialized, we can delete all constraints */
   if( reopt->activeconss == NULL )
      return TRUE;

   return !SCIPhashmapExists(reopt->activeconss, (void*)cons);
}
