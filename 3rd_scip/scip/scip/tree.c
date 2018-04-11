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

/**@file   tree.c
 * @brief  methods for branch and bound tree
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/visual.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/relax.h"
#include "scip/var.h"
#include "scip/implics.h"
#include "scip/primal.h"
#include "scip/tree.h"
#include "scip/reopt.h"
#include "scip/conflictstore.h"
#include "scip/solve.h"
#include "scip/cons.h"
#include "scip/nodesel.h"
#include "scip/prop.h"
#include "scip/debug.h"
#include "scip/prob.h"
#include "scip/scip.h"
#include "scip/pub_message.h"
#include "lpi/lpi.h"


#define MAXREPROPMARK       511  /**< maximal subtree repropagation marker; must correspond to node data structure */


/*
 * dynamic memory arrays
 */

/** resizes children arrays to be able to store at least num nodes */
static
SCIP_RETCODE treeEnsureChildrenMem(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of node slots in array */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   if( num > tree->childrensize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&tree->children, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&tree->childrenprio, newsize) );
      tree->childrensize = newsize;
   }
   assert(num <= tree->childrensize);

   return SCIP_OKAY;
}

/** resizes path array to be able to store at least num nodes */
static
SCIP_RETCODE treeEnsurePathMem(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of node slots in path */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   if( num > tree->pathsize )
   {
      int newsize;

      newsize = SCIPsetCalcPathGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&tree->path, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&tree->pathnlpcols, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&tree->pathnlprows, newsize) );
      tree->pathsize = newsize;
   }
   assert(num <= tree->pathsize);

   return SCIP_OKAY;
}

/** resizes pendingbdchgs array to be able to store at least num nodes */
static
SCIP_RETCODE treeEnsurePendingbdchgsMem(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of node slots in path */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   if( num > tree->pendingbdchgssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&tree->pendingbdchgs, newsize) );
      tree->pendingbdchgssize = newsize;
   }
   assert(num <= tree->pendingbdchgssize);

   return SCIP_OKAY;
}




/*
 * Node methods
 */

/** node comparator for best lower bound */
SCIP_DECL_SORTPTRCOMP(SCIPnodeCompLowerbound)
{  /*lint --e{715}*/
   assert(elem1 != NULL);
   assert(elem2 != NULL);

   if( ((SCIP_NODE*)elem1)->lowerbound < ((SCIP_NODE*)elem2)->lowerbound )
      return -1;
   else if( ((SCIP_NODE*)elem1)->lowerbound > ((SCIP_NODE*)elem2)->lowerbound )
      return +1;
   else
      return 0;
}

/** increases the reference counter of the LP state in the fork */
static
void forkCaptureLPIState(
   SCIP_FORK*            fork,               /**< fork data */
   int                   nuses               /**< number to add to the usage counter */
   )
{
   assert(fork != NULL);
   assert(fork->nlpistateref >= 0);
   assert(nuses > 0);

   fork->nlpistateref += nuses;
   SCIPdebugMessage("captured LPI state of fork %p %d times -> new nlpistateref=%d\n", (void*)fork, nuses, fork->nlpistateref);
}

/** decreases the reference counter of the LP state in the fork */
static
SCIP_RETCODE forkReleaseLPIState(
   SCIP_FORK*            fork,               /**< fork data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(fork != NULL);
   assert(fork->nlpistateref > 0);
   assert(blkmem != NULL);
   assert(lp != NULL);

   fork->nlpistateref--;
   if( fork->nlpistateref == 0 )
   {
      SCIP_CALL( SCIPlpFreeState(lp, blkmem, &(fork->lpistate)) );
   }

   SCIPdebugMessage("released LPI state of fork %p -> new nlpistateref=%d\n", (void*)fork, fork->nlpistateref);

   return SCIP_OKAY;
}

/** increases the reference counter of the LP state in the subroot */
static
void subrootCaptureLPIState(
   SCIP_SUBROOT*         subroot,            /**< subroot data */
   int                   nuses               /**< number to add to the usage counter */
   )
{
   assert(subroot != NULL);
   assert(subroot->nlpistateref >= 0);
   assert(nuses > 0);

   subroot->nlpistateref += nuses;
   SCIPdebugMessage("captured LPI state of subroot %p %d times -> new nlpistateref=%d\n", 
      (void*)subroot, nuses, subroot->nlpistateref);
}

/** decreases the reference counter of the LP state in the subroot */
static
SCIP_RETCODE subrootReleaseLPIState(
   SCIP_SUBROOT*         subroot,            /**< subroot data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(subroot != NULL);
   assert(subroot->nlpistateref > 0);
   assert(blkmem != NULL);
   assert(lp != NULL);

   subroot->nlpistateref--;
   if( subroot->nlpistateref == 0 )
   {
      SCIP_CALL( SCIPlpFreeState(lp, blkmem, &(subroot->lpistate)) );
   }

   SCIPdebugMessage("released LPI state of subroot %p -> new nlpistateref=%d\n", (void*)subroot, subroot->nlpistateref);

   return SCIP_OKAY;
}

/** increases the reference counter of the LP state in the fork or subroot node */
SCIP_RETCODE SCIPnodeCaptureLPIState(
   SCIP_NODE*            node,               /**< fork/subroot node */
   int                   nuses               /**< number to add to the usage counter */
   )
{
   assert(node != NULL);

   SCIPdebugMessage("capture %d times LPI state of node #%" SCIP_LONGINT_FORMAT " at depth %d (current: %d)\n",
      nuses, SCIPnodeGetNumber(node), SCIPnodeGetDepth(node),
      SCIPnodeGetType(node) == SCIP_NODETYPE_FORK ? node->data.fork->nlpistateref : node->data.subroot->nlpistateref);

   switch( SCIPnodeGetType(node) )
   {  
   case SCIP_NODETYPE_FORK:
      forkCaptureLPIState(node->data.fork, nuses);
      break;
   case SCIP_NODETYPE_SUBROOT:
      subrootCaptureLPIState(node->data.subroot, nuses);
      break;
   default:
      SCIPerrorMessage("node for capturing the LPI state is neither fork nor subroot\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;  /*lint !e527*/
   }  /*lint !e788*/
   return SCIP_OKAY;
}

/** decreases the reference counter of the LP state in the fork or subroot node */
SCIP_RETCODE SCIPnodeReleaseLPIState(
   SCIP_NODE*            node,               /**< fork/subroot node */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(node != NULL);

   SCIPdebugMessage("release LPI state of node #%" SCIP_LONGINT_FORMAT " at depth %d (current: %d)\n",
      SCIPnodeGetNumber(node), SCIPnodeGetDepth(node),
      SCIPnodeGetType(node) == SCIP_NODETYPE_FORK ? node->data.fork->nlpistateref : node->data.subroot->nlpistateref);
   switch( SCIPnodeGetType(node) )
   {  
   case SCIP_NODETYPE_FORK:
      return forkReleaseLPIState(node->data.fork, blkmem, lp);
   case SCIP_NODETYPE_SUBROOT:
      return subrootReleaseLPIState(node->data.subroot, blkmem, lp);
   default:
      SCIPerrorMessage("node for releasing the LPI state is neither fork nor subroot\n");
      return SCIP_INVALIDDATA;
   }  /*lint !e788*/
}

/** creates probingnode data without LP information */
static
SCIP_RETCODE probingnodeCreate(
   SCIP_PROBINGNODE**    probingnode,        /**< pointer to probingnode data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(probingnode != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, probingnode) );

   (*probingnode)->lpistate = NULL;
   (*probingnode)->lpinorms = NULL;
   (*probingnode)->ninitialcols = SCIPlpGetNCols(lp);
   (*probingnode)->ninitialrows = SCIPlpGetNRows(lp);
   (*probingnode)->ncols = (*probingnode)->ninitialcols;
   (*probingnode)->nrows = (*probingnode)->ninitialrows;
   (*probingnode)->origobjvars = NULL;
   (*probingnode)->origobjvals = NULL;
   (*probingnode)->nchgdobjs = 0;

   SCIPdebugMessage("created probingnode information (%d cols, %d rows)\n", (*probingnode)->ncols, (*probingnode)->nrows);

   return SCIP_OKAY;
}

/** updates LP information in probingnode data */
static
SCIP_RETCODE probingnodeUpdate(
   SCIP_PROBINGNODE*     probingnode,        /**< probingnode data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Bool storenorms = FALSE;

   assert(probingnode != NULL);
   assert(SCIPtreeIsPathComplete(tree));
   assert(lp != NULL);

   /* free old LP state */
   if( probingnode->lpistate != NULL )
   {
      SCIP_CALL( SCIPlpFreeState(lp, blkmem, &probingnode->lpistate) );
   }

   /* free old LP norms */
   if( probingnode->lpinorms != NULL )
   {
      SCIP_CALL( SCIPlpFreeNorms(lp, blkmem, &probingnode->lpinorms) );
      probingnode->lpinorms = NULL;
      storenorms = TRUE;
   }

   /* get current LP state */
   if( lp->flushed && lp->solved )
   {
      SCIP_CALL( SCIPlpGetState(lp, blkmem, &probingnode->lpistate) );

      /* if LP norms were stored at this node before, store the new ones */
      if( storenorms )
      {
         SCIP_CALL( SCIPlpGetNorms(lp, blkmem, &probingnode->lpinorms) );
      }
      probingnode->lpwasprimfeas = lp->primalfeasible;
      probingnode->lpwasprimchecked = lp->primalchecked;
      probingnode->lpwasdualfeas = lp->dualfeasible;
      probingnode->lpwasdualchecked = lp->dualchecked;
   }
   else
      probingnode->lpistate = NULL;

   probingnode->ncols = SCIPlpGetNCols(lp);
   probingnode->nrows = SCIPlpGetNRows(lp);

   SCIPdebugMessage("updated probingnode information (%d cols, %d rows)\n", probingnode->ncols, probingnode->nrows);

   return SCIP_OKAY;
}

/** frees probingnode data */
static
SCIP_RETCODE probingnodeFree(
   SCIP_PROBINGNODE**    probingnode,        /**< probingnode data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(probingnode != NULL);
   assert(*probingnode != NULL);

   /* free the associated LP state */
   if( (*probingnode)->lpistate != NULL )
   {
      SCIP_CALL( SCIPlpFreeState(lp, blkmem, &(*probingnode)->lpistate) );
   }
   /* free the associated LP norms */
   if( (*probingnode)->lpinorms != NULL )
   {
      SCIP_CALL( SCIPlpFreeNorms(lp, blkmem, &(*probingnode)->lpinorms) );
   }

   /* free objective information */
   if( (*probingnode)->nchgdobjs > 0 )
   {
      assert((*probingnode)->origobjvars != NULL);
      assert((*probingnode)->origobjvals != NULL);

      BMSfreeMemoryArray(&(*probingnode)->origobjvars);
      BMSfreeMemoryArray(&(*probingnode)->origobjvals);
   }

   BMSfreeBlockMemory(blkmem, probingnode);

   return SCIP_OKAY;
}

/** initializes junction data */
static
SCIP_RETCODE junctionInit(
   SCIP_JUNCTION*        junction,           /**< pointer to junction data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(junction != NULL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->focusnode != NULL);

   junction->nchildren = tree->nchildren;

   /* increase the LPI state usage counter of the current LP fork */
   if( tree->focuslpstatefork != NULL )
   {
      SCIP_CALL( SCIPnodeCaptureLPIState(tree->focuslpstatefork, tree->nchildren) );
   }

   return SCIP_OKAY;
}

/** creates pseudofork data */
static
SCIP_RETCODE pseudoforkCreate(
   SCIP_PSEUDOFORK**     pseudofork,         /**< pointer to pseudofork data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(pseudofork != NULL);
   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->focusnode != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, pseudofork) );

   (*pseudofork)->addedcols = NULL;
   (*pseudofork)->addedrows = NULL;
   (*pseudofork)->naddedcols = SCIPlpGetNNewcols(lp);
   (*pseudofork)->naddedrows = SCIPlpGetNNewrows(lp);
   (*pseudofork)->nchildren = tree->nchildren;

   SCIPdebugMessage("creating pseudofork information with %d children (%d new cols, %d new rows)\n",
      (*pseudofork)->nchildren, (*pseudofork)->naddedcols, (*pseudofork)->naddedrows);

   if( (*pseudofork)->naddedcols > 0 )
   {
      /* copy the newly created columns to the pseudofork's col array */
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*pseudofork)->addedcols, SCIPlpGetNewcols(lp), (*pseudofork)->naddedcols) ); /*lint !e666*/
   }
   if( (*pseudofork)->naddedrows > 0 )
   {
      int i;

      /* copy the newly created rows to the pseudofork's row array */
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*pseudofork)->addedrows, SCIPlpGetNewrows(lp), (*pseudofork)->naddedrows) ); /*lint !e666*/

      /* capture the added rows */
      for( i = 0; i < (*pseudofork)->naddedrows; ++i )
         SCIProwCapture((*pseudofork)->addedrows[i]);
   }

   /* increase the LPI state usage counter of the current LP fork */
   if( tree->focuslpstatefork != NULL )
   {
      SCIP_CALL( SCIPnodeCaptureLPIState(tree->focuslpstatefork, tree->nchildren) );
   }

   return SCIP_OKAY;
}

/** frees pseudofork data */
static
SCIP_RETCODE pseudoforkFree(
   SCIP_PSEUDOFORK**     pseudofork,         /**< pseudofork data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(pseudofork != NULL);
   assert(*pseudofork != NULL);
   assert((*pseudofork)->nchildren == 0);
   assert(blkmem != NULL);
   assert(set != NULL);

   /* release the added rows */
   for( i = 0; i < (*pseudofork)->naddedrows; ++i )
   {
      SCIP_CALL( SCIProwRelease(&(*pseudofork)->addedrows[i], blkmem, set, lp) );
   }

   BMSfreeBlockMemoryArrayNull(blkmem, &(*pseudofork)->addedcols, (*pseudofork)->naddedcols);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*pseudofork)->addedrows, (*pseudofork)->naddedrows);
   BMSfreeBlockMemory(blkmem, pseudofork);

   return SCIP_OKAY;
}

/** creates fork data */
static
SCIP_RETCODE forkCreate(
   SCIP_FORK**           fork,               /**< pointer to fork data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(fork != NULL);
   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);
   assert(tree->nchildren < (1 << 30));
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->focusnode != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, fork) );

   SCIP_CALL( SCIPlpGetState(lp, blkmem, &((*fork)->lpistate)) );
   (*fork)->lpwasprimfeas = lp->primalfeasible;
   (*fork)->lpwasprimchecked = lp->primalchecked;
   (*fork)->lpwasdualfeas = lp->dualfeasible;
   (*fork)->lpwasdualchecked = lp->dualchecked;
   (*fork)->lpobjval = SCIPlpGetObjval(lp, set, prob);
   (*fork)->nlpistateref = 0;
   (*fork)->addedcols = NULL;
   (*fork)->addedrows = NULL;
   (*fork)->naddedcols = SCIPlpGetNNewcols(lp);
   (*fork)->naddedrows = SCIPlpGetNNewrows(lp);
   (*fork)->nchildren = (unsigned int) tree->nchildren;

   SCIPsetDebugMsg(set, "creating fork information with %u children (%d new cols, %d new rows)\n", (*fork)->nchildren, (*fork)->naddedcols, (*fork)->naddedrows);

   if( (*fork)->naddedcols > 0 )
   {
      /* copy the newly created columns to the fork's col array */
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*fork)->addedcols, SCIPlpGetNewcols(lp), (*fork)->naddedcols) ); /*lint !e666*/
   }
   if( (*fork)->naddedrows > 0 )
   {
      int i;

      /* copy the newly created rows to the fork's row array */
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*fork)->addedrows, SCIPlpGetNewrows(lp), (*fork)->naddedrows) ); /*lint !e666*/

      /* capture the added rows */
      for( i = 0; i < (*fork)->naddedrows; ++i )
         SCIProwCapture((*fork)->addedrows[i]);
   }

   /* capture the LPI state for the children */
   forkCaptureLPIState(*fork, tree->nchildren);

   return SCIP_OKAY;
}

/** frees fork data */
static
SCIP_RETCODE forkFree(
   SCIP_FORK**           fork,               /**< fork data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(fork != NULL);
   assert(*fork != NULL);
   assert((*fork)->nchildren == 0);
   assert((*fork)->nlpistateref == 0);
   assert((*fork)->lpistate == NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   /* release the added rows */
   for( i = (*fork)->naddedrows - 1; i >= 0; --i )
   {
      SCIP_CALL( SCIProwRelease(&(*fork)->addedrows[i], blkmem, set, lp) );
   }

   BMSfreeBlockMemoryArrayNull(blkmem, &(*fork)->addedcols, (*fork)->naddedcols);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*fork)->addedrows, (*fork)->naddedrows);
   BMSfreeBlockMemory(blkmem, fork);

   return SCIP_OKAY;
}

#ifdef WITHSUBROOTS /** @todo test whether subroots should be created */
/** creates subroot data */
static
SCIP_RETCODE subrootCreate(
   SCIP_SUBROOT**        subroot,            /**< pointer to subroot data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(subroot != NULL);
   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->focusnode != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, subroot) );
   (*subroot)->lpobjval = SCIPlpGetObjval(lp, set, prob);
   (*subroot)->nlpistateref = 0;
   (*subroot)->ncols = SCIPlpGetNCols(lp);
   (*subroot)->nrows = SCIPlpGetNRows(lp);
   (*subroot)->nchildren = (unsigned int) tree->nchildren;
   SCIP_CALL( SCIPlpGetState(lp, blkmem, &((*subroot)->lpistate)) );
   (*subroot)->lpwasprimfeas = lp->primalfeasible;
   (*subroot)->lpwasprimchecked = lp->primalchecked;
   (*subroot)->lpwasdualfeas = lp->dualfeasible;
   (*subroot)->lpwasdualchecked = lp->dualchecked;

   if( (*subroot)->ncols != 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*subroot)->cols, SCIPlpGetCols(lp), (*subroot)->ncols) );
   }
   else
      (*subroot)->cols = NULL;
   if( (*subroot)->nrows != 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*subroot)->rows, SCIPlpGetRows(lp), (*subroot)->nrows) );
   }
   else
      (*subroot)->rows = NULL;

   /* capture the rows of the subroot */
   for( i = 0; i < (*subroot)->nrows; ++i )
      SCIProwCapture((*subroot)->rows[i]);

   /* capture the LPI state for the children */
   subrootCaptureLPIState(*subroot, tree->nchildren);

   return SCIP_OKAY;
}
#endif

/** frees subroot */
static
SCIP_RETCODE subrootFree(
   SCIP_SUBROOT**        subroot,            /**< subroot data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(subroot != NULL);
   assert(*subroot != NULL);
   assert((*subroot)->nchildren == 0);
   assert((*subroot)->nlpistateref == 0);
   assert((*subroot)->lpistate == NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   /* release the rows of the subroot */
   for( i = 0; i < (*subroot)->nrows; ++i )
   {
      SCIP_CALL( SCIProwRelease(&(*subroot)->rows[i], blkmem, set, lp) );
   }

   BMSfreeBlockMemoryArrayNull(blkmem, &(*subroot)->cols, (*subroot)->ncols);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*subroot)->rows, (*subroot)->nrows);
   BMSfreeBlockMemory(blkmem, subroot);

   return SCIP_OKAY;
}

/** removes given sibling node from the siblings array */
static
void treeRemoveSibling(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            sibling             /**< sibling node to remove */
   )
{
   int delpos;

   assert(tree != NULL);
   assert(sibling != NULL);
   assert(SCIPnodeGetType(sibling) == SCIP_NODETYPE_SIBLING);
   assert(sibling->data.sibling.arraypos >= 0 && sibling->data.sibling.arraypos < tree->nsiblings);
   assert(tree->siblings[sibling->data.sibling.arraypos] == sibling);
   assert(SCIPnodeGetType(tree->siblings[tree->nsiblings-1]) == SCIP_NODETYPE_SIBLING);

   delpos = sibling->data.sibling.arraypos;

   /* move last sibling in array to position of removed sibling */
   tree->siblings[delpos] = tree->siblings[tree->nsiblings-1];
   tree->siblingsprio[delpos] = tree->siblingsprio[tree->nsiblings-1];
   tree->siblings[delpos]->data.sibling.arraypos = delpos;
   sibling->data.sibling.arraypos = -1;
   tree->nsiblings--;
}

/** adds given child node to children array of focus node */
static
SCIP_RETCODE treeAddChild(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            child,              /**< child node to add */
   SCIP_Real             nodeselprio         /**< node selection priority of child node */
   )
{
   assert(tree != NULL);
   assert(child != NULL);
   assert(SCIPnodeGetType(child) == SCIP_NODETYPE_CHILD);
   assert(child->data.child.arraypos == -1);

   SCIP_CALL( treeEnsureChildrenMem(tree, set, tree->nchildren+1) );
   tree->children[tree->nchildren] = child;
   tree->childrenprio[tree->nchildren] = nodeselprio;
   child->data.child.arraypos = tree->nchildren;
   tree->nchildren++;

   return SCIP_OKAY;
}

/** removes given child node from the children array */
static
void treeRemoveChild(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            child               /**< child node to remove */
   )
{
   int delpos;

   assert(tree != NULL);
   assert(child != NULL);
   assert(SCIPnodeGetType(child) == SCIP_NODETYPE_CHILD);
   assert(child->data.child.arraypos >= 0 && child->data.child.arraypos < tree->nchildren);
   assert(tree->children[child->data.child.arraypos] == child);
   assert(SCIPnodeGetType(tree->children[tree->nchildren-1]) == SCIP_NODETYPE_CHILD);

   delpos = child->data.child.arraypos;

   /* move last child in array to position of removed child */
   tree->children[delpos] = tree->children[tree->nchildren-1];
   tree->childrenprio[delpos] = tree->childrenprio[tree->nchildren-1];
   tree->children[delpos]->data.child.arraypos = delpos;
   child->data.child.arraypos = -1;
   tree->nchildren--;
}

/** makes node a child of the given parent node, which must be the focus node; if the child is a probing node,
 *  the parent node can also be a refocused node or a probing node
 */
static
SCIP_RETCODE nodeAssignParent(
   SCIP_NODE*            node,               /**< child node */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            parent,             /**< parent (= focus) node (or NULL, if node is root) */
   SCIP_Real             nodeselprio         /**< node selection priority of child node */
   )
{
   assert(node != NULL);
   assert(node->parent == NULL);
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
   assert(node->conssetchg == NULL);
   assert(node->domchg == NULL);
   assert(SCIPsetIsInfinity(set, -node->lowerbound)); /* node was just created */
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] == parent);
   assert(parent == tree->focusnode || SCIPnodeGetType(parent) == SCIP_NODETYPE_PROBINGNODE);
   assert(parent == NULL || SCIPnodeGetType(parent) == SCIP_NODETYPE_FOCUSNODE
      || (SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE
         && (SCIPnodeGetType(parent) == SCIP_NODETYPE_REFOCUSNODE
            || SCIPnodeGetType(parent) == SCIP_NODETYPE_PROBINGNODE)));

   /* link node to parent */
   node->parent = parent;
   if( parent != NULL )
   {
      assert(parent->lowerbound <= parent->estimate);
      node->lowerbound = parent->lowerbound;
      node->estimate = parent->estimate;
      node->depth = parent->depth+1; /*lint !e732*/
      if( parent->depth >= SCIP_MAXTREEDEPTH )
      {
         SCIPerrorMessage("maximal depth level exceeded\n");
         return SCIP_MAXDEPTHLEVEL;
      }
   }
   SCIPsetDebugMsg(set, "assigning parent #%" SCIP_LONGINT_FORMAT " to node #%" SCIP_LONGINT_FORMAT " at depth %d\n",
      parent != NULL ? SCIPnodeGetNumber(parent) : -1, SCIPnodeGetNumber(node), SCIPnodeGetDepth(node));

   /* register node in the childlist of the focus (the parent) node */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD )
   {
      assert(parent == NULL || SCIPnodeGetType(parent) == SCIP_NODETYPE_FOCUSNODE);
      SCIP_CALL( treeAddChild(tree, set, node, nodeselprio) );
   }

   return SCIP_OKAY;
}

/** decreases number of children of the parent, frees it if no children are left */
static
SCIP_RETCODE nodeReleaseParent(
   SCIP_NODE*            node,               /**< child node */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_NODE* parent;

   assert(node != NULL);
   assert(blkmem != NULL);
   assert(tree != NULL);

   SCIPsetDebugMsg(set, "releasing parent-child relationship of node #%" SCIP_LONGINT_FORMAT " at depth %d of type %d with parent #%" SCIP_LONGINT_FORMAT " of type %d\n",
      SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), SCIPnodeGetType(node),
      node->parent != NULL ? SCIPnodeGetNumber(node->parent) : -1,
      node->parent != NULL ? (int)SCIPnodeGetType(node->parent) : -1);
   parent = node->parent;
   if( parent != NULL )
   {
      SCIP_Bool freeParent;
      SCIP_Bool singleChild;

      freeParent = FALSE;
      singleChild = FALSE;
      switch( SCIPnodeGetType(parent) )
      {
      case SCIP_NODETYPE_FOCUSNODE:
         assert(parent->active);
         assert(SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE
            || SCIPnodeGetType(node) == SCIP_NODETYPE_LEAF);
         if( SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD )
            treeRemoveChild(tree, node);
         /* don't kill the focus node at this point => freeParent = FALSE */
         break;
      case SCIP_NODETYPE_PROBINGNODE:
         assert(SCIPtreeProbing(tree));
         /* probing nodes have to be freed individually => freeParent = FALSE */
         break;
      case SCIP_NODETYPE_SIBLING:
         SCIPerrorMessage("sibling cannot be a parent node\n");
         return SCIP_INVALIDDATA;
      case SCIP_NODETYPE_CHILD:
         SCIPerrorMessage("child cannot be a parent node\n");
         return SCIP_INVALIDDATA;
      case SCIP_NODETYPE_LEAF:
         SCIPerrorMessage("leaf cannot be a parent node\n");
         return SCIP_INVALIDDATA;
      case SCIP_NODETYPE_DEADEND:
         SCIPerrorMessage("dead-end cannot be a parent node\n");
         return SCIP_INVALIDDATA;
      case SCIP_NODETYPE_JUNCTION:
         assert(parent->data.junction.nchildren > 0);
         parent->data.junction.nchildren--;
         freeParent = (parent->data.junction.nchildren == 0); /* free parent if it has no more children */
         singleChild = (parent->data.junction.nchildren == 1);
         break;
      case SCIP_NODETYPE_PSEUDOFORK:
         assert(parent->data.pseudofork != NULL);
         assert(parent->data.pseudofork->nchildren > 0);
         parent->data.pseudofork->nchildren--;
         freeParent = (parent->data.pseudofork->nchildren == 0); /* free parent if it has no more children */
         singleChild = (parent->data.pseudofork->nchildren == 1);
         break;
      case SCIP_NODETYPE_FORK:
         assert(parent->data.fork != NULL);
         assert(parent->data.fork->nchildren > 0);
         parent->data.fork->nchildren--;
         freeParent = (parent->data.fork->nchildren == 0); /* free parent if it has no more children */
         singleChild = (parent->data.fork->nchildren == 1);
         break;
      case SCIP_NODETYPE_SUBROOT:
         assert(parent->data.subroot != NULL);
         assert(parent->data.subroot->nchildren > 0);
         parent->data.subroot->nchildren--;
         freeParent = (parent->data.subroot->nchildren == 0); /* free parent if it has no more children */
         singleChild = (parent->data.subroot->nchildren == 1);
         break;
      case SCIP_NODETYPE_REFOCUSNODE:
         /* the only possible child a refocused node can have in its refocus state is the probing root node;
          * we don't want to free the refocused node, because we first have to convert it back to its original
          * type (where it possibly has children) => freeParent = FALSE
          */
         assert(SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
         assert(!SCIPtreeProbing(tree));
         break;
      default:
         SCIPerrorMessage("unknown node type %d\n", SCIPnodeGetType(parent));
         return SCIP_INVALIDDATA;
      }

      /* free parent, if it is not on the current active path */
      if( freeParent && !parent->active )
      {
         SCIP_CALL( SCIPnodeFree(&node->parent, blkmem, set, stat, eventqueue, tree, lp) );
      }

      /* update the effective root depth
       * in reoptimization we must not increase the effective root depth
       */
      assert(tree->effectiverootdepth >= 0);
      if( singleChild && SCIPnodeGetDepth(parent) == tree->effectiverootdepth && !set->reopt_enable )
      {
         tree->effectiverootdepth++;
         SCIPsetDebugMsg(set, "unlinked node #%" SCIP_LONGINT_FORMAT " in depth %d -> new effective root depth: %d\n",
            SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), tree->effectiverootdepth);
      }
   }

   return SCIP_OKAY;
}

/** creates a node data structure */
static
SCIP_RETCODE nodeCreate(
   SCIP_NODE**           node,               /**< pointer to node data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(node != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, node) );
   (*node)->parent = NULL;
   (*node)->conssetchg = NULL;
   (*node)->domchg = NULL;
   (*node)->number = 0;
   (*node)->lowerbound = -SCIPsetInfinity(set);
   (*node)->estimate = -SCIPsetInfinity(set);
   (*node)->reoptid = 0;
   (*node)->reopttype = (unsigned int) SCIP_REOPTTYPE_NONE;
   (*node)->depth = 0;
   (*node)->active = FALSE;
   (*node)->cutoff = FALSE;
   (*node)->reprop = FALSE;
   (*node)->repropsubtreemark = 0;

   return SCIP_OKAY;
}

/** creates a child node of the focus node */
SCIP_RETCODE SCIPnodeCreateChild(
   SCIP_NODE**           node,               /**< pointer to node data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Real             nodeselprio,        /**< node selection priority of new node */
   SCIP_Real             estimate            /**< estimate for (transformed) objective value of best feasible solution in subtree */
   )
{
   assert(node != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->pathlen == 0 || tree->path != NULL);
   assert((tree->pathlen == 0) == (tree->focusnode == NULL));
   assert(tree->focusnode == NULL || tree->focusnode == tree->path[tree->pathlen-1]);
   assert(tree->focusnode == NULL || SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);

   stat->ncreatednodes++;
   stat->ncreatednodesrun++;

   /* create the node data structure */
   SCIP_CALL( nodeCreate(node, blkmem, set) );
   (*node)->number = stat->ncreatednodesrun;

   /* mark node to be a child node */
   (*node)->nodetype = SCIP_NODETYPE_CHILD; /*lint !e641*/
   (*node)->data.child.arraypos = -1;

   /* make focus node the parent of the new child */
   SCIP_CALL( nodeAssignParent(*node, blkmem, set, tree, tree->focusnode, nodeselprio) );

   /* update the estimate of the child */
   SCIPnodeSetEstimate(*node, set, estimate);

   /* output node creation to visualization file */
   SCIP_CALL( SCIPvisualNewChild(stat->visual, set, stat, *node) );

   SCIPsetDebugMsg(set, "created child node #%" SCIP_LONGINT_FORMAT " at depth %u (prio: %g)\n", SCIPnodeGetNumber(*node), (*node)->depth, nodeselprio);

   return SCIP_OKAY;
}

/** frees node */
SCIP_RETCODE SCIPnodeFree(
   SCIP_NODE**           node,               /**< node data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Bool isroot;

   assert(node != NULL);
   assert(*node != NULL);
   assert(!(*node)->active);
   assert(blkmem != NULL);
   assert(tree != NULL);

   SCIPsetDebugMsg(set, "free node #%" SCIP_LONGINT_FORMAT " at depth %d of type %d\n", SCIPnodeGetNumber(*node), SCIPnodeGetDepth(*node), SCIPnodeGetType(*node));

   /* inform solution debugger, that the node has been freed */
   SCIP_CALL( SCIPdebugRemoveNode(blkmem, set, *node) );

   /* check, if the node to be freed is the root node */
   isroot = (SCIPnodeGetDepth(*node) == 0);

   /* free nodetype specific data, and release no longer needed LPI states */
   switch( SCIPnodeGetType(*node) )
   {
   case SCIP_NODETYPE_FOCUSNODE:
      assert(tree->focusnode == *node);
      assert(!SCIPtreeProbing(tree));
      SCIPerrorMessage("cannot free focus node - has to be converted into a dead end first\n");
      return SCIP_INVALIDDATA;
   case SCIP_NODETYPE_PROBINGNODE:
      assert(SCIPtreeProbing(tree));
      assert(SCIPnodeGetDepth(tree->probingroot) <= SCIPnodeGetDepth(*node));
      assert(SCIPnodeGetDepth(*node) > 0);
      SCIP_CALL( probingnodeFree(&((*node)->data.probingnode), blkmem, lp) );
      break;
   case SCIP_NODETYPE_SIBLING:
      assert((*node)->data.sibling.arraypos >= 0);
      assert((*node)->data.sibling.arraypos < tree->nsiblings);
      assert(tree->siblings[(*node)->data.sibling.arraypos] == *node);
      if( tree->focuslpstatefork != NULL )
      {
         assert(SCIPnodeGetType(tree->focuslpstatefork) == SCIP_NODETYPE_FORK
            || SCIPnodeGetType(tree->focuslpstatefork) == SCIP_NODETYPE_SUBROOT);
         SCIP_CALL( SCIPnodeReleaseLPIState(tree->focuslpstatefork, blkmem, lp) );
      }
      treeRemoveSibling(tree, *node);
      break;
   case SCIP_NODETYPE_CHILD:
      assert((*node)->data.child.arraypos >= 0);
      assert((*node)->data.child.arraypos < tree->nchildren);
      assert(tree->children[(*node)->data.child.arraypos] == *node);
      /* The children capture the LPI state at the moment, where the focus node is
       * converted into a junction, pseudofork, fork, or subroot, and a new node is focused.
       * At the same time, they become siblings or leaves, such that freeing a child
       * of the focus node doesn't require to release the LPI state;
       * we don't need to call treeRemoveChild(), because this is done in nodeReleaseParent()
       */
      break;
   case SCIP_NODETYPE_LEAF:
      if( (*node)->data.leaf.lpstatefork != NULL )
      {
         SCIP_CALL( SCIPnodeReleaseLPIState((*node)->data.leaf.lpstatefork, blkmem, lp) );
      }
      break;
   case SCIP_NODETYPE_DEADEND:
   case SCIP_NODETYPE_JUNCTION:
      break;
   case SCIP_NODETYPE_PSEUDOFORK:
      SCIP_CALL( pseudoforkFree(&((*node)->data.pseudofork), blkmem, set, lp) );
      break;
   case SCIP_NODETYPE_FORK:

      /* release special root LPI state capture which is used to keep the root LPI state over the whole solving
       * process
       */
      if( isroot )
      {
         SCIP_CALL( SCIPnodeReleaseLPIState(*node, blkmem, lp) );
      }
      SCIP_CALL( forkFree(&((*node)->data.fork), blkmem, set, lp) );
      break;
   case SCIP_NODETYPE_SUBROOT:
      SCIP_CALL( subrootFree(&((*node)->data.subroot), blkmem, set, lp) );
      break;
   case SCIP_NODETYPE_REFOCUSNODE:
      SCIPerrorMessage("cannot free node as long it is refocused\n");
      return SCIP_INVALIDDATA;
   default:
      SCIPerrorMessage("unknown node type %d\n", SCIPnodeGetType(*node));
      return SCIP_INVALIDDATA;
   }

   /* free common data */
   SCIP_CALL( SCIPconssetchgFree(&(*node)->conssetchg, blkmem, set) );
   SCIP_CALL( SCIPdomchgFree(&(*node)->domchg, blkmem, set, eventqueue, lp) );
   SCIP_CALL( nodeReleaseParent(*node, blkmem, set, stat, eventqueue, tree, lp) );

   /* check, if the node is the current probing root */
   if( *node == tree->probingroot )
   {
      assert(SCIPnodeGetType(*node) == SCIP_NODETYPE_PROBINGNODE);
      tree->probingroot = NULL;
   }

   BMSfreeBlockMemory(blkmem, node);

   /* delete the tree's root node pointer, if the freed node was the root */
   if( isroot )
      tree->root = NULL;

   return SCIP_OKAY;
}

/** cuts off node and whole sub tree from branch and bound tree */
SCIP_RETCODE SCIPnodeCutoff(
   SCIP_NODE*            node,               /**< node that should be cut off */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_Real oldbound;

   assert(node != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);

   if( set->reopt_enable )
   {
      assert(reopt != NULL);
      /* check if the node should be stored for reoptimization */
      SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, node, SCIP_EVENTTYPE_NODEINFEASIBLE, lp, SCIPlpGetSolstat(lp),
            tree->root == node, tree->focusnode == node, node->lowerbound, tree->effectiverootdepth) );
   }

   oldbound = node->lowerbound;
   node->cutoff = TRUE;
   node->lowerbound = SCIPsetInfinity(set);
   node->estimate = SCIPsetInfinity(set);
   if( node->active )
      tree->cutoffdepth = MIN(tree->cutoffdepth, (int)node->depth);

   /* update primal integral */
   if( node->depth == 0 )
   {
      stat->rootlowerbound = SCIPsetInfinity(set);
      if( set->misc_calcintegral )
         SCIPstatUpdatePrimalDualIntegral(stat, set, transprob, origprob, SCIPsetInfinity(set), SCIPsetInfinity(set));
   }
   else if( set->misc_calcintegral && SCIPsetIsEQ(set, oldbound, stat->lastlowerbound) )
   {
      SCIP_Real lowerbound;
      lowerbound = SCIPtreeGetLowerbound(tree, set);

      /* updating the primal integral is only necessary if dual bound has increased since last evaluation */
      if( lowerbound > stat->lastlowerbound )
         SCIPstatUpdatePrimalDualIntegral(stat, set, transprob, origprob, SCIPsetInfinity(set), SCIPsetInfinity(set));
   }

   SCIPvisualCutoffNode(stat->visual, set, stat, node, TRUE);

   SCIPsetDebugMsg(set, "cutting off %s node #%" SCIP_LONGINT_FORMAT " at depth %d (cutoffdepth: %d)\n",
      node->active ? "active" : "inactive", SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), tree->cutoffdepth);

   return SCIP_OKAY;
}

/** marks node, that propagation should be applied again the next time, a node of its subtree is focused */
void SCIPnodePropagateAgain(
   SCIP_NODE*            node,               /**< node that should be propagated again */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(node != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);

   if( !node->reprop )
   {
      node->reprop = TRUE;
      if( node->active )
         tree->repropdepth = MIN(tree->repropdepth, (int)node->depth);

      SCIPvisualMarkedRepropagateNode(stat->visual, stat, node);

      SCIPsetDebugMsg(set, "marked %s node #%" SCIP_LONGINT_FORMAT " at depth %d to be propagated again (repropdepth: %d)\n",
         node->active ? "active" : "inactive", SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), tree->repropdepth);
   }
}

/** marks node, that it is completely propagated in the current repropagation subtree level */
void SCIPnodeMarkPropagated(
   SCIP_NODE*            node,               /**< node that should be marked to be propagated */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(node != NULL);
   assert(tree != NULL);

   if( node->parent != NULL )
      node->repropsubtreemark = node->parent->repropsubtreemark; /*lint !e732*/
   node->reprop = FALSE;

   /* if the node was the highest repropagation node in the path, update the repropdepth in the tree data */
   if( node->active && node->depth == tree->repropdepth )
   {
      do
      {
         assert(tree->repropdepth < tree->pathlen);
         assert(tree->path[tree->repropdepth]->active);
         assert(!tree->path[tree->repropdepth]->reprop);
         tree->repropdepth++;
      }
      while( tree->repropdepth < tree->pathlen && !tree->path[tree->repropdepth]->reprop );
      if( tree->repropdepth == tree->pathlen )
         tree->repropdepth = INT_MAX;
   }
}

/** moves the subtree repropagation counter to the next value */
static
void treeNextRepropsubtreecount(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   tree->repropsubtreecount++;
   tree->repropsubtreecount %= (MAXREPROPMARK+1);
}

/** applies propagation on the node, that was marked to be propagated again */
static
SCIP_RETCODE nodeRepropagate(
   SCIP_NODE*            node,               /**< node to apply propagation on */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_NODETYPE oldtype;
   SCIP_NODE* oldfocusnode;
   SCIP_NODE* oldfocuslpfork;
   SCIP_NODE* oldfocuslpstatefork;
   SCIP_NODE* oldfocussubroot;
   SCIP_Longint oldfocuslpstateforklpcount;
   int oldnchildren;
   int oldnsiblings;
   SCIP_Bool oldfocusnodehaslp;
   SCIP_Longint oldnboundchgs;
   SCIP_Bool initialreprop;
   SCIP_Bool clockisrunning;

   assert(node != NULL);
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_FOCUSNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_JUNCTION
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_PSEUDOFORK
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_FORK
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_SUBROOT);
   assert(node->active);
   assert(node->reprop || node->repropsubtreemark != node->parent->repropsubtreemark);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(SCIPeventqueueIsDelayed(eventqueue));
   assert(cutoff != NULL);

   SCIPsetDebugMsg(set, "propagating again node #%" SCIP_LONGINT_FORMAT " at depth %d\n", SCIPnodeGetNumber(node), SCIPnodeGetDepth(node));
   initialreprop = node->reprop;

   SCIPvisualRepropagatedNode(stat->visual, stat, node);

   /* process the delayed events in order to flush the problem changes */
   SCIP_CALL( SCIPeventqueueProcess(eventqueue, blkmem, set, primal, lp, branchcand, eventfilter) );

   /* stop node activation timer */
   clockisrunning = SCIPclockIsRunning(stat->nodeactivationtime);
   if( clockisrunning )
      SCIPclockStop(stat->nodeactivationtime, set);

   /* mark the node refocused and temporarily install it as focus node */
   oldtype = (SCIP_NODETYPE)node->nodetype;
   oldfocusnode = tree->focusnode;
   oldfocuslpfork = tree->focuslpfork;
   oldfocuslpstatefork = tree->focuslpstatefork;
   oldfocussubroot = tree->focussubroot;
   oldfocuslpstateforklpcount = tree->focuslpstateforklpcount;
   oldnchildren = tree->nchildren;
   oldnsiblings = tree->nsiblings;
   oldfocusnodehaslp = tree->focusnodehaslp;
   node->nodetype = SCIP_NODETYPE_REFOCUSNODE; /*lint !e641*/
   tree->focusnode = node;
   tree->focuslpfork = NULL;
   tree->focuslpstatefork = NULL;
   tree->focussubroot = NULL;
   tree->focuslpstateforklpcount = -1;
   tree->nchildren = 0;
   tree->nsiblings = 0;
   tree->focusnodehaslp = FALSE;

   /* propagate the domains again */
   oldnboundchgs = stat->nboundchgs;
   SCIP_CALL( SCIPpropagateDomains(blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand,
         eventqueue, conflict, cliquetable, SCIPnodeGetDepth(node), 0, SCIP_PROPTIMING_ALWAYS, cutoff) );
   assert(!node->reprop || *cutoff);
   assert(node->parent == NULL || node->repropsubtreemark == node->parent->repropsubtreemark);
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_REFOCUSNODE);
   assert(tree->focusnode == node);
   assert(tree->focuslpfork == NULL);
   assert(tree->focuslpstatefork == NULL);
   assert(tree->focussubroot == NULL);
   assert(tree->focuslpstateforklpcount == -1);
   assert(tree->nchildren == 0);
   assert(tree->nsiblings == 0);
   assert(tree->focusnodehaslp == FALSE);
   assert(stat->nboundchgs >= oldnboundchgs);
   stat->nreprops++;
   stat->nrepropboundchgs += stat->nboundchgs - oldnboundchgs;
   if( *cutoff )
      stat->nrepropcutoffs++;

   SCIPsetDebugMsg(set, "repropagation %" SCIP_LONGINT_FORMAT " at depth %u changed %" SCIP_LONGINT_FORMAT " bounds (total reprop bound changes: %" SCIP_LONGINT_FORMAT "), cutoff: %u\n",
      stat->nreprops, node->depth, stat->nboundchgs - oldnboundchgs, stat->nrepropboundchgs, *cutoff);

   /* if a propagation marked with the reprop flag was successful, we want to repropagate the whole subtree */
   /**@todo because repropsubtree is only a bit flag, we cannot mark a whole subtree a second time for
    *       repropagation; use a (small) part of the node's bits to be able to store larger numbers,
    *       and update tree->repropsubtreelevel with this number
    */
   if( initialreprop && !(*cutoff) && stat->nboundchgs > oldnboundchgs )
   {
      treeNextRepropsubtreecount(tree);
      node->repropsubtreemark = tree->repropsubtreecount; /*lint !e732*/
      SCIPsetDebugMsg(set, "initial repropagation at depth %u changed %" SCIP_LONGINT_FORMAT " bounds -> repropagating subtree (new mark: %d)\n",
         node->depth, stat->nboundchgs - oldnboundchgs, tree->repropsubtreecount);
      assert((int)(node->repropsubtreemark) == tree->repropsubtreecount); /* bitfield must be large enough */
   }

   /* reset the node's type and reinstall the old focus node */
   node->nodetype = oldtype; /*lint !e641*/
   tree->focusnode = oldfocusnode;
   tree->focuslpfork = oldfocuslpfork;
   tree->focuslpstatefork = oldfocuslpstatefork;
   tree->focussubroot = oldfocussubroot;
   tree->focuslpstateforklpcount = oldfocuslpstateforklpcount;
   tree->nchildren = oldnchildren;
   tree->nsiblings = oldnsiblings;
   tree->focusnodehaslp = oldfocusnodehaslp;

   /* make the domain change data static again to save memory */
   if( (SCIP_NODETYPE)node->nodetype != SCIP_NODETYPE_FOCUSNODE )
   {
      SCIP_CALL( SCIPdomchgMakeStatic(&node->domchg, blkmem, set, eventqueue, lp) );
   }

   /* start node activation timer again */
   if( clockisrunning )
      SCIPclockStart(stat->nodeactivationtime, set);

   /* delay events in path switching */
   SCIP_CALL( SCIPeventqueueDelay(eventqueue) );

   /* mark the node to be cut off if a cutoff was detected */
   if( *cutoff )
   {
      SCIP_CALL( SCIPnodeCutoff(node, set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
   }

   return SCIP_OKAY;
}

/** informs node, that it is now on the active path and applies any domain and constraint set changes */
static
SCIP_RETCODE nodeActivate(
   SCIP_NODE*            node,               /**< node to activate */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reotimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   assert(node != NULL);
   assert(!node->active);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(cutoff != NULL);

   SCIPsetDebugMsg(set, "activate node #%" SCIP_LONGINT_FORMAT " at depth %d of type %d (reprop subtree mark: %u)\n",
      SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), SCIPnodeGetType(node), node->repropsubtreemark);

   /* apply domain and constraint set changes */
   SCIP_CALL( SCIPconssetchgApply(node->conssetchg, blkmem, set, stat, node->depth,
         (SCIPnodeGetType(node) == SCIP_NODETYPE_FOCUSNODE)) );
   SCIP_CALL( SCIPdomchgApply(node->domchg, blkmem, set, stat, lp, branchcand, eventqueue, node->depth, cutoff) );

   /* mark node active */
   node->active = TRUE;
   stat->nactivatednodes++;

   /* check if the domain change produced a cutoff */
   if( *cutoff )
   {
      /* try to repropagate the node to see, if the propagation also leads to a conflict and a conflict constraint
       * could be generated; if propagation conflict analysis is turned off, repropagating the node makes no
       * sense, since it is already cut off
       */
      node->reprop = set->conf_enable && set->conf_useprop;

      /* mark the node to be cut off */
      SCIP_CALL( SCIPnodeCutoff(node, set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
   }

   /* propagate node again, if the reprop flag is set; in the new focus node, no repropagation is necessary, because
    * the focus node is propagated anyways
    */
   if( SCIPnodeGetType(node) != SCIP_NODETYPE_FOCUSNODE
      && (node->reprop || (node->parent != NULL && node->repropsubtreemark != node->parent->repropsubtreemark)) )
   {
      SCIP_Bool propcutoff;

      SCIP_CALL( nodeRepropagate(node, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, conflict,
            eventfilter, eventqueue, cliquetable, &propcutoff) );
      *cutoff = *cutoff || propcutoff;
   }

   return SCIP_OKAY;
}

/** informs node, that it is no longer on the active path and undoes any domain and constraint set changes */
static
SCIP_RETCODE nodeDeactivate(
   SCIP_NODE*            node,               /**< node to deactivate */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_Bool freeNode;

   assert(node != NULL);
   assert(node->active);
   assert(tree != NULL);
   assert(SCIPnodeGetType(node) != SCIP_NODETYPE_FOCUSNODE);

   SCIPsetDebugMsg(set, "deactivate node #%" SCIP_LONGINT_FORMAT " at depth %d of type %d (reprop subtree mark: %u)\n",
      SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), SCIPnodeGetType(node), node->repropsubtreemark);

   /* undo domain and constraint set changes */
   SCIP_CALL( SCIPdomchgUndo(node->domchg, blkmem, set, stat, lp, branchcand, eventqueue) );
   SCIP_CALL( SCIPconssetchgUndo(node->conssetchg, blkmem, set, stat) );

   /* mark node inactive */
   node->active = FALSE;

   /* count number of deactivated nodes (ignoring probing switches) */
   if( !SCIPtreeProbing(tree) )
      stat->ndeactivatednodes++;

   /* free node if it is a dead-end node, i.e., has no children */
   switch( SCIPnodeGetType(node) )   
   {
   case SCIP_NODETYPE_FOCUSNODE:
   case SCIP_NODETYPE_PROBINGNODE:
   case SCIP_NODETYPE_SIBLING:
   case SCIP_NODETYPE_CHILD:
   case SCIP_NODETYPE_LEAF:
   case SCIP_NODETYPE_DEADEND:
   case SCIP_NODETYPE_REFOCUSNODE:
      freeNode = FALSE;
      break;
   case SCIP_NODETYPE_JUNCTION:
      freeNode = (node->data.junction.nchildren == 0); 
      break;
   case SCIP_NODETYPE_PSEUDOFORK:
      freeNode = (node->data.pseudofork->nchildren == 0); 
      break;
   case SCIP_NODETYPE_FORK:
      freeNode = (node->data.fork->nchildren == 0); 
      break;
   case SCIP_NODETYPE_SUBROOT:
      freeNode = (node->data.subroot->nchildren == 0); 
      break;
   default:
      SCIPerrorMessage("unknown node type %d\n", SCIPnodeGetType(node));
      return SCIP_INVALIDDATA;
   }
   if( freeNode ) 
   {
      SCIP_CALL( SCIPnodeFree(&node, blkmem, set, stat, eventqueue, tree, lp) );
   }

   return SCIP_OKAY;
}

/** adds constraint locally to the node and captures it; activates constraint, if node is active;
 *  if a local constraint is added to the root node, it is automatically upgraded into a global constraint
 */
SCIP_RETCODE SCIPnodeAddCons(
   SCIP_NODE*            node,               /**< node to add constraint to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   assert(node != NULL);
   assert(cons != NULL);
   assert(cons->validdepth <= SCIPnodeGetDepth(node));
   assert(tree != NULL);
   assert(tree->effectiverootdepth >= 0);
   assert(tree->root != NULL);
   assert(SCIPconsIsGlobal(cons) || SCIPnodeGetDepth(node) > tree->effectiverootdepth);

#ifndef NDEBUG
   /* check if we add this constraint to the same scip, where we create the constraint */
   if( cons->scip != set->scip )
   {
      SCIPerrorMessage("try to add a constraint of another scip instance\n");
      return SCIP_INVALIDDATA;
   }
#endif

   /* add constraint addition to the node's constraint set change data, and activate constraint if node is active */
   SCIP_CALL( SCIPconssetchgAddAddedCons(&node->conssetchg, blkmem, set, stat, cons, node->depth,
         (SCIPnodeGetType(node) == SCIP_NODETYPE_FOCUSNODE), node->active) );
   assert(node->conssetchg != NULL);
   assert(node->conssetchg->addedconss != NULL);
   assert(!node->active || SCIPconsIsActive(cons));

   /* if the constraint is added to an active node which is not a probing node, increment the corresponding counter */
   if( node->active && SCIPnodeGetType(node) != SCIP_NODETYPE_PROBINGNODE )
      stat->nactiveconssadded++;

   return SCIP_OKAY;
}

/** locally deletes constraint at the given node by disabling its separation, enforcing, and propagation capabilities
 *  at the node; captures constraint; disables constraint, if node is active
 */
SCIP_RETCODE SCIPnodeDelCons(
   SCIP_NODE*            node,               /**< node to add constraint to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_CONS*            cons                /**< constraint to locally delete */
   )
{
   assert(node != NULL);
   assert(tree != NULL);
   assert(cons != NULL);

   SCIPsetDebugMsg(set, "disabling constraint <%s> at node at depth %u\n", cons->name, node->depth);

   /* add constraint disabling to the node's constraint set change data */
   SCIP_CALL( SCIPconssetchgAddDisabledCons(&node->conssetchg, blkmem, set, cons) );
   assert(node->conssetchg != NULL);
   assert(node->conssetchg->disabledconss != NULL);

   /* disable constraint, if node is active */
   if( node->active && cons->enabled && !cons->updatedisable )
   {
      SCIP_CALL( SCIPconsDisable(cons, set, stat) );
   }

   return SCIP_OKAY;
}

/** returns all constraints added to a given node */
void SCIPnodeGetAddedConss(
   SCIP_NODE*            node,               /**< node */
   SCIP_CONS**           addedconss,         /**< array to store the constraints */
   int*                  naddedconss,        /**< number of added constraints */
   int                   addedconsssize      /**< size of the constraint array */
   )
{
   int cons;

   assert(node != NULL );
   assert(node->conssetchg != NULL);
   assert(node->conssetchg->addedconss != NULL);
   assert(node->conssetchg->naddedconss >= 1);

   *naddedconss = node->conssetchg->naddedconss;

   /* check the size and return if the array is not large enough */
   if( addedconsssize < *naddedconss )
      return;

   /* fill the array */
   for( cons = 0; cons < *naddedconss; cons++ )
   {
      addedconss[cons] = node->conssetchg->addedconss[cons];
   }

   return;
}

/** returns the number of added constraints to the given node */
int SCIPnodeGetNAddedConss(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   if( node->conssetchg == NULL )
      return 0;
   else
      return node->conssetchg->naddedconss;
}

/** adds the given bound change to the list of pending bound changes */
static
SCIP_RETCODE treeAddPendingBdchg(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node to add bound change to */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change, or NULL */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             probingchange       /**< is the bound change a temporary setting due to probing? */
   )
{
   assert(tree != NULL);

   /* make sure that enough memory is allocated for the pendingbdchgs array */
   SCIP_CALL( treeEnsurePendingbdchgsMem(tree, set, tree->npendingbdchgs+1) );

   /* capture the variable */
   SCIPvarCapture(var);

   /* add the bound change to the pending list */
   tree->pendingbdchgs[tree->npendingbdchgs].node = node;
   tree->pendingbdchgs[tree->npendingbdchgs].var = var;
   tree->pendingbdchgs[tree->npendingbdchgs].newbound = newbound;
   tree->pendingbdchgs[tree->npendingbdchgs].boundtype = boundtype;
   tree->pendingbdchgs[tree->npendingbdchgs].infercons = infercons;
   tree->pendingbdchgs[tree->npendingbdchgs].inferprop = inferprop;
   tree->pendingbdchgs[tree->npendingbdchgs].inferinfo = inferinfo;
   tree->pendingbdchgs[tree->npendingbdchgs].probingchange = probingchange;
   tree->npendingbdchgs++;

   /* check global pending boundchanges against debug solution */
   if( node->depth == 0 )
   {
#ifndef NDEBUG
      SCIP_Real bound = newbound;

      /* get bound adjusted for integrality(, this should already be done) */
      SCIPvarAdjustBd(var, set, boundtype, &bound);

      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
	 /* check that the bound is feasible */
	 if( bound > SCIPvarGetUbGlobal(var) )
	 {
	    /* due to numerics we only want to be feasible in feasibility tolerance */
	    assert(SCIPsetIsFeasLE(set, bound, SCIPvarGetUbGlobal(var)));
	    bound = SCIPvarGetUbGlobal(var);
	 }
      }
      else
      {
	 assert(boundtype == SCIP_BOUNDTYPE_UPPER);

	 /* check that the bound is feasible */
	 if( bound < SCIPvarGetLbGlobal(var) )
	 {
	    /* due to numerics we only want to be feasible in feasibility tolerance */
	    assert(SCIPsetIsFeasGE(set, bound, SCIPvarGetLbGlobal(var)));
	    bound = SCIPvarGetLbGlobal(var);
	 }
      }
      /* check that the given bound was already adjusted for integrality */
      assert(SCIPsetIsEQ(set, newbound, bound));
#endif
      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
	 /* check bound on debugging solution */
	 SCIP_CALL( SCIPdebugCheckLbGlobal(set->scip, var, newbound) ); /*lint !e506 !e774*/
      }
      else
      {
	 assert(boundtype == SCIP_BOUNDTYPE_UPPER);

	 /* check bound on debugging solution */
	 SCIP_CALL( SCIPdebugCheckUbGlobal(set->scip, var, newbound) ); /*lint !e506 !e774*/
      }
   }

   return SCIP_OKAY;
}

/** adds bound change with inference information to focus node, child of focus node, or probing node;
 *  if possible, adjusts bound to integral value;
 *  at most one of infercons and inferprop may be non-NULL
 */
SCIP_RETCODE SCIPnodeAddBoundinfer(
   SCIP_NODE*            node,               /**< node to add bound change to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change, or NULL */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             probingchange       /**< is the bound change a temporary setting due to probing? */
   )
{
   SCIP_VAR* infervar;
   SCIP_BOUNDTYPE inferboundtype;
   SCIP_Real oldlb;
   SCIP_Real oldub;
   SCIP_Real oldbound;

   assert(node != NULL);
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_FOCUSNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_PROBINGNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_CHILD
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_REFOCUSNODE
      || node->depth == 0);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->effectiverootdepth >= 0);
   assert(tree->root != NULL);
   assert(var != NULL);
   assert(node->active || (infercons == NULL && inferprop == NULL));
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_PROBINGNODE || !probingchange);

   SCIPsetDebugMsg(set, "adding boundchange at node %llu at depth %u to variable <%s>: old bounds=[%g,%g], new %s bound: %g (infer%s=<%s>, inferinfo=%d)\n",
      node->number, node->depth, SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", newbound, infercons != NULL ? "cons" : "prop",
      infercons != NULL ? SCIPconsGetName(infercons) : (inferprop != NULL ? SCIPpropGetName(inferprop) : "-"), inferinfo);

   /* remember variable as inference variable, and get corresponding active variable, bound and bound type */
   infervar = var;
   inferboundtype = boundtype;

   SCIP_CALL( SCIPvarGetProbvarBound(&var, &newbound, &boundtype) );

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPerrorMessage("cannot change bounds of multi-aggregated variable <%s>\n", SCIPvarGetName(var));
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   if( (int) node->depth <= tree->effectiverootdepth )
   {
      oldlb = SCIPvarGetLbGlobal(var);
      oldub = SCIPvarGetUbGlobal(var);
   }
   else
   {
      oldlb = SCIPvarGetLbLocal(var);
      oldub = SCIPvarGetUbLocal(var);
   }
   assert(SCIPsetIsLE(set, oldlb, oldub));

   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      /* adjust lower bound w.r.t. to integrality */
      SCIPvarAdjustLb(var, set, &newbound);
      assert(SCIPsetIsGT(set, newbound, oldlb));
      assert(SCIPsetIsFeasLE(set, newbound, oldub));
      oldbound = oldlb;
      newbound = MIN(newbound, oldub);

      if ( set->stage == SCIP_STAGE_SOLVING && SCIPsetIsInfinity(set, newbound) )
      {
         SCIPerrorMessage("cannot change lower bound of variable <%s> to infinity.\n", SCIPvarGetName(var));
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);

      /* adjust the new upper bound */
      SCIPvarAdjustUb(var, set, &newbound);
      assert(SCIPsetIsLT(set, newbound, oldub));
      assert(SCIPsetIsFeasGE(set, newbound, oldlb));
      oldbound = oldub;
      newbound = MAX(newbound, oldlb);

      if ( set->stage == SCIP_STAGE_SOLVING && SCIPsetIsInfinity(set, -newbound) )
      {
         SCIPerrorMessage("cannot change upper bound of variable <%s> to minus infinity.\n", SCIPvarGetName(var));
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   SCIPsetDebugMsg(set, " -> transformed to active variable <%s>: old bounds=[%g,%g], new %s bound: %g, obj: %g\n",
      SCIPvarGetName(var), oldlb, oldub, boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", newbound,
      SCIPvarGetObj(var));

   /* if the bound change takes place at an active node but is conflicting with the current local bounds,
    * we cannot apply it immediately because this would introduce inconsistencies to the bound change data structures
    * in the tree and to the bound change information data in the variable;
    * instead we have to remember the bound change as a pending bound change and mark the affected nodes on the active
    * path to be infeasible
    */
   if( node->active )
   {
      int conflictingdepth;

      conflictingdepth = SCIPvarGetConflictingBdchgDepth(var, set, boundtype, newbound);

      if( conflictingdepth >= 0 )
      {
         /* 0 would mean the bound change conflicts with a global bound */
         assert(conflictingdepth > 0);
         assert(conflictingdepth < tree->pathlen);

         SCIPsetDebugMsg(set, " -> bound change <%s> %s %g violates current local bounds [%g,%g] since depth %d: remember for later application\n",
            SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", newbound,
            SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), conflictingdepth);

         /* remember the pending bound change */
         SCIP_CALL( treeAddPendingBdchg(tree, set, node, var, newbound, boundtype, infercons, inferprop, inferinfo, 
               probingchange) );

         /* mark the node with the conflicting bound change to be cut off */
         SCIP_CALL( SCIPnodeCutoff(tree->path[conflictingdepth], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );

         return SCIP_OKAY;
      }
   }

   SCIPstatIncrement(stat, set, nboundchgs);

   /* if we are in probing mode we have to additionally count the bound changes for the probing statistic */
   if( tree->probingroot != NULL )
      SCIPstatIncrement(stat, set, nprobboundchgs);

   /* if the node is the root node: change local and global bound immediately */
   if( SCIPnodeGetDepth(node) <= tree->effectiverootdepth )
   {
      assert(node->active || tree->focusnode == NULL );
      assert(SCIPnodeGetType(node) != SCIP_NODETYPE_PROBINGNODE);
      assert(!probingchange);

      SCIPsetDebugMsg(set, " -> bound change in root node: perform global bound change\n");
      SCIP_CALL( SCIPvarChgBdGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newbound, boundtype) );

      if( set->stage == SCIP_STAGE_SOLVING )
      {
         /* the root should be repropagated due to the bound change */
         SCIPnodePropagateAgain(tree->root, set, stat, tree);
         SCIPsetDebugMsg(set, "marked root node to be repropagated due to global bound change <%s>:[%g,%g] -> [%g,%g] found in depth %u\n",
            SCIPvarGetName(var), oldlb, oldub, boundtype == SCIP_BOUNDTYPE_LOWER ? newbound : oldlb,
            boundtype == SCIP_BOUNDTYPE_LOWER ? oldub : newbound, node->depth);
      }

      return SCIP_OKAY;
   }

   /* if the node is a child, or the bound is a temporary probing bound
    *  - the bound change is a branching decision
    *  - the child's lower bound can be updated due to the changed pseudo solution
    * otherwise:
    *  - the bound change is an inference
    */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD || probingchange )
   {
      SCIP_Real newpseudoobjval;
      SCIP_Real lpsolval;

      assert(!node->active || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);

      /* get the solution value of variable in last solved LP on the active path:
       *  - if the LP was solved at the current node, the LP values of the columns are valid
       *  - if the last solved LP was the one in the current lpstatefork, the LP value in the columns are still valid
       *  - otherwise, the LP values are invalid
       */
      if( SCIPtreeHasCurrentNodeLP(tree)
         || (tree->focuslpstateforklpcount == stat->lpcount && SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN) )
      {
         lpsolval = SCIPvarGetLPSol(var);
      }
      else
         lpsolval = SCIP_INVALID;

      /* remember the bound change as branching decision (infervar/infercons/inferprop are not important: use NULL) */
      SCIP_CALL( SCIPdomchgAddBoundchg(&node->domchg, blkmem, set, var, newbound, boundtype, SCIP_BOUNDCHGTYPE_BRANCHING, 
            lpsolval, NULL, NULL, NULL, 0, inferboundtype) );

      /* update the child's lower bound */
      if( set->misc_exactsolve )
         newpseudoobjval = SCIPlpGetModifiedProvedPseudoObjval(lp, set, var, oldbound, newbound, boundtype);
      else
         newpseudoobjval = SCIPlpGetModifiedPseudoObjval(lp, set, transprob, var, oldbound, newbound, boundtype);
      SCIPnodeUpdateLowerbound(node, stat, set, tree, transprob, origprob, newpseudoobjval);
   }
   else
   {
      /* check the infered bound change on the debugging solution */
      SCIP_CALL( SCIPdebugCheckInference(blkmem, set, node, var, newbound, boundtype) ); /*lint !e506 !e774*/

      /* remember the bound change as inference (lpsolval is not important: use 0.0) */
      SCIP_CALL( SCIPdomchgAddBoundchg(&node->domchg, blkmem, set, var, newbound, boundtype,
            infercons != NULL ? SCIP_BOUNDCHGTYPE_CONSINFER : SCIP_BOUNDCHGTYPE_PROPINFER, 
            0.0, infervar, infercons, inferprop, inferinfo, inferboundtype) );
   }

   assert(node->domchg != NULL);
   assert(node->domchg->domchgdyn.domchgtype == SCIP_DOMCHGTYPE_DYNAMIC); /*lint !e641*/
   assert(node->domchg->domchgdyn.boundchgs != NULL);
   assert(node->domchg->domchgdyn.nboundchgs > 0);
   assert(node->domchg->domchgdyn.boundchgs[node->domchg->domchgdyn.nboundchgs-1].var == var);
   assert(node->domchg->domchgdyn.boundchgs[node->domchg->domchgdyn.nboundchgs-1].newbound == newbound); /*lint !e777*/

   /* if node is active, apply the bound change immediately */
   if( node->active )
   {
      SCIP_Bool cutoff;

      /**@todo if the node is active, it currently must either be the effective root (see above) or the current node;
       *       if a bound change to an intermediate active node should be added, we must make sure, the bound change
       *       information array of the variable stays sorted (new info must be sorted in instead of putting it to
       *       the end of the array), and we should identify now redundant bound changes that are applied at a
       *       later node on the active path
       */
      assert(SCIPtreeGetCurrentNode(tree) == node); 
      SCIP_CALL( SCIPboundchgApply(&node->domchg->domchgdyn.boundchgs[node->domchg->domchgdyn.nboundchgs-1],
            blkmem, set, stat, lp, branchcand, eventqueue, node->depth, node->domchg->domchgdyn.nboundchgs-1, &cutoff) );
      assert(node->domchg->domchgdyn.boundchgs[node->domchg->domchgdyn.nboundchgs-1].var == var);
      assert(!cutoff);
   }

   return SCIP_OKAY;
}

/** adds bound change to focus node, or child of focus node, or probing node;
 *  if possible, adjusts bound to integral value
 */
SCIP_RETCODE SCIPnodeAddBoundchg(
   SCIP_NODE*            node,               /**< node to add bound change to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   SCIP_Bool             probingchange       /**< is the bound change a temporary setting due to probing? */
   )
{
   SCIP_CALL( SCIPnodeAddBoundinfer(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
         cliquetable, var, newbound, boundtype, NULL, NULL, 0, probingchange) );

   return SCIP_OKAY;
}

/** adds hole with inference information to focus node, child of focus node, or probing node;
 *  if possible, adjusts bound to integral value;
 *  at most one of infercons and inferprop may be non-NULL
 */
SCIP_RETCODE SCIPnodeAddHoleinfer(
   SCIP_NODE*            node,               /**< node to add bound change to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             left,               /**< left bound of open interval defining the hole (left,right) */
   SCIP_Real             right,              /**< right bound of open interval defining the hole (left,right) */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change, or NULL */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             probingchange,      /**< is the bound change a temporary setting due to probing? */
   SCIP_Bool*            added               /**< pointer to store whether the hole was added, or NULL */
   )
{
#if 0
   SCIP_VAR* infervar;
#endif

   assert(node != NULL);
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_FOCUSNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_PROBINGNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_CHILD
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_REFOCUSNODE
      || node->depth == 0);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->effectiverootdepth >= 0);
   assert(tree->root != NULL);
   assert(var != NULL);
   assert(node->active || (infercons == NULL && inferprop == NULL));
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_PROBINGNODE || !probingchange);

   /* the interval should not be empty */
   assert(SCIPsetIsLT(set, left, right));

#ifndef NDEBUG
   {
      SCIP_Real adjustedleft;
      SCIP_Real adjustedright;

      adjustedleft = left;
      adjustedright = right;

      SCIPvarAdjustUb(var, set, &adjustedleft);
      SCIPvarAdjustLb(var, set, &adjustedright);

      assert(SCIPsetIsEQ(set, left, adjustedleft));
      assert(SCIPsetIsEQ(set, right, adjustedright));
   }
#endif

   /* the hole should lay within the lower and upper bounds */
   assert(SCIPsetIsGE(set, left, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsLE(set, right, SCIPvarGetUbLocal(var)));

   SCIPsetDebugMsg(set, "adding hole (%g,%g) at node at depth %u to variable <%s>: bounds=[%g,%g], (infer%s=<%s>, inferinfo=%d)\n",
      left, right, node->depth, SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), infercons != NULL ? "cons" : "prop",
      infercons != NULL ? SCIPconsGetName(infercons) : (inferprop != NULL ? SCIPpropGetName(inferprop) : "-"), inferinfo);

#if 0
   /* remember variable as inference variable, and get corresponding active variable, bound and bound type */
   infervar = var;
#endif
   SCIP_CALL( SCIPvarGetProbvarHole(&var, &left, &right) );

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPerrorMessage("cannot change bounds of multi-aggregated variable <%s>\n", SCIPvarGetName(var));
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   SCIPsetDebugMsg(set, " -> transformed to active variable <%s>: hole (%g,%g), obj: %g\n", SCIPvarGetName(var), left, right, SCIPvarGetObj(var));

   stat->nholechgs++;

   /* if we are in probing mode we have to additionally count the bound changes for the probing statistic */
   if( tree->probingroot != NULL )
      stat->nprobholechgs++;

   /* if the node is the root node: change local and global bound immediately */
   if( SCIPnodeGetDepth(node) <= tree->effectiverootdepth )
   {
      assert(node->active || tree->focusnode == NULL );
      assert(SCIPnodeGetType(node) != SCIP_NODETYPE_PROBINGNODE);
      assert(!probingchange);

      SCIPsetDebugMsg(set, " -> hole added in root node: perform global domain change\n");
      SCIP_CALL( SCIPvarAddHoleGlobal(var, blkmem, set, stat, eventqueue, left, right, added) );

      if( set->stage == SCIP_STAGE_SOLVING && (*added) )
      {
         /* the root should be repropagated due to the bound change */
         SCIPnodePropagateAgain(tree->root, set, stat, tree);
         SCIPsetDebugMsg(set, "marked root node to be repropagated due to global added hole <%s>: (%g,%g) found in depth %u\n",
            SCIPvarGetName(var), left, right, node->depth);
      }

      return SCIP_OKAY;
   }

   /**@todo add adding of local domain holes */

   (*added) = FALSE;
   SCIPerrorMessage("WARNING: currently domain holes can only be handled globally!\n");

   stat->nholechgs--;

   /* if we are in probing mode we have to additionally count the bound changes for the probing statistic */
   if( tree->probingroot != NULL )
      stat->nprobholechgs--;

   return SCIP_OKAY;
}

/** adds hole change to focus node, or child of focus node */
SCIP_RETCODE SCIPnodeAddHolechg(
   SCIP_NODE*            node,               /**< node to add bound change to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             left,               /**< left bound of open interval defining the hole (left,right) */
   SCIP_Real             right,              /**< right bound of open interval defining the hole (left,right) */
   SCIP_Bool             probingchange,      /**< is the bound change a temporary setting due to probing? */
   SCIP_Bool*            added               /**< pointer to store whether the hole was added, or NULL */
   )
{
   assert(node != NULL);
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_FOCUSNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_PROBINGNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_CHILD);
   assert(blkmem != NULL);

   SCIPsetDebugMsg(set, "adding hole (%g,%g) at node at depth %u of variable <%s>\n",
      left, right, node->depth, SCIPvarGetName(var));

   SCIP_CALL( SCIPnodeAddHoleinfer(node, blkmem, set, stat, tree, eventqueue, var, left, right,
         NULL, NULL, 0, probingchange, added) );

   /**@todo apply hole change on active nodes and issue event */

   return SCIP_OKAY;
}

/** applies the pending bound changes */
static
SCIP_RETCODE treeApplyPendingBdchgs(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   SCIP_VAR* var;
   int npendingbdchgs;
   int conflictdepth;
   int i;

   assert(tree != NULL);

   npendingbdchgs = tree->npendingbdchgs;
   for( i = 0; i < npendingbdchgs; ++i )
   {
      var = tree->pendingbdchgs[i].var;
      assert(SCIPnodeGetDepth(tree->pendingbdchgs[i].node) < tree->cutoffdepth);

      conflictdepth = SCIPvarGetConflictingBdchgDepth(var, set, tree->pendingbdchgs[i].boundtype,
         tree->pendingbdchgs[i].newbound);

      /* It can happen, that a pending bound change conflicts with the global bounds, because when it was collected, it
       * just conflicted with the local bounds, but a conflicting global bound change was applied afterwards. In this
       * case, we can cut off the node where the pending bound change should be applied.
       */
      if( conflictdepth == 0 )
      {
         SCIP_CALL( SCIPnodeCutoff(tree->pendingbdchgs[i].node, set, stat, tree, transprob, origprob, reopt, lp, blkmem) );

         if( ((int) tree->pendingbdchgs[i].node->depth) <= tree->effectiverootdepth )
            return SCIP_OKAY;
         else
            continue;
      }

      assert(conflictdepth == -1);

      SCIPsetDebugMsg(set, "applying pending bound change <%s>[%g,%g] %s %g\n", SCIPvarGetName(var),
         SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
         tree->pendingbdchgs[i].boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
         tree->pendingbdchgs[i].newbound);

      /* ignore bounds that are now redundant (for example, multiple entries in the pendingbdchgs for the same
       * variable)
       */
      if( tree->pendingbdchgs[i].boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_Real lb;

         lb = SCIPvarGetLbLocal(var);
         if( !SCIPsetIsGT(set, tree->pendingbdchgs[i].newbound, lb) )
         {
            /* release the variable */
            SCIP_CALL( SCIPvarRelease(&var, blkmem, set, eventqueue, lp) );
            continue;
         }
      }
      else
      {
         SCIP_Real ub;

         assert(tree->pendingbdchgs[i].boundtype == SCIP_BOUNDTYPE_UPPER);
         ub = SCIPvarGetUbLocal(var);
         if( !SCIPsetIsLT(set, tree->pendingbdchgs[i].newbound, ub) )
         {
            /* release the variable */
            SCIP_CALL( SCIPvarRelease(&var, blkmem, set, eventqueue, lp) );
            continue;
         }
      }

      SCIP_CALL( SCIPnodeAddBoundinfer(tree->pendingbdchgs[i].node, blkmem, set, stat, transprob, origprob, tree, reopt,
            lp, branchcand, eventqueue, cliquetable, var, tree->pendingbdchgs[i].newbound, tree->pendingbdchgs[i].boundtype,
            tree->pendingbdchgs[i].infercons, tree->pendingbdchgs[i].inferprop, tree->pendingbdchgs[i].inferinfo,
            tree->pendingbdchgs[i].probingchange) );
      assert(tree->npendingbdchgs == npendingbdchgs); /* this time, the bound change can be applied! */

      /* release the variable */
      SCIP_CALL( SCIPvarRelease(&var, blkmem, set, eventqueue, lp) );
   }
   tree->npendingbdchgs = 0;

   return SCIP_OKAY;
}

/** if given value is larger than the node's lower bound, sets the node's lower bound to the new value */
void SCIPnodeUpdateLowerbound(
   SCIP_NODE*            node,               /**< node to update lower bound for */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   assert(node != NULL);
   assert(stat != NULL);

   if( newbound > node->lowerbound )
   {
      SCIP_Real oldbound;

      oldbound = node->lowerbound;
      node->lowerbound = newbound;
      node->estimate = MAX(node->estimate, newbound);
      if( node->depth == 0 )
      {
         stat->rootlowerbound = newbound;
         if( set->misc_calcintegral )
            SCIPstatUpdatePrimalDualIntegral(stat, set, transprob, origprob, SCIPsetInfinity(set), newbound);
      }
      else if( set->misc_calcintegral && SCIPsetIsEQ(set, oldbound, stat->lastlowerbound) )
      {
         SCIP_Real lowerbound;
         lowerbound = SCIPtreeGetLowerbound(tree, set);
         assert(newbound >= lowerbound);

         /* updating the primal integral is only necessary if dual bound has increased since last evaluation */
         if( lowerbound > stat->lastlowerbound )
            SCIPstatUpdatePrimalDualIntegral(stat, set, transprob, origprob, SCIPsetInfinity(set), lowerbound);
      }
   }
}

/** updates lower bound of node using lower bound of LP */
SCIP_RETCODE SCIPnodeUpdateLowerboundLP(
   SCIP_NODE*            node,               /**< node to set lower bound for */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp                  /**< LP data */
   )
{
   SCIP_Real lpobjval;

   assert(set != NULL);
   assert(lp->flushed);

   /* in case of iteration or time limit, the LP value may not be a valid dual bound */
   /* @todo check for dual feasibility of LP solution and use sub-optimal solution if they are dual feasible */
   if( lp->lpsolstat == SCIP_LPSOLSTAT_ITERLIMIT || lp->lpsolstat == SCIP_LPSOLSTAT_TIMELIMIT )
      return SCIP_OKAY;

   if( set->misc_exactsolve )
   {
      SCIP_CALL( SCIPlpGetProvedLowerbound(lp, set, &lpobjval) );
   }
   else
      lpobjval = SCIPlpGetObjval(lp, set, transprob);

   SCIPnodeUpdateLowerbound(node, stat, set, tree, transprob, origprob, lpobjval);

   return SCIP_OKAY;
}


/** change the node selection priority of the given child */
void SCIPchildChgNodeselPrio(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            child,              /**< child to update the node selection priority */
   SCIP_Real             priority            /**< node selection priority value */
   )
{
   int pos;

   assert( SCIPnodeGetType(child) == SCIP_NODETYPE_CHILD );

   pos = child->data.child.arraypos;
   assert( pos >= 0 );

   tree->childrenprio[pos] = priority;
}


/** sets the node's estimated bound to the new value */
void SCIPnodeSetEstimate(
   SCIP_NODE*            node,               /**< node to update lower bound for */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newestimate         /**< new estimated bound for the node */
   )
{
   assert(node != NULL);
   assert(set != NULL);
   assert(SCIPsetIsRelGE(set, newestimate, node->lowerbound));

   node->estimate = newestimate;
}

/** propagates implications of binary fixings at the given node triggered by the implication graph and the clique table */
SCIP_RETCODE SCIPnodePropagateImplics(
   SCIP_NODE*            node,               /**< node to propagate implications on */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   int nboundchgs;
   int i;

   assert(node != NULL);
   assert(SCIPnodeIsActive(node));
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_FOCUSNODE
      || SCIPnodeGetType(node) == SCIP_NODETYPE_REFOCUSNODE
      || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
   assert(cutoff != NULL);

   SCIPsetDebugMsg(set, "implication graph propagation of node #%" SCIP_LONGINT_FORMAT " in depth %d\n",
      SCIPnodeGetNumber(node), SCIPnodeGetDepth(node));

   *cutoff = FALSE;

   /* propagate all fixings of binary variables performed at this node */
   nboundchgs = SCIPdomchgGetNBoundchgs(node->domchg);
   for( i = 0; i < nboundchgs && !(*cutoff); ++i )
   {
      SCIP_BOUNDCHG* boundchg;
      SCIP_VAR* var;

      boundchg = SCIPdomchgGetBoundchg(node->domchg, i);

      /* ignore redundant bound changes */
      if( SCIPboundchgIsRedundant(boundchg) )
         continue;

      var = SCIPboundchgGetVar(boundchg);
      if( SCIPvarIsBinary(var) )
      {
         SCIP_Bool varfixing;
         int nimpls;
         SCIP_VAR** implvars;
         SCIP_BOUNDTYPE* impltypes;
         SCIP_Real* implbounds;
         SCIP_CLIQUE** cliques;
         int ncliques;
         int j;

         varfixing = (SCIPboundchgGetBoundtype(boundchg) == SCIP_BOUNDTYPE_LOWER);
         nimpls = SCIPvarGetNImpls(var, varfixing);
         implvars = SCIPvarGetImplVars(var, varfixing);
         impltypes = SCIPvarGetImplTypes(var, varfixing);
         implbounds = SCIPvarGetImplBounds(var, varfixing);

         /* apply implications */
         for( j = 0; j < nimpls; ++j )
         {
            SCIP_Real lb;
            SCIP_Real ub;

            /* @note should this be checked here (because SCIPnodeAddBoundinfer fails for multi-aggregated variables)
             *       or should SCIPnodeAddBoundinfer() just return for multi-aggregated variables?
             */
            if( SCIPvarGetStatus(SCIPvarGetProbvar(implvars[j])) == SCIP_VARSTATUS_MULTAGGR )
               continue;

            /* check for infeasibility */
            lb = SCIPvarGetLbLocal(implvars[j]);
            ub = SCIPvarGetUbLocal(implvars[j]);
            if( impltypes[j] == SCIP_BOUNDTYPE_LOWER )
            {
               if( SCIPsetIsFeasGT(set, implbounds[j], ub) )
               {
                  *cutoff = TRUE;
                  return SCIP_OKAY;
               }
               if( SCIPsetIsFeasLE(set, implbounds[j], lb) )
                  continue;
            }
            else
            {
               if( SCIPsetIsFeasLT(set, implbounds[j], lb) )
               {
                  *cutoff = TRUE;
                  return SCIP_OKAY;
               }
               if( SCIPsetIsFeasGE(set, implbounds[j], ub) )
                  continue;
            }

            /* @note the implication might affect a fixed variable (after resolving (multi-)aggregations);
             *       normally, the implication should have been deleted in that case, but this is only possible
             *       if the implied variable has the reverse implication stored as a variable bound;
             *       due to numerics, the variable bound may not be present and so the implication is not deleted
             */
            if( SCIPvarGetStatus(SCIPvarGetProbvar(implvars[j])) == SCIP_VARSTATUS_FIXED )
               continue;

            /* apply the implication */
            SCIP_CALL( SCIPnodeAddBoundinfer(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
                  eventqueue, cliquetable, implvars[j], implbounds[j], impltypes[j], NULL, NULL, 0, FALSE) );
         }

         /* apply cliques */
         ncliques = SCIPvarGetNCliques(var, varfixing);
         cliques = SCIPvarGetCliques(var, varfixing);
         for( j = 0; j < ncliques; ++j )
         {
            SCIP_VAR** vars;
            SCIP_Bool* values;
            int nvars;
            int k;

            nvars = SCIPcliqueGetNVars(cliques[j]);
            vars = SCIPcliqueGetVars(cliques[j]);
            values = SCIPcliqueGetValues(cliques[j]);
            for( k = 0; k < nvars; ++k )
            {
               SCIP_Real lb;
               SCIP_Real ub;

               assert(SCIPvarIsBinary(vars[k]));

               if( SCIPvarGetStatus(vars[k]) == SCIP_VARSTATUS_MULTAGGR )
                  continue;

               if( vars[k] == var && values[k] == varfixing )
                  continue;

               /* check for infeasibility */
               lb = SCIPvarGetLbLocal(vars[k]);
               ub = SCIPvarGetUbLocal(vars[k]);
               if( values[k] == FALSE )
               {
                  if( ub < 0.5 )
                  {
                     *cutoff = TRUE;
                     return SCIP_OKAY;
                  }
                  if( lb > 0.5 )
                     continue;
               }
               else
               {
                  if( lb > 0.5 )
                  {
                     *cutoff = TRUE;
                     return SCIP_OKAY;
                  }
                  if( ub < 0.5 )
                     continue;
               }

               if( SCIPvarGetStatus(SCIPvarGetProbvar(vars[k])) == SCIP_VARSTATUS_FIXED )
                  continue;

               /* apply the clique implication */
               SCIP_CALL( SCIPnodeAddBoundinfer(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
                     eventqueue, cliquetable, vars[k], (SCIP_Real)(!values[k]), values[k] ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER,
                     NULL, NULL, 0, FALSE) );
            }
         }
      }
   }

   return SCIP_OKAY;
}




/*
 * Path Switching
 */

/** updates the LP sizes of the active path starting at the given depth */
static
SCIP_RETCODE treeUpdatePathLPSize(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   startdepth          /**< depth to start counting */
   )
{
   SCIP_NODE* node;
   int ncols;
   int nrows;
   int i;

   assert(tree != NULL);
   assert(startdepth >= 0);
   assert(startdepth <= tree->pathlen);

   if( startdepth == 0 )
   {
      ncols = 0;
      nrows = 0;
   }
   else
   {
      ncols = tree->pathnlpcols[startdepth-1];
      nrows = tree->pathnlprows[startdepth-1];
   }

   for( i = startdepth; i < tree->pathlen; ++i )
   {
      node = tree->path[i];
      assert(node != NULL);
      assert(node->active);
      assert((int)(node->depth) == i);

      switch( SCIPnodeGetType(node) )
      {
      case SCIP_NODETYPE_FOCUSNODE:
         assert(i == tree->pathlen-1 || SCIPtreeProbing(tree));
         break;
      case SCIP_NODETYPE_PROBINGNODE:
         assert(SCIPtreeProbing(tree));
         assert(i >= 1);
         assert(SCIPnodeGetType(tree->path[i-1]) == SCIP_NODETYPE_FOCUSNODE
            || (ncols == node->data.probingnode->ninitialcols && nrows == node->data.probingnode->ninitialrows));
         assert(ncols <= node->data.probingnode->ncols || !tree->focuslpconstructed);
         assert(nrows <= node->data.probingnode->nrows || !tree->focuslpconstructed);
         if( i < tree->pathlen-1 )
         {
            ncols = node->data.probingnode->ncols;
            nrows = node->data.probingnode->nrows;
         }
         else
         {
            /* for the current probing node, the initial LP size is stored in the path */
            ncols = node->data.probingnode->ninitialcols;
            nrows = node->data.probingnode->ninitialrows;
         }
         break;
      case SCIP_NODETYPE_SIBLING:
         SCIPerrorMessage("sibling cannot be in the active path\n");
         SCIPABORT();
         return SCIP_INVALIDDATA;  /*lint !e527*/
      case SCIP_NODETYPE_CHILD:
         SCIPerrorMessage("child cannot be in the active path\n");
         SCIPABORT();
         return SCIP_INVALIDDATA;  /*lint !e527*/
      case SCIP_NODETYPE_LEAF:
         SCIPerrorMessage("leaf cannot be in the active path\n");
         SCIPABORT();
         return SCIP_INVALIDDATA;  /*lint !e527*/
      case SCIP_NODETYPE_DEADEND:
         SCIPerrorMessage("dead-end cannot be in the active path\n");
         SCIPABORT();
         return SCIP_INVALIDDATA;  /*lint !e527*/
      case SCIP_NODETYPE_JUNCTION:
         break;
      case SCIP_NODETYPE_PSEUDOFORK:
         assert(node->data.pseudofork != NULL);
         ncols += node->data.pseudofork->naddedcols;
         nrows += node->data.pseudofork->naddedrows;
         break;
      case SCIP_NODETYPE_FORK:
         assert(node->data.fork != NULL);
         ncols += node->data.fork->naddedcols;
         nrows += node->data.fork->naddedrows;
         break;
      case SCIP_NODETYPE_SUBROOT:
         assert(node->data.subroot != NULL);
         ncols = node->data.subroot->ncols;
         nrows = node->data.subroot->nrows;
         break;
      case SCIP_NODETYPE_REFOCUSNODE:
         SCIPerrorMessage("node cannot be of type REFOCUSNODE at this point\n");
         SCIPABORT();
         return SCIP_INVALIDDATA;  /*lint !e527*/
      default:
         SCIPerrorMessage("unknown node type %d\n", SCIPnodeGetType(node));
         SCIPABORT();
         return SCIP_INVALIDDATA;  /*lint !e527*/
      }
      tree->pathnlpcols[i] = ncols;
      tree->pathnlprows[i] = nrows;
   }
   return SCIP_OKAY;
}

/** finds the common fork node, the new LP state defining fork, and the new focus subroot, if the path is switched to
 *  the given node
 */
static
void treeFindSwitchForks(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            node,               /**< new focus node, or NULL */
   SCIP_NODE**           commonfork,         /**< pointer to store common fork node of old and new focus node */
   SCIP_NODE**           newlpfork,          /**< pointer to store the new LP defining fork node */
   SCIP_NODE**           newlpstatefork,     /**< pointer to store the new LP state defining fork node */
   SCIP_NODE**           newsubroot,         /**< pointer to store the new subroot node */
   SCIP_Bool*            cutoff              /**< pointer to store whether the given node can be cut off and no path switching
                                              *   should be performed */
   )
{
   SCIP_NODE* fork;
   SCIP_NODE* lpfork;
   SCIP_NODE* lpstatefork;
   SCIP_NODE* subroot;

   assert(tree != NULL);
   assert(tree->root != NULL);
   assert((tree->focusnode == NULL) == !tree->root->active);
   assert(tree->focuslpfork == NULL || tree->focusnode != NULL);
   assert(tree->focuslpfork == NULL || tree->focuslpfork->depth < tree->focusnode->depth);
   assert(tree->focuslpstatefork == NULL || tree->focuslpfork != NULL);
   assert(tree->focuslpstatefork == NULL || tree->focuslpstatefork->depth <= tree->focuslpfork->depth);
   assert(tree->focussubroot == NULL || tree->focuslpstatefork != NULL);
   assert(tree->focussubroot == NULL || tree->focussubroot->depth <= tree->focuslpstatefork->depth);
   assert(tree->cutoffdepth >= 0);
   assert(tree->cutoffdepth == INT_MAX || tree->cutoffdepth < tree->pathlen);
   assert(tree->cutoffdepth == INT_MAX || tree->path[tree->cutoffdepth]->cutoff);
   assert(tree->repropdepth >= 0);
   assert(tree->repropdepth == INT_MAX || tree->repropdepth < tree->pathlen);
   assert(tree->repropdepth == INT_MAX || tree->path[tree->repropdepth]->reprop);
   assert(commonfork != NULL);
   assert(newlpfork != NULL);
   assert(newlpstatefork != NULL);
   assert(newsubroot != NULL);
   assert(cutoff != NULL);

   *commonfork = NULL;
   *newlpfork = NULL;
   *newlpstatefork = NULL;
   *newsubroot = NULL;
   *cutoff = FALSE;

   /* if the new focus node is NULL, there is no common fork node, and the new LP fork, LP state fork, and subroot
    * are NULL
    */
   if( node == NULL )
   {
      tree->cutoffdepth = INT_MAX;
      tree->repropdepth = INT_MAX;
      return;
   }

   /* check if the new node is marked to be cut off */
   if( node->cutoff )
   {
      *cutoff = TRUE;
      return;
   }

   /* if the old focus node is NULL, there is no common fork node, and we have to search the new LP fork, LP state fork
    * and subroot
    */
   if( tree->focusnode == NULL )
   {
      assert(!tree->root->active);
      assert(tree->pathlen == 0);
      assert(tree->cutoffdepth == INT_MAX);
      assert(tree->repropdepth == INT_MAX);

      lpfork = node;
      while( SCIPnodeGetType(lpfork) != SCIP_NODETYPE_PSEUDOFORK
         && SCIPnodeGetType(lpfork) != SCIP_NODETYPE_FORK && SCIPnodeGetType(lpfork) != SCIP_NODETYPE_SUBROOT )
      {
         lpfork = lpfork->parent;
         if( lpfork == NULL )
            return;
         if( lpfork->cutoff )
         {
            *cutoff = TRUE;
            return;
         }
      }
      *newlpfork = lpfork;

      lpstatefork = lpfork;
      while( SCIPnodeGetType(lpstatefork) != SCIP_NODETYPE_FORK && SCIPnodeGetType(lpstatefork) != SCIP_NODETYPE_SUBROOT )
      {
         lpstatefork = lpstatefork->parent;
         if( lpstatefork == NULL )
            return;
         if( lpstatefork->cutoff )
         {
            *cutoff = TRUE;
            return;
         }
      }
      *newlpstatefork = lpstatefork;

      subroot = lpstatefork;
      while( SCIPnodeGetType(subroot) != SCIP_NODETYPE_SUBROOT )
      {
         subroot = subroot->parent;
         if( subroot == NULL )
            return;
         if( subroot->cutoff )
         {
            *cutoff = TRUE;
            return;
         }
      }
      *newsubroot = subroot;

      fork = subroot;
      while( fork->parent != NULL )
      {
         fork = fork->parent;
         if( fork->cutoff )
         {
            *cutoff = TRUE;
            return;
         }
      }
      return;
   }

   /* find the common fork node, the new LP defining fork, the new LP state defining fork, and the new focus subroot */
   fork = node;
   lpfork = NULL;
   lpstatefork = NULL;
   subroot = NULL;
   assert(fork != NULL);

   while( !fork->active )
   {
      fork = fork->parent;
      assert(fork != NULL); /* because the root is active, there must be a common fork node */

      if( fork->cutoff )
      {
         *cutoff = TRUE;
         return;
      }
      if( lpfork == NULL
         && (SCIPnodeGetType(fork) == SCIP_NODETYPE_PSEUDOFORK
            || SCIPnodeGetType(fork) == SCIP_NODETYPE_FORK || SCIPnodeGetType(fork) == SCIP_NODETYPE_SUBROOT) )
         lpfork = fork;
      if( lpstatefork == NULL
         && (SCIPnodeGetType(fork) == SCIP_NODETYPE_FORK || SCIPnodeGetType(fork) == SCIP_NODETYPE_SUBROOT) )
         lpstatefork = fork;
      if( subroot == NULL && SCIPnodeGetType(fork) == SCIP_NODETYPE_SUBROOT )
         subroot = fork;
   }
   assert(lpfork == NULL || !lpfork->active || lpfork == fork);
   assert(lpstatefork == NULL || !lpstatefork->active || lpstatefork == fork);
   assert(subroot == NULL || !subroot->active || subroot == fork);
   SCIPdebugMessage("find switch forks: forkdepth=%u\n", fork->depth);

   /* if the common fork node is below the current cutoff depth, the cutoff node is an ancestor of the common fork
    * and thus an ancestor of the new focus node, s.t. the new node can also be cut off
    */
   assert((int)fork->depth != tree->cutoffdepth);
   if( (int)fork->depth > tree->cutoffdepth )
   {
#ifndef NDEBUG
      while( !fork->cutoff )
      {
         fork = fork->parent;
         assert(fork != NULL);
      }
      assert((int)fork->depth >= tree->cutoffdepth);
#endif
      *cutoff = TRUE;
      return;
   }
   tree->cutoffdepth = INT_MAX;

   /* if not already found, continue searching the LP defining fork; it cannot be deeper than the common fork */
   if( lpfork == NULL )
   {
      if( tree->focuslpfork != NULL && (int)(tree->focuslpfork->depth) > fork->depth )
      {
         /* focuslpfork is not on the same active path as the new node: we have to continue searching */
         lpfork = fork;
         while( lpfork != NULL
            && SCIPnodeGetType(lpfork) != SCIP_NODETYPE_PSEUDOFORK
            && SCIPnodeGetType(lpfork) != SCIP_NODETYPE_FORK
            && SCIPnodeGetType(lpfork) != SCIP_NODETYPE_SUBROOT )
         {
            assert(lpfork->active);
            lpfork = lpfork->parent;
         }
      }
      else
      {
         /* focuslpfork is on the same active path as the new node: old and new node have the same lpfork */
         lpfork = tree->focuslpfork;
      }
      assert(lpfork == NULL || (int)(lpfork->depth) <= fork->depth);
      assert(lpfork == NULL || lpfork->active);
   }
   assert(lpfork == NULL
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_PSEUDOFORK
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT);
   SCIPdebugMessage("find switch forks: lpforkdepth=%d\n", lpfork == NULL ? -1 : (int)(lpfork->depth));

   /* if not already found, continue searching the LP state defining fork; it cannot be deeper than the
    * LP defining fork and the common fork
    */
   if( lpstatefork == NULL )
   {
      if( tree->focuslpstatefork != NULL && (int)(tree->focuslpstatefork->depth) > fork->depth )
      {
         /* focuslpstatefork is not on the same active path as the new node: we have to continue searching */
         if( lpfork != NULL && lpfork->depth < fork->depth )
            lpstatefork = lpfork;
         else
            lpstatefork = fork;
         while( lpstatefork != NULL
            && SCIPnodeGetType(lpstatefork) != SCIP_NODETYPE_FORK
            && SCIPnodeGetType(lpstatefork) != SCIP_NODETYPE_SUBROOT )
         {
            assert(lpstatefork->active);
            lpstatefork = lpstatefork->parent;
         }
      }
      else
      {
         /* focuslpstatefork is on the same active path as the new node: old and new node have the same lpstatefork */
         lpstatefork = tree->focuslpstatefork;
      }
      assert(lpstatefork == NULL || (int)(lpstatefork->depth) <= fork->depth);
      assert(lpstatefork == NULL || lpstatefork->active);
   }
   assert(lpstatefork == NULL
      || SCIPnodeGetType(lpstatefork) == SCIP_NODETYPE_FORK
      || SCIPnodeGetType(lpstatefork) == SCIP_NODETYPE_SUBROOT);
   assert(lpstatefork == NULL || (lpfork != NULL && lpstatefork->depth <= lpfork->depth));
   SCIPdebugMessage("find switch forks: lpstateforkdepth=%d\n", lpstatefork == NULL ? -1 : (int)(lpstatefork->depth));

   /* if not already found, continue searching the subroot; it cannot be deeper than the LP defining fork, the
    * LP state fork and the common fork
    */
   if( subroot == NULL )
   {
      if( tree->focussubroot != NULL && (int)(tree->focussubroot->depth) > fork->depth )
      {
         /* focussubroot is not on the same active path as the new node: we have to continue searching */
         if( lpstatefork != NULL && lpstatefork->depth < fork->depth )
            subroot = lpstatefork;
         else if( lpfork != NULL && lpfork->depth < fork->depth )
            subroot = lpfork;
         else
            subroot = fork;
         while( subroot != NULL && SCIPnodeGetType(subroot) != SCIP_NODETYPE_SUBROOT )
         {
            assert(subroot->active);
            subroot = subroot->parent;
         }
      }
      else
         subroot = tree->focussubroot;
      assert(subroot == NULL || subroot->depth <= fork->depth);
      assert(subroot == NULL || subroot->active);
   }
   assert(subroot == NULL || SCIPnodeGetType(subroot) == SCIP_NODETYPE_SUBROOT);
   assert(subroot == NULL || (lpstatefork != NULL && subroot->depth <= lpstatefork->depth));
   SCIPdebugMessage("find switch forks: subrootdepth=%d\n", subroot == NULL ? -1 : (int)(subroot->depth));

   /* if a node prior to the common fork should be repropagated, we select the node to be repropagated as common
    * fork in order to undo all bound changes up to this node, repropagate the node, and redo the bound changes
    * afterwards
    */
   if( (int)fork->depth > tree->repropdepth )
   {
      fork = tree->path[tree->repropdepth];
      assert(fork->active);
      assert(fork->reprop);
   }

   *commonfork = fork;
   *newlpfork = lpfork;
   *newlpstatefork = lpstatefork;
   *newsubroot = subroot;

#ifndef NDEBUG
   while( fork != NULL )
   {
      assert(fork->active);
      assert(!fork->cutoff);
      assert(fork->parent == NULL || !fork->parent->reprop);
      fork = fork->parent;
   }
#endif
   tree->repropdepth = INT_MAX;
}

/** switches the active path to the new focus node, applies domain and constraint set changes */
static
SCIP_RETCODE treeSwitchPath(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_NODE*            fork,               /**< common fork node of old and new focus node, or NULL */
   SCIP_NODE*            focusnode,          /**< new focus node, or NULL */
   SCIP_Bool*            cutoff              /**< pointer to store whether the new focus node can be cut off */
   )
{
   int focusnodedepth;  /* depth of the new focus node, or -1 if focusnode == NULL */
   int forkdepth;       /* depth of the common subroot/fork/pseudofork/junction node, or -1 if no common fork exists */
   int i;

   assert(tree != NULL);
   assert(fork == NULL || (fork->active && !fork->cutoff));
   assert(fork == NULL || focusnode != NULL);
   assert(focusnode == NULL || (!focusnode->active && !focusnode->cutoff));
   assert(focusnode == NULL || SCIPnodeGetType(focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   SCIPsetDebugMsg(set, "switch path: old pathlen=%d\n", tree->pathlen);

   /* get the nodes' depths */
   focusnodedepth = (focusnode != NULL ? (int)focusnode->depth : -1);
   forkdepth = (fork != NULL ? (int)fork->depth : -1);
   assert(forkdepth <= focusnodedepth);
   assert(forkdepth < tree->pathlen);

   /* delay events in path switching */
   SCIP_CALL( SCIPeventqueueDelay(eventqueue) );

   /* undo the domain and constraint set changes of the old active path by deactivating the path's nodes */
   for( i = tree->pathlen-1; i > forkdepth; --i )
   {
      SCIP_CALL( nodeDeactivate(tree->path[i], blkmem, set, stat, tree, lp, branchcand, eventqueue) );
   }
   tree->pathlen = forkdepth+1;

   /* apply the pending bound changes */
   SCIP_CALL( treeApplyPendingBdchgs(tree, reopt, blkmem, set, stat, transprob, origprob, lp, branchcand, eventqueue, cliquetable) );

   /* create the new active path */
   SCIP_CALL( treeEnsurePathMem(tree, set, focusnodedepth+1) );
   while( focusnode != fork )
   {
      assert(focusnode != NULL);
      assert(!focusnode->active);
      assert(!focusnode->cutoff);
      tree->path[focusnode->depth] = focusnode;
      focusnode = focusnode->parent;
   }

   /* fork might be cut off when applying the pending bound changes */
   if( fork != NULL && fork->cutoff )
      *cutoff = TRUE;
   else if( fork != NULL && fork->reprop )
   {
     /* propagate common fork again, if the reprop flag is set */
      assert(tree->path[forkdepth] == fork);
      assert(fork->active);
      assert(!fork->cutoff);

      SCIP_CALL( nodeRepropagate(fork, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, conflict,
            eventfilter, eventqueue, cliquetable, cutoff) );
   }
   assert(fork != NULL || !(*cutoff));

   /* Apply domain and constraint set changes of the new path by activating the path's nodes;
    * on the way, domain propagation might be applied again to the path's nodes, which can result in the cutoff of
    * the node (and its subtree).
    * We only activate all nodes down to the parent of the new focus node, because the events in this process are
    * delayed, which means that multiple changes of a bound of a variable are merged (and might even be cancelled out,
    * if the bound is first relaxed when deactivating a node on the old path and then tightened to the same value
    * when activating a node on the new path).
    * This is valid for all nodes down to the parent of the new focus node, since they have already been propagated.
    * Bound change events on the new focus node, however, must not be cancelled out, since they need to be propagated
    * and thus, the event must be thrown and catched by the constraint handlers to mark constraints for propagation.
    */
   for( i = forkdepth+1; i < focusnodedepth && !(*cutoff); ++i )
   {
      assert(!tree->path[i]->cutoff);
      assert(tree->pathlen == i);

      /* activate the node, and apply domain propagation if the reprop flag is set */
      tree->pathlen++;
      SCIP_CALL( nodeActivate(tree->path[i], blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand,
            conflict, eventfilter, eventqueue, cliquetable, cutoff) );
   }

   /* process the delayed events */
   SCIP_CALL( SCIPeventqueueProcess(eventqueue, blkmem, set, primal, lp, branchcand, eventfilter) );

   /* activate the new focus node; there is no need to delay these events */
   if( !(*cutoff) && (i == focusnodedepth) )
   {
      assert(!tree->path[focusnodedepth]->cutoff);
      assert(tree->pathlen == focusnodedepth);

      /* activate the node, and apply domain propagation if the reprop flag is set */
      tree->pathlen++;
      SCIP_CALL( nodeActivate(tree->path[focusnodedepth], blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand,
            conflict, eventfilter, eventqueue, cliquetable, cutoff) );
   }

   /* mark last node of path to be cut off, if a cutoff was found */
   if( *cutoff )
   {
      assert(tree->pathlen > 0);
      assert(tree->path[tree->pathlen-1]->active);
      SCIP_CALL( SCIPnodeCutoff(tree->path[tree->pathlen-1], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
   }

   /* count the new LP sizes of the path */
   SCIP_CALL( treeUpdatePathLPSize(tree, forkdepth+1) );

   SCIPsetDebugMsg(set, "switch path: new pathlen=%d\n", tree->pathlen);

   return SCIP_OKAY;
}

/** loads the subroot's LP data */
static
SCIP_RETCODE subrootConstructLP(
   SCIP_NODE*            subroot,            /**< subroot node to construct LP for */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int ncols;
   int nrows;
   int c;
   int r;

   assert(subroot != NULL);
   assert(SCIPnodeGetType(subroot) == SCIP_NODETYPE_SUBROOT);
   assert(subroot->data.subroot != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   cols = subroot->data.subroot->cols;
   rows = subroot->data.subroot->rows;
   ncols = subroot->data.subroot->ncols;
   nrows = subroot->data.subroot->nrows;

   assert(ncols == 0 || cols != NULL);
   assert(nrows == 0 || rows != NULL);

   for( c = 0; c < ncols; ++c )
   {
      SCIP_CALL( SCIPlpAddCol(lp, set, cols[c], subroot->depth) );
   }
   for( r = 0; r < nrows; ++r )
   {
      SCIP_CALL( SCIPlpAddRow(lp, blkmem, set, eventqueue, eventfilter, rows[r], subroot->depth) );
   }

   return SCIP_OKAY;
}

/** loads the fork's additional LP data */
static
SCIP_RETCODE forkAddLP(
   SCIP_NODE*            fork,               /**< fork node to construct additional LP for */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int ncols;
   int nrows;
   int c;
   int r;

   assert(fork != NULL);
   assert(SCIPnodeGetType(fork) == SCIP_NODETYPE_FORK);
   assert(fork->data.fork != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   cols = fork->data.fork->addedcols;
   rows = fork->data.fork->addedrows;
   ncols = fork->data.fork->naddedcols;
   nrows = fork->data.fork->naddedrows;

   assert(ncols == 0 || cols != NULL);
   assert(nrows == 0 || rows != NULL);

   for( c = 0; c < ncols; ++c )
   {
      SCIP_CALL( SCIPlpAddCol(lp, set, cols[c], fork->depth) );
   }
   for( r = 0; r < nrows; ++r )
   {
      SCIP_CALL( SCIPlpAddRow(lp, blkmem, set, eventqueue, eventfilter, rows[r], fork->depth) );
   }

   return SCIP_OKAY;
}

/** loads the pseudofork's additional LP data */
static
SCIP_RETCODE pseudoforkAddLP(
   SCIP_NODE*            pseudofork,         /**< pseudofork node to construct additional LP for */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int ncols;
   int nrows;
   int c;
   int r;

   assert(pseudofork != NULL);
   assert(SCIPnodeGetType(pseudofork) == SCIP_NODETYPE_PSEUDOFORK);
   assert(pseudofork->data.pseudofork != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   cols = pseudofork->data.pseudofork->addedcols;
   rows = pseudofork->data.pseudofork->addedrows;
   ncols = pseudofork->data.pseudofork->naddedcols;
   nrows = pseudofork->data.pseudofork->naddedrows;

   assert(ncols == 0 || cols != NULL);
   assert(nrows == 0 || rows != NULL);

   for( c = 0; c < ncols; ++c )
   {
      SCIP_CALL( SCIPlpAddCol(lp, set, cols[c], pseudofork->depth) );
   }
   for( r = 0; r < nrows; ++r )
   {
      SCIP_CALL( SCIPlpAddRow(lp, blkmem, set, eventqueue, eventfilter, rows[r], pseudofork->depth) );
   }

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** checks validity of active path */
static
void treeCheckPath(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   SCIP_NODE* node;
   int ncols;
   int nrows;
   int d;

   assert(tree != NULL);
   assert(tree->path != NULL);

   ncols = 0;
   nrows = 0;
   for( d = 0; d < tree->pathlen; ++d )
   {
      node = tree->path[d];
      assert(node != NULL);
      assert((int)(node->depth) == d);
      switch( SCIPnodeGetType(node) )
      {  
      case SCIP_NODETYPE_PROBINGNODE:
         assert(SCIPtreeProbing(tree));
         assert(d >= 1);
         assert(SCIPnodeGetType(tree->path[d-1]) == SCIP_NODETYPE_FOCUSNODE
            || (ncols == node->data.probingnode->ninitialcols && nrows == node->data.probingnode->ninitialrows));
         assert(ncols <= node->data.probingnode->ncols || !tree->focuslpconstructed);
         assert(nrows <= node->data.probingnode->nrows || !tree->focuslpconstructed);
         if( d < tree->pathlen-1 )
         {
            ncols = node->data.probingnode->ncols;
            nrows = node->data.probingnode->nrows;
         }
         else
         {
            /* for the current probing node, the initial LP size is stored in the path */
            ncols = node->data.probingnode->ninitialcols;
            nrows = node->data.probingnode->ninitialrows;
         }
         break;
      case SCIP_NODETYPE_JUNCTION:
         break;
      case SCIP_NODETYPE_PSEUDOFORK:
         ncols += node->data.pseudofork->naddedcols;
         nrows += node->data.pseudofork->naddedrows;
         break;
      case SCIP_NODETYPE_FORK:
         ncols += node->data.fork->naddedcols;
         nrows += node->data.fork->naddedrows;
         break;
      case SCIP_NODETYPE_SUBROOT:
         ncols = node->data.subroot->ncols;
         nrows = node->data.subroot->nrows;
         break;
      case SCIP_NODETYPE_FOCUSNODE:
      case SCIP_NODETYPE_REFOCUSNODE:
         assert(d == tree->pathlen-1 || SCIPtreeProbing(tree));
         break;
      default:
         SCIPerrorMessage("node at depth %d on active path has to be of type JUNCTION, PSEUDOFORK, FORK, SUBROOT, FOCUSNODE, REFOCUSNODE, or PROBINGNODE, but is %d\n",
            d, SCIPnodeGetType(node));
         SCIPABORT();
      }  /*lint !e788*/
      assert(tree->pathnlpcols[d] == ncols);
      assert(tree->pathnlprows[d] == nrows);
   }
}
#else
#define treeCheckPath(tree) /**/
#endif

/** constructs the LP relaxation of the focus node */
SCIP_RETCODE SCIPtreeLoadLP(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool*            initroot            /**< pointer to store whether the root LP relaxation has to be initialized */
   )
{
   SCIP_NODE* lpfork;
   int lpforkdepth;
   int d;

   assert(tree != NULL);
   assert(!tree->focuslpconstructed);
   assert(tree->path != NULL);
   assert(tree->pathlen > 0);
   assert(tree->focusnode != NULL);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(SCIPnodeGetDepth(tree->focusnode) == tree->pathlen-1);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode == tree->path[tree->pathlen-1]);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(initroot != NULL);

   SCIPsetDebugMsg(set, "load LP for current fork node #%" SCIP_LONGINT_FORMAT " at depth %d\n",
      tree->focuslpfork == NULL ? -1 : SCIPnodeGetNumber(tree->focuslpfork),
      tree->focuslpfork == NULL ? -1 : SCIPnodeGetDepth(tree->focuslpfork));
   SCIPsetDebugMsg(set, "-> old LP has %d cols and %d rows\n", SCIPlpGetNCols(lp), SCIPlpGetNRows(lp));
   SCIPsetDebugMsg(set, "-> correct LP has %d cols and %d rows\n",
      tree->correctlpdepth >= 0 ? tree->pathnlpcols[tree->correctlpdepth] : 0,
      tree->correctlpdepth >= 0 ? tree->pathnlprows[tree->correctlpdepth] : 0);
   SCIPsetDebugMsg(set, "-> old correctlpdepth: %d\n", tree->correctlpdepth);

   treeCheckPath(tree);

   lpfork = tree->focuslpfork;

   /* find out the lpfork's depth (or -1, if lpfork is NULL) */
   if( lpfork == NULL )
   {
      assert(tree->correctlpdepth == -1 || tree->pathnlpcols[tree->correctlpdepth] == 0);
      assert(tree->correctlpdepth == -1 || tree->pathnlprows[tree->correctlpdepth] == 0);
      assert(tree->focuslpstatefork == NULL);
      assert(tree->focussubroot == NULL);
      lpforkdepth = -1;
   }
   else
   {
      assert(SCIPnodeGetType(lpfork) == SCIP_NODETYPE_PSEUDOFORK
         || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT);
      assert(lpfork->active);
      assert(tree->path[lpfork->depth] == lpfork);
      lpforkdepth = lpfork->depth;
   }
   assert(lpforkdepth < tree->pathlen-1); /* lpfork must not be the last (the focus) node of the active path */

   /* find out, if we are in the same subtree */
   if( tree->correctlpdepth >= 0 )
   {
      /* same subtree: shrink LP to the deepest node with correct LP */
      assert(lpforkdepth == -1 || tree->pathnlpcols[tree->correctlpdepth] <= tree->pathnlpcols[lpforkdepth]);
      assert(lpforkdepth == -1 || tree->pathnlprows[tree->correctlpdepth] <= tree->pathnlprows[lpforkdepth]);
      assert(lpforkdepth >= 0 || tree->pathnlpcols[tree->correctlpdepth] == 0);
      assert(lpforkdepth >= 0 || tree->pathnlprows[tree->correctlpdepth] == 0);
      SCIP_CALL( SCIPlpShrinkCols(lp, set, tree->pathnlpcols[tree->correctlpdepth]) );
      SCIP_CALL( SCIPlpShrinkRows(lp, blkmem, set, eventqueue, eventfilter, tree->pathnlprows[tree->correctlpdepth]) );
   }
   else
   {
      /* other subtree: fill LP with the subroot LP data */
      SCIP_CALL( SCIPlpClear(lp, blkmem, set, eventqueue, eventfilter) );
      if( tree->focussubroot != NULL )
      {
         SCIP_CALL( subrootConstructLP(tree->focussubroot, blkmem, set, eventqueue, eventfilter, lp) );
         tree->correctlpdepth = tree->focussubroot->depth; 
      }
   }

   assert(lpforkdepth < tree->pathlen);

   /* add the missing columns and rows */
   for( d = tree->correctlpdepth+1; d <= lpforkdepth; ++d )
   {
      SCIP_NODE* pathnode;

      pathnode = tree->path[d];
      assert(pathnode != NULL);
      assert((int)(pathnode->depth) == d);
      assert(SCIPnodeGetType(pathnode) == SCIP_NODETYPE_JUNCTION
         || SCIPnodeGetType(pathnode) == SCIP_NODETYPE_PSEUDOFORK
         || SCIPnodeGetType(pathnode) == SCIP_NODETYPE_FORK);
      if( SCIPnodeGetType(pathnode) == SCIP_NODETYPE_FORK )
      {
         SCIP_CALL( forkAddLP(pathnode, blkmem, set, eventqueue, eventfilter, lp) );
      }
      else if( SCIPnodeGetType(pathnode) == SCIP_NODETYPE_PSEUDOFORK )
      {
         SCIP_CALL( pseudoforkAddLP(pathnode, blkmem, set, eventqueue, eventfilter, lp) );
      }
   }
   tree->correctlpdepth = MAX(tree->correctlpdepth, lpforkdepth);
   assert(lpforkdepth == -1 || tree->pathnlpcols[tree->correctlpdepth] == tree->pathnlpcols[lpforkdepth]);
   assert(lpforkdepth == -1 || tree->pathnlprows[tree->correctlpdepth] == tree->pathnlprows[lpforkdepth]);
   assert(lpforkdepth == -1 || SCIPlpGetNCols(lp) == tree->pathnlpcols[lpforkdepth]);
   assert(lpforkdepth == -1 || SCIPlpGetNRows(lp) == tree->pathnlprows[lpforkdepth]);
   assert(lpforkdepth >= 0 || SCIPlpGetNCols(lp) == 0);
   assert(lpforkdepth >= 0 || SCIPlpGetNRows(lp) == 0);

   /* mark the LP's size, such that we know which rows and columns were added in the new node */
   SCIPlpMarkSize(lp);

   SCIPsetDebugMsg(set, "-> new correctlpdepth: %d\n", tree->correctlpdepth);
   SCIPsetDebugMsg(set, "-> new LP has %d cols and %d rows\n", SCIPlpGetNCols(lp), SCIPlpGetNRows(lp));

   /* if the correct LP depth is still -1, the root LP relaxation has to be initialized */
   *initroot = (tree->correctlpdepth == -1);

   /* mark the LP of the focus node constructed */
   tree->focuslpconstructed = TRUE;

   return SCIP_OKAY;
}

/** loads LP state for fork/subroot of the focus node */
SCIP_RETCODE SCIPtreeLoadLPState(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_NODE* lpstatefork;
   SCIP_Bool updatefeas;
   SCIP_Bool checkbdchgs;
   int lpstateforkdepth;
   int d;

   assert(tree != NULL);
   assert(tree->focuslpconstructed);
   assert(tree->path != NULL);
   assert(tree->pathlen > 0);
   assert(tree->focusnode != NULL);
   assert(tree->correctlpdepth < tree->pathlen);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(SCIPnodeGetDepth(tree->focusnode) == tree->pathlen-1);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode == tree->path[tree->pathlen-1]);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   SCIPsetDebugMsg(set, "load LP state for current fork node #%" SCIP_LONGINT_FORMAT " at depth %d\n",
      tree->focuslpstatefork == NULL ? -1 : SCIPnodeGetNumber(tree->focuslpstatefork),
      tree->focuslpstatefork == NULL ? -1 : SCIPnodeGetDepth(tree->focuslpstatefork));

   lpstatefork = tree->focuslpstatefork;

   /* if there is no LP state defining fork, nothing can be done */
   if( lpstatefork == NULL )
      return SCIP_OKAY;

   /* get the lpstatefork's depth */
   assert(SCIPnodeGetType(lpstatefork) == SCIP_NODETYPE_FORK || SCIPnodeGetType(lpstatefork) == SCIP_NODETYPE_SUBROOT);
   assert(lpstatefork->active);
   assert(tree->path[lpstatefork->depth] == lpstatefork);
   lpstateforkdepth = lpstatefork->depth;
   assert(lpstateforkdepth < tree->pathlen-1); /* lpstatefork must not be the last (the focus) node of the active path */
   assert(lpstateforkdepth <= tree->correctlpdepth); /* LP must have been constructed at least up to the fork depth */
   assert(tree->pathnlpcols[tree->correctlpdepth] >= tree->pathnlpcols[lpstateforkdepth]); /* LP can only grow */
   assert(tree->pathnlprows[tree->correctlpdepth] >= tree->pathnlprows[lpstateforkdepth]); /* LP can only grow */

   /* load LP state */
   if( tree->focuslpstateforklpcount != stat->lpcount )
   {
      if( SCIPnodeGetType(lpstatefork) == SCIP_NODETYPE_FORK )
      {
         assert(lpstatefork->data.fork != NULL);
         SCIP_CALL( SCIPlpSetState(lp, blkmem, set, eventqueue, lpstatefork->data.fork->lpistate,
               lpstatefork->data.fork->lpwasprimfeas, lpstatefork->data.fork->lpwasprimchecked,
               lpstatefork->data.fork->lpwasdualfeas, lpstatefork->data.fork->lpwasdualchecked) );
      }
      else
      {
         assert(SCIPnodeGetType(lpstatefork) == SCIP_NODETYPE_SUBROOT);
         assert(lpstatefork->data.subroot != NULL);
         SCIP_CALL( SCIPlpSetState(lp, blkmem, set, eventqueue, lpstatefork->data.subroot->lpistate,
               lpstatefork->data.subroot->lpwasprimfeas, lpstatefork->data.subroot->lpwasprimchecked,
               lpstatefork->data.subroot->lpwasdualfeas, lpstatefork->data.subroot->lpwasdualchecked) );
      }
      updatefeas = !lp->solved || !lp->solisbasic;
      checkbdchgs = TRUE;
   }
   else
   {
      updatefeas = TRUE;

      /* we do not need to check the bounds, since primalfeasible is updated anyway when flushing the LP */
      checkbdchgs = FALSE;
   }

   if( updatefeas )
   {
      /* check whether the size of the LP increased (destroying primal/dual feasibility) */
      lp->primalfeasible = lp->primalfeasible
         && (tree->pathnlprows[tree->correctlpdepth] == tree->pathnlprows[lpstateforkdepth]);
      lp->primalchecked = lp->primalchecked
         && (tree->pathnlprows[tree->correctlpdepth] == tree->pathnlprows[lpstateforkdepth]);
      lp->dualfeasible = lp->dualfeasible
         && (tree->pathnlpcols[tree->correctlpdepth] == tree->pathnlpcols[lpstateforkdepth]);
      lp->dualchecked = lp->dualchecked
         && (tree->pathnlpcols[tree->correctlpdepth] == tree->pathnlpcols[lpstateforkdepth]);

      /* check the path from LP fork to focus node for domain changes (destroying primal feasibility of LP basis) */
      if( checkbdchgs )
      {
         for( d = lpstateforkdepth; d < (int)(tree->focusnode->depth) && lp->primalfeasible; ++d )
         {
            assert(d < tree->pathlen);
            lp->primalfeasible = (tree->path[d]->domchg == NULL || tree->path[d]->domchg->domchgbound.nboundchgs == 0);
            lp->primalchecked = lp->primalfeasible;
         }
      }
   }

   SCIPsetDebugMsg(set, "-> primalfeasible=%u, dualfeasible=%u\n", lp->primalfeasible, lp->dualfeasible);

   return SCIP_OKAY;
}




/*
 * Node Conversion
 */

/** converts node into LEAF and moves it into the array of the node queue
 *  if node's lower bound is greater or equal than the given upper bound, the node is deleted;
 *  otherwise, it is moved to the node queue; anyways, the given pointer is NULL after the call
 */
static
SCIP_RETCODE nodeToLeaf(
   SCIP_NODE**           node,               /**< pointer to child or sibling node to convert */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_NODE*            lpstatefork,        /**< LP state defining fork of the node */
   SCIP_Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   )
{
   assert(SCIPnodeGetType(*node) == SCIP_NODETYPE_SIBLING || SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD
      || SCIPnodeGetType(*node) == SCIP_NODETYPE_FOCUSNODE);
   assert(stat != NULL);
   assert(lpstatefork == NULL || lpstatefork->depth < (*node)->depth);
   assert(lpstatefork == NULL || lpstatefork->active || SCIPsetIsGE(set, (*node)->lowerbound, cutoffbound));
   assert(lpstatefork == NULL
      || SCIPnodeGetType(lpstatefork) == SCIP_NODETYPE_FORK
      || SCIPnodeGetType(lpstatefork) == SCIP_NODETYPE_SUBROOT);

   /* convert node into leaf */
   SCIPsetDebugMsg(set, "convert node #%" SCIP_LONGINT_FORMAT " at depth %d to leaf with lpstatefork #%" SCIP_LONGINT_FORMAT " at depth %d\n",
      SCIPnodeGetNumber(*node), SCIPnodeGetDepth(*node),
      lpstatefork == NULL ? -1 : SCIPnodeGetNumber(lpstatefork),
      lpstatefork == NULL ? -1 : SCIPnodeGetDepth(lpstatefork));
   (*node)->nodetype = SCIP_NODETYPE_LEAF; /*lint !e641*/
   (*node)->data.leaf.lpstatefork = lpstatefork;

#ifndef NDEBUG
   /* check, if the LP state fork is the first node with LP state information on the path back to the root */
   if( !SCIPsetIsInfinity(set, -cutoffbound) ) /* if the node was cut off in SCIPnodeFocus(), the lpstatefork is invalid */
   {
      SCIP_NODE* pathnode;
      pathnode = (*node)->parent;
      while( pathnode != NULL && pathnode != lpstatefork )
      {
         assert(SCIPnodeGetType(pathnode) == SCIP_NODETYPE_JUNCTION
            || SCIPnodeGetType(pathnode) == SCIP_NODETYPE_PSEUDOFORK);
         pathnode = pathnode->parent;
      }
      assert(pathnode == lpstatefork);
   }
#endif

   /* if node is good enough to keep, put it on the node queue */
   if( SCIPsetIsLT(set, (*node)->lowerbound, cutoffbound) )
   {
      /* insert leaf in node queue */
      SCIP_CALL( SCIPnodepqInsert(tree->leaves, set, *node) );

      /* make the domain change data static to save memory */
      SCIP_CALL( SCIPdomchgMakeStatic(&(*node)->domchg, blkmem, set, eventqueue, lp) );

      /* node is now member of the node queue: delete the pointer to forbid further access */
      *node = NULL;
   }
   else
   {
      if( set->reopt_enable )
      {
         assert(reopt != NULL);
         /* check if the node should be stored for reoptimization */
         SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, *node, SCIP_EVENTTYPE_NODEINFEASIBLE, lp, SCIPlpGetSolstat(lp),
               tree->root == *node, tree->focusnode == *node, (*node)->lowerbound, tree->effectiverootdepth) );
      }

      /* delete node due to bound cut off */
      SCIPvisualCutoffNode(stat->visual, set, stat, *node, FALSE);
      SCIP_CALL( SCIPnodeFree(node, blkmem, set, stat, eventqueue, tree, lp) );
   }
   assert(*node == NULL);

   return SCIP_OKAY;
}

/** removes variables from the problem, that are marked to be deletable, and were created at the focusnode;
 *  only removes variables that were created at the focusnode, unless inlp is TRUE (e.g., when the node is cut off, anyway)
 */
static
SCIP_RETCODE focusnodeCleanupVars(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             inlp                /**< should variables in the LP be deleted, too?*/
   )
{
   SCIP_VAR* var;
   int i;
   int ndelvars;
   SCIP_Bool needdel;
   SCIP_Bool deleted;

   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode != NULL);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(lp != NULL);

   /* check the settings, whether variables should be deleted */
   needdel = (tree->focusnode == tree->root ? set->price_delvarsroot : set->price_delvars);

   if( !needdel )
      return SCIP_OKAY;

   ndelvars = 0;

   /* also delete variables currently in the LP, thus remove all new variables from the LP, first */
   if( inlp )
   {
      /* remove all additions to the LP at this node */
      SCIP_CALL( SCIPlpShrinkCols(lp, set, SCIPlpGetNCols(lp) - SCIPlpGetNNewcols(lp)) );

      SCIP_CALL( SCIPlpFlush(lp, blkmem, set, eventqueue) );
   }

   /* mark variables as deleted */
   for( i = 0; i < transprob->nvars; i++ )
   {
      var = transprob->vars[i];
      assert(var != NULL);

      /* check whether variable is deletable */
      if( SCIPvarIsDeletable(var) )
      {
         if( !SCIPvarIsInLP(var) )
         {
            /* fix the variable to 0, first */
            assert(!SCIPsetIsFeasPositive(set, SCIPvarGetLbGlobal(var)));
            assert(!SCIPsetIsFeasNegative(set, SCIPvarGetUbGlobal(var)));

            if( !SCIPsetIsFeasZero(set, SCIPvarGetLbGlobal(var)) )
            {
               SCIP_CALL( SCIPnodeAddBoundchg(tree->root, blkmem, set, stat, transprob, origprob,
                     tree, reopt, lp, branchcand, eventqueue, cliquetable, var, 0.0, SCIP_BOUNDTYPE_LOWER, FALSE) );
            }
            if( !SCIPsetIsFeasZero(set, SCIPvarGetUbGlobal(var)) )
            {
               SCIP_CALL( SCIPnodeAddBoundchg(tree->root, blkmem, set, stat, transprob, origprob,
                     tree, reopt, lp, branchcand, eventqueue, cliquetable, var, 0.0, SCIP_BOUNDTYPE_UPPER, FALSE) );
            }

            SCIP_CALL( SCIPprobDelVar(transprob, blkmem, set, eventqueue, var, &deleted) );

            if( deleted )
               ndelvars++;
         }
         else
         {
            /* mark variable to be non-deletable, because it will be contained in the basis information
             * at this node and must not be deleted from now on
             */
            SCIPvarMarkNotDeletable(var);
         }
      }
   }

   SCIPsetDebugMsg(set, "delvars at node %" SCIP_LONGINT_FORMAT ", deleted %d vars\n", stat->nnodes, ndelvars);

   if( ndelvars > 0 )
   {
      /* perform the variable deletions from the problem */
      SCIP_CALL( SCIPprobPerformVarDeletions(transprob, blkmem, set, stat, eventqueue, cliquetable, lp, branchcand) );
   }

   return SCIP_OKAY;
}

/** converts the focus node into a dead-end node */
static
SCIP_RETCODE focusnodeToDeadend(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode != NULL);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(tree->nchildren == 0);

   SCIPsetDebugMsg(set, "focusnode #%" SCIP_LONGINT_FORMAT " to dead-end at depth %d\n",
      SCIPnodeGetNumber(tree->focusnode), SCIPnodeGetDepth(tree->focusnode));

   /* remove variables from the problem that are marked as deletable and were created at this node */
   SCIP_CALL( focusnodeCleanupVars(blkmem, set, stat, eventqueue, transprob, origprob, tree, reopt, lp, branchcand, cliquetable, TRUE) );

   tree->focusnode->nodetype = SCIP_NODETYPE_DEADEND; /*lint !e641*/

   /* release LPI state */
   if( tree->focuslpstatefork != NULL )
   {
      SCIP_CALL( SCIPnodeReleaseLPIState(tree->focuslpstatefork, blkmem, lp) );
   }

   return SCIP_OKAY;
}

/** converts the focus node into a leaf node (if it was postponed) */
static
SCIP_RETCODE focusnodeToLeaf(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_NODE*            lpstatefork,        /**< LP state defining fork of the node */
   SCIP_Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */

   )
{
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode != NULL);
   assert(tree->focusnode->active);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);

   SCIPsetDebugMsg(set, "focusnode #%" SCIP_LONGINT_FORMAT " to leaf at depth %d\n",
      SCIPnodeGetNumber(tree->focusnode), SCIPnodeGetDepth(tree->focusnode));

   SCIP_CALL( nodeToLeaf(&tree->focusnode, blkmem, set, stat, eventqueue, tree, reopt, lp, lpstatefork, cutoffbound));

   return SCIP_OKAY;
}

/** converts the focus node into a junction node */
static
SCIP_RETCODE focusnodeToJunction(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode != NULL);
   assert(tree->focusnode->active); /* otherwise, no children could be created at the focus node */
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(SCIPlpGetNNewcols(lp) == 0);

   SCIPsetDebugMsg(set, "focusnode #%" SCIP_LONGINT_FORMAT " to junction at depth %d\n",
      SCIPnodeGetNumber(tree->focusnode), SCIPnodeGetDepth(tree->focusnode));

   /* convert node into junction */
   tree->focusnode->nodetype = SCIP_NODETYPE_JUNCTION; /*lint !e641*/

   SCIP_CALL( junctionInit(&tree->focusnode->data.junction, tree) );

   /* release LPI state */
   if( tree->focuslpstatefork != NULL )
   {
      SCIP_CALL( SCIPnodeReleaseLPIState(tree->focuslpstatefork, blkmem, lp) );
   }

   /* make the domain change data static to save memory */
   SCIP_CALL( SCIPdomchgMakeStatic(&tree->focusnode->domchg, blkmem, set, eventqueue, lp) );

   return SCIP_OKAY;
}

/** converts the focus node into a pseudofork node */
static
SCIP_RETCODE focusnodeToPseudofork(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   SCIP_PSEUDOFORK* pseudofork;

   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode != NULL);
   assert(tree->focusnode->active); /* otherwise, no children could be created at the focus node */
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(tree->nchildren > 0);
   assert(lp != NULL);

   SCIPsetDebugMsg(set, "focusnode #%" SCIP_LONGINT_FORMAT " to pseudofork at depth %d\n",
      SCIPnodeGetNumber(tree->focusnode), SCIPnodeGetDepth(tree->focusnode));

   /* remove variables from the problem that are marked as deletable and were created at this node */
   SCIP_CALL( focusnodeCleanupVars(blkmem, set, stat, eventqueue, transprob, origprob, tree, reopt, lp, branchcand, cliquetable, FALSE) );

   /* create pseudofork data */
   SCIP_CALL( pseudoforkCreate(&pseudofork, blkmem, tree, lp) );

   tree->focusnode->nodetype = SCIP_NODETYPE_PSEUDOFORK; /*lint !e641*/
   tree->focusnode->data.pseudofork = pseudofork;

   /* release LPI state */
   if( tree->focuslpstatefork != NULL )
   {
      SCIP_CALL( SCIPnodeReleaseLPIState(tree->focuslpstatefork, blkmem, lp) );
   }

   /* make the domain change data static to save memory */
   SCIP_CALL( SCIPdomchgMakeStatic(&tree->focusnode->domchg, blkmem, set, eventqueue, lp) );

   return SCIP_OKAY;
}

/** converts the focus node into a fork node */
static
SCIP_RETCODE focusnodeToFork(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   SCIP_FORK* fork;
   SCIP_Bool lperror;

   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode != NULL);
   assert(tree->focusnode->active); /* otherwise, no children could be created at the focus node */
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(tree->nchildren > 0);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved || lp->resolvelperror);

   SCIPsetDebugMsg(set, "focusnode #%" SCIP_LONGINT_FORMAT " to fork at depth %d\n",
      SCIPnodeGetNumber(tree->focusnode), SCIPnodeGetDepth(tree->focusnode));

   /* usually, the LP should be solved to optimality; otherwise, numerical troubles occured,
    * and we have to forget about the LP and transform the node into a junction (see below)
    */
   lperror = FALSE;
   if( !lp->resolvelperror && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      /* clean up newly created part of LP to keep only necessary columns and rows */
      SCIP_CALL( SCIPlpCleanupNew(lp, blkmem, set, stat, eventqueue, eventfilter, (tree->focusnode->depth == 0)) );

      /* resolve LP after cleaning up */
      if( !lp->solved || !lp->flushed )
      {
         SCIPsetDebugMsg(set, "resolving LP after cleanup\n");
         SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, transprob, -1LL, FALSE, FALSE, TRUE, &lperror) );
      }
   }
   assert(lp->flushed);
   assert(lp->solved || lperror || lp->resolvelperror);

   /* There are two reasons, that the (reduced) LP is not solved to optimality:
    *  - The primal heuristics (called after the current node's LP was solved) found a new 
    *    solution, that is better than the current node's lower bound.
    *    (But in this case, all children should be cut off and the node should be converted
    *    into a dead-end instead of a fork.)
    *  - Something numerically weird happened after cleaning up or after resolving a diving or probing LP.
    * The only thing we can do, is to completely forget about the LP and treat the node as
    * if it was only a pseudo-solution node. Therefore we have to remove all additional
    * columns and rows from the LP and convert the node into a junction.
    * However, the node's lower bound is kept, thus automatically throwing away nodes that
    * were cut off due to a primal solution.
    */
   if( lperror || lp->resolvelperror || SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %" SCIP_LONGINT_FORMAT ") numerical troubles: LP %" SCIP_LONGINT_FORMAT " not optimal -- convert node into junction instead of fork\n",
         stat->nnodes, stat->nlps);

      /* remove all additions to the LP at this node */
      SCIP_CALL( SCIPlpShrinkCols(lp, set, SCIPlpGetNCols(lp) - SCIPlpGetNNewcols(lp)) );
      SCIP_CALL( SCIPlpShrinkRows(lp, blkmem, set, eventqueue, eventfilter, SCIPlpGetNRows(lp) - SCIPlpGetNNewrows(lp)) );

      /* convert node into a junction */
      SCIP_CALL( focusnodeToJunction(blkmem, set, eventqueue, tree, lp) );

      return SCIP_OKAY;
   }
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);

   /* remove variables from the problem that are marked as deletable, were created at this node and are not contained in the LP */
   SCIP_CALL( focusnodeCleanupVars(blkmem, set, stat, eventqueue, transprob, origprob, tree, reopt, lp, branchcand, cliquetable, FALSE) );

   assert(lp->flushed);
   assert(lp->solved);

   /* create fork data */
   SCIP_CALL( forkCreate(&fork, blkmem, set, transprob, tree, lp) );

   tree->focusnode->nodetype = SCIP_NODETYPE_FORK; /*lint !e641*/
   tree->focusnode->data.fork = fork;

   /* capture the LPI state of the root node to ensure that the LPI state of the root stays for the whole solving
    * process
    */
   if( tree->focusnode == tree->root )
      forkCaptureLPIState(fork, 1);

   /* release LPI state */
   if( tree->focuslpstatefork != NULL )
   {
      SCIP_CALL( SCIPnodeReleaseLPIState(tree->focuslpstatefork, blkmem, lp) );
   }

   /* make the domain change data static to save memory */
   SCIP_CALL( SCIPdomchgMakeStatic(&tree->focusnode->domchg, blkmem, set, eventqueue, lp) );

   return SCIP_OKAY;
}

#ifdef WITHSUBROOTS  /** @todo test whether subroots should be created */
/** converts the focus node into a subroot node */
static
SCIP_RETCODE focusnodeToSubroot(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   SCIP_SUBROOT* subroot;
   SCIP_Bool lperror;

   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode != NULL);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(tree->focusnode->active); /* otherwise, no children could be created at the focus node */
   assert(tree->nchildren > 0);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   SCIPsetDebugMsg(set, "focusnode #%" SCIP_LONGINT_FORMAT " to subroot at depth %d\n",
      SCIPnodeGetNumber(tree->focusnode), SCIPnodeGetDepth(tree->focusnode));

   /* usually, the LP should be solved to optimality; otherwise, numerical troubles occured,
    * and we have to forget about the LP and transform the node into a junction (see below)
    */
   lperror = FALSE;
   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      /* clean up whole LP to keep only necessary columns and rows */
#ifdef SCIP_DISABLED_CODE
      if( tree->focusnode->depth == 0 )
      {
         SCIP_CALL( SCIPlpCleanupAll(lp, blkmem, set, stat, eventqueue, eventfilter, (tree->focusnode->depth == 0)) );
      }
      else
#endif
      {
         SCIP_CALL( SCIPlpRemoveAllObsoletes(lp, blkmem, set, stat, eventqueue, eventfilter) );
      }

      /* resolve LP after cleaning up */
      if( !lp->solved || !lp->flushed )
      {
         SCIPsetDebugMsg(set, "resolving LP after cleanup\n");
         SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, transprob, -1LL, FALSE, FALSE, TRUE, &lperror) );
      }
   }
   assert(lp->flushed);
   assert(lp->solved || lperror);

   /* There are two reasons, that the (reduced) LP is not solved to optimality:
    *  - The primal heuristics (called after the current node's LP was solved) found a new
    *    solution, that is better than the current node's lower bound.
    *    (But in this case, all children should be cut off and the node should be converted
    *    into a dead-end instead of a subroot.)
    *  - Something numerically weird happened after cleaning up.
    * The only thing we can do, is to completely forget about the LP and treat the node as
    * if it was only a pseudo-solution node. Therefore we have to remove all additional
    * columns and rows from the LP and convert the node into a junction.
    * However, the node's lower bound is kept, thus automatically throwing away nodes that
    * were cut off due to a primal solution.
    */
   if( lperror || SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %" SCIP_LONGINT_FORMAT ") numerical troubles: LP %" SCIP_LONGINT_FORMAT " not optimal -- convert node into junction instead of subroot\n",
         stat->nnodes, stat->nlps);

      /* remove all additions to the LP at this node */
      SCIP_CALL( SCIPlpShrinkCols(lp, set, SCIPlpGetNCols(lp) - SCIPlpGetNNewcols(lp)) );
      SCIP_CALL( SCIPlpShrinkRows(lp, blkmem, set, eventqueue, eventfilter, SCIPlpGetNRows(lp) - SCIPlpGetNNewrows(lp)) );

      /* convert node into a junction */
      SCIP_CALL( focusnodeToJunction(blkmem, set, eventqueue, tree, lp) );

      return SCIP_OKAY;
   }
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);

   /* remove variables from the problem that are marked as deletable, were created at this node and are not contained in the LP */
   SCIP_CALL( focusnodeCleanupVars(blkmem, set, stat, eventqueue, transprob, origprob, tree, lp, branchcand, cliquetable, FALSE) );

   assert(lp->flushed);
   assert(lp->solved);


   /* create subroot data */
   SCIP_CALL( subrootCreate(&subroot, blkmem, set, transprob, tree, lp) );

   tree->focusnode->nodetype = SCIP_NODETYPE_SUBROOT; /*lint !e641*/
   tree->focusnode->data.subroot = subroot;

   /* update the LP column and row counter for the converted node */
   SCIP_CALL( treeUpdatePathLPSize(tree, tree->focusnode->depth) );

   /* release LPI state */
   if( tree->focuslpstatefork != NULL )
   {
      SCIP_CALL( SCIPnodeReleaseLPIState(tree->focuslpstatefork, blkmem, lp) );
   }

   /* make the domain change data static to save memory */
   SCIP_CALL( SCIPdomchgMakeStatic(&tree->focusnode->domchg, blkmem, set, eventqueue, lp) );

   return SCIP_OKAY;
}
#endif

/** puts all nodes in the array on the node queue and makes them LEAFs */
static
SCIP_RETCODE treeNodesToQueue(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_NODE**           nodes,              /**< array of nodes to put on the queue */
   int*                  nnodes,             /**< pointer to number of nodes in the array */
   SCIP_NODE*            lpstatefork,        /**< LP state defining fork of the nodes */
   SCIP_Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   )
{
   int i;

   assert(tree != NULL);
   assert(set != NULL);
   assert(nnodes != NULL);
   assert(*nnodes == 0 || nodes != NULL);

   for( i = 0; i < *nnodes; ++i )
   {
      /* convert node to LEAF and put it into leaves queue, or delete it if it's lower bound exceeds the cutoff bound */
      SCIP_CALL( nodeToLeaf(&nodes[i], blkmem, set, stat, eventqueue, tree, reopt, lp, lpstatefork, cutoffbound) );
      assert(nodes[i] == NULL);
   }
   *nnodes = 0;

   return SCIP_OKAY;
}

/** converts children into siblings, clears children array */
static
void treeChildrenToSiblings(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   SCIP_NODE** tmpnodes;
   SCIP_Real* tmpprios;
   int tmpnodessize;
   int i;

   assert(tree != NULL);
   assert(tree->nsiblings == 0);

   tmpnodes = tree->siblings;
   tmpprios = tree->siblingsprio;
   tmpnodessize = tree->siblingssize;

   tree->siblings = tree->children;
   tree->siblingsprio = tree->childrenprio;
   tree->nsiblings = tree->nchildren;
   tree->siblingssize = tree->childrensize;

   tree->children = tmpnodes;
   tree->childrenprio = tmpprios;
   tree->nchildren = 0;
   tree->childrensize = tmpnodessize;

   for( i = 0; i < tree->nsiblings; ++i )
   {
      assert(SCIPnodeGetType(tree->siblings[i]) == SCIP_NODETYPE_CHILD);
      tree->siblings[i]->nodetype = SCIP_NODETYPE_SIBLING; /*lint !e641*/

      /* because CHILD and SIBLING structs contain the same data in the same order, we do not have to copy it */
      assert(&(tree->siblings[i]->data.sibling.arraypos) == &(tree->siblings[i]->data.child.arraypos));
   }
}

/** installs a child, a sibling, or a leaf node as the new focus node */
SCIP_RETCODE SCIPnodeFocus(
   SCIP_NODE**           node,               /**< pointer to node to focus (or NULL to remove focus); the node
                                              *   is freed, if it was cut off due to a cut off subtree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the given node can be cut off */
   SCIP_Bool             postponed,          /**< was the current focus node postponed? */
   SCIP_Bool             exitsolve           /**< are we in exitsolve stage, so we only need to loose the children */
   )
{  /*lint --e{715}*/
   SCIP_NODE* oldfocusnode;
   SCIP_NODE* fork;
   SCIP_NODE* lpfork;
   SCIP_NODE* lpstatefork;
   SCIP_NODE* subroot;
   SCIP_NODE* childrenlpstatefork;
   int oldcutoffdepth;

   assert(node != NULL);
   assert(*node == NULL
      || SCIPnodeGetType(*node) == SCIP_NODETYPE_SIBLING
      || SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD
      || SCIPnodeGetType(*node) == SCIP_NODETYPE_LEAF);
   assert(*node == NULL || !(*node)->active);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(lp != NULL);
   assert(cutoff != NULL);

   SCIPsetDebugMsg(set, "focusing node #%" SCIP_LONGINT_FORMAT " of type %d in depth %d\n",
      *node != NULL ? SCIPnodeGetNumber(*node) : -1, *node != NULL ? (int)SCIPnodeGetType(*node) : 0,
      *node != NULL ? SCIPnodeGetDepth(*node) : -1);

   /* remember old cutoff depth in order to know, whether the children and siblings can be deleted */
   oldcutoffdepth = tree->cutoffdepth;

   /* find the common fork node, the new LP defining fork, and the new focus subroot,
    * thereby checking, if the new node can be cut off
    */
   treeFindSwitchForks(tree, *node, &fork, &lpfork, &lpstatefork, &subroot, cutoff);
   SCIPsetDebugMsg(set, "focus node: focusnodedepth=%d, forkdepth=%d, lpforkdepth=%d, lpstateforkdepth=%d, subrootdepth=%d, cutoff=%u\n",
      *node != NULL ? (*node)->depth : -1, fork != NULL ? fork->depth : -1, /*lint !e705 */
      lpfork != NULL ? lpfork->depth : -1, lpstatefork != NULL ? lpstatefork->depth : -1, /*lint !e705 */
      subroot != NULL ? subroot->depth : -1, *cutoff); /*lint !e705 */

   /* free the new node, if it is located in a cut off subtree */
   if( *cutoff )
   {
      assert(*node != NULL);
      assert(tree->cutoffdepth == oldcutoffdepth);
      if( SCIPnodeGetType(*node) == SCIP_NODETYPE_LEAF )
      {
         SCIP_CALL( SCIPnodepqRemove(tree->leaves, set, *node) );
      }
      SCIPvisualCutoffNode(stat->visual, set, stat, *node, FALSE);

      if( set->reopt_enable )
      {
         assert(reopt != NULL);
         /* check if the node should be stored for reoptimization */
         SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, *node, SCIP_EVENTTYPE_NODEINFEASIBLE, lp, SCIPlpGetSolstat(lp),
               tree->root == (*node), tree->focusnode == (*node), (*node)->lowerbound, tree->effectiverootdepth) );
      }

      SCIP_CALL( SCIPnodeFree(node, blkmem, set, stat, eventqueue, tree, lp) );

      return SCIP_OKAY;
   }

   assert(tree->cutoffdepth == INT_MAX);
   assert(fork == NULL || fork->active);
   assert(lpstatefork == NULL || lpfork != NULL);
   assert(subroot == NULL || lpstatefork != NULL);

   /* remember the depth of the common fork node for LP updates */
   SCIPsetDebugMsg(set, "focus node: old correctlpdepth=%d\n", tree->correctlpdepth);
   if( subroot == tree->focussubroot && fork != NULL && lpfork != NULL )
   {
      /* we are in the same subtree with valid LP fork: the LP is correct at most upto the common fork depth */
      assert(subroot == NULL || subroot->active);
      tree->correctlpdepth = MIN(tree->correctlpdepth, (int)fork->depth);
   }
   else
   {
      /* we are in a different subtree, or no valid LP fork exists: the LP is completely incorrect */
      assert(subroot == NULL || !subroot->active
         || (tree->focussubroot != NULL && (int)(tree->focussubroot->depth) > subroot->depth));
      tree->correctlpdepth = -1;
   }

   /* if the LP state fork changed, the lpcount information for the new LP state fork is unknown */
   if( lpstatefork != tree->focuslpstatefork )
      tree->focuslpstateforklpcount = -1;

   /* in exitsolve we only need to take care of open children
    *
    * @note because we might do a 'newstart' and converted cuts to constraints might have rendered the LP in the current
    *       focusnode unsolved the latter code would have resolved the LP unnecessarily
    */
   if( exitsolve && tree->nchildren > 0 )
   {
      SCIPsetDebugMsg(set, " -> deleting the %d children (in exitsolve) of the old focus node\n", tree->nchildren);
      SCIP_CALL( treeNodesToQueue(tree, reopt, blkmem, set, stat, eventqueue, lp, tree->children, &tree->nchildren, NULL, -SCIPsetInfinity(set)) );
      assert(tree->nchildren == 0);
   }

   /* if the old focus node was cut off, we can delete its children;
    * if the old focus node's parent was cut off, we can also delete the focus node's siblings
    */
   if( tree->focusnode != NULL && oldcutoffdepth <= (int)tree->focusnode->depth )
   {
      SCIPsetDebugMsg(set, "path to old focus node of depth %u was cut off at depth %d\n", tree->focusnode->depth, oldcutoffdepth);

      /* delete the focus node's children by converting them to leaves with a cutoffbound of -SCIPsetInfinity(set);
       * we cannot delete them directly, because in SCIPnodeFree(), the children array is changed, which is the
       * same array we would have to iterate over here;
       * the children don't have an LP fork, because the old focus node is not yet converted into a fork or subroot
       */
      SCIPsetDebugMsg(set, " -> deleting the %d children of the old focus node\n", tree->nchildren);
      SCIP_CALL( treeNodesToQueue(tree, reopt, blkmem, set, stat, eventqueue, lp, tree->children, &tree->nchildren, NULL, -SCIPsetInfinity(set)) );
      assert(tree->nchildren == 0);

      if( oldcutoffdepth < (int)tree->focusnode->depth )
      {
         /* delete the focus node's siblings by converting them to leaves with a cutoffbound of -SCIPsetInfinity(set);
          * we cannot delete them directly, because in SCIPnodeFree(), the siblings array is changed, which is the
          * same array we would have to iterate over here;
          * the siblings have the same LP state fork as the old focus node
          */
         SCIPsetDebugMsg(set, " -> deleting the %d siblings of the old focus node\n", tree->nsiblings);
         SCIP_CALL( treeNodesToQueue(tree, reopt, blkmem, set, stat, eventqueue, lp, tree->siblings, &tree->nsiblings, tree->focuslpstatefork,
               -SCIPsetInfinity(set)) );
         assert(tree->nsiblings == 0);
      }
   }

   /* convert the old focus node into a fork or subroot node, if it has children;
    * otherwise, convert it into a dead-end, which will be freed later in treeSwitchPath();
    * if the node was postponed, make it a leaf.
    */
   childrenlpstatefork = tree->focuslpstatefork;

   assert(!postponed || *node == NULL);
   assert(!postponed || tree->focusnode != NULL);

   if( postponed )
   {
      assert(tree->nchildren == 0);
      assert(*node == NULL);

      /* if the node is infeasible, convert it into a deadend; otherwise, put it into the LEAF queue */
      if( SCIPsetIsGE(set, tree->focusnode->lowerbound, primal->cutoffbound) )
      {
         /* in case the LP was not constructed (due to the parameter settings for example) we have the finally remember the
          * old size of the LP (if it was constructed in an earlier node) before we change the current node into a deadend
          */
         if( !tree->focuslpconstructed )
            SCIPlpMarkSize(lp);

         /* convert old focus node into deadend */
         SCIP_CALL( focusnodeToDeadend(blkmem, set, stat, eventqueue, transprob, origprob, tree, reopt, lp, branchcand,
               cliquetable) );
      }
      else
      {
         SCIP_CALL( focusnodeToLeaf(blkmem, set, stat, eventqueue, tree, reopt, lp, tree->focuslpstatefork,
               SCIPsetInfinity(set)) );
      }
   }
   else if( tree->nchildren > 0 )
   {
      SCIP_Bool selectedchild;

      assert(tree->focusnode != NULL);
      assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
      assert(oldcutoffdepth == INT_MAX);

      /* check whether the next focus node is a child of the old focus node */
      selectedchild = (*node != NULL && SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD);

      if( tree->focusnodehaslp && lp->isrelax )
      {
         assert(tree->focuslpconstructed);

#ifdef WITHSUBROOTS  /** @todo test whether subroots should be created, decide: old focus node becomes fork or subroot */
         if( tree->focusnode->depth > 0 && tree->focusnode->depth % 25 == 0 )
         {
            /* convert old focus node into a subroot node */
            SCIP_CALL( focusnodeToSubroot(blkmem, set, messagehdlr, stat, eventqueue, eventfilter, transprob, origprob, tree, lp, branchcand) );
            if( *node != NULL && SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD
               && SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_SUBROOT )
               subroot = tree->focusnode;
         }
         else
#endif
         {
            /* convert old focus node into a fork node */
            SCIP_CALL( focusnodeToFork(blkmem, set, messagehdlr, stat, eventqueue, eventfilter, transprob, origprob, tree,
                  reopt, lp, branchcand, cliquetable) );
         }

         /* check, if the conversion into a subroot or fork was successful */
         if( SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FORK
            || SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_SUBROOT )
         {
            childrenlpstatefork = tree->focusnode;

            /* if a child of the old focus node was selected as new focus node, the old node becomes the new focus
             * LP fork and LP state fork
             */
            if( selectedchild )
            {
               lpfork = tree->focusnode;
               tree->correctlpdepth = tree->focusnode->depth;
               lpstatefork = tree->focusnode;
               tree->focuslpstateforklpcount = stat->lpcount;
            }
         }

         /* update the path's LP size */
         tree->pathnlpcols[tree->focusnode->depth] = SCIPlpGetNCols(lp);
         tree->pathnlprows[tree->focusnode->depth] = SCIPlpGetNRows(lp);
      }
      else if( tree->focuslpconstructed && (SCIPlpGetNNewcols(lp) > 0 || SCIPlpGetNNewrows(lp) > 0) )
      {
         /* convert old focus node into pseudofork */
         SCIP_CALL( focusnodeToPseudofork(blkmem, set, stat, eventqueue, transprob, origprob, tree, reopt, lp,
               branchcand, cliquetable) );
         assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_PSEUDOFORK);

         /* update the path's LP size */
         tree->pathnlpcols[tree->focusnode->depth] = SCIPlpGetNCols(lp);
         tree->pathnlprows[tree->focusnode->depth] = SCIPlpGetNRows(lp);

         /* if a child of the old focus node was selected as new focus node, the old node becomes the new focus LP fork */
         if( selectedchild )
         {
            lpfork = tree->focusnode;
            tree->correctlpdepth = tree->focusnode->depth;
         }
      }
      else
      {
         /* in case the LP was not constructed (due to the parameter settings for example) we have the finally remember the
          * old size of the LP (if it was constructed in an earlier node) before we change the current node into a junction
          */
         SCIPlpMarkSize(lp);

         /* convert old focus node into junction */
         SCIP_CALL( focusnodeToJunction(blkmem, set, eventqueue, tree, lp) );
      }
   }
   else if( tree->focusnode != NULL )
   {
      /* in case the LP was not constructed (due to the parameter settings for example) we have the finally remember the
       * old size of the LP (if it was constructed in an earlier node) before we change the current node into a deadend
       */
      if( !tree->focuslpconstructed )
         SCIPlpMarkSize(lp);

      /* convert old focus node into deadend */
      SCIP_CALL( focusnodeToDeadend(blkmem, set, stat, eventqueue, transprob, origprob, tree, reopt, lp, branchcand, cliquetable) );
   }
   assert(subroot == NULL || SCIPnodeGetType(subroot) == SCIP_NODETYPE_SUBROOT);
   assert(lpstatefork == NULL
      || SCIPnodeGetType(lpstatefork) == SCIP_NODETYPE_SUBROOT
      || SCIPnodeGetType(lpstatefork) == SCIP_NODETYPE_FORK);
   assert(childrenlpstatefork == NULL
      || SCIPnodeGetType(childrenlpstatefork) == SCIP_NODETYPE_SUBROOT
      || SCIPnodeGetType(childrenlpstatefork) == SCIP_NODETYPE_FORK);
   assert(lpfork == NULL
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_PSEUDOFORK);
   SCIPsetDebugMsg(set, "focus node: new correctlpdepth=%d\n", tree->correctlpdepth);

   /* set up the new lists of siblings and children */
   oldfocusnode = tree->focusnode;
   if( *node == NULL )
   {
      /* move siblings to the queue, make them LEAFs */
      SCIP_CALL( treeNodesToQueue(tree, reopt, blkmem, set, stat, eventqueue, lp, tree->siblings, &tree->nsiblings, tree->focuslpstatefork,
            primal->cutoffbound) );

      /* move children to the queue, make them LEAFs */
      SCIP_CALL( treeNodesToQueue(tree, reopt, blkmem, set, stat, eventqueue, lp, tree->children, &tree->nchildren, childrenlpstatefork,
            primal->cutoffbound) );
   }
   else
   {
      SCIP_NODE* bestleaf;

      switch( SCIPnodeGetType(*node) )
      {  
      case SCIP_NODETYPE_SIBLING:
         /* reset plunging depth, if the selected node is better than all leaves */
         bestleaf = SCIPtreeGetBestLeaf(tree);
         if( bestleaf == NULL || SCIPnodepqCompare(tree->leaves, set, *node, bestleaf) <= 0 )
            stat->plungedepth = 0;

         /* move children to the queue, make them LEAFs */
         SCIP_CALL( treeNodesToQueue(tree, reopt, blkmem, set, stat, eventqueue, lp, tree->children, &tree->nchildren, childrenlpstatefork,
               primal->cutoffbound) );

         /* remove selected sibling from the siblings array */
         treeRemoveSibling(tree, *node);

         SCIPsetDebugMsg(set, "selected sibling node, lowerbound=%g, plungedepth=%d\n", (*node)->lowerbound, stat->plungedepth);
         break;

      case SCIP_NODETYPE_CHILD:
         /* reset plunging depth, if the selected node is better than all leaves; otherwise, increase plunging depth */
         bestleaf = SCIPtreeGetBestLeaf(tree);
         if( bestleaf == NULL || SCIPnodepqCompare(tree->leaves, set, *node, bestleaf) <= 0 )
            stat->plungedepth = 0;
         else
            stat->plungedepth++;

         /* move siblings to the queue, make them LEAFs */
         SCIP_CALL( treeNodesToQueue(tree, reopt, blkmem, set, stat, eventqueue, lp, tree->siblings, &tree->nsiblings, tree->focuslpstatefork,
               primal->cutoffbound) );

         /* remove selected child from the children array */      
         treeRemoveChild(tree, *node);

         /* move remaining children to the siblings array, make them SIBLINGs */
         treeChildrenToSiblings(tree);

         SCIPsetDebugMsg(set, "selected child node, lowerbound=%g, plungedepth=%d\n", (*node)->lowerbound, stat->plungedepth);
         break;

      case SCIP_NODETYPE_LEAF:
         /* move siblings to the queue, make them LEAFs */
         SCIP_CALL( treeNodesToQueue(tree, reopt, blkmem, set, stat, eventqueue, lp, tree->siblings, &tree->nsiblings, tree->focuslpstatefork,
               primal->cutoffbound) );

         /* encounter an early backtrack if there is a child which does not exceed given reference bound */
         if( !SCIPsetIsInfinity(set, stat->referencebound) )
         {
            int c;

            /* loop over children and stop if we find a child with a lower bound below given reference bound */
            for( c = 0; c < tree->nchildren; ++c )
            {
               if( SCIPsetIsLT(set, SCIPnodeGetLowerbound(tree->children[c]), stat->referencebound) )
               {
                  ++stat->nearlybacktracks;
                  break;
               }
            }
         }
         /* move children to the queue, make them LEAFs */
         SCIP_CALL( treeNodesToQueue(tree, reopt, blkmem, set, stat, eventqueue, lp, tree->children, &tree->nchildren, childrenlpstatefork,
               primal->cutoffbound) );

         /* remove node from the queue */
         SCIP_CALL( SCIPnodepqRemove(tree->leaves, set, *node) );

         stat->plungedepth = 0;
         if( SCIPnodeGetDepth(*node) > 0 )
            stat->nbacktracks++;
         SCIPsetDebugMsg(set, "selected leaf node, lowerbound=%g, plungedepth=%d\n", (*node)->lowerbound, stat->plungedepth);
         break;

      default:
         SCIPerrorMessage("selected node is neither sibling, child, nor leaf (nodetype=%d)\n", SCIPnodeGetType(*node));
         return SCIP_INVALIDDATA;
      }  /*lint !e788*/

      /* convert node into the focus node */
      (*node)->nodetype = SCIP_NODETYPE_FOCUSNODE; /*lint !e641*/
   }
   assert(tree->nchildren == 0);

   /* set new focus node, LP fork, LP state fork, and subroot */
   assert(subroot == NULL || (lpstatefork != NULL && subroot->depth <= lpstatefork->depth));
   assert(lpstatefork == NULL || (lpfork != NULL && lpstatefork->depth <= lpfork->depth));
   assert(lpfork == NULL || (*node != NULL && lpfork->depth < (*node)->depth));
   tree->focusnode = *node;
   tree->focuslpfork = lpfork;
   tree->focuslpstatefork = lpstatefork;
   tree->focussubroot = subroot;
   tree->focuslpconstructed = FALSE;
   lp->resolvelperror = FALSE;

   /* track the path from the old focus node to the new node, and perform domain and constraint set changes */
   SCIP_CALL( treeSwitchPath(tree, reopt, blkmem, set, stat, transprob, origprob, primal, lp, branchcand, conflict,
         eventfilter, eventqueue, cliquetable, fork, *node, cutoff) );
   assert(tree->pathlen >= 0);
   assert(*node != NULL || tree->pathlen == 0);
   assert(*node == NULL || tree->pathlen-1 <= (int)(*node)->depth);

   /* if the old focus node is a dead end (has no children), delete it */
   if( oldfocusnode != NULL && SCIPnodeGetType(oldfocusnode) == SCIP_NODETYPE_DEADEND )
   {
      int appliedeffectiverootdepth;

      appliedeffectiverootdepth = tree->appliedeffectiverootdepth;
      assert(appliedeffectiverootdepth <= tree->effectiverootdepth);

      SCIP_CALL( SCIPnodeFree(&oldfocusnode, blkmem, set, stat, eventqueue, tree, lp) );
      assert(tree->effectiverootdepth < tree->pathlen || *node == NULL || *cutoff);

      if( tree->effectiverootdepth > appliedeffectiverootdepth && *node != NULL && !(*cutoff) )
      {
         int d;

         /* promote the constraint set and bound changes up to the new effective root to be global changes */
         SCIPsetDebugMsg(set, "effective root is now at depth %d: applying constraint set and bound changes to global problem\n",
            tree->effectiverootdepth);

         for( d = appliedeffectiverootdepth + 1; d <= tree->effectiverootdepth; ++d )
         {
            SCIP_Bool nodecutoff;

            SCIPsetDebugMsg(set, " -> applying constraint set changes of depth %d\n", d);
            SCIP_CALL( SCIPconssetchgMakeGlobal(&tree->path[d]->conssetchg, blkmem, set, stat, transprob, reopt) );
            SCIPsetDebugMsg(set, " -> applying bound changes of depth %d\n", d);
            SCIP_CALL( SCIPdomchgApplyGlobal(tree->path[d]->domchg, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable,
                  &nodecutoff) );

            if( nodecutoff )
            {
               SCIP_CALL( SCIPnodeCutoff(tree->path[d], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
               *cutoff = TRUE;
            }
         }

         tree->appliedeffectiverootdepth = tree->effectiverootdepth;
      }
   }
   assert(*cutoff || SCIPtreeIsPathComplete(tree));

   return SCIP_OKAY;
}




/*
 * Tree methods
 */

/** creates an initialized tree data structure */
SCIP_RETCODE SCIPtreeCreate(
   SCIP_TREE**           tree,               /**< pointer to tree data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting leaves in the priority queue */
   )
{
   int p;

   assert(tree != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocMemory(tree) );

   (*tree)->root = NULL;

   SCIP_CALL( SCIPnodepqCreate(&(*tree)->leaves, set, nodesel) );

   /* allocate one slot for the prioritized and the unprioritized bound change */
   for( p = 0; p <= 1; ++p )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*tree)->divebdchgdirs[p], 1) ); /*lint !e866*/
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*tree)->divebdchgvars[p], 1) ); /*lint !e866*/
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*tree)->divebdchgvals[p], 1) ); /*lint !e866*/
      (*tree)->ndivebdchanges[p] = 0;
      (*tree)->divebdchgsize[p] = 1;
   }

   (*tree)->path = NULL;
   (*tree)->focusnode = NULL;
   (*tree)->focuslpfork = NULL;
   (*tree)->focuslpstatefork = NULL;
   (*tree)->focussubroot = NULL;
   (*tree)->children = NULL;
   (*tree)->siblings = NULL;
   (*tree)->probingroot = NULL;
   (*tree)->childrenprio = NULL;
   (*tree)->siblingsprio = NULL;
   (*tree)->pathnlpcols = NULL;
   (*tree)->pathnlprows = NULL;
   (*tree)->probinglpistate = NULL;
   (*tree)->probinglpinorms = NULL;
   (*tree)->pendingbdchgs = NULL;
   (*tree)->probdiverelaxsol = NULL;
   (*tree)->pendingbdchgssize = 0;
   (*tree)->npendingbdchgs = 0;
   (*tree)->focuslpstateforklpcount = -1;
   (*tree)->childrensize = 0;
   (*tree)->nchildren = 0;
   (*tree)->siblingssize = 0;
   (*tree)->nsiblings = 0;
   (*tree)->pathlen = 0;
   (*tree)->pathsize = 0;
   (*tree)->effectiverootdepth = 0;
   (*tree)->appliedeffectiverootdepth = 0;
   (*tree)->correctlpdepth = -1;
   (*tree)->cutoffdepth = INT_MAX;
   (*tree)->repropdepth = INT_MAX;
   (*tree)->repropsubtreecount = 0;
   (*tree)->focusnodehaslp = FALSE;
   (*tree)->probingnodehaslp = FALSE;
   (*tree)->focuslpconstructed = FALSE;
   (*tree)->cutoffdelayed = FALSE;
   (*tree)->probinglpwasflushed = FALSE;
   (*tree)->probinglpwassolved = FALSE;
   (*tree)->probingloadlpistate = FALSE;
   (*tree)->probinglpwasrelax = FALSE;
   (*tree)->probingsolvedlp = FALSE;
   (*tree)->forcinglpmessage = FALSE;
   (*tree)->sbprobing = FALSE;
   (*tree)->probinglpwasprimfeas = TRUE;
   (*tree)->probinglpwasdualfeas = TRUE;
   (*tree)->probdiverelaxstored = FALSE;
   (*tree)->probdiverelaxincludeslp = FALSE;

   return SCIP_OKAY;
}

/** frees tree data structure */
SCIP_RETCODE SCIPtreeFree(
   SCIP_TREE**           tree,               /**< pointer to tree data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int p;

   assert(tree != NULL);
   assert(*tree != NULL);
   assert((*tree)->nchildren == 0);
   assert((*tree)->nsiblings == 0);
   assert((*tree)->focusnode == NULL);
   assert(!SCIPtreeProbing(*tree));

   SCIPsetDebugMsg(set, "free tree\n");

   /* free node queue */
   SCIP_CALL( SCIPnodepqFree(&(*tree)->leaves, blkmem, set, stat, eventqueue, *tree, lp) );

   /* free diving bound change storage */
   for( p = 0; p <= 1; ++p )
   {
      BMSfreeBlockMemoryArray(blkmem, &(*tree)->divebdchgdirs[p], (*tree)->divebdchgsize[p]); /*lint !e866*/
      BMSfreeBlockMemoryArray(blkmem, &(*tree)->divebdchgvals[p], (*tree)->divebdchgsize[p]); /*lint !e866*/
      BMSfreeBlockMemoryArray(blkmem, &(*tree)->divebdchgvars[p], (*tree)->divebdchgsize[p]); /*lint !e866*/
   }

   /* free pointer arrays */
   BMSfreeMemoryArrayNull(&(*tree)->path);
   BMSfreeMemoryArrayNull(&(*tree)->children);
   BMSfreeMemoryArrayNull(&(*tree)->siblings);
   BMSfreeMemoryArrayNull(&(*tree)->childrenprio);
   BMSfreeMemoryArrayNull(&(*tree)->siblingsprio);
   BMSfreeMemoryArrayNull(&(*tree)->pathnlpcols);
   BMSfreeMemoryArrayNull(&(*tree)->pathnlprows);
   BMSfreeMemoryArrayNull(&(*tree)->probdiverelaxsol);
   BMSfreeMemoryArrayNull(&(*tree)->pendingbdchgs);

   BMSfreeMemory(tree);

   return SCIP_OKAY;
}

/** clears and resets tree data structure and deletes all nodes */
SCIP_RETCODE SCIPtreeClear(
   SCIP_TREE*            tree,               /**< tree data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int v;

   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(tree->nsiblings == 0);
   assert(tree->focusnode == NULL);
   assert(!SCIPtreeProbing(tree));

   SCIPsetDebugMsg(set, "clearing tree\n");

   /* clear node queue */
   SCIP_CALL( SCIPnodepqClear(tree->leaves, blkmem, set, stat, eventqueue, tree, lp) );
   assert(tree->root == NULL);

   /* we have to remove the captures of the variables within the pending bound change data structure */
   for( v = tree->npendingbdchgs-1; v >= 0; --v )
   {
      SCIP_VAR* var;

      var = tree->pendingbdchgs[v].var;
      assert(var != NULL);

      /* release the variable */
      SCIP_CALL( SCIPvarRelease(&var, blkmem, set, eventqueue, lp) );
   }

   /* mark working arrays to be empty and reset data */
   tree->focuslpstateforklpcount = -1;
   tree->nchildren = 0;
   tree->nsiblings = 0;
   tree->pathlen = 0;
   tree->effectiverootdepth = 0;
   tree->appliedeffectiverootdepth = 0;
   tree->correctlpdepth = -1;
   tree->cutoffdepth = INT_MAX;
   tree->repropdepth = INT_MAX;
   tree->repropsubtreecount = 0;
   tree->npendingbdchgs = 0;
   tree->focusnodehaslp = FALSE;
   tree->probingnodehaslp = FALSE;
   tree->cutoffdelayed = FALSE;
   tree->probinglpwasflushed = FALSE;
   tree->probinglpwassolved = FALSE;
   tree->probingloadlpistate = FALSE;
   tree->probinglpwasrelax = FALSE;
   tree->probingsolvedlp = FALSE;

   return SCIP_OKAY;
}

/** creates the root node of the tree and puts it into the leaves queue */
SCIP_RETCODE SCIPtreeCreateRoot(
   SCIP_TREE*            tree,               /**< tree data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(tree->nsiblings == 0);
   assert(tree->root == NULL);
   assert(tree->focusnode == NULL);
   assert(!SCIPtreeProbing(tree));

   /* create root node */
   SCIP_CALL( SCIPnodeCreateChild(&tree->root, blkmem, set, stat, tree, 0.0, -SCIPsetInfinity(set)) );
   assert(tree->nchildren == 1);

#ifndef NDEBUG
   /* check, if the sizes in the data structures match the maximal numbers defined here */
   tree->root->depth = SCIP_MAXTREEDEPTH + 1;
   tree->root->repropsubtreemark = MAXREPROPMARK;
   assert(tree->root->depth - 1 == SCIP_MAXTREEDEPTH); /*lint !e650*/
   assert(tree->root->repropsubtreemark == MAXREPROPMARK);
   tree->root->depth++;             /* this should produce an overflow and reset the value to 0 */
   tree->root->repropsubtreemark++; /* this should produce an overflow and reset the value to 0 */
   assert(tree->root->depth == 0);
   assert((SCIP_NODETYPE)tree->root->nodetype == SCIP_NODETYPE_CHILD);
   assert(!tree->root->active);
   assert(!tree->root->cutoff);
   assert(!tree->root->reprop);
   assert(tree->root->repropsubtreemark == 0);
#endif

   /* move root to the queue, convert it to LEAF */
   SCIP_CALL( treeNodesToQueue(tree, reopt, blkmem, set, stat, eventqueue, lp, tree->children, &tree->nchildren, NULL,
         SCIPsetInfinity(set)) );

   return SCIP_OKAY;
}

/** creates a temporary presolving root node of the tree and installs it as focus node */
SCIP_RETCODE SCIPtreeCreatePresolvingRoot(
   SCIP_TREE*            tree,               /**< tree data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   SCIP_Bool cutoff;

   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(tree->nsiblings == 0);
   assert(tree->root == NULL);
   assert(tree->focusnode == NULL);
   assert(!SCIPtreeProbing(tree));

   /* create temporary presolving root node */
   SCIP_CALL( SCIPtreeCreateRoot(tree, reopt, blkmem, set, stat, eventqueue, lp) );
   assert(tree->root != NULL);

   /* install the temporary root node as focus node */
   SCIP_CALL( SCIPnodeFocus(&tree->root, blkmem, set, messagehdlr, stat, transprob, origprob, primal, tree, reopt, lp, branchcand,
         conflict, conflictstore, eventfilter, eventqueue, cliquetable, &cutoff, FALSE, FALSE) );
   assert(!cutoff);

   return SCIP_OKAY;
}

/** frees the temporary presolving root and resets tree data structure */
SCIP_RETCODE SCIPtreeFreePresolvingRoot(
   SCIP_TREE*            tree,               /**< tree data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   SCIP_NODE* node;
   SCIP_Bool cutoff;

   assert(tree != NULL);
   assert(tree->root != NULL);
   assert(tree->focusnode == tree->root);
   assert(tree->pathlen == 1);

   /* unfocus the temporary root node */
   node = NULL;
   SCIP_CALL( SCIPnodeFocus(&node, blkmem, set, messagehdlr, stat, transprob, origprob, primal, tree, reopt, lp, branchcand,
         conflict, conflictstore, eventfilter, eventqueue, cliquetable, &cutoff, FALSE, FALSE) );
   assert(!cutoff);
   assert(tree->root == NULL);
   assert(tree->focusnode == NULL);
   assert(tree->pathlen == 0);

   /* reset tree data structure */
   SCIP_CALL( SCIPtreeClear(tree, blkmem, set, stat, eventqueue, lp) );

   return SCIP_OKAY;
}

/** returns the node selector associated with the given node priority queue */
SCIP_NODESEL* SCIPtreeGetNodesel(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqGetNodesel(tree->leaves);
}

/** sets the node selector used for sorting the nodes in the priority queue, and resorts the queue if necessary */
SCIP_RETCODE SCIPtreeSetNodesel(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   )
{
   assert(tree != NULL);
   assert(stat != NULL);

   if( SCIPnodepqGetNodesel(tree->leaves) != nodesel )
   {
      /* change the node selector used in the priority queue and resort the queue */
      SCIP_CALL( SCIPnodepqSetNodesel(&tree->leaves, set, nodesel) );

      /* issue message */
      if( stat->nnodes > 0 )
      {
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %" SCIP_LONGINT_FORMAT ") switching to node selector <%s>\n", stat->nnodes, SCIPnodeselGetName(nodesel));
      }
   }

   return SCIP_OKAY;
}

/** cuts off nodes with lower bound not better than given cutoff bound */
SCIP_RETCODE SCIPtreeCutoff(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   )
{
   SCIP_NODE* node;
   int i;

   assert(tree != NULL);
   assert(stat != NULL);
   assert(lp != NULL);

   /* if we are in diving mode, it is not allowed to cut off nodes, because this can lead to deleting LP rows which
    * would modify the currently unavailable (due to diving modifications) SCIP_LP
    *  -> the cutoff must be delayed and executed after the diving ends
    */
   if( SCIPlpDiving(lp) )
   {
      tree->cutoffdelayed = TRUE;
      return SCIP_OKAY;
   }

   tree->cutoffdelayed = FALSE;

   /* cut off leaf nodes in the queue */
   SCIP_CALL( SCIPnodepqBound(tree->leaves, blkmem, set, stat, eventqueue, tree, reopt, lp, cutoffbound) );

   /* cut off siblings: we have to loop backwards, because a removal leads to moving the last node in empty slot */
   for( i = tree->nsiblings-1; i >= 0; --i )
   {
      node = tree->siblings[i];
      if( SCIPsetIsGE(set, node->lowerbound, cutoffbound) )
      {
         SCIPsetDebugMsg(set, "cut off sibling #%" SCIP_LONGINT_FORMAT " at depth %d with lowerbound=%g at position %d\n",
            SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), node->lowerbound, i);

         if( set->reopt_enable )
         {
            assert(reopt != NULL);
            /* check if the node should be stored for reoptimization */
            SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, node, SCIP_EVENTTYPE_NODEINFEASIBLE, lp, SCIPlpGetSolstat(lp),
                  tree->root == node, tree->focusnode == node, node->lowerbound, tree->effectiverootdepth) );
         }

         SCIPvisualCutoffNode(stat->visual, set, stat, node, FALSE);

         SCIP_CALL( SCIPnodeFree(&node, blkmem, set, stat, eventqueue, tree, lp) );
      }
   }

   /* cut off children: we have to loop backwards, because a removal leads to moving the last node in empty slot */
   for( i = tree->nchildren-1; i >= 0; --i )
   {
      node = tree->children[i];
      if( SCIPsetIsGE(set, node->lowerbound, cutoffbound) )
      {
         SCIPsetDebugMsg(set, "cut off child #%" SCIP_LONGINT_FORMAT " at depth %d with lowerbound=%g at position %d\n",
            SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), node->lowerbound, i);

         if( set->reopt_enable )
         {
            assert(reopt != NULL);
            /* check if the node should be stored for reoptimization */
            SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, node, SCIP_EVENTTYPE_NODEINFEASIBLE, lp, SCIPlpGetSolstat(lp),
                  tree->root == node, tree->focusnode == node, node->lowerbound, tree->effectiverootdepth) );
         }

         SCIPvisualCutoffNode(stat->visual, set, stat, node, FALSE);

         SCIP_CALL( SCIPnodeFree(&node, blkmem, set, stat, eventqueue, tree, lp) );
      }
   }

   return SCIP_OKAY;
}

/** calculates the node selection priority for moving the given variable's LP value to the given target value;
 *  this node selection priority can be given to the SCIPcreateChild() call
 */
SCIP_Real SCIPtreeCalcNodeselPriority(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   SCIP_BRANCHDIR        branchdir,          /**< type of branching that was performed: upwards, downwards, or fixed 
                                              * fixed should only be used, when both bounds changed 
                                              */
   SCIP_Real             targetvalue         /**< new value of the variable in the child node */
   )
{
   SCIP_Real prio;
   SCIP_Real varsol;
   SCIP_Real varrootsol;
   SCIP_Real downinfs;
   SCIP_Real upinfs;
   SCIP_Bool isroot;
   SCIP_Bool haslp;

   assert(set != NULL);

   /* extract necessary information */
   isroot = (SCIPtreeGetCurrentDepth(tree) == 0);
   haslp = SCIPtreeHasFocusNodeLP(tree);
   varsol = SCIPvarGetSol(var, haslp);
   varrootsol = SCIPvarGetRootSol(var);
   downinfs = SCIPvarGetAvgInferences(var, stat, SCIP_BRANCHDIR_DOWNWARDS);
   upinfs = SCIPvarGetAvgInferences(var, stat, SCIP_BRANCHDIR_UPWARDS);

   switch( branchdir )
   {
   case SCIP_BRANCHDIR_DOWNWARDS:
      switch( SCIPvarGetBranchDirection(var) )
      {
      case SCIP_BRANCHDIR_DOWNWARDS:
         prio = +1.0;
         break;
      case SCIP_BRANCHDIR_UPWARDS:
         prio = -1.0;
         break;
      case SCIP_BRANCHDIR_AUTO:
         switch( set->nodesel_childsel )
         {
         case 'd':
            prio = +1.0;
            break;
         case 'u':
            prio = -1.0;
            break;
         case 'p':
            prio = -SCIPvarGetPseudocost(var, stat, targetvalue - varsol);
            break;
         case 'i':
            prio = downinfs;
            break;
         case 'l':
            prio = targetvalue - varsol;
            break;
         case 'r':
            prio = varrootsol - varsol;
            break;
         case 'h':
            prio = downinfs + SCIPsetEpsilon(set);
            if( !isroot && haslp )
               prio *= (varrootsol - varsol + 1.0);
            break;
         default:
            SCIPerrorMessage("invalid child selection rule <%c>\n", set->nodesel_childsel);
            prio = 0.0;
            break;
         }
         break;
      default:
         SCIPerrorMessage("invalid preferred branching direction <%d> of variable <%s>\n", 
            SCIPvarGetBranchDirection(var), SCIPvarGetName(var));
         prio = 0.0;
         break;
      }
      break;
   case SCIP_BRANCHDIR_UPWARDS:
      /* the branch is directed upwards */
      switch( SCIPvarGetBranchDirection(var) )
      {
      case SCIP_BRANCHDIR_DOWNWARDS:
         prio = -1.0;
         break;
      case SCIP_BRANCHDIR_UPWARDS:
         prio = +1.0;
         break;
      case SCIP_BRANCHDIR_AUTO:
         switch( set->nodesel_childsel )
         {
         case 'd':
            prio = -1.0;
            break;
         case 'u':
            prio = +1.0;
            break;
         case 'p':
            prio = -SCIPvarGetPseudocost(var, stat, targetvalue - varsol);
            break;
         case 'i':
            prio = upinfs;
            break;
         case 'l':
            prio = varsol - targetvalue;
            break;
         case 'r':
            prio = varsol - varrootsol;
            break;
         case 'h':
            prio = upinfs  + SCIPsetEpsilon(set);
            if( !isroot && haslp )
               prio *= (varsol - varrootsol + 1.0);
            break;
         default:
            SCIPerrorMessage("invalid child selection rule <%c>\n", set->nodesel_childsel);
            prio = 0.0;
            break;
         }
         /* since choosing the upwards direction is usually superior than the downwards direction (see results of
          * Achterberg's thesis (2007)), we break ties towards upwards branching
          */
         prio += SCIPsetEpsilon(set);
         break;

      default:
         SCIPerrorMessage("invalid preferred branching direction <%d> of variable <%s>\n", 
            SCIPvarGetBranchDirection(var), SCIPvarGetName(var));
         prio = 0.0;
         break;
      }
      break;
   case SCIP_BRANCHDIR_FIXED:
      prio = SCIPsetInfinity(set);
      break;
   case SCIP_BRANCHDIR_AUTO:
   default:
      SCIPerrorMessage("invalid branching direction <%d> of variable <%s>\n", 
         SCIPvarGetBranchDirection(var), SCIPvarGetName(var));
      prio = 0.0;
      break;
   }

   return prio;
}

/** calculates an estimate for the objective of the best feasible solution contained in the subtree after applying the given 
 *  branching; this estimate can be given to the SCIPcreateChild() call
 */
SCIP_Real SCIPtreeCalcChildEstimate(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   SCIP_Real             targetvalue         /**< new value of the variable in the child node */
   )
{
   SCIP_Real estimateinc;
   SCIP_Real estimate;
   SCIP_Real varsol;

   assert(tree != NULL);
   assert(var != NULL);

   estimate = SCIPnodeGetEstimate(tree->focusnode);
   varsol = SCIPvarGetSol(var, SCIPtreeHasFocusNodeLP(tree));

   /* compute increase above parent node's (i.e., focus node's) estimate value */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      estimateinc = SCIPvarGetPseudocost(var, stat, targetvalue - varsol);
   else
   {
      SCIP_Real pscdown;
      SCIP_Real pscup;

      /* calculate estimate based on pseudo costs:
       *   estimate = lowerbound + sum(min{f_j * pscdown_j, (1-f_j) * pscup_j})
       *            = parentestimate - min{f_b * pscdown_b, (1-f_b) * pscup_b} + (targetvalue-oldvalue)*{pscdown_b or pscup_b}
       */
      pscdown = SCIPvarGetPseudocost(var, stat, SCIPsetFeasFloor(set, varsol) - varsol);
      pscup = SCIPvarGetPseudocost(var, stat, SCIPsetFeasCeil(set, varsol) - varsol);
      estimateinc = SCIPvarGetPseudocost(var, stat, targetvalue - varsol) - MIN(pscdown, pscup);
   }

   /* due to rounding errors estimateinc might be slightly negative; in this case return the parent node's estimate */
   if( estimateinc > 0.0 )
      estimate += estimateinc;

   return estimate;
}

/** branches on a variable x
 *  if x is a continuous variable, then two child nodes will be created
 *  (x <= x', x >= x')
 *  but if the bounds of x are such that their relative difference is smaller than epsilon,
 *  the variable is fixed to val (if not SCIP_INVALID) or a well chosen alternative in the current node,
 *  i.e., no children are created
 *  if x is not a continuous variable, then:
 *  if solution value x' is fractional, two child nodes will be created
 *  (x <= floor(x'), x >= ceil(x')),
 *  if solution value is integral, the x' is equal to lower or upper bound of the branching
 *  variable and the bounds of x are finite, then two child nodes will be created
 *  (x <= x", x >= x"+1 with x" = floor((lb + ub)/2)),
 *  otherwise (up to) three child nodes will be created
 *  (x <= x'-1, x == x', x >= x'+1)
 *  if solution value is equal to one of the bounds and the other bound is infinite, only two child nodes
 *  will be created (the third one would be infeasible anyway)
 */
SCIP_RETCODE SCIPtreeBranchVar(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             val,                /**< value to branch on or SCIP_INVALID for branching on current LP/pseudo solution. 
                                              *   A branching value is required for branching on continuous variables */
   SCIP_NODE**           downchild,          /**< pointer to return the left child with variable rounded down, or NULL */
   SCIP_NODE**           eqchild,            /**< pointer to return the middle child with variable fixed, or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child with variable rounded up, or NULL */
   )
{
   SCIP_NODE* node;
   SCIP_Real priority;
   SCIP_Real estimate;

   SCIP_Real downub;
   SCIP_Real fixval;
   SCIP_Real uplb;
   SCIP_Real lpval;

   SCIP_Bool validval;

   assert(tree != NULL);
   assert(set != NULL);
   assert(var != NULL);

   /* initialize children pointer */
   if( downchild != NULL )
      *downchild = NULL;
   if( eqchild != NULL )
      *eqchild = NULL;
   if( upchild != NULL )
      *upchild = NULL;

   /* store whether a valid value was given for branching */
   validval = (val != SCIP_INVALID);  /*lint !e777 */

   /* get the corresponding active problem variable
    * if branching value is given, then transform it to the value of the active variable */
   if( validval )
   {
      SCIP_Real scalar;
      SCIP_Real constant;

      scalar   = 1.0;
      constant = 0.0;

      SCIP_CALL( SCIPvarGetProbvarSum(&var, set, &scalar, &constant) );

      if( scalar == 0.0 )
      {
         SCIPerrorMessage("cannot branch on fixed variable <%s>\n", SCIPvarGetName(var));
         return SCIP_INVALIDDATA;
      }

      /* we should have givenvariable = scalar * activevariable + constant */
      val = (val - constant) / scalar;
   }   
   else
      var = SCIPvarGetProbvar(var);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPerrorMessage("cannot branch on fixed or multi-aggregated variable <%s>\n", SCIPvarGetName(var));
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   /* ensure, that branching on continuous variables will only be performed when a branching point is given. */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS && !validval )
   {
      SCIPerrorMessage("Cannot branch on continuous variables without a given branching value.\n", SCIPvarGetName(var));
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   assert(SCIPvarIsActive(var));
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, SCIPvarGetLbLocal(var)));
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, SCIPvarGetUbLocal(var)));
   assert(SCIPsetIsLT(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

   /* update the information for the focus node before creating children */
   SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, tree->focusnode) );

   /* get value of variable in current LP or pseudo solution */
   lpval = SCIPvarGetSol(var, tree->focusnodehaslp);

   /* if there was no explicit value given for branching, branch on current LP or pseudo solution value */
   if( !validval )
   {
      val = lpval;

      /* avoid branching on infinite values in pseudo solution */
      if( SCIPsetIsInfinity(set, -val) || SCIPsetIsInfinity(set, val) )
      {
         val = SCIPvarGetWorstBoundLocal(var);

         /* if both bounds are infinite, choose zero as branching point */
         if( SCIPsetIsInfinity(set, -val) || SCIPsetIsInfinity(set, val) )
         {
            assert(SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)));
            assert(SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var)));         
            val = 0.0;
         }
      }
   }

   assert(SCIPsetIsFeasGE(set, val, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsFeasLE(set, val, SCIPvarGetUbLocal(var)));
   /* see comment in SCIPbranchVarVal */
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS ||
      SCIPrelDiff(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var)) <= 2.02 * SCIPsetEpsilon(set) ||
      SCIPsetIsInfinity(set, -2.1*SCIPvarGetLbLocal(var)) || SCIPsetIsInfinity(set, 2.1*SCIPvarGetUbLocal(var)) ||
      (SCIPsetIsLT(set, 2.1*SCIPvarGetLbLocal(var), 2.1*val) && SCIPsetIsLT(set, 2.1*val, 2.1*SCIPvarGetUbLocal(var))) );

   downub = SCIP_INVALID;
   fixval = SCIP_INVALID;
   uplb = SCIP_INVALID;

   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      if( SCIPsetIsRelEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
      {
         SCIPsetDebugMsg(set, "fixing continuous variable <%s> with value %g and bounds [%.15g, %.15g], priority %d (current lower bound: %g)\n",
            SCIPvarGetName(var), val, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetBranchPriority(var), SCIPnodeGetLowerbound(tree->focusnode));

         /* if val is at least epsilon away from both bounds, then we change both bounds to this value
          * otherwise, we fix the variable to its worst bound
          */
         if( SCIPsetIsGT(set, val, SCIPvarGetLbLocal(var)) && SCIPsetIsLT(set, val, SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(tree->focusnode, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                  branchcand, eventqueue, NULL, var, val, SCIP_BOUNDTYPE_LOWER, FALSE) );
            SCIP_CALL( SCIPnodeAddBoundchg(tree->focusnode, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                  branchcand, eventqueue, NULL, var, val, SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
         else if( SCIPvarGetObj(var) >= 0.0 )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(tree), blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, NULL, var, SCIPvarGetUbLocal(var), SCIP_BOUNDTYPE_LOWER, FALSE) );
         }
         else
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(tree), blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, NULL, var, SCIPvarGetLbLocal(var), SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
      }
      else if( SCIPrelDiff(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var)) <= 2.02 * SCIPsetEpsilon(set) )
      {
         /* if the only way to branch is such that in both sides the relative domain width becomes smaller epsilon,
          * then fix the variable in both branches right away
          *
          * however, if one of the bounds is at infinity (and thus the other bound is at most 2eps away from the same infinity (in relative sense),
          * then fix the variable to the non-infinite value, as we cannot fix a variable to infinity
          */
         SCIPsetDebugMsg(set, "continuous branch on variable <%s> with bounds [%.15g, %.15g], priority %d (current lower bound: %g), node %p\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetBranchPriority(var), SCIPnodeGetLowerbound(tree->focusnode), (void*)tree->focusnode);
         if( SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)) )
         {
            assert(!SCIPsetIsInfinity(set, -SCIPvarGetUbLocal(var)));
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(tree), blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, NULL, var, SCIPvarGetUbLocal(var), SCIP_BOUNDTYPE_LOWER, FALSE) );
         }
         else if( SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var)) )
         {
            assert(!SCIPsetIsInfinity(set, SCIPvarGetLbLocal(var)));
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(tree), blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, NULL, var, SCIPvarGetLbLocal(var), SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
         else
         {
            downub = SCIPvarGetLbLocal(var);
            uplb = SCIPvarGetUbLocal(var);
         }
      }
      else
      {
         /* in the general case, there is enough space for two branches
          * a sophisticated user should have also chosen the branching value such that it is not very close to the bounds
          * so here we only ensure that it is at least epsilon away from both bounds
          */
         SCIPsetDebugMsg(set, "continuous branch on variable <%s> with value %g, priority %d (current lower bound: %g)\n",
            SCIPvarGetName(var), val, SCIPvarGetBranchPriority(var), SCIPnodeGetLowerbound(tree->focusnode));
         downub = MIN(val, SCIPvarGetUbLocal(var) - SCIPsetEpsilon(set)); /*lint !e666*/
         uplb   = MAX(val, SCIPvarGetLbLocal(var) + SCIPsetEpsilon(set)); /*lint !e666*/
      }
   }
   else if( SCIPsetIsFeasIntegral(set, val) )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      /* if there was no explicit value given for branching, the variable has a finite domain and the current LP/pseudo
       * solution is one of the bounds, we branch in the center of the domain */
      if( !validval && !SCIPsetIsInfinity(set, -lb) && !SCIPsetIsInfinity(set, ub) 
         && (SCIPsetIsFeasEQ(set, val, lb) || SCIPsetIsFeasEQ(set, val, ub)) )
      {
         SCIP_Real center;

         /* create child nodes with x <= x", and x >= x"+1 with x" = floor((lb + ub)/2);
          * if x" is integral, make the interval smaller in the child in which the current solution x'
          * is still feasible
          */
         center = (ub + lb) / 2.0;
         if( val <= center )
         {
            downub = SCIPsetFeasFloor(set, center);
            uplb = downub + 1.0;
         }
         else
         {
            uplb = SCIPsetFeasCeil(set, center);
            downub = uplb - 1.0;
         }
      }
      else
      {
         /* create child nodes with x <= x'-1, x = x', and x >= x'+1 */
         assert(SCIPsetIsEQ(set, SCIPsetFeasCeil(set, val), SCIPsetFeasFloor(set, val)));

         fixval = SCIPsetFeasCeil(set, val); /* get rid of numerical issues */

         /* create child node with x <= x'-1, if this would be feasible */
         if( SCIPsetIsFeasGE(set, fixval-1.0, lb) )
            downub = fixval - 1.0;

         /* create child node with x >= x'+1, if this would be feasible */
         if( SCIPsetIsFeasLE(set, fixval+1.0, ub) )
            uplb = fixval + 1.0;
      }
      SCIPsetDebugMsg(set, "integral branch on variable <%s> with value %g, priority %d (current lower bound: %g)\n",
         SCIPvarGetName(var), val, SCIPvarGetBranchPriority(var), SCIPnodeGetLowerbound(tree->focusnode));
   }
   else
   {
      /* create child nodes with x <= floor(x'), and x >= ceil(x') */
      downub = SCIPsetFeasFloor(set, val);
      uplb = downub + 1.0;
      assert( SCIPsetIsRelEQ(set, SCIPsetCeil(set, val), uplb) );
      SCIPsetDebugMsg(set, "fractional branch on variable <%s> with value %g, root value %g, priority %d (current lower bound: %g)\n",
         SCIPvarGetName(var), val, SCIPvarGetRootSol(var), SCIPvarGetBranchPriority(var), SCIPnodeGetLowerbound(tree->focusnode));
   }

   /* perform the branching;
    * set the node selection priority in a way, s.t. a node is preferred whose branching goes in the same direction
    * as the deviation from the variable's root solution
    */
   if( downub != SCIP_INVALID )    /*lint !e777*/
   {
      /* create child node x <= downub */
      priority = SCIPtreeCalcNodeselPriority(tree, set, stat, var, SCIP_BRANCHDIR_DOWNWARDS, downub);
      /* if LP solution is cutoff in child, compute a new estimate
       * otherwise we cannot expect a direct change in the best solution, so we keep the estimate of the parent node */
      if( SCIPsetIsGT(set, lpval, downub) )
         estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, downub);
      else
         estimate = SCIPnodeGetEstimate(tree->focusnode);
      SCIPsetDebugMsg(set, " -> creating child: <%s> <= %g (priority: %g, estimate: %g)\n",
         SCIPvarGetName(var), downub, priority, estimate);
      SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, priority, estimate) );
      SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
            NULL, var, downub, SCIP_BOUNDTYPE_UPPER, FALSE) );
      /* output branching bound change to visualization file */
      SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, node) );

      if( downchild != NULL )
         *downchild = node;
   }

   if( fixval != SCIP_INVALID )    /*lint !e777*/
   {
      /* create child node with x = fixval */
      priority = SCIPtreeCalcNodeselPriority(tree, set, stat, var, SCIP_BRANCHDIR_FIXED, fixval);
      estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, fixval);
      SCIPsetDebugMsg(set, " -> creating child: <%s> == %g (priority: %g, estimate: %g)\n",
         SCIPvarGetName(var), fixval, priority, estimate);
      SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, priority, estimate) );
      if( !SCIPsetIsFeasEQ(set, SCIPvarGetLbLocal(var), fixval) )
      {
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               NULL, var, fixval, SCIP_BOUNDTYPE_LOWER, FALSE) );
      }
      if( !SCIPsetIsFeasEQ(set, SCIPvarGetUbLocal(var), fixval) )
      {
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               NULL, var, fixval, SCIP_BOUNDTYPE_UPPER, FALSE) );
      }
      /* output branching bound change to visualization file */
      SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, node) );

      if( eqchild != NULL )
         *eqchild = node;
   }

   if( uplb != SCIP_INVALID )    /*lint !e777*/
   {
      /* create child node with x >= uplb */
      priority = SCIPtreeCalcNodeselPriority(tree, set, stat, var, SCIP_BRANCHDIR_UPWARDS, uplb);
      if( SCIPsetIsLT(set, lpval, uplb) )
         estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, uplb);
      else
         estimate = SCIPnodeGetEstimate(tree->focusnode);
      SCIPsetDebugMsg(set, " -> creating child: <%s> >= %g (priority: %g, estimate: %g)\n",
         SCIPvarGetName(var), uplb, priority, estimate);
      SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, priority, estimate) );
      SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
            NULL, var, uplb, SCIP_BOUNDTYPE_LOWER, FALSE) );
      /* output branching bound change to visualization file */
      SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, node) );

      if( upchild != NULL )
         *upchild = node;
   }


   return SCIP_OKAY;
}

/** branches a variable x using the given domain hole; two child nodes will be created (x <= left, x >= right) */
SCIP_RETCODE SCIPtreeBranchVarHole(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             left,               /**< left side of the domain hole */
   SCIP_Real             right,              /**< right side of the domain hole */
   SCIP_NODE**           downchild,          /**< pointer to return the left child with variable rounded down, or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child with variable rounded up, or NULL */
   )
{
   SCIP_NODE* node;
   SCIP_Real priority;
   SCIP_Real estimate;
   SCIP_Real lpval;

   assert(tree != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPsetIsLT(set, left, SCIPvarGetUbLocal(var)));
   assert(SCIPsetIsGE(set, left, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsGT(set, right, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsLE(set, right, SCIPvarGetUbLocal(var)));
   assert(SCIPsetIsLE(set, left, right));

   /* initialize children pointer */
   if( downchild != NULL )
      *downchild = NULL;
   if( upchild != NULL )
      *upchild = NULL;

   /* get the corresponding active problem variable */
   SCIP_CALL( SCIPvarGetProbvarHole(&var, &left, &right) );

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPerrorMessage("cannot branch on fixed or multi-aggregated variable <%s>\n", SCIPvarGetName(var));
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   assert(SCIPvarIsActive(var));
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, SCIPvarGetLbLocal(var)));
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, SCIPvarGetUbLocal(var)));
   assert(SCIPsetIsLT(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

   assert(SCIPsetIsFeasGE(set, left, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsFeasLE(set, right, SCIPvarGetUbLocal(var)));

   /* adjust left and right side of the domain hole if the variable is integral */
   if( SCIPvarIsIntegral(var) )
   {
      left = SCIPsetFeasFloor(set, left);
      right = SCIPsetFeasCeil(set, right);
   }

   assert(SCIPsetIsLT(set, left, SCIPvarGetUbLocal(var)));
   assert(SCIPsetIsGE(set, left, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsGT(set, right, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsLE(set, right, SCIPvarGetUbLocal(var)));
   assert(SCIPsetIsLE(set, left, right));

   /* get value of variable in current LP or pseudo solution */
   lpval = SCIPvarGetSol(var, tree->focusnodehaslp);

   /* perform the branching;
    * set the node selection priority in a way, s.t. a node is preferred whose branching goes in the same direction
    * as the deviation from the variable's root solution
    */

   /* create child node x <= left */
   priority = SCIPtreeCalcNodeselPriority(tree, set, stat, var, SCIP_BRANCHDIR_DOWNWARDS, left);

   /* if LP solution is cutoff in child, compute a new estimate
    * otherwise we cannot expect a direct change in the best solution, so we keep the estimate of the parent node
    */
   if( SCIPsetIsGT(set, lpval, left) )
      estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, left);
   else
      estimate = SCIPnodeGetEstimate(tree->focusnode);

   SCIPsetDebugMsg(set, " -> creating child: <%s> <= %g (priority: %g, estimate: %g)\n",
      SCIPvarGetName(var), left, priority, estimate);

   SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, priority, estimate) );
   SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, NULL,
         var, left, SCIP_BOUNDTYPE_UPPER, FALSE) );
   /* output branching bound change to visualization file */
   SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, node) );

   if( downchild != NULL )
      *downchild = node;

   /* create child node with x >= right */
   priority = SCIPtreeCalcNodeselPriority(tree, set, stat, var, SCIP_BRANCHDIR_UPWARDS, right);

   if( SCIPsetIsLT(set, lpval, right) )
      estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, right);
   else
      estimate = SCIPnodeGetEstimate(tree->focusnode);

   SCIPsetDebugMsg(set, " -> creating child: <%s> >= %g (priority: %g, estimate: %g)\n",
      SCIPvarGetName(var), right, priority, estimate);

   SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, priority, estimate) );
   SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
         NULL, var, right, SCIP_BOUNDTYPE_LOWER, FALSE) );
   /* output branching bound change to visualization file */
   SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, node) );

   if( upchild != NULL )
      *upchild = node;

   return SCIP_OKAY;
}

/** n-ary branching on a variable x
 * Branches on variable x such that up to n/2 children are created on each side of the usual branching value.
 * The branching value is selected as in SCIPtreeBranchVar().
 * If n is 2 or the variables local domain is too small for a branching into n pieces, SCIPtreeBranchVar() is called.
 * The parameters minwidth and widthfactor determine the domain width of the branching variable in the child nodes.
 * If n is odd, one child with domain width 'width' and having the branching value in the middle is created.
 * Otherwise, two children with domain width 'width' and being left and right of the branching value are created.
 * Next further nodes to the left and right are created, where width is multiplied by widthfactor with increasing distance from the first nodes.
 * The initial width is calculated such that n/2 nodes are created to the left and to the right of the branching value.
 * If this value is below minwidth, the initial width is set to minwidth, which may result in creating less than n nodes.
 *
 * Giving a large value for widthfactor results in creating children with small domain when close to the branching value
 * and large domain when closer to the current variable bounds. That is, setting widthfactor to a very large value and n to 3
 * results in a ternary branching where the branching variable is mostly fixed in the middle child.
 * Setting widthfactor to 1.0 results in children where the branching variable always has the same domain width
 * (except for one child if the branching value is not in the middle).
 */
SCIP_RETCODE SCIPtreeBranchVarNary(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             val,                /**< value to branch on or SCIP_INVALID for branching on current LP/pseudo solution.
                                              *   A branching value is required for branching on continuous variables */
   int                   n,                  /**< attempted number of children to be created, must be >= 2 */
   SCIP_Real             minwidth,           /**< minimal domain width in children */
   SCIP_Real             widthfactor,        /**< multiplier for children domain width with increasing distance from val, must be >= 1.0 */
   int*                  nchildren           /**< buffer to store number of created children, or NULL */
   )
{
   SCIP_NODE* node;
   SCIP_Real priority;
   SCIP_Real estimate;
   SCIP_Real lpval;
   SCIP_Real width;
   SCIP_Bool validval;
   SCIP_Real left;
   SCIP_Real right;
   SCIP_Real bnd;
   int i;

   assert(tree != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(n >= 2);
   assert(minwidth >= 0.0);

   /* if binary branching is requested or we have not enough space for n children, delegate to SCIPtreeBranchVar */
   if( n == 2 ||
      2.0 * minwidth >= SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var) ||
      SCIPrelDiff(SCIPvarGetUbLocal(SCIPvarGetProbvar(var)), SCIPvarGetLbLocal(SCIPvarGetProbvar(var))) <= n * SCIPsetEpsilon(set) )
   {
      SCIP_NODE* downchild;
      SCIP_NODE* fixchild;
      SCIP_NODE* upchild;

      SCIP_CALL( SCIPtreeBranchVar(tree, reopt, blkmem, set, stat, transprob, origprob, lp, branchcand, eventqueue, var, val,
            &downchild, &fixchild, &upchild) );

      if( nchildren != NULL )
         *nchildren = (downchild != NULL ? 1 : 0) + (fixchild != NULL ? 1 : 0) + (upchild != NULL ? 1 : 0);

      return SCIP_OKAY;
   }

   /* store whether a valid value was given for branching */
   validval = (val != SCIP_INVALID);  /*lint !e777 */

   /* get the corresponding active problem variable
    * if branching value is given, then transform it to the value of the active variable */
   if( validval )
   {
      SCIP_Real scalar;
      SCIP_Real constant;

      scalar   = 1.0;
      constant = 0.0;

      SCIP_CALL( SCIPvarGetProbvarSum(&var, set, &scalar, &constant) );

      if( scalar == 0.0 )
      {
         SCIPerrorMessage("cannot branch on fixed variable <%s>\n", SCIPvarGetName(var));
         return SCIP_INVALIDDATA;
      }

      /* we should have givenvariable = scalar * activevariable + constant */
      val = (val - constant) / scalar;
   }
   else
      var = SCIPvarGetProbvar(var);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPerrorMessage("cannot branch on fixed or multi-aggregated variable <%s>\n", SCIPvarGetName(var));
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   /* ensure, that branching on continuous variables will only be performed when a branching point is given. */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS && !validval )
   {
      SCIPerrorMessage("Cannot branch on continuous variables without a given branching value.\n", SCIPvarGetName(var));
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   assert(SCIPvarIsActive(var));
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, SCIPvarGetLbLocal(var)));
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, SCIPvarGetUbLocal(var)));
   assert(SCIPsetIsLT(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

   /* get value of variable in current LP or pseudo solution */
   lpval = SCIPvarGetSol(var, tree->focusnodehaslp);

   /* if there was no explicit value given for branching, branch on current LP or pseudo solution value */
   if( !validval )
   {
      val = lpval;

      /* avoid branching on infinite values in pseudo solution */
      if( SCIPsetIsInfinity(set, -val) || SCIPsetIsInfinity(set, val) )
      {
         val = SCIPvarGetWorstBoundLocal(var);

         /* if both bounds are infinite, choose zero as branching point */
         if( SCIPsetIsInfinity(set, -val) || SCIPsetIsInfinity(set, val) )
         {
            assert(SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)));
            assert(SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var)));
            val = 0.0;
         }
      }
   }

   assert(SCIPsetIsFeasGE(set, val, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsFeasLE(set, val, SCIPvarGetUbLocal(var)));
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS ||
      SCIPsetIsRelEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) ||
      (SCIPsetIsLT(set, 2.1*SCIPvarGetLbLocal(var), 2.1*val) && SCIPsetIsLT(set, 2.1*val, 2.1*SCIPvarGetUbLocal(var))) );  /* see comment in SCIPbranchVarVal */

   /* calculate minimal distance of val from bounds */
   width = SCIP_REAL_MAX;
   if( !SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)) )
   {
      width = val - SCIPvarGetLbLocal(var);
   }
   if( !SCIPsetIsInfinity(set,  SCIPvarGetUbLocal(var)) )
   {
      width = MIN(width, SCIPvarGetUbLocal(var) - val); /*lint !e666*/
   }
   /* calculate initial domain width of child nodes
    * if we have at least one finite bound, choose width such that we have roughly the same number of nodes left and right of val
    */
   if( width == SCIP_REAL_MAX ) /*lint !e777*/
   {
      /* unbounded variable, let's create a child with a small domain */
      width = 1.0;
   }
   else if( widthfactor == 1.0 )
   {
      /* most domains get same size */
      width /= n/2; /*lint !e653*/ /* rounding is ok at this point */
   }
   else
   {
      /* width is increased by widthfactor for each child
       * if n is even, compute width such that we can create n/2 nodes with width
       * width, widthfactor*width, ..., widthfactor^(n/2)*width on each side, i.e.,
       *      sum(width * widthfactor^(i-1), i = 1..n/2) = min(ub-val, val-lb)
       *  <-> width * (widthfactor^(n/2) - 1) / (widthfactor - 1) = min(ub-val, val-lb)
       *
       * if n is odd, compute width such that we can create one middle node with width width
       * and n/2 nodes with width widthfactor*width, ..., widthfactor^(n/2)*width on each side, i.e.,
       *      width/2 + sum(width * widthfactor^i, i = 1..n/2) = min(ub-val, val-lb)
       *  <-> width * (1/2 + widthfactor * (widthfactor^(n/2) - 1) / (widthfactor - 1) = min(ub-val, val-lb)
       */
      assert(widthfactor > 1.0);
      if( n % 2 == 0 )
         width *= (widthfactor - 1.0) / (pow(widthfactor, (SCIP_Real)(n/2)) - 1.0); /*lint !e653*/
      else
         width /= 0.5 + widthfactor * (pow(widthfactor, (SCIP_Real)(n/2)) - 1.0) / (widthfactor - 1.0); /*lint !e653*/
   }
   if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      minwidth = MAX(1.0, minwidth);
   if( width < minwidth )
      width = minwidth;
   assert(SCIPsetIsPositive(set, width));

   SCIPsetDebugMsg(set, "%d-ary branching on variable <%s> [%g, %g] around %g, initial width = %g\n",
      n, SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), val, width);

   if( nchildren != NULL )
      *nchildren = 0;

   /* initialize upper bound on children left of val and children right of val
    * if we are supposed to create an odd number of children, then create a child that has val in the middle of its domain */
   if( n % 2 == 1 )
   {
      left  = val - width/2.0;
      right = val + width/2.0;
      SCIPvarAdjustLb(var, set, &left);
      SCIPvarAdjustUb(var, set, &right);

      /* create child node left <= x <= right, if left <= right */
      if( left <= right )
      {
         priority = SCIPtreeCalcNodeselPriority(tree, set, stat, var, SCIP_BRANCHDIR_FIXED, val); /* ????????????? how to compute priority for such a child? */
         /* if LP solution is cutoff in child, compute a new estimate
          * otherwise we cannot expect a direct change in the best solution, so we keep the estimate of the parent node */
         if( SCIPsetIsLT(set, lpval, left) )
            estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, left);
         else if( SCIPsetIsGT(set, lpval, right) )
            estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, right);
         else
            estimate = SCIPnodeGetEstimate(tree->focusnode);

         SCIPsetDebugMsg(set, " -> creating middle child: %g <= <%s> <= %g (priority: %g, estimate: %g, width: %g)\n",
            left, SCIPvarGetName(var), right, priority, estimate, right - left);

         SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, priority, estimate) );
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
               eventqueue, NULL, var, left , SCIP_BOUNDTYPE_LOWER, FALSE) );
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               NULL, var, right, SCIP_BOUNDTYPE_UPPER, FALSE) );
         /* output branching bound change to visualization file */
         SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, node) );

         if( nchildren != NULL )
            ++*nchildren;
      }
      --n;

      if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      {
         /* if it's a discrete variable, we can use left-1 and right+1 as upper and lower bounds for following nodes on the left and right, resp. */
         left  -= 1.0;
         right += 1.0;
      }

      width *= widthfactor;
   }
   else
   {
      if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      {
         left  = SCIPsetFloor(set, val);
         right = SCIPsetCeil(set, val);
         if( right - left < 0.5 )
            left -= 1.0;
      }
      else if( SCIPsetIsZero(set, val) )
      {
         left  = 0.0;
         right = 0.0;
      }
      else
      {
         left  = val;
         right = val;
      }
   }

   assert(n % 2 == 0);
   n /= 2;
   for( i = 0; i < n; ++i )
   {
      /* create child node left - width <= x <= left, if left > lb(x) or x is discrete */
      if( SCIPsetIsRelLT(set, SCIPvarGetLbLocal(var), left) || SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      {
         /* new lower bound should be variables lower bound, if we are in the last round or left - width is very close to lower bound
          * otherwise we take left - width
          */
         if( i == n-1 || SCIPsetIsRelEQ(set, SCIPvarGetLbLocal(var), left - width))
         {
            bnd = SCIPvarGetLbLocal(var);
         }
         else
         {
            bnd = left - width;
            SCIPvarAdjustLb(var, set, &bnd);
            bnd = MAX(SCIPvarGetLbLocal(var), bnd); /*lint !e666*/
         }
         assert(SCIPsetIsRelLT(set, bnd, left));

         /* the nodeselection priority of nodes is decreased as more as they are away from val */
         priority = SCIPtreeCalcNodeselPriority(tree, set, stat, var, SCIP_BRANCHDIR_DOWNWARDS, bnd) / (i+1);
         /* if LP solution is cutoff in child, compute a new estimate
          * otherwise we cannot expect a direct change in the best solution, so we keep the estimate of the parent node */
         if( SCIPsetIsLT(set, lpval, bnd) )
            estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, bnd);
         else if( SCIPsetIsGT(set, lpval, left) )
            estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, left);
         else
            estimate = SCIPnodeGetEstimate(tree->focusnode);

         SCIPsetDebugMsg(set, " -> creating left  child: %g <= <%s> <= %g (priority: %g, estimate: %g, width: %g)\n",
            bnd, SCIPvarGetName(var), left, priority, estimate, left - bnd);

         SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, priority, estimate) );
         if( SCIPsetIsGT(set, bnd, SCIPvarGetLbLocal(var)) )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               NULL, var, bnd, SCIP_BOUNDTYPE_LOWER, FALSE) );
         }
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
            NULL, var, left, SCIP_BOUNDTYPE_UPPER, FALSE) );
         /* output branching bound change to visualization file */
         SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, node) );

         if( nchildren != NULL )
            ++*nchildren;

         left = bnd;
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
            left -= 1.0;
      }

      /* create child node right <= x <= right + width, if right < ub(x) */
      if( SCIPsetIsRelGT(set, SCIPvarGetUbLocal(var), right) || SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      {
         /* new upper bound should be variables upper bound, if we are in the last round or right + width is very close to upper bound
          * otherwise we take right + width
          */
         if( i == n-1 || SCIPsetIsRelEQ(set, SCIPvarGetUbLocal(var), right + width))
         {
            bnd = SCIPvarGetUbLocal(var);
         }
         else
         {
            bnd = right + width;
            SCIPvarAdjustUb(var, set, &bnd);
            bnd = MIN(SCIPvarGetUbLocal(var), bnd); /*lint !e666*/
         }
         assert(SCIPsetIsRelGT(set, bnd, right));

         /* the nodeselection priority of nodes is decreased as more as they are away from val */
         priority = SCIPtreeCalcNodeselPriority(tree, set, stat, var, SCIP_BRANCHDIR_UPWARDS, bnd) / (i+1);
         /* if LP solution is cutoff in child, compute a new estimate
          * otherwise we cannot expect a direct change in the best solution, so we keep the estimate of the parent node */
         if( SCIPsetIsLT(set, lpval, right) )
            estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, right);
         else if( SCIPsetIsGT(set, lpval, bnd) )
            estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, bnd);
         else
            estimate = SCIPnodeGetEstimate(tree->focusnode);

         SCIPsetDebugMsg(set, " -> creating right child: %g <= <%s> <= %g (priority: %g, estimate: %g, width: %g)\n",
            right, SCIPvarGetName(var), bnd, priority, estimate, bnd - right);

         SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, priority, estimate) );
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
            NULL, var, right, SCIP_BOUNDTYPE_LOWER, FALSE) );
         if( SCIPsetIsLT(set, bnd, SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               NULL, var, bnd, SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
         /* output branching bound change to visualization file */
         SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, node) );

         if( nchildren != NULL )
            ++*nchildren;

         right = bnd;
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
            right += 1.0;
      }

      width *= widthfactor;
   }

   return SCIP_OKAY;
}

/** adds a diving bound change to the tree together with the information if this is a bound change
 *  for the preferred direction or not
 */
#define ARRAYGROWTH 5
SCIP_RETCODE SCIPtreeAddDiveBoundChange(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_VAR*             var,                /**< variable to apply the bound change to */
   SCIP_BRANCHDIR        dir,                /**< direction of the bound change */
   SCIP_Real             value,              /**< value to adjust this variable bound to */
   SCIP_Bool             preferred           /**< is this a bound change for the preferred child? */
   )
{
   int idx = preferred ? 0 : 1;
   int pos = tree->ndivebdchanges[idx];

   assert(pos < tree->divebdchgsize[idx]);

   if( pos == tree->divebdchgsize[idx] - 1 )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &tree->divebdchgdirs[idx], tree->divebdchgsize[idx], tree->divebdchgsize[idx] + ARRAYGROWTH) ); /*lint !e866*/
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &tree->divebdchgvars[idx], tree->divebdchgsize[idx], tree->divebdchgsize[idx] + ARRAYGROWTH) ); /*lint !e866*/
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &tree->divebdchgvals[idx], tree->divebdchgsize[idx], tree->divebdchgsize[idx] + ARRAYGROWTH) ); /*lint !e866*/
      tree->divebdchgsize[idx] += ARRAYGROWTH;
   }

   tree->divebdchgvars[idx][pos] = var;
   tree->divebdchgdirs[idx][pos] = dir;
   tree->divebdchgvals[idx][pos] = value;

   ++tree->ndivebdchanges[idx];

   return SCIP_OKAY;
}

/** get the dive bound change data for the preferred or the alternative direction */
void SCIPtreeGetDiveBoundChangeData(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_VAR***           variables,          /**< pointer to store variables for the specified direction */
   SCIP_BRANCHDIR**      directions,         /**< pointer to store the branching directions */
   SCIP_Real**           values,             /**< pointer to store bound change values */
   int*                  ndivebdchgs,        /**< pointer to store the number of dive bound changes */
   SCIP_Bool             preferred           /**< should the dive bound changes for the preferred child be output? */
   )
{
   int idx = preferred ? 0 : 1;

   assert(variables != NULL);
   assert(directions != NULL);
   assert(values != NULL);
   assert(ndivebdchgs != NULL);

   *variables = tree->divebdchgvars[idx];
   *directions = tree->divebdchgdirs[idx];
   *values = tree->divebdchgvals[idx];
   *ndivebdchgs = tree->ndivebdchanges[idx];
}

/** clear the tree bound change data structure */
void SCIPtreeClearDiveBoundChanges(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   int p;

   for( p = 0; p < 2; ++p )
      tree->ndivebdchanges[p] = 0;
}

/** creates a probing child node of the current node, which must be the focus node, the current refocused node,
 *  or another probing node; if the current node is the focus or a refocused node, the created probing node is
 *  installed as probing root node
 */
static
SCIP_RETCODE treeCreateProbingNode(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_NODE* currentnode;
   SCIP_NODE* node;
   SCIP_RETCODE retcode;

   assert(tree != NULL);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->pathlen > 0);
   assert(blkmem != NULL);
   assert(set != NULL);

   /* get the current node */
   currentnode = SCIPtreeGetCurrentNode(tree);
   assert(SCIPnodeGetType(currentnode) == SCIP_NODETYPE_FOCUSNODE
      || SCIPnodeGetType(currentnode) == SCIP_NODETYPE_REFOCUSNODE
      || SCIPnodeGetType(currentnode) == SCIP_NODETYPE_PROBINGNODE);
   assert((SCIPnodeGetType(currentnode) == SCIP_NODETYPE_PROBINGNODE) == SCIPtreeProbing(tree));

   /* create the node data structure */
   SCIP_CALL( nodeCreate(&node, blkmem, set) );
   assert(node != NULL);

   /* mark node to be a probing node */
   node->nodetype = SCIP_NODETYPE_PROBINGNODE; /*lint !e641*/

   /* create the probingnode data */
   SCIP_CALL( probingnodeCreate(&node->data.probingnode, blkmem, lp) );

   /* make the current node the parent of the new probing node */
   retcode = nodeAssignParent(node, blkmem, set, tree, currentnode, 0.0);

   /* if we reached the maximal depth level we clean up the allocated memory and stop */
   if( retcode == SCIP_MAXDEPTHLEVEL )
   {
      SCIP_CALL( probingnodeFree(&(node->data.probingnode), blkmem, lp) );
      BMSfreeBlockMemory(blkmem, &node);
   }
   SCIP_CALL( retcode );
   assert(SCIPnodeGetDepth(node) == tree->pathlen);

   /* check, if the node is the probing root node */
   if( tree->probingroot == NULL )
   {
      tree->probingroot = node;
      SCIPsetDebugMsg(set, "created probing root node #%" SCIP_LONGINT_FORMAT " at depth %d\n",
         SCIPnodeGetNumber(node), SCIPnodeGetDepth(node));
   }
   else
   {
      assert(SCIPnodeGetType(tree->probingroot) == SCIP_NODETYPE_PROBINGNODE);
      assert(SCIPnodeGetDepth(tree->probingroot) < SCIPnodeGetDepth(node));

      SCIPsetDebugMsg(set, "created probing child node #%" SCIP_LONGINT_FORMAT " at depth %d, probing depth %d\n",
         SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), SCIPnodeGetDepth(node) - SCIPnodeGetDepth(tree->probingroot));
   }

   /* create the new active path */
   SCIP_CALL( treeEnsurePathMem(tree, set, tree->pathlen+1) );
   node->active = TRUE;
   tree->path[tree->pathlen] = node;
   tree->pathlen++;

   /* update the path LP size for the previous node and set the (initial) path LP size for the newly created node */
   SCIP_CALL( treeUpdatePathLPSize(tree, tree->pathlen-2) );

   /* mark the LP's size */
   SCIPlpMarkSize(lp);
   assert(tree->pathlen >= 2);
   assert(lp->firstnewrow == tree->pathnlprows[tree->pathlen-1]); /* marked LP size should be initial size of new node */
   assert(lp->firstnewcol == tree->pathnlpcols[tree->pathlen-1]);

   /* the current probing node does not yet have a solved LP */
   tree->probingnodehaslp = FALSE;

   return SCIP_OKAY;
}

/** switches to probing mode and creates a probing root */
SCIP_RETCODE SCIPtreeStartProbing(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_Bool             strongbranching     /**< is the probing mode used for strongbranching? */
   )
{
   assert(tree != NULL);
   assert(tree->probinglpistate == NULL);
   assert(tree->probinglpinorms == NULL);
   assert(!SCIPtreeProbing(tree));
   assert(lp != NULL);

   SCIPsetDebugMsg(set, "probing started in depth %d (LP flushed: %u, LP solved: %u, solstat: %d), probing root in depth %d\n",
      tree->pathlen-1, lp->flushed, lp->solved, SCIPlpGetSolstat(lp), tree->pathlen);

   /* store all marked constraints for propagation */
   SCIP_CALL( SCIPconshdlrsStorePropagationStatus(set, set->conshdlrs, set->nconshdlrs) );

   /* inform LP about probing mode */
   SCIP_CALL( SCIPlpStartProbing(lp) );

   assert(!lp->divingobjchg);

   /* remember, whether the LP was flushed and solved */
   tree->probinglpwasflushed = lp->flushed;
   tree->probinglpwassolved = lp->solved;
   tree->probingloadlpistate = FALSE;
   tree->probinglpwasrelax = lp->isrelax;
   lp->isrelax = TRUE;
   tree->probingsolvedlp = FALSE;
   tree->probingobjchanged = FALSE;
   lp->divingobjchg = FALSE;
   tree->probingsumchgdobjs = 0;
   tree->sbprobing = strongbranching;

   /* remember the LP state in order to restore the LP solution quickly after probing */
   /**@todo could the lp state be worth storing if the LP is not flushed (and hence not solved)? */
   if( lp->flushed && lp->solved )
   {
      SCIP_CALL( SCIPlpGetState(lp, blkmem, &tree->probinglpistate) );
      SCIP_CALL( SCIPlpGetNorms(lp, blkmem, &tree->probinglpinorms) );
      tree->probinglpwasprimfeas = lp->primalfeasible;
      tree->probinglpwasprimchecked = lp->primalchecked;
      tree->probinglpwasdualfeas = lp->dualfeasible;
      tree->probinglpwasdualchecked = lp->dualchecked;
   }

   /* remember the relaxation solution to reset it later */
   if( SCIPrelaxationIsSolValid(relaxation) )
   {
      SCIP_CALL( SCIPtreeStoreRelaxSol(tree, set, relaxation, transprob) );
   }

   /* create temporary probing root node */
   SCIP_CALL( treeCreateProbingNode(tree, blkmem, set, lp) );
   assert(SCIPtreeProbing(tree));

   return SCIP_OKAY;
}

/** creates a new probing child node in the probing path */
SCIP_RETCODE SCIPtreeCreateProbingNode(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(SCIPtreeProbing(tree));

   SCIPsetDebugMsg(set, "new probing child in depth %d (probing depth: %d)\n", tree->pathlen, tree->pathlen-1 - SCIPnodeGetDepth(tree->probingroot));

   /* create temporary probing root node */
   SCIP_CALL( treeCreateProbingNode(tree, blkmem, set, lp) );

   return SCIP_OKAY;
}

/** sets the LP state for the current probing node
 *
 *  @note state and norms are stored at the node and later released by SCIP; therefore, the pointers are set
 *        to NULL by the method
 *
 *  @note the pointers to state and norms must not be NULL; however, they may point to a NULL pointer if the
 *        respective information should not be set
 */
SCIP_RETCODE SCIPtreeSetProbingLPState(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPISTATE**       lpistate,           /**< pointer to LP state information (like basis information) */
   SCIP_LPINORMS**       lpinorms,           /**< pointer to LP pricing norms information */
   SCIP_Bool             primalfeas,         /**< primal feasibility when LP state information was stored */
   SCIP_Bool             dualfeas            /**< dual feasibility when LP state information was stored */
   )
{
   SCIP_NODE* node;

   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));
   assert(lpistate != NULL);
   assert(lpinorms != NULL);

   /* get the current probing node */
   node = SCIPtreeGetCurrentNode(tree);

   /* this check is necessary to avoid cppcheck warnings */
   if( node == NULL )
      return SCIP_INVALIDDATA;

   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
   assert(node->data.probingnode != NULL);

   /* free already present LP state */
   if( node->data.probingnode->lpistate != NULL )
   {
      SCIP_CALL( SCIPlpFreeState(lp, blkmem, &(node->data.probingnode->lpistate)) );
   }

   /* free already present LP pricing norms */
   if( node->data.probingnode->lpinorms != NULL )
   {
      SCIP_CALL( SCIPlpFreeNorms(lp, blkmem, &(node->data.probingnode->lpinorms)) );
   }

   node->data.probingnode->lpistate = *lpistate;
   node->data.probingnode->lpinorms = *lpinorms;
   node->data.probingnode->lpwasprimfeas = primalfeas;
   node->data.probingnode->lpwasdualfeas = dualfeas;

   /* set the pointers to NULL to avoid that they are still used and modified by the caller */
   *lpistate = NULL;
   *lpinorms = NULL;

   tree->probingloadlpistate = TRUE;

   return SCIP_OKAY;
}

/** loads the LP state for the current probing node */
SCIP_RETCODE SCIPtreeLoadProbingLPState(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));

   /* loading the LP state is only necessary if we backtracked */
   if( tree->probingloadlpistate )
   {
      SCIP_NODE* node;
      SCIP_LPISTATE* lpistate;
      SCIP_LPINORMS* lpinorms;
      SCIP_Bool lpwasprimfeas = FALSE;
      SCIP_Bool lpwasprimchecked = FALSE;
      SCIP_Bool lpwasdualfeas = FALSE;
      SCIP_Bool lpwasdualchecked = FALSE;

      /* get the current probing node */
      node = SCIPtreeGetCurrentNode(tree);
      assert(node != NULL);
      assert(SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);

      /* search the last node where an LP state information was attached */
      lpistate = NULL;
      lpinorms = NULL;
      do
      {
         assert(SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
         assert(node->data.probingnode != NULL);
         if( node->data.probingnode->lpistate != NULL )
         {
            lpistate = node->data.probingnode->lpistate;
            lpinorms = node->data.probingnode->lpinorms;
            lpwasprimfeas = node->data.probingnode->lpwasprimfeas;
            lpwasprimchecked = node->data.probingnode->lpwasprimchecked;
            lpwasdualfeas = node->data.probingnode->lpwasdualfeas;
            lpwasdualchecked = node->data.probingnode->lpwasdualchecked;
            break;
         }
         node = node->parent;
         assert(node != NULL); /* the root node cannot be a probing node! */
      }
      while( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE );

      /* if there was no LP information stored in the probing nodes, use the one stored before probing started */
      if( lpistate == NULL )
      {
         lpistate = tree->probinglpistate;
         lpinorms = tree->probinglpinorms;
         lpwasprimfeas = tree->probinglpwasprimfeas;
         lpwasprimchecked = tree->probinglpwasprimchecked;
         lpwasdualfeas = tree->probinglpwasdualfeas;
         lpwasdualchecked = tree->probinglpwasdualchecked;
      }

      /* set the LP state */
      if( lpistate != NULL )
      {
         SCIP_CALL( SCIPlpSetState(lp, blkmem, set, eventqueue, lpistate,
               lpwasprimfeas, lpwasprimchecked, lpwasdualfeas, lpwasdualchecked) );
      }

      /* set the LP pricing norms */
      if( lpinorms != NULL )
      {
         SCIP_CALL( SCIPlpSetNorms(lp, blkmem, lpinorms) );
      }

      /* now we don't need to load the LP state again until the next backtracking */
      tree->probingloadlpistate = FALSE;
   }

   return SCIP_OKAY;
}

/** marks the probing node to have a solved LP relaxation */
SCIP_RETCODE SCIPtreeMarkProbingNodeHasLP(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_NODE* node;

   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));

   /* mark the probing node to have an LP */
   tree->probingnodehaslp = TRUE;

   /* get current probing node */
   node = SCIPtreeGetCurrentNode(tree);
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
   assert(node != NULL && node->data.probingnode != NULL);

   /* update LP information in probingnode data */
   /* cppcheck-suppress nullPointer */
   SCIP_CALL( probingnodeUpdate(node->data.probingnode, blkmem, tree, lp) );

   return SCIP_OKAY;
}

/** undoes all changes to the problem applied in probing up to the given probing depth */
static
SCIP_RETCODE treeBacktrackProbing(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_PRIMAL*          primal,             /**< primal data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   int                   probingdepth        /**< probing depth of the node in the probing path that should be reactivated,
                                              *   -1 to even deactivate the probing root, thus exiting probing mode */
   )
{
   int newpathlen;
   int i;

   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));
   assert(tree->probingroot != NULL);
   assert(tree->focusnode != NULL);
   assert(SCIPnodeGetType(tree->probingroot) == SCIP_NODETYPE_PROBINGNODE);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE
      || SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_REFOCUSNODE);
   assert(tree->probingroot->parent == tree->focusnode);
   assert(SCIPnodeGetDepth(tree->probingroot) == SCIPnodeGetDepth(tree->focusnode)+1);
   assert(tree->pathlen >= 2);
   assert(SCIPnodeGetType(tree->path[tree->pathlen-1]) == SCIP_NODETYPE_PROBINGNODE);
   assert(-1 <= probingdepth && probingdepth <= SCIPtreeGetProbingDepth(tree));

   treeCheckPath(tree);

   newpathlen = SCIPnodeGetDepth(tree->probingroot) + probingdepth + 1;
   assert(newpathlen >= 1); /* at least root node of the tree remains active */

   /* check if we have to do any backtracking */
   if( newpathlen < tree->pathlen )
   {
      int ncols;
      int nrows;

      /* the correct LP size of the node to which we backtracked is stored as initial LP size for its child */
      assert(SCIPnodeGetType(tree->path[newpathlen]) == SCIP_NODETYPE_PROBINGNODE);
      ncols = tree->path[newpathlen]->data.probingnode->ninitialcols;
      nrows = tree->path[newpathlen]->data.probingnode->ninitialrows;
      assert(ncols >= tree->pathnlpcols[newpathlen-1] || !tree->focuslpconstructed);
      assert(nrows >= tree->pathnlprows[newpathlen-1] || !tree->focuslpconstructed);

      while( tree->pathlen > newpathlen )
      {
         SCIP_NODE* node;

         node = tree->path[tree->pathlen-1];

         assert(SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
         assert(tree->pathlen-1 == SCIPnodeGetDepth(node));
         assert(tree->pathlen-1 >= SCIPnodeGetDepth(tree->probingroot));

         if( node->data.probingnode->nchgdobjs > 0 )
         {
            /* @todo only do this if we don't backtrack to the root node - in that case, we can just restore the unchanged
             *       objective values
             */
            for( i = node->data.probingnode->nchgdobjs - 1; i >= 0; --i )
            {
               assert(tree->probingobjchanged);

               SCIP_CALL( SCIPvarChgObj(node->data.probingnode->origobjvars[i], blkmem, set, transprob, primal, lp,
                     eventqueue, node->data.probingnode->origobjvals[i]) );
            }
            tree->probingsumchgdobjs -= node->data.probingnode->nchgdobjs;
            assert(tree->probingsumchgdobjs >= 0);

            /* reset probingobjchanged flag and cutoff bound */
            if( tree->probingsumchgdobjs == 0 )
            {
               SCIPlpUnmarkDivingObjChanged(lp);
               tree->probingobjchanged = FALSE;

               SCIP_CALL( SCIPlpSetCutoffbound(lp, set, transprob, primal->cutoffbound) );
            }

            /* recompute global and local pseudo objective values */
            SCIPlpRecomputeLocalAndGlobalPseudoObjval(lp, set, transprob);
         }

         /* undo bound changes by deactivating the probing node */
         SCIP_CALL( nodeDeactivate(node, blkmem, set, stat, tree, lp, branchcand, eventqueue) );

         /* free the probing node */
         SCIP_CALL( SCIPnodeFree(&tree->path[tree->pathlen-1], blkmem, set, stat, eventqueue, tree, lp) );
         tree->pathlen--;
      }
      assert(tree->pathlen == newpathlen);

      /* reset the path LP size to the initial size of the probing node */
      if( SCIPnodeGetType(tree->path[tree->pathlen-1]) == SCIP_NODETYPE_PROBINGNODE )
      {
         tree->pathnlpcols[tree->pathlen-1] = tree->path[tree->pathlen-1]->data.probingnode->ninitialcols;
         tree->pathnlprows[tree->pathlen-1] = tree->path[tree->pathlen-1]->data.probingnode->ninitialrows;
      }
      else
         assert(SCIPnodeGetType(tree->path[tree->pathlen-1]) == SCIP_NODETYPE_FOCUSNODE);
      treeCheckPath(tree);

      /* undo LP extensions */
      SCIP_CALL( SCIPlpShrinkCols(lp, set, ncols) );
      SCIP_CALL( SCIPlpShrinkRows(lp, blkmem, set, eventqueue, eventfilter, nrows) );
      tree->probingloadlpistate = TRUE; /* LP state must be reloaded if the next LP is solved */

      /* reset the LP's marked size to the initial size of the LP at the node stored in the path */
      assert(lp->nrows >= tree->pathnlprows[tree->pathlen-1] || !tree->focuslpconstructed);
      assert(lp->ncols >= tree->pathnlpcols[tree->pathlen-1] || !tree->focuslpconstructed);
      SCIPlpSetSizeMark(lp, tree->pathnlprows[tree->pathlen-1], tree->pathnlpcols[tree->pathlen-1]);

      /* if the highest cutoff or repropagation depth is inside the deleted part of the probing path,
       * reset them to infinity
       */
      if( tree->cutoffdepth >= tree->pathlen )
      {
         /* apply the pending bound changes */
         SCIP_CALL( treeApplyPendingBdchgs(tree, reopt, blkmem, set, stat, transprob, origprob, lp, branchcand, eventqueue, cliquetable) );

         /* applying the pending bound changes might have changed the cutoff depth; so the highest cutoff depth might
          * be outside of the deleted part of the probing path now
          */
         if( tree->cutoffdepth >= tree->pathlen )
            tree->cutoffdepth = INT_MAX;
      }
      if( tree->repropdepth >= tree->pathlen )
         tree->repropdepth = INT_MAX;
   }

   SCIPsetDebugMsg(set, "probing backtracked to depth %d (%d cols, %d rows)\n", tree->pathlen-1, SCIPlpGetNCols(lp), SCIPlpGetNRows(lp));

   return SCIP_OKAY;
}

/** undoes all changes to the problem applied in probing up to the given probing depth;
 *  the changes of the probing node of the given probing depth are the last ones that remain active;
 *  changes that were applied before calling SCIPtreeCreateProbingNode() cannot be undone
 */
SCIP_RETCODE SCIPtreeBacktrackProbing(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_PRIMAL*          primal,             /**< primal data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   int                   probingdepth        /**< probing depth of the node in the probing path that should be reactivated */
   )
{
   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));
   assert(0 <= probingdepth && probingdepth <= SCIPtreeGetProbingDepth(tree));

   /* undo the domain and constraint set changes and free the temporary probing nodes below the given probing depth */
   SCIP_CALL( treeBacktrackProbing(tree, reopt, blkmem, set, stat, transprob, origprob, lp, primal, branchcand,
         eventqueue, eventfilter, cliquetable, probingdepth) );

   assert(SCIPtreeProbing(tree));
   assert(SCIPnodeGetType(SCIPtreeGetCurrentNode(tree)) == SCIP_NODETYPE_PROBINGNODE);

   return SCIP_OKAY;
}

/** switches back from probing to normal operation mode, frees all nodes on the probing path, restores bounds of all
 *  variables and restores active constraints arrays of focus node
 */
SCIP_RETCODE SCIPtreeEndProbing(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_PRIMAL*          primal,             /**< Primal LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));
   assert(tree->probingroot != NULL);
   assert(tree->focusnode != NULL);
   assert(SCIPnodeGetType(tree->probingroot) == SCIP_NODETYPE_PROBINGNODE);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE
      || SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_REFOCUSNODE);
   assert(tree->probingroot->parent == tree->focusnode);
   assert(SCIPnodeGetDepth(tree->probingroot) == SCIPnodeGetDepth(tree->focusnode)+1);
   assert(tree->pathlen >= 2);
   assert(SCIPnodeGetType(tree->path[tree->pathlen-1]) == SCIP_NODETYPE_PROBINGNODE);
   assert(set != NULL);

   /* undo the domain and constraint set changes of the temporary probing nodes and free the probing nodes */
   SCIP_CALL( treeBacktrackProbing(tree, reopt, blkmem, set, stat, transprob, origprob, lp, primal, branchcand,
         eventqueue, eventfilter, cliquetable, -1) );
   assert(tree->probingsumchgdobjs == 0);
   assert(!tree->probingobjchanged);
   assert(!lp->divingobjchg);
   assert(lp->cutoffbound == primal->cutoffbound); /*lint !e777*/
   assert(SCIPtreeGetCurrentNode(tree) == tree->focusnode);
   assert(!SCIPtreeProbing(tree));

   /* if the LP was flushed before probing starts, flush it again */
   if( tree->probinglpwasflushed )
   {
      SCIP_CALL( SCIPlpFlush(lp, blkmem, set, eventqueue) );

      /* if the LP was solved before probing starts, solve it again to restore the LP solution */
      if( tree->probinglpwassolved )
      {
         SCIP_Bool lperror;

         /* reset the LP state before probing started */
         if( tree->probinglpistate == NULL )
         {
            assert(tree->probinglpinorms == NULL);
            SCIP_CALL( SCIPlpiClearState(lp->lpi) );
            lp->primalfeasible = (lp->nlpicols == 0 && lp->nlpirows == 0);
            lp->primalchecked = (lp->nlpicols == 0 && lp->nlpirows == 0);
            lp->dualfeasible = (lp->nlpicols == 0 && lp->nlpirows == 0);
            lp->dualchecked = (lp->nlpicols == 0 && lp->nlpirows == 0);
            lp->solisbasic = FALSE;
         }
         else
         {
            SCIP_CALL( SCIPlpSetState(lp, blkmem, set, eventqueue, tree->probinglpistate,
                  tree->probinglpwasprimfeas, tree->probinglpwasprimchecked, tree->probinglpwasdualfeas,
                  tree->probinglpwasdualchecked) );
            SCIP_CALL( SCIPlpFreeState(lp, blkmem, &tree->probinglpistate) );

            if( tree->probinglpinorms != NULL )
            {
               SCIP_CALL( SCIPlpSetNorms(lp, blkmem, tree->probinglpinorms) );
               SCIP_CALL( SCIPlpFreeNorms(lp, blkmem, &tree->probinglpinorms) );
               tree->probinglpinorms = NULL;
            }
         }
         SCIPlpSetIsRelax(lp, tree->probinglpwasrelax);

         /* resolve LP to reset solution */
         SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, transprob, -1LL, FALSE, FALSE, FALSE, &lperror) );
         if( lperror )
         {
            SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "(node %" SCIP_LONGINT_FORMAT ") unresolved numerical troubles while resolving LP %" SCIP_LONGINT_FORMAT " after probing\n",
               stat->nnodes, stat->nlps);
            lp->resolvelperror = TRUE;
            tree->focusnodehaslp = FALSE;
         }
         else if( SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OPTIMAL 
            && SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_INFEASIBLE
            && SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY
            && SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OBJLIMIT )
         {
            SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "LP was not resolved to a sufficient status after probing\n");
            lp->resolvelperror = TRUE;
            tree->focusnodehaslp = FALSE;
         }
         else if( tree->focuslpconstructed && SCIPlpIsRelax(lp) && SCIPprobAllColsInLP(transprob, set, lp))
         {
            SCIP_CALL( SCIPnodeUpdateLowerboundLP(tree->focusnode, set, stat, tree, transprob, origprob, lp) );
         }
      }
   }
   else
      lp->flushed = FALSE;

   assert(tree->probinglpistate == NULL);

   /* if no LP was solved during probing and the LP before probing was not solved, then it should not be solved now */
   assert(tree->probingsolvedlp || tree->probinglpwassolved || !lp->solved);

   /* if the LP was solved (and hence flushed) before probing, then lp->solved should be TRUE unless we occured an error
    * during resolving right above
    */
   assert(!tree->probinglpwassolved || lp->solved || lp->resolvelperror);

   /* if the LP was not solved before probing it should be marked unsolved now; this can occur if a probing LP was
    * solved in between
    */
   if( !tree->probinglpwassolved )
   {
      lp->solved = FALSE;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   /* if the LP was solved during probing, but had been unsolved before probing started, we discard the LP state */
   if( set->lp_clearinitialprobinglp && tree->probingsolvedlp && !tree->probinglpwassolved )
   {
      SCIPsetDebugMsg(set, "clearing lp state at end of probing mode because LP was initially unsolved\n");
      SCIP_CALL( SCIPlpiClearState(lp->lpi) );
   }

   /* if a relaxation was stored before probing, restore it now */
   if( tree->probdiverelaxstored )
   {
      SCIP_CALL( SCIPtreeRestoreRelaxSol(tree, set, relaxation, transprob) );
   }

   assert(tree->probingobjchanged == SCIPlpDivingObjChanged(lp));

   /* reset flags */
   tree->probinglpwasflushed = FALSE;
   tree->probinglpwassolved = FALSE;
   tree->probingloadlpistate = FALSE;
   tree->probinglpwasrelax = FALSE;
   tree->probingsolvedlp = FALSE;
   tree->sbprobing = FALSE;

   /* inform LP about end of probing mode */
   SCIP_CALL( SCIPlpEndProbing(lp) );

   /* reset all marked constraints for propagation */
   SCIP_CALL( SCIPconshdlrsResetPropagationStatus(set, blkmem, set->conshdlrs, set->nconshdlrs) );

   SCIPsetDebugMsg(set, "probing ended in depth %d (LP flushed: %u, solstat: %d)\n", tree->pathlen-1, lp->flushed, SCIPlpGetSolstat(lp));

   return SCIP_OKAY;
}

/** stores relaxation solution before diving or probing */
SCIP_RETCODE SCIPtreeStoreRelaxSol(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_PROB*            transprob           /**< transformed problem after presolve */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(tree != NULL);
   assert(set != NULL);
   assert(relaxation != NULL);
   assert(transprob != NULL);
   assert(SCIPrelaxationIsSolValid(relaxation));

   nvars = transprob->nvars;
   vars = transprob->vars;

   /* check if memory still needs to be allocated */
   if( tree->probdiverelaxsol == NULL )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&(tree->probdiverelaxsol), nvars) );
   }

   /* iterate over all variables to save the relaxation solution */
   for( v = 0; v < nvars; ++v )
      tree->probdiverelaxsol[v] = SCIPvarGetRelaxSol(vars[v], set);

   tree->probdiverelaxstored = TRUE;
   tree->probdiverelaxincludeslp = SCIPrelaxationIsLpIncludedForSol(relaxation);

   return SCIP_OKAY;
}

/** restores relaxation solution after diving or probing */
SCIP_RETCODE SCIPtreeRestoreRelaxSol(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_PROB*            transprob           /**< transformed problem after presolve */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(tree != NULL);
   assert(set != NULL);
   assert(tree->probdiverelaxstored);
   assert(tree->probdiverelaxsol != NULL);

   nvars = transprob->nvars;
   vars = transprob->vars;

   /* iterate over all variables to restore the relaxation solution */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPvarSetRelaxSol(vars[v], set, relaxation, tree->probdiverelaxsol[v], TRUE) );
   }

   tree->probdiverelaxstored = FALSE;
   SCIPrelaxationSetSolValid(relaxation, TRUE, tree->probdiverelaxincludeslp);

   return SCIP_OKAY;
}

/** gets the best child of the focus node w.r.t. the node selection priority assigned by the branching rule */
SCIP_NODE* SCIPtreeGetPrioChild(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   SCIP_NODE* bestnode;
   SCIP_Real bestprio;
   int i;

   assert(tree != NULL);

   bestnode = NULL;
   bestprio = SCIP_REAL_MIN;
   for( i = 0; i < tree->nchildren; ++i )
   {
      if( tree->childrenprio[i] > bestprio )
      {
         bestnode = tree->children[i];
         bestprio = tree->childrenprio[i];
      }
   }
   assert((tree->nchildren == 0) == (bestnode == NULL));

   return bestnode;
}

/** gets the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule */
SCIP_NODE* SCIPtreeGetPrioSibling(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   SCIP_NODE* bestnode;
   SCIP_Real bestprio;
   int i;

   assert(tree != NULL);

   bestnode = NULL;
   bestprio = SCIP_REAL_MIN;
   for( i = 0; i < tree->nsiblings; ++i )
   {
      if( tree->siblingsprio[i] > bestprio )
      {
         bestnode = tree->siblings[i];
         bestprio = tree->siblingsprio[i];
      }
   }
   assert((tree->nsiblings == 0) == (bestnode == NULL));

   return bestnode;
}

/** gets the best child of the focus node w.r.t. the node selection strategy */
SCIP_NODE* SCIPtreeGetBestChild(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODE* bestnode;
   int i;

   assert(tree != NULL);

   nodesel = SCIPnodepqGetNodesel(tree->leaves);
   assert(nodesel != NULL);

   bestnode = NULL;
   for( i = 0; i < tree->nchildren; ++i )
   {
      if( bestnode == NULL || SCIPnodeselCompare(nodesel, set, tree->children[i], bestnode) < 0 )
      {
         bestnode = tree->children[i];
      }
   }

   return bestnode;
}

/** gets the best sibling of the focus node w.r.t. the node selection strategy */
SCIP_NODE* SCIPtreeGetBestSibling(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODE* bestnode;
   int i;

   assert(tree != NULL);

   nodesel = SCIPnodepqGetNodesel(tree->leaves);
   assert(nodesel != NULL);

   bestnode = NULL;
   for( i = 0; i < tree->nsiblings; ++i )
   {
      if( bestnode == NULL || SCIPnodeselCompare(nodesel, set, tree->siblings[i], bestnode) < 0 )
      {
         bestnode = tree->siblings[i];
      }
   }

   return bestnode;
}

/** gets the best leaf from the node queue w.r.t. the node selection strategy */
SCIP_NODE* SCIPtreeGetBestLeaf(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqFirst(tree->leaves);
}

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy */
SCIP_NODE* SCIPtreeGetBestNode(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODE* bestchild;
   SCIP_NODE* bestsibling;
   SCIP_NODE* bestleaf;
   SCIP_NODE* bestnode;

   assert(tree != NULL);

   nodesel = SCIPnodepqGetNodesel(tree->leaves);
   assert(nodesel != NULL);

   /* get the best child, sibling, and leaf */
   bestchild = SCIPtreeGetBestChild(tree, set);
   bestsibling = SCIPtreeGetBestSibling(tree, set);
   bestleaf = SCIPtreeGetBestLeaf(tree);

   /* return the best of the three */
   bestnode = bestchild;
   if( bestsibling != NULL && (bestnode == NULL || SCIPnodeselCompare(nodesel, set, bestsibling, bestnode) < 0) )
      bestnode = bestsibling;
   if( bestleaf != NULL && (bestnode == NULL || SCIPnodeselCompare(nodesel, set, bestleaf, bestnode) < 0) )
      bestnode = bestleaf;

   assert(SCIPtreeGetNLeaves(tree) == 0 || bestnode != NULL);

   return bestnode;
}

/** gets the minimal lower bound of all nodes in the tree */
SCIP_Real SCIPtreeGetLowerbound(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_Real lowerbound;
   int i;

   assert(tree != NULL);
   assert(set != NULL);

   /* get the lower bound from the queue */
   lowerbound = SCIPnodepqGetLowerbound(tree->leaves, set);

   /* compare lower bound with children */
   for( i = 0; i < tree->nchildren; ++i )
   {
      assert(tree->children[i] != NULL);
      lowerbound = MIN(lowerbound, tree->children[i]->lowerbound); 
   }

   /* compare lower bound with siblings */
   for( i = 0; i < tree->nsiblings; ++i )
   {
      assert(tree->siblings[i] != NULL);
      lowerbound = MIN(lowerbound, tree->siblings[i]->lowerbound); 
   }

   /* compare lower bound with focus node */
   if( tree->focusnode != NULL )
   {
      lowerbound = MIN(lowerbound, tree->focusnode->lowerbound);
   }

   return lowerbound;
}

/** gets the node with minimal lower bound of all nodes in the tree (child, sibling, or leaf) */
SCIP_NODE* SCIPtreeGetLowerboundNode(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_NODE* lowerboundnode;
   SCIP_Real lowerbound;
   SCIP_Real bestprio;
   int i;

   assert(tree != NULL);
   assert(set != NULL);

   /* get the lower bound from the queue */
   lowerboundnode = SCIPnodepqGetLowerboundNode(tree->leaves, set);
   lowerbound = lowerboundnode != NULL ? lowerboundnode->lowerbound : SCIPsetInfinity(set);
   bestprio = -SCIPsetInfinity(set);

   /* compare lower bound with children */
   for( i = 0; i < tree->nchildren; ++i )
   {
      assert(tree->children[i] != NULL);
      if( SCIPsetIsLE(set, tree->children[i]->lowerbound, lowerbound) )
      {
         if( SCIPsetIsLT(set, tree->children[i]->lowerbound, lowerbound) || tree->childrenprio[i] > bestprio )
         {
            lowerboundnode = tree->children[i]; 
            lowerbound = lowerboundnode->lowerbound; 
            bestprio = tree->childrenprio[i];
         }
      }
   }

   /* compare lower bound with siblings */
   for( i = 0; i < tree->nsiblings; ++i )
   {
      assert(tree->siblings[i] != NULL);
      if( SCIPsetIsLE(set, tree->siblings[i]->lowerbound, lowerbound) )
      {
         if( SCIPsetIsLT(set, tree->siblings[i]->lowerbound, lowerbound) || tree->siblingsprio[i] > bestprio )
         {
            lowerboundnode = tree->siblings[i]; 
            lowerbound = lowerboundnode->lowerbound; 
            bestprio = tree->siblingsprio[i];
         }
      }
   }

   return lowerboundnode;
}

/** gets the average lower bound of all nodes in the tree */
SCIP_Real SCIPtreeGetAvgLowerbound(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Real             cutoffbound         /**< global cutoff bound */
   )
{
   SCIP_Real lowerboundsum;
   int nnodes;
   int i;

   assert(tree != NULL);

   /* get sum of lower bounds from nodes in the queue */
   lowerboundsum = SCIPnodepqGetLowerboundSum(tree->leaves);
   nnodes = SCIPtreeGetNLeaves(tree);

   /* add lower bound of focus node */
   if( tree->focusnode != NULL && tree->focusnode->lowerbound < cutoffbound )
   {
      lowerboundsum += tree->focusnode->lowerbound;
      nnodes++;
   }

   /* add lower bounds of siblings */
   for( i = 0; i < tree->nsiblings; ++i )
   {
      assert(tree->siblings[i] != NULL);
      lowerboundsum += tree->siblings[i]->lowerbound;
   }
   nnodes += tree->nsiblings;

   /* add lower bounds of children */
   for( i = 0; i < tree->nchildren; ++i )
   {
      assert(tree->children[i] != NULL);
      lowerboundsum += tree->children[i]->lowerbound;
   }
   nnodes += tree->nchildren;

   return nnodes == 0 ? 0.0 : lowerboundsum/nnodes;
}




/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPnodeGetType
#undef SCIPnodeGetNumber
#undef SCIPnodeGetDepth
#undef SCIPnodeGetLowerbound
#undef SCIPnodeGetEstimate
#undef SCIPnodeGetDomchg
#undef SCIPnodeGetParent
#undef SCIPnodeGetConssetchg
#undef SCIPnodeIsActive
#undef SCIPnodeIsPropagatedAgain
#undef SCIPtreeGetNLeaves
#undef SCIPtreeGetNChildren
#undef SCIPtreeGetNSiblings
#undef SCIPtreeGetNNodes
#undef SCIPtreeIsPathComplete
#undef SCIPtreeProbing
#undef SCIPtreeGetProbingRoot
#undef SCIPtreeGetProbingDepth
#undef SCIPtreeGetFocusNode
#undef SCIPtreeGetFocusDepth
#undef SCIPtreeHasFocusNodeLP
#undef SCIPtreeSetFocusNodeLP 
#undef SCIPtreeIsFocusNodeLPConstructed
#undef SCIPtreeInRepropagation
#undef SCIPtreeGetCurrentNode
#undef SCIPtreeGetCurrentDepth
#undef SCIPtreeHasCurrentNodeLP
#undef SCIPtreeGetEffectiveRootDepth
#undef SCIPtreeGetRootNode
#undef SCIPtreeProbingObjChanged
#undef SCIPtreeMarkProbingObjChanged

/** gets the type of the node */
SCIP_NODETYPE SCIPnodeGetType(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return (SCIP_NODETYPE)(node->nodetype);
}

/** gets successively assigned number of the node */
SCIP_Longint SCIPnodeGetNumber(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->number;
}

/** gets the depth of the node */
int SCIPnodeGetDepth(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return (int) node->depth;
}

/** gets the lower bound of the node */
SCIP_Real SCIPnodeGetLowerbound(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->lowerbound;
}

/** gets the estimated value of the best feasible solution in subtree of the node */
SCIP_Real SCIPnodeGetEstimate(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->estimate;
}

/** gets the reoptimization type of this node */
SCIP_REOPTTYPE SCIPnodeGetReopttype(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return (SCIP_REOPTTYPE)node->reopttype;
}

/** sets the reoptimization type of this node */
void SCIPnodeSetReopttype(
   SCIP_NODE*            node,               /**< node */
   SCIP_REOPTTYPE        reopttype           /**< reoptimization type */
   )
{
   assert(node != NULL);
   assert(reopttype == SCIP_REOPTTYPE_NONE
       || reopttype == SCIP_REOPTTYPE_TRANSIT
       || reopttype == SCIP_REOPTTYPE_INFSUBTREE
       || reopttype == SCIP_REOPTTYPE_STRBRANCHED
       || reopttype == SCIP_REOPTTYPE_LOGICORNODE
       || reopttype == SCIP_REOPTTYPE_LEAF
       || reopttype == SCIP_REOPTTYPE_PRUNED
       || reopttype == SCIP_REOPTTYPE_FEASIBLE);

   node->reopttype = (unsigned int) reopttype;
}

/** gets the unique id to identify the node during reoptimization; the id is 0 if the node is the root or not part of
 * the reoptimization tree
 */
unsigned int SCIPnodeGetReoptID(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->reoptid; /*lint !e732*/
}

/** set a unique id to identify the node during reoptimization */
void SCIPnodeSetReoptID(
   SCIP_NODE*            node,               /**< node */
   unsigned int          id                  /**< unique id */
   )
{
   assert(node != NULL);
   assert(id <= 536870911); /* id has only 29 bits and needs to be smaller than 2^29 */

   node->reoptid = id;
}

/** gets the domain change information of the node, i.e., the information about the differences in the
 *  variables domains to the parent node
 */
SCIP_DOMCHG* SCIPnodeGetDomchg(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->domchg;
}

/** counts the number of bound changes due to branching, constraint propagation, and propagation */
void SCIPnodeGetNDomchg(
   SCIP_NODE*            node,               /**< node */
   int*                  nbranchings,        /**< pointer to store number of branchings (or NULL if not needed) */
   int*                  nconsprop,          /**< pointer to store number of constraint propagations (or NULL if not needed) */
   int*                  nprop               /**< pointer to store number of propagations (or NULL if not needed) */
   )
{  /*lint --e{641}*/
   SCIP_Bool count_branchings;
   SCIP_Bool count_consprop;
   SCIP_Bool count_prop;
   int i;

   assert(node != NULL);

   count_branchings = (nbranchings != NULL);
   count_consprop = (nconsprop != NULL);
   count_prop = (nprop != NULL);

   /* set counter to zero */
   if( count_branchings )
      *nbranchings = 0;
   if( count_consprop )
      *nconsprop = 0;
   if( count_prop )
      *nprop = 0;

   if( node->domchg != NULL )
   {
      for( i = 0; i < (int) node->domchg->domchgbound.nboundchgs; i++ )
      {
         if( count_branchings && node->domchg->domchgbound.boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
            (*nbranchings)++;
         else if( count_consprop && node->domchg->domchgbound.boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER )
            (*nconsprop)++;
         else if( count_prop && node->domchg->domchgbound.boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER )
            (*nprop)++;
      }
   }
}

/* return the number of bound changes based on dual information.
 *
 * currently, this methods works only for bound changes made by strong branching on binary variables. we need this
 * method to ensure optimality within reoptimization.
 *
 * since the bound changes made by strong branching are stored as SCIP_BOUNDCHGTYPE_CONSINFER or SCIP_BOUNDCHGTYPE_PROPINFER
 * with no constraint or propagator, resp., we are are interested in bound changes with these attributes.
 *
 * all bound changes of type SCIP_BOUNDCHGTYPE_BRANCHING are stored in the beginning of the bound change array, afterwards,
 * we can find the other two types. thus, we start the search at the end of the list and stop when reaching the first
 * bound change of type SCIP_BOUNDCHGTYPE_BRANCHING.
 */
int SCIPnodeGetNDualBndchgs(
   SCIP_NODE*            node                /**< node */
   )
{  /*lint --e{641}*/
   SCIP_BOUNDCHG* boundchgs;
   int i;
   int nboundchgs;
   int npseudobranchvars;

   assert(node != NULL);

   if( node->domchg == NULL )
      return 0;

   nboundchgs = (int)node->domchg->domchgbound.nboundchgs;
   boundchgs = node->domchg->domchgbound.boundchgs;

   npseudobranchvars = 0;

   assert(boundchgs != NULL);
   assert(nboundchgs >= 0);

   /* count the number of pseudo-branching decisions; pseudo-branching decisions have to be in the ending of the bound change
    * array
    */
   for( i = nboundchgs-1; i >= 0; i--)
   {
      SCIP_Bool isint;

      isint = boundchgs[i].var->vartype == SCIP_VARTYPE_BINARY || boundchgs[i].var->vartype == SCIP_VARTYPE_INTEGER
           || boundchgs[i].var->vartype == SCIP_VARTYPE_IMPLINT;

      if( isint && ((boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER
            && boundchgs[i].data.inferencedata.reason.cons == NULL)
        || (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER
            && boundchgs[i].data.inferencedata.reason.prop == NULL)) )
         npseudobranchvars++;
      else if( boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
         break;
   }

   return npseudobranchvars;
}

/** returns the set of variable branchings that were performed in the parent node to create this node */
void SCIPnodeGetDualBoundchgs(
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR**            vars,               /**< array of variables on which the bound change is based on dual information */
   SCIP_Real*            bounds,             /**< array of bounds which are based on dual information */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of boundtypes which are based on dual information */
   int*                  nvars,              /**< number of variables on which the bound change is based on dual information
                                              *   if this is larger than the array size, arrays should be reallocated and method
                                              *   should be called again */
   int                   varssize            /**< available slots in arrays */
   )
{  /*lint --e{641}*/
   SCIP_BOUNDCHG* boundchgs;
   int nboundchgs;
   int i;

   assert(node != NULL);
   assert(vars != NULL);
   assert(bounds != NULL);
   assert(boundtypes != NULL);
   assert(nvars != NULL);
   assert(varssize >= 0);

   (*nvars) = 0;

   if( SCIPnodeGetDepth(node) == 0 || node->domchg == NULL )
      return;

   nboundchgs = (int)node->domchg->domchgbound.nboundchgs;
   boundchgs = node->domchg->domchgbound.boundchgs;

   assert(boundchgs != NULL);
   assert(nboundchgs >= 0);

   /* count the number of pseudo-branching decisions; pseudo-branching decisions have to be in the ending of the bound change
    * array
    */
   for( i = nboundchgs-1; i >= 0; i--)
   {
      if( boundchgs[i].var->vartype == SCIP_VARTYPE_BINARY || boundchgs[i].var->vartype == SCIP_VARTYPE_INTEGER
       || boundchgs[i].var->vartype == SCIP_VARTYPE_IMPLINT )
      {
         if( (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER
               && boundchgs[i].data.inferencedata.reason.cons == NULL)
          || (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER
                && boundchgs[i].data.inferencedata.reason.prop == NULL) )
            (*nvars)++;
         else if( boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
            break;
      }
   }

   /* if the arrays have enough space store the branching decisions */
   if( varssize >= *nvars )
   {
      int j;
      j = 0;
      for( i = i+1; i < nboundchgs; i++)
      {
         assert( boundchgs[i].boundchgtype != SCIP_BOUNDCHGTYPE_BRANCHING );
         if( boundchgs[i].var->vartype == SCIP_VARTYPE_BINARY || boundchgs[i].var->vartype == SCIP_VARTYPE_INTEGER
          || boundchgs[i].var->vartype == SCIP_VARTYPE_IMPLINT )
         {
            if( (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER
                  && boundchgs[i].data.inferencedata.reason.cons == NULL)
             || (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER
                  && boundchgs[i].data.inferencedata.reason.prop == NULL) )
            {
               vars[j] = boundchgs[i].var;
               bounds[j] = boundchgs[i].newbound;
               boundtypes[j] = (SCIP_BOUNDTYPE) boundchgs[i].boundtype;
               j++;
            }
         }
      }
   }
}

/** gets the parent node of a node in the branch-and-bound tree, if any */
SCIP_NODE* SCIPnodeGetParent(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->parent;
}

/** returns the set of variable branchings that were performed in the parent node to create this node */
void SCIPnodeGetParentBranchings(
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR**            branchvars,         /**< array of variables on which the branching has been performed in the parent node */
   SCIP_Real*            branchbounds,       /**< array of bounds which the branching in the parent node set */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of boundtypes which the branching in the parent node set */
   int*                  nbranchvars,        /**< number of variables on which branching has been performed in the parent node
                                              *   if this is larger than the array size, arrays should be reallocated and method
                                              *   should be called again */
   int                   branchvarssize      /**< available slots in arrays */
   )
{
   SCIP_BOUNDCHG* boundchgs;
   int nboundchgs;
   int i;

   assert(node != NULL);
   assert(branchvars != NULL);
   assert(branchbounds != NULL);
   assert(boundtypes != NULL);
   assert(nbranchvars != NULL);
   assert(branchvarssize >= 0);

   (*nbranchvars) = 0;

   if( SCIPnodeGetDepth(node) == 0 || node->domchg == NULL )
      return;

   nboundchgs = (int)node->domchg->domchgbound.nboundchgs;
   boundchgs = node->domchg->domchgbound.boundchgs;

   assert(boundchgs != NULL);
   assert(nboundchgs >= 0);

   /* count the number of branching decisions; branching decisions have to be in the beginning of the bound change
    * array
    */
   for( i = 0; i < nboundchgs; i++)
   {
      if( boundchgs[i].boundchgtype != SCIP_BOUNDCHGTYPE_BRANCHING ) /*lint !e641*/
         break;

      (*nbranchvars)++;
   }

#ifndef NDEBUG
   /* check that the remaining bound change are no branching decisions */
   for( ; i < nboundchgs; i++)
      assert(boundchgs[i].boundchgtype != SCIP_BOUNDCHGTYPE_BRANCHING); /*lint !e641*/
#endif

   /* if the arrays have enough space store the branching decisions */
   if( branchvarssize >= *nbranchvars )
   {
      for( i = 0; i < *nbranchvars; i++)
      {
         assert( boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING ); /*lint !e641*/
         branchvars[i] = boundchgs[i].var;
         boundtypes[i] = (SCIP_BOUNDTYPE) boundchgs[i].boundtype;
         branchbounds[i] = boundchgs[i].newbound;
      }
   }
}

/** returns the set of variable branchings that were performed in all ancestor nodes (nodes on the path to the root) to create this node */
void SCIPnodeGetAncestorBranchings(
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR**            branchvars,         /**< array of variables on which the branchings has been performed in all ancestors */
   SCIP_Real*            branchbounds,       /**< array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of boundtypes which the branchings in all ancestors set */
   int*                  nbranchvars,        /**< number of variables on which branchings have been performed in all ancestors
                                              *   if this is larger than the array size, arrays should be reallocated and method
                                              *   should be called again */
   int                   branchvarssize      /**< available slots in arrays */
   )
{
   assert(node != NULL);
   assert(branchvars != NULL);
   assert(branchbounds != NULL);
   assert(boundtypes != NULL);
   assert(nbranchvars != NULL);
   assert(branchvarssize >= 0);

   (*nbranchvars) = 0;

   while( SCIPnodeGetDepth(node) != 0 )
   {
      int nodenbranchvars;
      int start;
      int size;

      start = *nbranchvars < branchvarssize - 1 ? *nbranchvars : branchvarssize - 1;
      size = *nbranchvars > branchvarssize ? 0 : branchvarssize-(*nbranchvars);

      SCIPnodeGetParentBranchings(node, &branchvars[start], &branchbounds[start], &boundtypes[start], &nodenbranchvars, size);
      *nbranchvars += nodenbranchvars;

      node = node->parent;
   }
}

/** returns the set of variable branchings that were performed between the given @p node and the given @p parent node. */
void SCIPnodeGetAncestorBranchingsPart(
   SCIP_NODE*            node,               /**< node data */
   SCIP_NODE*            parent,             /**< node data of the last ancestor node */
   SCIP_VAR**            branchvars,         /**< array of variables on which the branchings has been performed in all ancestors */
   SCIP_Real*            branchbounds,       /**< array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of boundtypes which the branchings in all ancestors set */
   int*                  nbranchvars,        /**< number of variables on which branchings have been performed in all ancestors
                                              *   if this is larger than the array size, arrays should be reallocated and method
                                              *   should be called again */
   int                   branchvarssize      /**< available slots in arrays */
   )
{
   assert(node != NULL);
   assert(parent != NULL);
   assert(branchvars != NULL);
   assert(branchbounds != NULL);
   assert(boundtypes != NULL);
   assert(nbranchvars != NULL);
   assert(branchvarssize >= 0);

   (*nbranchvars) = 0;

   while( node != parent )
   {
      int nodenbranchvars;
      int start;
      int size;

      start = *nbranchvars < branchvarssize - 1 ? *nbranchvars : branchvarssize - 1;
      size = *nbranchvars > branchvarssize ? 0 : branchvarssize-(*nbranchvars);

      SCIPnodeGetParentBranchings(node, &branchvars[start], &branchbounds[start], &boundtypes[start], &nodenbranchvars, size);
      *nbranchvars += nodenbranchvars;

      node = node->parent;
   }
}

/** return all bound changes based on constraint propagation; stop saving the bound changes if we reach a branching
 *  decision based on a dual information
 */
void SCIPnodeGetConsProps(
   SCIP_NODE*            node,               /**< node */
   SCIP_VAR**            vars,               /**< array of variables on which constraint propagation triggers a bound change */
   SCIP_Real*            varbounds,          /**< array of bounds set by constraint propagation */
   SCIP_BOUNDTYPE*       varboundtypes,      /**< array of boundtypes set by constraint propagation */
   int*                  nconspropvars,      /**< number of variables on which constraint propagation triggers a bound change
                                              *   if this is larger than the array size, arrays should be reallocated and method
                                              *   should be called again */
   int                   conspropvarssize    /**< available slots in arrays */
   )
{  /*lint --e{641}*/
   SCIP_BOUNDCHG* boundchgs;
   int nboundchgs;
   int first_dual;
   int nskip;
   int i;

   assert(node != NULL);
   assert(vars != NULL);
   assert(varbounds != NULL);
   assert(varboundtypes != NULL);
   assert(nconspropvars != NULL);
   assert(conspropvarssize >= 0);

   (*nconspropvars) = 0;

   if( SCIPnodeGetDepth(node) == 0 || node->domchg == NULL )
      return;

   nboundchgs = (int)node->domchg->domchgbound.nboundchgs;
   boundchgs = node->domchg->domchgbound.boundchgs;

   assert(boundchgs != NULL);
   assert(nboundchgs >= 0);

   SCIPnodeGetNDomchg(node, &nskip, NULL, NULL);
   i = nskip;

   while( i < nboundchgs
       && !(boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER && boundchgs[i].data.inferencedata.reason.cons == NULL)
       && !(boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER && boundchgs[i].data.inferencedata.reason.prop == NULL) )
      i++;

   first_dual = i;

   /* count the number of bound changes because of constraint propagation and propagation */
   for(i = nskip; i < first_dual; i++)
   {
      assert(boundchgs[i].boundchgtype != SCIP_BOUNDCHGTYPE_BRANCHING);

      if( (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER && boundchgs[i].data.inferencedata.reason.cons != NULL)
       || (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER && boundchgs[i].data.inferencedata.reason.prop != NULL) )
      {
         if( boundchgs[i].var->vartype != SCIP_VARTYPE_CONTINUOUS )
            (*nconspropvars)++;
      }
      else if( (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER && boundchgs[i].data.inferencedata.reason.cons == NULL)
            || (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER && boundchgs[i].data.inferencedata.reason.prop == NULL))
         break;
   }

   /* if the arrays have enough space store the branching decisions */
   if( conspropvarssize >= *nconspropvars )
   {
      int pos;

      for(i = nskip, pos = 0; i < first_dual; i++)
      {
         if( boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER && boundchgs[i].data.inferencedata.reason.cons != NULL )
         {
            if( boundchgs[i].var->vartype != SCIP_VARTYPE_CONTINUOUS )
            {
               vars[pos] = boundchgs[i].var;
               varboundtypes[pos] = (SCIP_BOUNDTYPE) boundchgs[i].boundtype;
               varbounds[pos] = boundchgs[i].newbound;
               pos++;
            }
         }
      }
   }

   return;
}

/** gets all bound changes applied after the first bound change based on dual information.
 *
 *  @note: currently, we can only detect bound changes based in dual information if they arise from strong branching.
 */
void SCIPnodeGetBdChgsAfterDual(
   SCIP_NODE*            node,               /**< node */
   SCIP_VAR**            vars,               /**< array of variables on which the branching has been performed in the parent node */
   SCIP_Real*            varbounds,          /**< array of bounds which the branching in the parent node set */
   SCIP_BOUNDTYPE*       varboundtypes,      /**< array of boundtypes which the branching in the parent node set */
   int                   start,              /**< first free slot in the arrays */
   int*                  nbranchvars,        /**< number of variables on which branching has been performed in the parent node
                                              *   if this is larger than the array size, arrays should be reallocated and method
                                              *   should be called again */
   int                   branchvarssize      /**< available slots in arrays */
   )
{  /*lint --e{641}*/
   SCIP_BOUNDCHG* boundchgs;
   int nboundchgs;
   int first_dual;
   int i;

   assert(node != NULL);
   assert(vars != NULL);
   assert(varbounds != NULL);
   assert(varboundtypes != NULL);
   assert(nbranchvars != NULL);
   assert(branchvarssize >= 0);

   (*nbranchvars) = 0;

   if( SCIPnodeGetDepth(node) == 0 || node->domchg == NULL )
      return;

   nboundchgs = (int)node->domchg->domchgbound.nboundchgs;
   boundchgs = node->domchg->domchgbound.boundchgs;

   assert(boundchgs != NULL);
   assert(nboundchgs >= 0);

   /* find the first based on dual information */
   i = 0;
   while( i < nboundchgs
       && !(boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER && boundchgs[i].data.inferencedata.reason.cons == NULL)
       && !(boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER && boundchgs[i].data.inferencedata.reason.prop == NULL) )
      i++;

   first_dual = i;

   /* count the number of branching decisions; branching decisions have to be in the beginning of the bound change array */
   for( ; i < nboundchgs; i++)
   {
      assert(boundchgs[i].boundchgtype != SCIP_BOUNDCHGTYPE_BRANCHING);

      if( (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER && boundchgs[i].data.inferencedata.reason.cons != NULL)
       || (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER && boundchgs[i].data.inferencedata.reason.prop != NULL) )
      {
         if( boundchgs[i].var->vartype != SCIP_VARTYPE_CONTINUOUS )
            (*nbranchvars)++;
      }
   }

   /* if the arrays have enough space store the branching decisions */
   if( branchvarssize >= *nbranchvars )
   {
      int p;
      for(i = first_dual, p = start; i < nboundchgs; i++)
      {
         if( (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER && boundchgs[i].data.inferencedata.reason.cons != NULL)
          || (boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER && boundchgs[i].data.inferencedata.reason.prop != NULL)  )
         {
            if( boundchgs[i].var->vartype != SCIP_VARTYPE_CONTINUOUS )
            {
               vars[p] = boundchgs[i].var;
               varboundtypes[p] = (SCIP_BOUNDTYPE) boundchgs[i].boundtype;
               varbounds[p] = boundchgs[i].newbound;
               p++;
            }
         }
      }
   }
}

/** outputs the path into given file stream in GML format */
SCIP_RETCODE SCIPnodePrintAncestorBranchings(
   SCIP_NODE*            node,               /**< node data */
   FILE*                 file                /**< file to output the path */
   )
{
   int nbranchings;

   nbranchings = 0;

   /* print opening in GML format */
   SCIPgmlWriteOpening(file, TRUE);

   while( SCIPnodeGetDepth(node) != 0 )
   {
      SCIP_BOUNDCHG* boundchgs;
      char label[SCIP_MAXSTRLEN];
      int nboundchgs;
      int i;

      nboundchgs = (int)node->domchg->domchgbound.nboundchgs;
      boundchgs = node->domchg->domchgbound.boundchgs;

      for( i = 0; i < nboundchgs; i++)
      {
         if( boundchgs[i].boundchgtype != SCIP_BOUNDCHGTYPE_BRANCHING ) /*lint !e641*/
            break;

         (void) SCIPsnprintf(label, SCIP_MAXSTRLEN, "%s %s %g", SCIPvarGetName(boundchgs[i].var),
            (SCIP_BOUNDTYPE) boundchgs[i].boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", boundchgs[i].newbound);

         SCIPgmlWriteNode(file, (unsigned int)nbranchings, label, "circle", NULL, NULL);

         if( nbranchings > 0 )
         {
            SCIPgmlWriteArc(file, (unsigned int)nbranchings, (unsigned int)(nbranchings-1), NULL, NULL);
         }

         nbranchings++;
      }

      node = node->parent;
   }

   /* print closing in GML format */
   SCIPgmlWriteClosing(file);

   return SCIP_OKAY;
}

/**  returns the set of variable branchings that were performed in all ancestor nodes (nodes on the path to the root) to create this node
 *  sorted by the nodes, starting from the current node going up to the root
 */
void SCIPnodeGetAncestorBranchingPath(
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR**            branchvars,         /**< array of variables on which the branchings has been performed in all ancestors */
   SCIP_Real*            branchbounds,       /**< array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of boundtypes which the branchings in all ancestors set */
   int*                  nbranchvars,        /**< number of variables on which branchings have been performed in all ancestors
                                              *   if this is larger than the array size, arrays should be reallocated and method
                                              *   should be called again */
   int                   branchvarssize,     /**< available slots in arrays */
   int*                  nodeswitches,       /**< marks, where in the arrays the branching decisions of the next node on the path
                                              *   start branchings performed at the parent of node always start at position 0.
                                              *   For single variable branching, nodeswitches[i] = i holds */
   int*                  nnodes,             /**< number of nodes in the nodeswitch array */
   int                   nodeswitchsize      /**< available slots in node switch array */
   )
{
   assert(node != NULL);
   assert(branchvars != NULL);
   assert(branchbounds != NULL);
   assert(boundtypes != NULL);
   assert(nbranchvars != NULL);
   assert(branchvarssize >= 0);

   (*nbranchvars) = 0;
   (*nnodes) = 0;

   /* go up to the root, in the root no domains were changed due to branching */
   while( SCIPnodeGetDepth(node) != 0 )
   {
      int nodenbranchvars;
      int start;
      int size;

      /* calculate the start position for the current node and the maximum remaining slots in the arrays */
      start = *nbranchvars < branchvarssize - 1 ? *nbranchvars : branchvarssize - 1;
      size = *nbranchvars > branchvarssize ? 0 : branchvarssize-(*nbranchvars);
      if( *nnodes < nodeswitchsize )         
         nodeswitches[*nnodes] = start;

      /* get branchings for a single node */
      SCIPnodeGetParentBranchings(node, &branchvars[start], &branchbounds[start], &boundtypes[start], &nodenbranchvars, size);
      *nbranchvars += nodenbranchvars;
      (*nnodes)++;

      node = node->parent;
   }
}

/** checks for two nodes whether they share the same root path, i.e., whether one is an ancestor of the other */
SCIP_Bool SCIPnodesSharePath(
   SCIP_NODE*            node1,              /**< node data */
   SCIP_NODE*            node2               /**< node data */
   )
{
   assert(node1 != NULL);
   assert(node2 != NULL);
   assert(SCIPnodeGetDepth(node1) >= 0);
   assert(SCIPnodeGetDepth(node2) >= 0);

   /* if node2 is deeper than node1, follow the path until the level of node2 */
   while( SCIPnodeGetDepth(node1) < SCIPnodeGetDepth(node2) )
      node2 = node2->parent;

   /* if node1 is deeper than node2, follow the path until the level of node1 */
   while( SCIPnodeGetDepth(node2) < SCIPnodeGetDepth(node1) )
      node1 = node1->parent;

   assert(SCIPnodeGetDepth(node2) == SCIPnodeGetDepth(node1));

   return (node1 == node2);
}

/** finds the common ancestor node of two given nodes */
SCIP_NODE* SCIPnodesGetCommonAncestor(
   SCIP_NODE*            node1,              /**< node data */
   SCIP_NODE*            node2               /**< node data */
   )
{
   assert(node1 != NULL);
   assert(node2 != NULL);
   assert(SCIPnodeGetDepth(node1) >= 0);
   assert(SCIPnodeGetDepth(node2) >= 0);

   /* if node2 is deeper than node1, follow the path until the level of node2 */
   while( SCIPnodeGetDepth(node1) < SCIPnodeGetDepth(node2) )
      node2 = node2->parent;

   /* if node1 is deeper than node2, follow the path until the level of node1 */
   while( SCIPnodeGetDepth(node2) < SCIPnodeGetDepth(node1) )
      node1 = node1->parent;

   /* move up level by level until you found a common ancestor */
   while( node1 != node2 )
   {
      node1 = node1->parent;
      node2 = node2->parent;
      assert(SCIPnodeGetDepth(node1) == SCIPnodeGetDepth(node2));
   }
   assert(SCIPnodeGetDepth(node1) >= 0);

   return node1;
}

/** returns whether node is in the path to the current node */
SCIP_Bool SCIPnodeIsActive(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->active;
}

/** returns whether the node is marked to be propagated again */
SCIP_Bool SCIPnodeIsPropagatedAgain(
   SCIP_NODE*            node                /**< node data */
   )
{
   assert(node != NULL);

   return node->reprop;
}

/* returns the set of changed constraints for a particular node */
SCIP_CONSSETCHG* SCIPnodeGetConssetchg(
   SCIP_NODE*            node                /**< node data */
   )
{
   assert(node != NULL);

   return node->conssetchg;
}

/** gets number of children of the focus node */
int SCIPtreeGetNChildren(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->nchildren;
}

/** gets number of siblings of the focus node  */
int SCIPtreeGetNSiblings(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->nsiblings;
}

/** gets number of leaves in the tree (excluding children and siblings of focus nodes) */
int SCIPtreeGetNLeaves(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqLen(tree->leaves);
}

/** gets number of open nodes in the tree (children + siblings + leaves) */
int SCIPtreeGetNNodes(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->nchildren + tree->nsiblings + SCIPtreeGetNLeaves(tree);
}

/** returns whether the active path goes completely down to the focus node */
SCIP_Bool SCIPtreeIsPathComplete(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->focusnode != NULL || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->focusnode != NULL);
   assert(tree->pathlen >= 2 || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] != NULL);
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1]->depth == tree->pathlen-1);
   assert(tree->focusnode == NULL || (int)tree->focusnode->depth >= tree->pathlen
      || tree->path[tree->focusnode->depth] == tree->focusnode);

   return (tree->focusnode == NULL || (int)tree->focusnode->depth < tree->pathlen);
}

/** returns whether the current node is a temporary probing node */
SCIP_Bool SCIPtreeProbing(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->probingroot == NULL || (SCIP_NODETYPE)tree->probingroot->nodetype == SCIP_NODETYPE_PROBINGNODE);
   assert(tree->probingroot == NULL || tree->pathlen > SCIPnodeGetDepth(tree->probingroot));
   assert(tree->probingroot == NULL || tree->path[SCIPnodeGetDepth(tree->probingroot)] == tree->probingroot);

   return (tree->probingroot != NULL);
}

/** returns the temporary probing root node, or NULL if the we are not in probing mode */
SCIP_NODE* SCIPtreeGetProbingRoot(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->probingroot == NULL || (SCIP_NODETYPE)tree->probingroot->nodetype == SCIP_NODETYPE_PROBINGNODE);
   assert(tree->probingroot == NULL || tree->pathlen > SCIPnodeGetDepth(tree->probingroot));
   assert(tree->probingroot == NULL || tree->path[SCIPnodeGetDepth(tree->probingroot)] == tree->probingroot);

   return tree->probingroot;
}

/** gets focus node of the tree */
SCIP_NODE* SCIPtreeGetFocusNode(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->focusnode != NULL || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->focusnode != NULL);
   assert(tree->pathlen >= 2 || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] != NULL);
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1]->depth == tree->pathlen-1);
   assert(tree->focusnode == NULL || (int)tree->focusnode->depth >= tree->pathlen
      || tree->path[tree->focusnode->depth] == tree->focusnode);

   return tree->focusnode;
}

/** gets depth of focus node in the tree */
int SCIPtreeGetFocusDepth(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->focusnode != NULL || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->focusnode != NULL);
   assert(tree->pathlen >= 2 || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] != NULL);
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1]->depth == tree->pathlen-1);
   assert(tree->focusnode == NULL || (int)tree->focusnode->depth >= tree->pathlen
      || tree->path[tree->focusnode->depth] == tree->focusnode);

   return tree->focusnode != NULL ? (int)tree->focusnode->depth : -1;
}

/** returns, whether the LP was or is to be solved in the focus node */
SCIP_Bool SCIPtreeHasFocusNodeLP(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->focusnodehaslp;
}

/** sets mark to solve or to ignore the LP while processing the focus node */
void SCIPtreeSetFocusNodeLP(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             solvelp             /**< should the LP be solved in focus node? */
   )
{
   assert(tree != NULL);

   tree->focusnodehaslp = solvelp;
}

/** returns whether the LP of the focus node is already constructed */
SCIP_Bool SCIPtreeIsFocusNodeLPConstructed(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->focuslpconstructed;
}

/** returns whether the focus node is already solved and only propagated again */
SCIP_Bool SCIPtreeInRepropagation(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return (tree->focusnode != NULL && SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_REFOCUSNODE);
}

/** gets current node of the tree, i.e. the last node in the active path, or NULL if no current node exists */
SCIP_NODE* SCIPtreeGetCurrentNode(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->focusnode != NULL || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->focusnode != NULL);
   assert(tree->pathlen >= 2 || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] != NULL);
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1]->depth == tree->pathlen-1);
   assert(tree->focusnode == NULL || (int)tree->focusnode->depth >= tree->pathlen
      || tree->path[tree->focusnode->depth] == tree->focusnode);

   return (tree->pathlen > 0 ? tree->path[tree->pathlen-1] : NULL);
}

/** gets depth of current node in the tree, i.e. the length of the active path minus 1, or -1 if no current node exists */
int SCIPtreeGetCurrentDepth(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->focusnode != NULL || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->focusnode != NULL);
   assert(tree->pathlen >= 2 || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] != NULL);
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1]->depth == tree->pathlen-1);
   assert(tree->focusnode == NULL || (int)tree->focusnode->depth >= tree->pathlen
      || tree->path[tree->focusnode->depth] == tree->focusnode);

   return tree->pathlen-1;
}

/** returns, whether the LP was or is to be solved in the current node */
SCIP_Bool SCIPtreeHasCurrentNodeLP(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(SCIPtreeIsPathComplete(tree));

   return SCIPtreeProbing(tree) ? tree->probingnodehaslp : SCIPtreeHasFocusNodeLP(tree);
}

/** returns the current probing depth, i.e. the number of probing sub nodes existing in the probing path */
int SCIPtreeGetProbingDepth(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));

   return SCIPtreeGetCurrentDepth(tree) - SCIPnodeGetDepth(tree->probingroot);
}

/** returns the depth of the effective root node (i.e. the first depth level of a node with at least two children) */
int SCIPtreeGetEffectiveRootDepth(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->effectiverootdepth >= 0);

   return tree->effectiverootdepth;
}

/** gets the root node of the tree */
SCIP_NODE* SCIPtreeGetRootNode(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->root;
}

/** returns whether we are in probing and the objective value of at least one column was changed */

SCIP_Bool SCIPtreeProbingObjChanged(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(SCIPtreeProbing(tree) || !tree->probingobjchanged);

   return tree->probingobjchanged;
}

/** marks the current probing node to have a changed objective function */
void SCIPtreeMarkProbingObjChanged(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));

   tree->probingobjchanged = TRUE;
}
