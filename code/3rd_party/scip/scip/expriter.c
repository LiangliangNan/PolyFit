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

/**@file   expriter.c
 * @ingroup OTHER_CFILES
 * @brief  functions for iterating over algebraic expressions
 * @author Benjamin Mueller
 * @author Stefan Vigerske
 */

/* enable this to record where active iterators were initialized
 * (not thread-safe; problematic when using several SCIP instances concurrently)
 */
/* #define SCIP_DEBUG_EXPRITER */

#include <assert.h>

#include "scip/expr.h"
#include "scip/pub_misc.h"
#include "scip/struct_expr.h"
#include "scip/struct_stat.h"

#define MINDFSSIZE  16 /**< minimum stack size for DFS*/
#define MINBFSSIZE  16 /**< minimum queue size for BFS */

#ifdef SCIP_DEBUG_EXPRITER
#include <execinfo.h>
#include <string.h>
#include <stdlib.h>

#define MAXSUBSCIPDEPTH   10   /**< minimal subscip-depth to no longer store backtrace */
#define MAXBACKTRACE      20   /**< maximal length of backtrace to store */

/** backtrace when iterator was initialized
 * - store per subscip-depth (easier than storing per SCIP instance)
 * - store per iterator position that can be active concurrently
 * - one string for each entry in backtrace
 * - each entry up to 200 characters
 */
char iterinitbacktrace[MAXSUBSCIPDEPTH][SCIP_EXPRITER_MAXNACTIVE][MAXBACKTRACE][200];
#endif

/*
 * local functions
 */

#ifdef SCIP_DEBUG_EXPRITER
/** obtain current backtrace and store it in iterinitbacktrace */
static
void storeBacktrace(
   int                   subscipdepth,       /**< current subscip depth */
   int                   iterpos             /**< iterator position where to store backtrace */
   )
{
   void* array[MAXBACKTRACE];
   char** strings;
   int size;
   int i;

   assert(subscipdepth >= 0);
   assert(iterpos >= 0);
   assert(iterpos < SCIP_EXPRITER_MAXNACTIVE);

   if( subscipdepth > MAXSUBSCIPDEPTH )
      return;

   size = backtrace(array, MAXBACKTRACE);
   strings = backtrace_symbols(array, size);
   if( strings == NULL )
      size = 0;

   for( i = 0; i < size; i++ )
      strncpy(iterinitbacktrace[subscipdepth][iterpos][i], strings[i], sizeof(iterinitbacktrace[0][0][0]));

   /* '\0' for remining backtrace entries */
   while( size < MAXBACKTRACE )
      iterinitbacktrace[subscipdepth][iterpos][size++][0] = '\0';

   free(strings);
}

static
void printBacktraces(
   int                   subscipdepth        /**< current subscip depth */
   )
{
   int i, j;

   assert(subscipdepth >= 0);
   if( subscipdepth >= MAXSUBSCIPDEPTH )
   {
      SCIPerrorMessage("subscip depth %d too high to report active iterators", subscipdepth);
      return;
   }

   for( i = 0; i < SCIP_EXPRITER_MAXNACTIVE-1; ++i )
   {
      SCIPerrorMessage("Active iterator %d created at:\n", i);
      for( j = 0; j < MAXBACKTRACE; ++j )
      {
         if( iterinitbacktrace[subscipdepth][i][j][0] == '\0' )
            break;
         SCIPerrorMessage("  %s\n", iterinitbacktrace[subscipdepth][i][j]);
      }
   }
}
#else
#define storeBacktrace(subscipdepth, iterpos)

/*lint -e{715}*/
static
void printBacktraces(
   int                   subscipdepth        /**< current subscip depth */
   )
{  /*lint --e{715}*/
   SCIPerrorMessage("Rebuild with SCIP_DEBUG_EXPRITER defined in src/scip/expriter.c to see where currently "
         "active iterators were initialized.\n");
}
#endif

static
void deinit(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL );

   if( !iterator->initialized )
      return;

   if( iterator->iterindex >= 0 )
   {
      /* the iterindex must be the one of the last initialized iterator */
      assert(iterator->iterindex == iterator->stat->nactiveexpriter-1);

      /* tell core that this iterator is no longer active */
      --iterator->stat->nactiveexpriter;

      iterator->iterindex = -1;
   }

   switch( iterator->itertype )
   {
      case SCIP_EXPRITER_BFS :
      {
         assert(iterator->queue != NULL);

         SCIPqueueFree(&iterator->queue);

         break;
      }

      case SCIP_EXPRITER_RTOPOLOGIC :
      {
         assert(iterator->dfsnvisited != NULL);
         assert(iterator->dfsexprs != NULL);

         /* free dfs arrays */
         BMSfreeBlockMemoryArray(iterator->blkmem, &iterator->dfsnvisited, iterator->dfssize);
         BMSfreeBlockMemoryArray(iterator->blkmem, &iterator->dfsexprs, iterator->dfssize);
         iterator->dfssize = 0;

         break;
      }

      case SCIP_EXPRITER_DFS :
      default: break;
   }
}

/** ensures minimum stack size of iterator's data */
static
SCIP_RETCODE ensureStackSize(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   int                   size                /**< minimum requires size */
   )
{
   assert(iterator != NULL);
   assert(iterator->blkmem != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_RTOPOLOGIC);
   assert(size >= 0);

   if( size > iterator->dfssize )
   {
      int newsize = size * 2;

      SCIP_ALLOC( BMSreallocBlockMemoryArray(iterator->blkmem, &iterator->dfsexprs, iterator->dfssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(iterator->blkmem, &iterator->dfsnvisited, iterator->dfssize, newsize) );
      iterator->dfssize = newsize;
   }

   return SCIP_OKAY;
}

/** adds an expression to the DFS stack */
static
void reverseTopologicalInsert(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(iterator != NULL);
   assert(expr != NULL);

   SCIP_CALL_ABORT( ensureStackSize(iterator, iterator->dfsnexprs + 1) );
   iterator->dfsexprs[iterator->dfsnexprs] = expr;
   iterator->dfsnvisited[iterator->dfsnexprs] = 0;
   ++(iterator->dfsnexprs);
}

/** moves to the next expression according to a reverse topological order */
static
SCIP_EXPR* doReverseTopologicalNext(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   SCIP_EXPR* expr;
   int childidx;

   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_RTOPOLOGIC);

   /* no expression left */
   if( iterator->dfsnexprs == 0 )
      return NULL;

   /* get expression on the top of the stack */
   expr = iterator->dfsexprs[iterator->dfsnexprs - 1];
   childidx = iterator->dfsnvisited[iterator->dfsnexprs - 1];

   /* remove the expression if all children have been visited */
   if( childidx >= SCIPexprGetNChildren(expr) )
   {
      --(iterator->dfsnexprs);
      return expr;
   }
   /* go to the next child */
   else
   {
      SCIP_EXPR* child = SCIPexprGetChildren(expr)[childidx];
      assert(child != NULL);

      /* mark that the child has been visited */
      ++(iterator->dfsnvisited[iterator->dfsnexprs-1]);

      /* do left-most step */
      while( SCIPexprGetNChildren(child) > 0 )
      {
         /* add child to the DFS stack */
         reverseTopologicalInsert(iterator, child);

         /* mark that the child has been visited; note that child is on top of the DFS stack */
         ++(iterator->dfsnvisited[iterator->dfsnexprs-1]);

         child = SCIPexprGetChildren(child)[0];
      }

      /* return last child; NOTE this child is not been added to the stack */
      return child;
   }
}

/** moves to the next expression according to the BFS rule */
static
SCIP_EXPR* doBfsNext(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   SCIP_EXPR* expr;
   int i;

   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_BFS);
   assert(iterator->queue != NULL);

   /* no expression left */
   if( SCIPqueueIsEmpty(iterator->queue) )
      return NULL;

   expr = (SCIP_EXPR*) SCIPqueueRemove(iterator->queue);
   assert(expr != NULL);

   assert(iterator->visitedtag == 0 || iterator->iterindex >= 0);
   assert(iterator->visitedtag == 0 || iterator->iterindex < SCIP_EXPRITER_MAXNACTIVE);
   /* we should have set the visitedtag when adding the expression to the queue */
   assert(iterator->visitedtag == 0 || expr->iterdata[iterator->iterindex].visitedtag == iterator->visitedtag);

   /* add all (possibly non-visited) children to the queue */
   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      SCIP_EXPR* child = SCIPexprGetChildren(expr)[i];
      assert(child != NULL);

      if( iterator->visitedtag != 0 )
      {
         assert(iterator->iterindex >= 0);
         assert(iterator->iterindex < SCIP_EXPRITER_MAXNACTIVE);

         /* skip children that have already been visited or have already been added to the queue */
         if( child->iterdata[iterator->iterindex].visitedtag == iterator->visitedtag )
            continue;

         /* mark child as being in the queue (will be inserted next) */
         child->iterdata[iterator->iterindex].visitedtag = iterator->visitedtag;
      }

      /* add child to the queue */
      SCIP_CALL_ABORT( SCIPqueueInsert(iterator->queue, child) );
   }

   return expr;
}

/** moves to the next expression according to the DFS rule */
static
SCIP_EXPR* doDfsNext(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   SCIP_EXPRITERDATA* iterdata;

   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->iterindex >= 0);

   if( iterator->curr == NULL )
      return NULL;

   iterdata = &iterator->curr->iterdata[iterator->iterindex];

   switch( iterator->dfsstage )
   {
      case SCIP_EXPRITER_VISITEDCHILD:
         /* consider next child */
         ++iterdata->currentchild;
         /* fall through */ /* no break */ /*lint -fallthrough*/

      case SCIP_EXPRITER_ENTEREXPR:
      {
         /* if there is an unvisited child (left), then go into visitingchild stage, otherwise go to leave stage */
         iterator->dfsstage = SCIP_EXPRITER_LEAVEEXPR;  /* expect that we will leave expr, and change mind to visitingchild below */
         while( iterdata->currentchild < iterator->curr->nchildren )
         {
            if( iterator->visitedtag == 0 || iterator->visitedtag != iterator->curr->children[iterdata->currentchild]->iterdata[iterator->iterindex].visitedtag )
            {
               /* if visitedtag is not used or child "currentchild" has not been visited yet, then go into visitingchild stage for this child */
               iterator->dfsstage = SCIP_EXPRITER_VISITINGCHILD;
               break;
            }
            ++iterdata->currentchild;
         }
         /* if leaving expr, then currentchild should be at nchildren */
         assert(iterator->dfsstage == SCIP_EXPRITER_VISITINGCHILD || iterdata->currentchild == iterator->curr->nchildren);
         /* if visiting child, then currentchild should be a valid index */
         assert(iterator->dfsstage == SCIP_EXPRITER_LEAVEEXPR || iterdata->currentchild < iterator->curr->nchildren);
         /* if visiting child, then either we don't care whether we visited it already or it has not been visited yet */
         assert(iterator->dfsstage == SCIP_EXPRITER_LEAVEEXPR || iterator->visitedtag == 0
         || iterator->visitedtag != iterator->curr->children[iterdata->currentchild]->iterdata[iterator->iterindex].visitedtag);

         return iterator->curr;
      }

      case SCIP_EXPRITER_VISITINGCHILD:
      {
         SCIP_EXPR* child;

         assert(iterdata->currentchild < iterator->curr->nchildren);

         /* remember the parent and set the first child that should be visited of the new root */
         child = iterator->curr->children[iterdata->currentchild];
         child->iterdata[iterator->iterindex].parent = iterator->curr;
         child->iterdata[iterator->iterindex].currentchild = 0;

         /* visit child */
         iterator->dfsstage = SCIP_EXPRITER_ENTEREXPR;

         return child;
      }

      case SCIP_EXPRITER_LEAVEEXPR:
      {
         /* go back to parent expression */

         /* remember that this expression has been visited */
         iterdata->visitedtag = iterator->visitedtag;

         /* be in visitedchild stage for the parent */
         iterator->dfsstage = SCIP_EXPRITER_VISITEDCHILD;

         return iterdata->parent;
      }

      default:
         /* unknown stage */
         SCIPABORT();
         return NULL;
   }
}

/*
 * private functions (expr.h)
 */

/** creates an expression iterator */
SCIP_RETCODE SCIPexpriterCreate(
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRITER**       iterator            /**< buffer to store expression iterator */
   )
{
   assert(stat != NULL);
   assert(blkmem  != NULL);
   assert(iterator != NULL);

   SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, iterator) );

   (*iterator)->stat = stat;
   (*iterator)->blkmem = blkmem;

   return SCIP_OKAY;
}

/** frees an expression iterator */
void SCIPexpriterFree(
   SCIP_EXPRITER**       iterator            /**< pointer to the expression iterator */
   )
{
   assert(iterator != NULL);
   assert(*iterator != NULL);
   assert((*iterator)->blkmem != NULL);

   deinit(*iterator);

   assert((*iterator)->queue == NULL);
   assert((*iterator)->dfsnvisited == NULL);
   assert((*iterator)->dfsexprs == NULL);

   /* free iterator */
   BMSfreeBlockMemory((*iterator)->blkmem, iterator);
}

/*
 * public functions (pub_expr.h)
 */

#ifdef NDEBUG
#undef SCIPexpriterIsInit
#undef SCIPexpriterGetCurrent
#undef SCIPexpriterGetStageDFS
#undef SCIPexpriterGetChildIdxDFS
#undef SCIPexpriterGetChildExprDFS
#undef SCIPexpriterGetParentDFS
#undef SCIPexpriterGetCurrentUserData
#undef SCIPexpriterGetChildUserDataDFS
#undef SCIPexpriterGetExprUserData
#undef SCIPexpriterSetCurrentUserData
#undef SCIPexpriterSetExprUserData
#undef SCIPexpriterSetChildUserData
#undef SCIPexpriterIsEnd
#endif

/** returns whether expression iterator is currently initialized */
SCIP_Bool SCIPexpriterIsInit(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);

   return iterator->initialized;
}

/** initializes an expression iterator
 *
 * @note If `expr` is NULL, then iterator will be set into ended-state (SCIPexpriterIsEnd() is TRUE). Useful if following with SCIPexpriterRestartDFS().
 *
 * If type is DFS, then `stopstages` will be set to \ref SCIP_EXPRITER_ENTEREXPR.
 * Use `SCIPexpriterSetStagesDFS` to change this.
 */
SCIP_RETCODE SCIPexpriterInit(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr,               /**< expression of the iterator, can be NULL */
   SCIP_EXPRITER_TYPE    type,               /**< type of expression iterator */
   SCIP_Bool             allowrevisit        /**< whether expression are allowed to be visited more than once */
   )
{
   assert(iterator != NULL);

   deinit(iterator);

   /* store the new type of the iterator */
   iterator->itertype = type;

   /* get iterindex, if necessary */
   if( !allowrevisit || type == SCIP_EXPRITER_DFS )
   {
      if( iterator->stat->nactiveexpriter + 1 >= SCIP_EXPRITER_MAXNACTIVE )
      {
         SCIPerrorMessage("Maximal number of active expression iterators reached at subscip-depth %d.\n",
               iterator->stat->subscipdepth);
         printBacktraces(iterator->stat->subscipdepth);
         return SCIP_MAXDEPTHLEVEL;
      }

      iterator->iterindex = iterator->stat->nactiveexpriter++;

      storeBacktrace(iterator->stat->subscipdepth, iterator->iterindex);
   }
   else
   {
      iterator->iterindex = -1;
   }

   /* get new tag to recognize visited expressions */
   if( !allowrevisit )
   {
      iterator->visitedtag = ++iterator->stat->exprlastvisitedtag;
   }
   else
   {
      iterator->visitedtag = 0L;
   }

   switch( iterator->itertype )
   {
      case SCIP_EXPRITER_BFS:
      {
         SCIP_CALL( SCIPqueueCreate(&iterator->queue, MINBFSSIZE, 2.0) );

         assert(iterator->queue != NULL);
         SCIPqueueClear(iterator->queue);

         if( expr == NULL )
         {
            iterator->curr = NULL;
            break;
         }

         SCIP_CALL( SCIPqueueInsert(iterator->queue, expr) );

         if( iterator->visitedtag != 0 )
         {
            assert(iterator->iterindex >= 0);
            assert(iterator->iterindex < SCIP_EXPRITER_MAXNACTIVE);
            assert(expr->iterdata[iterator->iterindex].visitedtag != iterator->visitedtag);

            /* mark expression as being in the queue */
            expr->iterdata[iterator->iterindex].visitedtag = iterator->visitedtag;
         }

         iterator->curr = SCIPexpriterGetNext(iterator);
         break;
      }

      case SCIP_EXPRITER_RTOPOLOGIC :
      {
         SCIP_CALL( ensureStackSize(iterator, MINDFSSIZE) );

         if( expr != NULL )
         {
            reverseTopologicalInsert(iterator, expr);
            iterator->curr = SCIPexpriterGetNext(iterator);
         }
         else
         {
            iterator->curr = NULL;
         }

         break;
      }

      case SCIP_EXPRITER_DFS :
      {
         assert(iterator->iterindex >= 0);

         iterator->stopstages = SCIP_EXPRITER_ENTEREXPR;
         iterator->curr = expr;

         if( expr == NULL )
            break;

         expr->iterdata[iterator->iterindex].currentchild = 0;
         expr->iterdata[iterator->iterindex].parent = NULL;
         iterator->dfsstage = SCIP_EXPRITER_ENTEREXPR;

         break;
      }
   }

   iterator->initialized = TRUE;

   return SCIP_OKAY;
}

/** restarts an already initialized expression iterator in DFS mode
 *
 * The expression iterator will continue from the given expression, not revisiting expressions that
 * this iterator has already been visited (if initialized with `allowrevisit=FALSE`) and giving access
 * to the same iterator specified expression data that may have been set already.
 * Also the stop-stages are not reset.
 *
 * If revisiting is forbidden and given expr has already been visited, then the iterator will behave
 * as on the end of iteration (SCIPexpriterIsEnd() is TRUE).
 * If the enterexpr stage is not one of the stop stages, then the iterator will be moved forward
 * (SCIPexpriterGetNext() is called).
 *
 * @return The current expression.
 */
SCIP_EXPR* SCIPexpriterRestartDFS(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr                /**< expression of the iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->initialized);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);

   /* if we forbid revisiting and root expr has already been visited, then set curr to NULL, that is, be at end of iterator */
   if( iterator->visitedtag > 0 && expr->iterdata[iterator->iterindex].visitedtag == iterator->visitedtag )
   {
      iterator->curr = NULL;
      return NULL;
   }

   /* set current to given expr, make it the root, and set stage to enterexpr */
   iterator->curr = expr;
   expr->iterdata[iterator->iterindex].currentchild = 0;
   expr->iterdata[iterator->iterindex].parent = NULL;
   iterator->dfsstage = SCIP_EXPRITER_ENTEREXPR;

   if( (iterator->stopstages & SCIP_EXPRITER_ENTEREXPR) == 0 )
      return SCIPexpriterGetNext(iterator);

   return iterator->curr;
}

/** specifies in which stages to stop a DFS iterator
 *
 * Parameter `stopstages` should be a bitwise OR of different \ref SCIP_EXPRITER_STAGE values
 *
 * If the current stage is not one of the `stopstages`, then the iterator will be moved on.
 */
void SCIPexpriterSetStagesDFS(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPRITER_STAGE   stopstages          /**< the stages in which to stop when iterating via DFS */
   )
{
   assert(iterator != NULL);

   if( (iterator->dfsstage & stopstages) == 0 )
   {
      iterator->stopstages = stopstages;
      (void) SCIPexpriterGetNext(iterator);
   }
   else
   {
      iterator->stopstages = stopstages;
   }
}

/** gets the current expression that the expression iterator points to */
SCIP_EXPR* SCIPexpriterGetCurrent(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);

   return iterator->curr;
}

/** gets the current stage that the expression iterator is in when using DFS
 *
 * If the iterator has finished (SCIPexpriterIsEnd() is TRUE), then the stage is undefined.
 */
SCIP_EXPRITER_STAGE SCIPexpriterGetStageDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);

   return iterator->dfsstage;
}

/** gets the index of the child that the expression iterator considers when in DFS mode and stage \ref SCIP_EXPRITER_VISITINGCHILD or \ref SCIP_EXPRITER_VISITEDCHILD */
int SCIPexpriterGetChildIdxDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->dfsstage == SCIP_EXPRITER_VISITINGCHILD || iterator->dfsstage == SCIP_EXPRITER_VISITEDCHILD);

   return iterator->curr->iterdata[iterator->iterindex].currentchild;
}

/** gets the child expression that the expression iterator considers when in DFS mode and stage \ref SCIP_EXPRITER_VISITINGCHILD or \ref SCIP_EXPRITER_VISITEDCHILD */
SCIP_EXPR* SCIPexpriterGetChildExprDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->dfsstage == SCIP_EXPRITER_VISITINGCHILD || iterator->dfsstage == SCIP_EXPRITER_VISITEDCHILD);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild >= 0);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild < iterator->curr->nchildren);

   return iterator->curr->children[iterator->curr->iterdata[iterator->iterindex].currentchild];
}

/** gives the parent of the current expression of an expression iteration if in DFS mode
 *
 * @return the expression from which the current expression has been accessed
 */
SCIP_EXPR* SCIPexpriterGetParentDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);

   return iterator->curr->iterdata[iterator->iterindex].parent;
}

/** gives the iterator specific user data of the current expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPRITER_USERDATA SCIPexpriterGetCurrentUserData(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);

   return iterator->curr->iterdata[iterator->iterindex].userdata;
}

/** gives the iterator specific user data of the current expressions current child
 *
 * @note The expression iterator mode must be in DFS mode and stage \ref SCIP_EXPRITER_VISITINGCHILD or \ref SCIP_EXPRITER_VISITEDCHILD
 */
SCIP_EXPRITER_USERDATA SCIPexpriterGetChildUserDataDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->dfsstage == SCIP_EXPRITER_VISITINGCHILD || iterator->dfsstage == SCIP_EXPRITER_VISITEDCHILD);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild >= 0);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild < iterator->curr->nchildren);

   return iterator->curr->children[iterator->curr->iterdata[iterator->iterindex].currentchild]->iterdata[iterator->iterindex].userdata;
}

/** gives the iterator specific user data of a given expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPRITER_USERDATA SCIPexpriterGetExprUserData(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr                /**< expression for which to get the userdata of this iterator */
   )
{
   assert(iterator != NULL);
   assert(expr != NULL);
   assert(iterator->iterindex >= 0);

   return expr->iterdata[iterator->iterindex].userdata;
}

/** sets the iterator specific user data of the current expression for an expression iteration if in DFS mode
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
void SCIPexpriterSetCurrentUserData(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPRITER_USERDATA userdata           /**< data to be stored */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);

   iterator->curr->iterdata[iterator->iterindex].userdata = userdata;
}

/** sets the iterator specific user data of a given expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
void SCIPexpriterSetExprUserData(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr,               /**< expression where to set iterator data */
   SCIP_EXPRITER_USERDATA userdata           /**< data to be stored in current child */
   )
{
   assert(iterator != NULL);
   assert(iterator->iterindex >= 0);

   expr->iterdata[iterator->iterindex].userdata = userdata;
}

/** sets the iterator specific user data of the current expressions current child
 *
 * @note The expression iterator mode must be in DFS mode and stage \ref SCIP_EXPRITER_VISITINGCHILD or \ref SCIP_EXPRITER_VISITEDCHILD
 */
void SCIPexpriterSetChildUserData(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPRITER_USERDATA userdata           /**< data to be stored in current child */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->dfsstage == SCIP_EXPRITER_VISITINGCHILD || iterator->dfsstage == SCIP_EXPRITER_VISITEDCHILD);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild >= 0);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild < iterator->curr->nchildren);

   iterator->curr->children[iterator->curr->iterdata[iterator->iterindex].currentchild]->iterdata[iterator->iterindex].userdata = userdata;
}

/** moves the iterator to the next expression according to the mode of the expression iterator
 *
 * @return the next expression, if any, and NULL otherwise
 */
SCIP_EXPR* SCIPexpriterGetNext(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   /* move to the next expression according to iterator type */
   switch( iterator->itertype )
   {
      case SCIP_EXPRITER_BFS:
      {
         iterator->curr = doBfsNext(iterator);
         break;
      }

      case SCIP_EXPRITER_RTOPOLOGIC :
      {
         iterator->curr = doReverseTopologicalNext(iterator);
         if( iterator->visitedtag != 0 )
         {
            assert(iterator->iterindex >= 0);
            assert(iterator->iterindex < SCIP_EXPRITER_MAXNACTIVE);

            /* skip already visited expressions */
            while( iterator->curr != NULL )
            {
               if( iterator->curr->iterdata[iterator->iterindex].visitedtag == iterator->visitedtag )
               {
                  /* if curr has already been visited, get next one
                   * TODO this isn't really efficient, since we still walk through already visited expressions
                   */
                  iterator->curr = doReverseTopologicalNext(iterator);
               }
               else
               {
                  /* curr has not been visited yet, so mark it as visited and interrupt loop */
                  iterator->curr->iterdata[iterator->iterindex].visitedtag = iterator->visitedtag;
                  break;
               }
            }
         }
         break;
      }

      case SCIP_EXPRITER_DFS :
      {
         assert(iterator->iterindex >= 0);

         /* get next until we are in a stopstage again
          * this might give expressions more than once, depending on what the stopstages are
          */
         do
         {
            iterator->curr = doDfsNext(iterator);
         }
         while( iterator->curr != NULL && (iterator->dfsstage & iterator->stopstages) == 0 );

         break;
      }
   }

   return iterator->curr;
}

/** moves a DFS iterator to one of the next expressions
 *
 * - If in \ref SCIP_EXPRITER_ENTEREXPR stage, then all children of that expression will be skipped.
 *   If \ref SCIP_EXPRITER_LEAVEEXPR is one of the `stopstages`, then it will be the next stage. Otherwise, the iterator will move further on (go to the parent, etc).
 * - If in \ref SCIP_EXPRITER_VISITINGCHILD stage, then the child that was going to be visited next will be skipped and the iterator will be moved on to the next child (if any).
 * - If in \ref SCIP_EXPRITER_VISITEDCHILD stage, then all remaining children will be skipped and we move on to the \ref SCIP_EXPRITER_LEAVEEXPR stage (if a stop stage, otherwise further on).
 * - It is not allowed to call this function when in \ref SCIP_EXPRITER_LEAVEEXPR stage.
 *
 * @return the next expression, if any, and NULL otherwise
 */
SCIP_EXPR* SCIPexpriterSkipDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->iterindex >= 0);

   switch( iterator->dfsstage )
   {
      case SCIP_EXPRITER_ENTEREXPR :
      case SCIP_EXPRITER_VISITEDCHILD :
      {
         /* move directly to leaveexpr */
         iterator->dfsstage = SCIP_EXPRITER_LEAVEEXPR;
         /* if leaveexpr is not a stopstage, then move on */
         while( iterator->curr != NULL && (iterator->dfsstage & iterator->stopstages) == 0 )
            iterator->curr = doDfsNext(iterator);
         return iterator->curr;
      }

      case SCIP_EXPRITER_VISITINGCHILD :
      {
         /* skip the child to be visited */
         /* pretend we just visited this child and get next */
         iterator->dfsstage = SCIP_EXPRITER_VISITEDCHILD;
         return SCIPexpriterGetNext(iterator);
      }

      case SCIP_EXPRITER_LEAVEEXPR :
      default :
         SCIPerrorMessage("SCIPexpriterSkipDFS called in invalid stage %u", iterator->dfsstage);
         SCIPABORT();
         return iterator->curr;
   }
}

/** returns whether the iterator visited all expressions already */
SCIP_Bool SCIPexpriterIsEnd(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);

   return iterator->curr == NULL;
}
