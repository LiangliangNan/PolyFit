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

/**@file   rbtree.c
 * @brief  intrusive red black tree datastructure
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/rbtree.h"

#define RED              ((uintptr_t)0x1u)
#define BLACK            ((uintptr_t)0x0u)
#define COLOR(node)      ((node)->parent & RED)
#define IS_RED(node)     ( (node) != NULL && COLOR(node) )
#define IS_BLACK(node)   ( (node) == NULL || !COLOR(node) )
#define MAKE_RED(node)   do { (node)->parent |= RED; } while(0)
#define MAKE_BLACK(node) do { (node)->parent &= ~RED; } while(0)
#define LEFT             0
#define RIGHT            1
#define OPPOSITE(dir)    ( 1 - (dir) )
#define PARENT(node)     ( (SCIP_RBTREENODE*)((node)->parent & ~RED) )
#define SET_PARENT(n, p) do { (n)->parent = (uintptr_t)(p) | COLOR(n); } while(0)
#define SET_COLOR(n, c)  do { if( c == RED ) { MAKE_RED(n); } else { MAKE_BLACK(n); } } while(0)

/* functions for red black tree management; see Cormen, Thomas H. Introduction to algorithms. MIT press, 2009. */

/** rotate the tree in the given direction */
static
void rbRotate(
   SCIP_RBTREENODE**     root,               /**< pointer to store new root of the tree */
   SCIP_RBTREENODE*      x,                  /**< node to perform rotation for */
   int                   dir                 /**< direction of rotation */
   )
{
   SCIP_RBTREENODE* p;
   SCIP_RBTREENODE* y = x->child[OPPOSITE(dir)];
   x->child[OPPOSITE(dir)] = y->child[dir];
   if( y->child[dir] != NULL )
   {
      SET_PARENT(y->child[dir], x);
   }

   p = PARENT(x);
   SET_PARENT(y, p);

   if( p == NULL )
      *root = y;
   else if( x == p->child[dir] )
      p->child[dir] = y;
   else
      p->child[OPPOSITE(dir)] = y;

   y->child[dir] = x;
   SET_PARENT(x, y);
}

/** restores the red-black tree property after an insert */
static
void rbInsertFixup(
   SCIP_RBTREENODE**     root,               /**< pointer to store new root of the tree */
   SCIP_RBTREENODE*      z                   /**< inserted node */
   )
{
   SCIP_RBTREENODE* p;
   p = PARENT(z);

   while( IS_RED(p) )
   {
      SCIP_RBTREENODE* pp;
      SCIP_RBTREENODE* y;
      int dir;

      pp = PARENT(p);
      dir = p == pp->child[LEFT] ? RIGHT : LEFT;

      y = pp->child[dir];
      if( IS_RED(y) )
      {
         MAKE_BLACK(p);
         MAKE_BLACK(y);
         MAKE_RED(pp);
         z = pp;
      }
      else
      {
         if( z == p->child[dir] )
         {
            z = p;
            rbRotate(root, z, OPPOSITE(dir));
            p = PARENT(z);
            pp = PARENT(p);
         }

         MAKE_BLACK(p);
         MAKE_RED(pp);
         rbRotate(root, pp, dir);
      }

      p = PARENT(z);
   }

   MAKE_BLACK(*root);
}

/** restores the red-black tree property after an insert */
static
void rbDeleteFixup(
   SCIP_RBTREENODE**     root,               /**< pointer to store new root of the tree */
   SCIP_RBTREENODE*      x,                  /**< start node for fixup */
   SCIP_RBTREENODE*      nil                 /**< fake node representing NULL to properly reassemble the tree */
   )
{
   while( x != *root && IS_BLACK(x) )
   {
      SCIP_RBTREENODE* p;
      SCIP_RBTREENODE* w;
      int dir;

      p = PARENT(x == NULL ? nil : x);
      dir = x == p->child[LEFT] ? RIGHT : LEFT;

      w = p->child[dir];
      assert(w != NULL);

      if( COLOR(w) == RED )
      {
         MAKE_BLACK(w);
         MAKE_RED(p);
         rbRotate(root, p, OPPOSITE(dir));
         assert(p == PARENT(x == NULL ? nil : x));
         w = p->child[dir];
         assert(w != NULL);
      }

      if( IS_BLACK(w->child[LEFT]) && IS_BLACK(w->child[RIGHT]) )
      {
         MAKE_RED(w);
         x = p;
      }
      else
      {
         if( IS_BLACK(w->child[dir]) )
         {
            MAKE_BLACK(w->child[OPPOSITE(dir)]);
            MAKE_RED(w);
            rbRotate(root, w, dir);
            assert(p == PARENT(x == NULL ? nil : x));
            w = p->child[dir];
         }
         SET_COLOR(w, COLOR(p));
         MAKE_BLACK(p);
         MAKE_BLACK(w->child[dir]);
         rbRotate(root, p, OPPOSITE(dir));
         x = *root;
      }
   }

   if( x != NULL )
   {
      MAKE_BLACK(x);
   }
}

/** replaces the subtree rooted at node u with the subtree rooted at node v */
static
void rbTransplant(
   SCIP_RBTREENODE**     root,               /**< pointer to store the new root */
   SCIP_RBTREENODE*      u,                  /**< node u */
   SCIP_RBTREENODE*      v,                  /**< node v */
   SCIP_RBTREENODE*      nil                 /**< fake node representing NULL to properly reassemble the tree */
   )
{
   SCIP_RBTREENODE* up;

   up = PARENT(u);

   if( up == NULL )
      *root = v;
   else if( u == up->child[LEFT] )
      up->child[LEFT] = v;
   else
      up->child[RIGHT] = v;

   if( v == NULL )
      v = nil;

   SET_PARENT(v, up);
}

/** get first element in tree with respect to sorting key */
SCIP_RBTREENODE* SCIPrbtreeFirst_call(
   SCIP_RBTREENODE*      root                /**< root of the tree */
   )
{
   if( root == NULL )
      return NULL;

   while(root->child[LEFT] != NULL)
      root = root->child[LEFT];

   return root;
}

/** get last element in tree with respect to sorting key */
SCIP_RBTREENODE* SCIPrbtreeLast_call(
   SCIP_RBTREENODE*      root                /**< root of the tree */
   )
{
   if( root == NULL )
      return NULL;

   while(root->child[RIGHT] != NULL)
      root = root->child[RIGHT];

   return root;
}

/** get successor of given element in the tree */
SCIP_RBTREENODE* SCIPrbtreeSuccessor_call(
   SCIP_RBTREENODE*      x                   /**< element to get successor for */
   )
{
   SCIP_RBTREENODE* y;
   if( x->child[RIGHT] != NULL )
      return SCIPrbtreeFirst_call(x->child[RIGHT]);

   y = PARENT(x);

   while( y != NULL && x == y->child[RIGHT] )
   {
      x = y;
      y = PARENT(y);
   }

   return y;
}

/** get predecessor of given element in the tree */
SCIP_RBTREENODE* SCIPrbtreePredecessor_call(
   SCIP_RBTREENODE*      x                   /**< element to get predecessor for */
   )
{
   SCIP_RBTREENODE* y;
   if( x->child[LEFT] != NULL )
      return SCIPrbtreeLast_call(x->child[LEFT]);

   y = PARENT(x);

   while( y != NULL && x == y->child[LEFT] )
   {
      x = y;
      y = PARENT(y);
   }

   return y;
}

/** delete the given node from the tree given by it's root node.
 *  The node must be contained in the tree rooted at root.
 */
void SCIPrbtreeDelete_call(
   SCIP_RBTREENODE**     root,               /**< root of the tree */
   SCIP_RBTREENODE*      node                /**< node to delete from tree */
   )
{
   SCIP_RBTREENODE nil;
   SCIP_RBTREENODE* y;
   SCIP_RBTREENODE* x;
   unsigned int yorigcolor;

   nil.parent = 0;

   y = node;
   yorigcolor = COLOR(y);

   if( node->child[LEFT] == NULL )
   {
      x = node->child[RIGHT];
      rbTransplant(root, node, x, &nil);
   }
   else if( node->child[RIGHT] == NULL )
   {
      x = node->child[LEFT];
      rbTransplant(root, node, x, &nil);
   }
   else
   {
      y = SCIPrbtreeFirst(node->child[RIGHT]);
      yorigcolor = COLOR(y);
      x = y->child[RIGHT];
      if( PARENT(y) == node )
      {
         SET_PARENT(x == NULL ? &nil : x, y);
      }
      else
      {
         rbTransplant(root, y, y->child[RIGHT], &nil);
         y->child[RIGHT] = node->child[RIGHT];
         SET_PARENT(y->child[RIGHT], y);
      }
      rbTransplant(root, node, y, &nil);
      y->child[LEFT] = node->child[LEFT];
      SET_PARENT(y->child[LEFT], y);
      SET_COLOR(y, COLOR(node));
   }

   if( yorigcolor == BLACK )
      rbDeleteFixup(root, x, &nil);
}

/** insert node into the tree given by it's root. Requires the
 *  future parent and the position of the parent as returned by
 *  the tree's find function defined using the SCIP_DEF_RBTREE_FIND
 *  macro.
 */
void SCIPrbtreeInsert_call(
   SCIP_RBTREENODE**     root,               /**< root of the tree */
   SCIP_RBTREENODE*      parent,             /**< future parent of node to be inserted */
   int                   pos,                /**< position of parent with respect to node, i.e.
                                              *   -1 if the parent's key comes before node and 1
                                              *   if the parent's key comes after the node
                                              */
   SCIP_RBTREENODE*      node                /**< node to insert into the tree */
   )
{
   SET_PARENT(node, parent);
   if( parent == NULL )
      *root = node;
   else if( pos > 0 )
      parent->child[LEFT] = node;
   else
      parent->child[RIGHT] = node;

   node->child[LEFT] = NULL;
   node->child[RIGHT] = NULL;
   MAKE_RED(node);
   rbInsertFixup(root, node);
}
