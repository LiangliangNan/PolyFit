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

/**@file   var.c
 * @brief  methods for problem variables
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Stefan Heinz
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 * @author Stefan Vigerske
 *
 * @todo Possibly implement the access of bounds of multi-aggregated variables by accessing the
 * corresponding linear constraint if it exists. This seems to require some work, since the linear
 * constraint has to be stored. Moreover, it has even to be created in case the original constraint
 * was deleted after multi-aggregation, but the bounds of the multi-aggregated variable should be
 * changed. This has to be done with care in order to not loose the performance gains of
 * multi-aggregation.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/prop.h"
#include "scip/relax.h"
#include "scip/var.h"
#include "scip/cons.h"
#include "scip/event.h"
#include "scip/history.h"
#include "scip/implics.h"
#include "scip/lp.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/stat.h"
#include "scip/tree.h"
#include "scip/reopt.h"

#include "scip/debug.h"

#include "scip/pub_message.h"
#include "scip/pub_history.h"

#define MAXIMPLSCLOSURE 100  /**< maximal number of descendants of implied variable for building closure
                              *   in implication graph */
#define MAXABSVBCOEF    1e+5 /**< maximal absolute coefficient in variable bounds added due to implications */

/*
 * hole, holelist, and domain methods
 */

/** creates a new holelist element */
static
SCIP_RETCODE holelistCreate(
   SCIP_HOLELIST**       holelist,           /**< pointer to holelist to create */
   BMS_BLKMEM*           blkmem,             /**< block memory for target holelist */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right               /**< right bound of open interval in new hole */
   )
{
   assert(holelist != NULL);
   assert(blkmem != NULL);
   assert(SCIPsetIsLT(set, left, right));

   SCIPsetDebugMsg(set, "create hole list element (%.15g,%.15g) in blkmem %p\n", left, right, (void*)blkmem);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, holelist) );
   (*holelist)->hole.left = left;
   (*holelist)->hole.right = right;
   (*holelist)->next = NULL;

   return SCIP_OKAY;
}

/** frees all elements in the holelist */
static
void holelistFree(
   SCIP_HOLELIST**       holelist,           /**< pointer to holelist to free */
   BMS_BLKMEM*           blkmem              /**< block memory for target holelist */
   )
{
   assert(holelist != NULL);
   assert(blkmem != NULL);

   while( *holelist != NULL )
   {
      SCIP_HOLELIST* next;

      SCIPdebugMessage("free hole list element (%.15g,%.15g) in blkmem %p\n", 
         (*holelist)->hole.left, (*holelist)->hole.right, (void*)blkmem);

      next = (*holelist)->next;
      BMSfreeBlockMemory(blkmem, holelist);
      assert(*holelist == NULL);

      *holelist = next;
   }
   assert(*holelist == NULL);
}

/** duplicates a list of holes */
static
SCIP_RETCODE holelistDuplicate(
   SCIP_HOLELIST**       target,             /**< pointer to target holelist */
   BMS_BLKMEM*           blkmem,             /**< block memory for target holelist */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_HOLELIST*        source              /**< holelist to duplicate */
   )
{
   assert(target != NULL);

   while( source != NULL )
   {
      assert(source->next == NULL || SCIPsetIsGE(set, source->next->hole.left, source->hole.right));
      SCIP_CALL( holelistCreate(target, blkmem, set, source->hole.left, source->hole.right) );
      source = source->next;
      target = &(*target)->next;
   }

   return SCIP_OKAY;
}

/** adds a hole to the domain */
static
SCIP_RETCODE domAddHole(
   SCIP_DOM*             dom,                /**< domain to add hole to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right,              /**< right bound of open interval in new hole */
   SCIP_Bool*            added               /**< pointer to store whether the hole was added (variable didn't had that hole before), or NULL */
   )
{
   SCIP_HOLELIST** insertpos;
   SCIP_HOLELIST* next;

   assert(dom != NULL);
   assert(added != NULL);

   /* search for the position of the new hole */
   insertpos = &dom->holelist;
   while( *insertpos != NULL && (*insertpos)->hole.left < left )
      insertpos = &(*insertpos)->next;

   /* check if new hole already exists in the hole list or is a sub hole of an existing one */
   if( *insertpos != NULL && (*insertpos)->hole.left == left && (*insertpos)->hole.right >= right )  /*lint !e777 */
   {
      SCIPsetDebugMsg(set, "new hole (%.15g,%.15g) is redundant through known hole (%.15g,%.15g)\n",
         left, right, (*insertpos)->hole.left, (*insertpos)->hole.right);
      *added = FALSE;
      return SCIP_OKAY;
   }

   /* add hole */
   *added = TRUE;

   next = *insertpos;
   SCIP_CALL( holelistCreate(insertpos, blkmem, set, left, right) );
   (*insertpos)->next = next;

   return SCIP_OKAY;
}

/** merges overlapping holes into single holes, computes and moves lower and upper bound, respectively */
/**@todo  the domMerge() method is currently called if a lower or an upper bound locally or globally changed; this could
 *        be more efficient if performed with the knowledge if it was a lower or an upper bound which triggered this
 *        merge */
static
void domMerge(
   SCIP_DOM*             dom,                /**< domain to merge */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            newlb,              /**< pointer to store new lower bound */
   SCIP_Real*            newub               /**< pointer to store new upper bound */
   )
{
   SCIP_HOLELIST** holelistptr;
   SCIP_HOLELIST** lastnextptr;
   SCIP_Real* lastrightptr;

   assert(dom != NULL);
   assert(SCIPsetIsLE(set, dom->lb, dom->ub));

#ifndef NDEBUG
   {
      /* check if the holelist is sorted w.r.t. to the left interval bounds */
      SCIP_Real lastleft;

      holelistptr = &dom->holelist;

      lastleft = -SCIPsetInfinity(set);

      while( *holelistptr != NULL )
      {
         if( (*holelistptr)->next != NULL )
         {
            assert( SCIPsetIsLE(set, lastleft, (*holelistptr)->hole.left) );
            lastleft = (*holelistptr)->hole.left;
         }

         holelistptr = &(*holelistptr)->next;
      }   
   }
#endif

   SCIPsetDebugMsg(set, "merge hole list\n");

   holelistptr = &dom->holelist;
   lastrightptr = &dom->lb;  /* lower bound is the right bound of the hole (-infinity,lb) */
   lastnextptr = holelistptr;

   while( *holelistptr != NULL )
   {
      SCIPsetDebugMsg(set, "check hole (%.15g,%.15g) last right interval was <%.15g>\n", (*holelistptr)->hole.left, (*holelistptr)->hole.right, *lastrightptr);

      /* check that the hole is not empty */
      assert(SCIPsetIsLT(set, (*holelistptr)->hole.left, (*holelistptr)->hole.right));

      if( SCIPsetIsGE(set, (*holelistptr)->hole.left, dom->ub) )
      {
         /* the remaining holes start behind the upper bound: remove them */
         SCIPsetDebugMsg(set, "remove remaining hole since upper bound <%.15g> is less then the left hand side of the current hole\n", dom->ub);
         holelistFree(holelistptr, blkmem);
         assert(*holelistptr == NULL);

         /* unlink this hole from the previous hole */
         *lastnextptr = NULL;
      }
      else if( SCIPsetIsGT(set, (*holelistptr)->hole.right, dom->ub) )
      {
         /* the hole overlaps the upper bound: decrease upper bound, remove this hole and all remaining holes */
         SCIPsetDebugMsg(set, "upper bound <%.15g> lays in current hole; store new upper bound and remove this and all remaining holes\n", dom->ub);

         assert(SCIPsetIsLT(set, (*holelistptr)->hole.left, dom->ub));

         /* adjust upper bound */
         dom->ub = (*holelistptr)->hole.left;

         if(newub != NULL )
            *newub = (*holelistptr)->hole.left;

         /* remove remaining hole list */
         holelistFree(holelistptr, blkmem);
         assert(*holelistptr == NULL);

         /* unlink this hole from the previous hole */
         *lastnextptr = NULL;
      }
      else if( SCIPsetIsGT(set, *lastrightptr, (*holelistptr)->hole.left) )
      {
         /* the right bound of the last hole is greater than the left bound of this hole: increase the right bound of
          * the last hole, delete this hole */
         SCIP_HOLELIST* nextholelist;

         if( SCIPsetIsEQ(set, *lastrightptr, dom->lb ) )
         {
            /* the reason for the overlap results from the lower bound hole (-infinity,lb); therefore, we can increase
             * the lower bound */
            SCIPsetDebugMsg(set, "lower bound <%.15g> lays in current hole; store new lower bound and remove hole\n", dom->lb);
            *lastrightptr = MAX(*lastrightptr, (*holelistptr)->hole.right);

            /* adjust lower bound */
            dom->lb = *lastrightptr;

            if(newlb != NULL )
               *newlb = *lastrightptr;
         }
         else
         {
            SCIPsetDebugMsg(set, "current hole overlaps with the previous one (...,%.15g); merge to (...,%.15g)\n",
               *lastrightptr, MAX(*lastrightptr, (*holelistptr)->hole.right) );
            *lastrightptr = MAX(*lastrightptr, (*holelistptr)->hole.right);
         }
         nextholelist = (*holelistptr)->next;
         (*holelistptr)->next = NULL;
         holelistFree(holelistptr, blkmem);

         /* connect the linked list after removing the hole */
         *lastnextptr = nextholelist;

         /* get next hole */
         *holelistptr = nextholelist;
      }
      else
      {
         /* the holes do not overlap: update lastholelist and lastrightptr */
         lastrightptr = &(*holelistptr)->hole.right;
         lastnextptr = &(*holelistptr)->next;

         /* get next hole */
         holelistptr = &(*holelistptr)->next;
      }
   }

#ifndef NDEBUG
   {
      /* check that holes are merged */
      SCIP_Real lastright;

      lastright = dom->lb; /* lower bound is the right bound of the hole (-infinity,lb) */
      holelistptr = &dom->holelist;

      while( *holelistptr != NULL )
      {
         /* check the the last right interval is smaller or equal to the current left interval (none overlapping) */
         assert( SCIPsetIsLE(set, lastright, (*holelistptr)->hole.left) );

         /* check the hole property (check that the hole is not empty) */
         assert( SCIPsetIsLT(set, (*holelistptr)->hole.left, (*holelistptr)->hole.right) );
         lastright = (*holelistptr)->hole.right;

         /* get next hole */
         holelistptr = &(*holelistptr)->next;
      }   

      /* check the the last right interval is smaller or equal to the upper bound (none overlapping) */     
      assert( SCIPsetIsLE(set, lastright, dom->ub) );
   }
#endif
}

/*
 * domain change methods
 */

/** ensures, that bound change info array for lower bound changes can store at least num entries */
static
SCIP_RETCODE varEnsureLbchginfosSize(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(var != NULL);
   assert(var->nlbchginfos <= var->lbchginfossize);
   assert(SCIPvarIsTransformed(var));

   if( num > var->lbchginfossize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &var->lbchginfos, var->lbchginfossize, newsize) );
      var->lbchginfossize = newsize;
   }
   assert(num <= var->lbchginfossize);

   return SCIP_OKAY;
}

/** ensures, that bound change info array for upper bound changes can store at least num entries */
static
SCIP_RETCODE varEnsureUbchginfosSize(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(var != NULL);
   assert(var->nubchginfos <= var->ubchginfossize);
   assert(SCIPvarIsTransformed(var));

   if( num > var->ubchginfossize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &var->ubchginfos, var->ubchginfossize, newsize) );
      var->ubchginfossize = newsize;
   }
   assert(num <= var->ubchginfossize);

   return SCIP_OKAY;
}

/** adds domain change info to the variable's lower bound change info array */
static
SCIP_RETCODE varAddLbchginfo(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound,           /**< new value for bound */
   int                   depth,              /**< depth in the tree, where the bound change takes place */
   int                   pos,                /**< position of the bound change in its bound change array */
   SCIP_VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   SCIP_CONS*            infercons,          /**< constraint that infered this bound change, or NULL */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_BOUNDTYPE        inferboundtype,     /**< type of bound for inference var: lower or upper bound */
   SCIP_BOUNDCHGTYPE     boundchgtype        /**< bound change type: branching decision or infered bound change */
   )
{
   assert(var != NULL);
   assert(SCIPsetIsLT(set, oldbound, newbound));
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, oldbound));
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));
   assert(!SCIPvarIsBinary(var) || SCIPsetIsEQ(set, oldbound, 0.0));
   assert(!SCIPvarIsBinary(var) || SCIPsetIsEQ(set, newbound, 1.0));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING || infervar != NULL);
   assert((boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER) == (infercons != NULL));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER || inferprop == NULL);

   SCIPsetDebugMsg(set, "adding lower bound change info to var <%s>[%g,%g]: depth=%d, pos=%d, infer%s=<%s>, inferinfo=%d, %g -> %g\n",
      SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, depth, pos, infercons != NULL ? "cons" : "prop",
      infercons != NULL ? SCIPconsGetName(infercons) : (inferprop != NULL ? SCIPpropGetName(inferprop) : "-"), inferinfo,
      oldbound, newbound);

   SCIP_CALL( varEnsureLbchginfosSize(var, blkmem, set, var->nlbchginfos+1) );
   var->lbchginfos[var->nlbchginfos].oldbound = oldbound;
   var->lbchginfos[var->nlbchginfos].newbound = newbound;
   var->lbchginfos[var->nlbchginfos].var = var;
   var->lbchginfos[var->nlbchginfos].bdchgidx.depth = depth;
   var->lbchginfos[var->nlbchginfos].bdchgidx.pos = pos;
   var->lbchginfos[var->nlbchginfos].pos = var->nlbchginfos; /*lint !e732*/
   var->lbchginfos[var->nlbchginfos].boundchgtype = boundchgtype; /*lint !e641*/
   var->lbchginfos[var->nlbchginfos].boundtype = SCIP_BOUNDTYPE_LOWER; /*lint !e641*/
   var->lbchginfos[var->nlbchginfos].redundant = FALSE;
   var->lbchginfos[var->nlbchginfos].inferboundtype = inferboundtype; /*lint !e641*/
   var->lbchginfos[var->nlbchginfos].inferencedata.var = infervar;
   var->lbchginfos[var->nlbchginfos].inferencedata.info = inferinfo;

   /**@note The "pos" data member of the bound change info has a size of 27 bits */
   assert(var->nlbchginfos < 1 << 27);

   switch( boundchgtype )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
      break;
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      assert(infercons != NULL);
      var->lbchginfos[var->nlbchginfos].inferencedata.reason.cons = infercons;
      break;
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      var->lbchginfos[var->nlbchginfos].inferencedata.reason.prop = inferprop;
      break;
   default:
      SCIPerrorMessage("invalid bound change type %d\n", boundchgtype);
      return SCIP_INVALIDDATA;
   }

   var->nlbchginfos++;

   assert(var->nlbchginfos < 2
      || SCIPbdchgidxIsEarlier(&var->lbchginfos[var->nlbchginfos-2].bdchgidx,
         &var->lbchginfos[var->nlbchginfos-1].bdchgidx));

   return SCIP_OKAY;
}

/** adds domain change info to the variable's upper bound change info array */
static
SCIP_RETCODE varAddUbchginfo(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound,           /**< new value for bound */
   int                   depth,              /**< depth in the tree, where the bound change takes place */
   int                   pos,                /**< position of the bound change in its bound change array */
   SCIP_VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   SCIP_CONS*            infercons,          /**< constraint that infered this bound change, or NULL */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_BOUNDTYPE        inferboundtype,     /**< type of bound for inference var: lower or upper bound */
   SCIP_BOUNDCHGTYPE     boundchgtype        /**< bound change type: branching decision or infered bound change */
   )
{
   assert(var != NULL);
   assert(SCIPsetIsGT(set, oldbound, newbound));
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, oldbound));
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));
   assert(!SCIPvarIsBinary(var) || SCIPsetIsEQ(set, oldbound, 1.0));
   assert(!SCIPvarIsBinary(var) || SCIPsetIsEQ(set, newbound, 0.0));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING || infervar != NULL);
   assert((boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER) == (infercons != NULL));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER || inferprop == NULL);

   SCIPsetDebugMsg(set, "adding upper bound change info to var <%s>[%g,%g]: depth=%d, pos=%d, infer%s=<%s>, inferinfo=%d, %g -> %g\n",
      SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, depth, pos, infercons != NULL ? "cons" : "prop",
      infercons != NULL ? SCIPconsGetName(infercons) : (inferprop != NULL ? SCIPpropGetName(inferprop) : "-"), inferinfo,
      oldbound, newbound);

   SCIP_CALL( varEnsureUbchginfosSize(var, blkmem, set, var->nubchginfos+1) );
   var->ubchginfos[var->nubchginfos].oldbound = oldbound;
   var->ubchginfos[var->nubchginfos].newbound = newbound;
   var->ubchginfos[var->nubchginfos].var = var;
   var->ubchginfos[var->nubchginfos].bdchgidx.depth = depth;
   var->ubchginfos[var->nubchginfos].bdchgidx.pos = pos;
   var->ubchginfos[var->nubchginfos].pos = var->nubchginfos; /*lint !e732*/
   var->ubchginfos[var->nubchginfos].boundchgtype = boundchgtype; /*lint !e641*/
   var->ubchginfos[var->nubchginfos].boundtype = SCIP_BOUNDTYPE_UPPER; /*lint !e641*/
   var->ubchginfos[var->nubchginfos].redundant = FALSE;
   var->ubchginfos[var->nubchginfos].inferboundtype = inferboundtype; /*lint !e641*/
   var->ubchginfos[var->nubchginfos].inferencedata.var = infervar;
   var->ubchginfos[var->nubchginfos].inferencedata.info = inferinfo;

   /**@note The "pos" data member of the bound change info has a size of 27 bits */
   assert(var->nubchginfos < 1 << 27);

   switch( boundchgtype )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
      break;
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      assert(infercons != NULL);
      var->ubchginfos[var->nubchginfos].inferencedata.reason.cons = infercons;
      break;
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      var->ubchginfos[var->nubchginfos].inferencedata.reason.prop = inferprop;
      break;
   default:
      SCIPerrorMessage("invalid bound change type %d\n", boundchgtype);
      return SCIP_INVALIDDATA;
   }

   var->nubchginfos++;

   assert(var->nubchginfos < 2
      || SCIPbdchgidxIsEarlier(&var->ubchginfos[var->nubchginfos-2].bdchgidx,
         &var->ubchginfos[var->nubchginfos-1].bdchgidx));

   return SCIP_OKAY;
}

/** applies single bound change */
SCIP_RETCODE SCIPboundchgApply(
   SCIP_BOUNDCHG*        boundchg,           /**< bound change to apply */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   depth,              /**< depth in the tree, where the bound change takes place */
   int                   pos,                /**< position of the bound change in its bound change array */
   SCIP_Bool*            cutoff              /**< pointer to store whether an infeasible bound change was detected */
   )
{
   SCIP_VAR* var;

   assert(boundchg != NULL);
   assert(stat != NULL);
   assert(depth > 0);
   assert(pos >= 0);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* ignore redundant bound changes */
   if( boundchg->redundant )
      return SCIP_OKAY;

   var = boundchg->var;
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(!SCIPvarIsIntegral(var) || SCIPsetIsIntegral(set, boundchg->newbound));

   /* apply bound change */
   switch( boundchg->boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      /* check, if the bound change is still active (could be replaced by inference due to repropagation of higher node) */
      if( SCIPsetIsGT(set, boundchg->newbound, var->locdom.lb) )
      {
         if( SCIPsetIsLE(set, boundchg->newbound, var->locdom.ub) )
         {
            /* add the bound change info to the variable's bound change info array */
            switch( boundchg->boundchgtype )
            {
            case SCIP_BOUNDCHGTYPE_BRANCHING:
               SCIPsetDebugMsg(set, " -> branching: new lower bound of <%s>[%g,%g]: %g\n",
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               SCIP_CALL( varAddLbchginfo(var, blkmem, set, var->locdom.lb, boundchg->newbound, depth, pos,
                     NULL, NULL, NULL, 0, SCIP_BOUNDTYPE_LOWER, SCIP_BOUNDCHGTYPE_BRANCHING) );
               stat->lastbranchvar = var;
               stat->lastbranchdir = SCIP_BRANCHDIR_UPWARDS;
               stat->lastbranchvalue = boundchg->newbound;
               break;

            case SCIP_BOUNDCHGTYPE_CONSINFER:
               assert(boundchg->data.inferencedata.reason.cons != NULL);
               SCIPsetDebugMsg(set, " -> constraint <%s> inference: new lower bound of <%s>[%g,%g]: %g\n",
                  SCIPconsGetName(boundchg->data.inferencedata.reason.cons),
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               SCIP_CALL( varAddLbchginfo(var, blkmem, set, var->locdom.lb, boundchg->newbound, depth, pos,
                     boundchg->data.inferencedata.var, boundchg->data.inferencedata.reason.cons, NULL,
                     boundchg->data.inferencedata.info,
                     (SCIP_BOUNDTYPE)(boundchg->inferboundtype), SCIP_BOUNDCHGTYPE_CONSINFER) );
               break;

            case SCIP_BOUNDCHGTYPE_PROPINFER:
               SCIPsetDebugMsg(set, " -> propagator <%s> inference: new lower bound of <%s>[%g,%g]: %g\n",
                  boundchg->data.inferencedata.reason.prop != NULL
                  ? SCIPpropGetName(boundchg->data.inferencedata.reason.prop) : "-",
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               SCIP_CALL( varAddLbchginfo(var, blkmem, set, var->locdom.lb, boundchg->newbound, depth, pos,
                     boundchg->data.inferencedata.var, NULL, boundchg->data.inferencedata.reason.prop,
                     boundchg->data.inferencedata.info,
                     (SCIP_BOUNDTYPE)(boundchg->inferboundtype), SCIP_BOUNDCHGTYPE_PROPINFER) );
               break;

            default:
               SCIPerrorMessage("invalid bound change type %d\n", boundchg->boundchgtype);
               return SCIP_INVALIDDATA;
            }

            /* change local bound of variable */
            SCIP_CALL( SCIPvarChgLbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, boundchg->newbound) );
         }
         else
         {
            SCIPsetDebugMsg(set, " -> cutoff: new lower bound of <%s>[%g,%g]: %g\n",
               SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
            *cutoff = TRUE;
            boundchg->redundant = TRUE; /* bound change has not entered the lbchginfos array of the variable! */
         }
      }
      else
      {
         /* mark bound change to be inactive */
         SCIPsetDebugMsg(set, " -> inactive %s: new lower bound of <%s>[%g,%g]: %g\n",
            (SCIP_BOUNDCHGTYPE)boundchg->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING ? "branching" : "inference",
            SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
         boundchg->redundant = TRUE;
      }
      break;

   case SCIP_BOUNDTYPE_UPPER:
      /* check, if the bound change is still active (could be replaced by inference due to repropagation of higher node) */
      if( SCIPsetIsLT(set, boundchg->newbound, var->locdom.ub) )
      {
         if( SCIPsetIsGE(set, boundchg->newbound, var->locdom.lb) )
         {
            /* add the bound change info to the variable's bound change info array */
            switch( boundchg->boundchgtype )
            {
            case SCIP_BOUNDCHGTYPE_BRANCHING:
               SCIPsetDebugMsg(set, " -> branching: new upper bound of <%s>[%g,%g]: %g\n",
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               SCIP_CALL( varAddUbchginfo(var, blkmem, set, var->locdom.ub, boundchg->newbound, depth, pos,
                     NULL, NULL, NULL, 0, SCIP_BOUNDTYPE_UPPER, SCIP_BOUNDCHGTYPE_BRANCHING) );
               stat->lastbranchvar = var;
               stat->lastbranchdir = SCIP_BRANCHDIR_DOWNWARDS;
               stat->lastbranchvalue = boundchg->newbound;
               break;

            case SCIP_BOUNDCHGTYPE_CONSINFER:
               assert(boundchg->data.inferencedata.reason.cons != NULL);
               SCIPsetDebugMsg(set, " -> constraint <%s> inference: new upper bound of <%s>[%g,%g]: %g\n",
                  SCIPconsGetName(boundchg->data.inferencedata.reason.cons),
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               SCIP_CALL( varAddUbchginfo(var, blkmem, set, var->locdom.ub, boundchg->newbound, depth, pos,
                     boundchg->data.inferencedata.var, boundchg->data.inferencedata.reason.cons, NULL,
                     boundchg->data.inferencedata.info,
                     (SCIP_BOUNDTYPE)(boundchg->inferboundtype), SCIP_BOUNDCHGTYPE_CONSINFER) );
               break;

            case SCIP_BOUNDCHGTYPE_PROPINFER:
               SCIPsetDebugMsg(set, " -> propagator <%s> inference: new upper bound of <%s>[%g,%g]: %g\n",
                  boundchg->data.inferencedata.reason.prop != NULL
                  ? SCIPpropGetName(boundchg->data.inferencedata.reason.prop) : "-",
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               SCIP_CALL( varAddUbchginfo(var, blkmem, set, var->locdom.ub, boundchg->newbound, depth, pos,
                     boundchg->data.inferencedata.var, NULL, boundchg->data.inferencedata.reason.prop,
                     boundchg->data.inferencedata.info,
                     (SCIP_BOUNDTYPE)(boundchg->inferboundtype), SCIP_BOUNDCHGTYPE_PROPINFER) );
               break;

            default:
               SCIPerrorMessage("invalid bound change type %d\n", boundchg->boundchgtype);
               return SCIP_INVALIDDATA;
            }

            /* change local bound of variable */
            SCIP_CALL( SCIPvarChgUbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, boundchg->newbound) );
         }
         else
         {
            SCIPsetDebugMsg(set, " -> cutoff: new upper bound of <%s>[%g,%g]: %g\n",
               SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
            *cutoff = TRUE;
            boundchg->redundant = TRUE; /* bound change has not entered the ubchginfos array of the variable! */
         }
      }
      else
      {
         /* mark bound change to be inactive */
         SCIPsetDebugMsg(set, " -> inactive %s: new upper bound of <%s>[%g,%g]: %g\n",
            (SCIP_BOUNDCHGTYPE)boundchg->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING ? "branching" : "inference",
            SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
         boundchg->redundant = TRUE;
      }
      break;

   default:
      SCIPerrorMessage("unknown bound type\n");
      return SCIP_INVALIDDATA;
   }

   /* update the branching and inference history */
   if( !boundchg->applied && !boundchg->redundant )
   {
      assert(var == boundchg->var);

      if( (SCIP_BOUNDCHGTYPE)boundchg->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
      {
         SCIP_CALL( SCIPvarIncNBranchings(var, blkmem, set, stat,
               (SCIP_BOUNDTYPE)boundchg->boundtype == SCIP_BOUNDTYPE_LOWER
               ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS, boundchg->newbound, depth) );
      }
      else if( stat->lastbranchvar != NULL )
      {
         /**@todo if last branching variable is unknown, retrieve it from the nodes' boundchg arrays */
         SCIP_CALL( SCIPvarIncInferenceSum(stat->lastbranchvar, blkmem, set, stat, stat->lastbranchdir, stat->lastbranchvalue, 1.0) );
      }
      boundchg->applied = TRUE;
   }

   return SCIP_OKAY;
}

/** undoes single bound change */
SCIP_RETCODE SCIPboundchgUndo(
   SCIP_BOUNDCHG*        boundchg,           /**< bound change to remove */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_VAR* var;

   assert(boundchg != NULL);
   assert(stat != NULL);

   /* ignore redundant bound changes */
   if( boundchg->redundant )
      return SCIP_OKAY;

   var = boundchg->var;
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   /* undo bound change: apply the previous bound change of variable */
   switch( boundchg->boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      var->nlbchginfos--;
      assert(var->nlbchginfos >= 0);
      assert(var->lbchginfos != NULL);
      assert( SCIPsetIsFeasEQ(set, var->lbchginfos[var->nlbchginfos].newbound, var->locdom.lb) ); /*lint !e777*/
      assert( SCIPsetIsFeasLE(set, boundchg->newbound, var->locdom.lb) ); /* current lb might be larger to intermediate global bound change */

      SCIPsetDebugMsg(set, "removed lower bound change info of var <%s>[%g,%g]: depth=%d, pos=%d, %g -> %g\n",
         SCIPvarGetName(var), var->locdom.lb, var->locdom.ub,
         var->lbchginfos[var->nlbchginfos].bdchgidx.depth, var->lbchginfos[var->nlbchginfos].bdchgidx.pos,
         var->lbchginfos[var->nlbchginfos].oldbound, var->lbchginfos[var->nlbchginfos].newbound);

      /* reinstall the previous local bound */
      SCIP_CALL( SCIPvarChgLbLocal(boundchg->var, blkmem, set, stat, lp, branchcand, eventqueue,
            var->lbchginfos[var->nlbchginfos].oldbound) );

      /* in case all bound changes are removed the local bound should match the global bound */
      assert(var->nlbchginfos > 0 || SCIPsetIsFeasEQ(set, var->locdom.lb, var->glbdom.lb));

      break;

   case SCIP_BOUNDTYPE_UPPER:
      var->nubchginfos--;
      assert(var->nubchginfos >= 0);
      assert(var->ubchginfos != NULL);
      assert( SCIPsetIsFeasEQ(set, var->ubchginfos[var->nubchginfos].newbound, var->locdom.ub) ); /*lint !e777*/
      assert( SCIPsetIsFeasGE(set, boundchg->newbound, var->locdom.ub) ); /* current ub might be smaller to intermediate global bound change */

      SCIPsetDebugMsg(set, "removed upper bound change info of var <%s>[%g,%g]: depth=%d, pos=%d, %g -> %g\n",
         SCIPvarGetName(var), var->locdom.lb, var->locdom.ub,
         var->ubchginfos[var->nubchginfos].bdchgidx.depth, var->ubchginfos[var->nubchginfos].bdchgidx.pos,
         var->ubchginfos[var->nubchginfos].oldbound, var->ubchginfos[var->nubchginfos].newbound);

      /* reinstall the previous local bound */
      SCIP_CALL( SCIPvarChgUbLocal(boundchg->var, blkmem, set, stat, lp, branchcand, eventqueue,
            var->ubchginfos[var->nubchginfos].oldbound) );

      /* in case all bound changes are removed the local bound should match the global bound */
      assert(var->nubchginfos > 0 || SCIPsetIsFeasEQ(set, var->locdom.ub, var->glbdom.ub));

      break;

   default:
      SCIPerrorMessage("unknown bound type\n");
      return SCIP_INVALIDDATA;
   }

   /* update last branching variable */
   if( (SCIP_BOUNDCHGTYPE)boundchg->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
   {
      stat->lastbranchvar = NULL;
      stat->lastbranchvalue = SCIP_UNKNOWN;
   }

   return SCIP_OKAY;
}

/** applies single bound change to the global problem by changing the global bound of the corresponding variable */
static
SCIP_RETCODE boundchgApplyGlobal(
   SCIP_BOUNDCHG*        boundchg,           /**< bound change to apply */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether an infeasible bound change was detected */
   )
{
   SCIP_VAR* var;
   SCIP_Real newbound;
   SCIP_BOUNDTYPE boundtype;

   assert(boundchg != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* ignore redundant bound changes */
   if( boundchg->redundant )
      return SCIP_OKAY;

   var = SCIPboundchgGetVar(boundchg);
   newbound = SCIPboundchgGetNewbound(boundchg);
   boundtype = SCIPboundchgGetBoundtype(boundchg);

   /* check if the bound change is redundant which can happen due to a (better) global bound change which was performed
    * after that bound change was applied
    *
    * @note a global bound change is not captured by the redundant member of the bound change data structure
    */
   if( (boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasLE(set, newbound, SCIPvarGetLbGlobal(var)))
      || (boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasGE(set, newbound, SCIPvarGetUbGlobal(var))) )
   {
      return SCIP_OKAY;
   }

   SCIPsetDebugMsg(set, "applying global bound change: <%s>[%g,%g] %s %g\n",
      SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
      boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", newbound);

   /* check for cutoff */
   if( (boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasGT(set, newbound, SCIPvarGetUbGlobal(var)))
      || (boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasLT(set, newbound, SCIPvarGetLbGlobal(var))) )
   {
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* apply bound change */
   SCIP_CALL( SCIPvarChgBdGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newbound, boundtype) );

   return SCIP_OKAY;
}

/** captures branching and inference data of bound change */
static
SCIP_RETCODE boundchgCaptureData(
   SCIP_BOUNDCHG*        boundchg            /**< bound change to remove */
   )
{
   assert(boundchg != NULL);

   /* capture variable associated with the bound change */
   assert(boundchg->var != NULL);
   SCIPvarCapture(boundchg->var);

   switch( boundchg->boundchgtype )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      break;

   case SCIP_BOUNDCHGTYPE_CONSINFER:
      assert(boundchg->data.inferencedata.var != NULL);
      assert(boundchg->data.inferencedata.reason.cons != NULL);
      SCIPconsCapture(boundchg->data.inferencedata.reason.cons);
      break;

   default:
      SCIPerrorMessage("invalid bound change type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** releases branching and inference data of bound change */
static
SCIP_RETCODE boundchgReleaseData(
   SCIP_BOUNDCHG*        boundchg,           /**< bound change to remove */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */

   )
{
   assert(boundchg != NULL);

   switch( boundchg->boundchgtype )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      break;

   case SCIP_BOUNDCHGTYPE_CONSINFER:
      assert(boundchg->data.inferencedata.var != NULL);
      assert(boundchg->data.inferencedata.reason.cons != NULL);
      SCIP_CALL( SCIPconsRelease(&boundchg->data.inferencedata.reason.cons, blkmem, set) );
      break;

   default:
      SCIPerrorMessage("invalid bound change type\n");
      return SCIP_INVALIDDATA;
   }

   /* release variable */
   assert(boundchg->var != NULL);
   SCIP_CALL( SCIPvarRelease(&boundchg->var, blkmem, set, eventqueue, lp) );


   return SCIP_OKAY;
}

/** creates empty domain change data with dynamic arrays */
static
SCIP_RETCODE domchgCreate(
   SCIP_DOMCHG**         domchg,             /**< pointer to domain change data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(domchg != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemorySize(blkmem, domchg, sizeof(SCIP_DOMCHGDYN)) );
   (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_DYNAMIC; /*lint !e641*/
   (*domchg)->domchgdyn.nboundchgs = 0;
   (*domchg)->domchgdyn.boundchgs = NULL;
   (*domchg)->domchgdyn.nholechgs = 0;
   (*domchg)->domchgdyn.holechgs = NULL;
   (*domchg)->domchgdyn.boundchgssize = 0;
   (*domchg)->domchgdyn.holechgssize = 0;

   return SCIP_OKAY;
}

/** frees domain change data */
SCIP_RETCODE SCIPdomchgFree(
   SCIP_DOMCHG**         domchg,             /**< pointer to domain change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(domchg != NULL);
   assert(blkmem != NULL);

   if( *domchg != NULL )
   {
      int i;

      /* release variables, branching and inference data associated with the bound changes */
      for( i = 0; i < (int)(*domchg)->domchgbound.nboundchgs; ++i )
      {
         SCIP_CALL( boundchgReleaseData(&(*domchg)->domchgbound.boundchgs[i], blkmem, set, eventqueue, lp) );
      }

      /* free memory for bound and hole changes */
      switch( (*domchg)->domchgdyn.domchgtype )
      {
      case SCIP_DOMCHGTYPE_BOUND:
         BMSfreeBlockMemoryArrayNull(blkmem, &(*domchg)->domchgbound.boundchgs, (*domchg)->domchgbound.nboundchgs);
         BMSfreeBlockMemorySize(blkmem, domchg, sizeof(SCIP_DOMCHGBOUND));
         break;
      case SCIP_DOMCHGTYPE_BOTH:
         BMSfreeBlockMemoryArrayNull(blkmem, &(*domchg)->domchgboth.boundchgs, (*domchg)->domchgboth.nboundchgs);
         BMSfreeBlockMemoryArrayNull(blkmem, &(*domchg)->domchgboth.holechgs, (*domchg)->domchgboth.nholechgs);
         BMSfreeBlockMemorySize(blkmem, domchg, sizeof(SCIP_DOMCHGBOTH));
         break;
      case SCIP_DOMCHGTYPE_DYNAMIC:
         BMSfreeBlockMemoryArrayNull(blkmem, &(*domchg)->domchgdyn.boundchgs, (*domchg)->domchgdyn.boundchgssize);
         BMSfreeBlockMemoryArrayNull(blkmem, &(*domchg)->domchgdyn.holechgs, (*domchg)->domchgdyn.holechgssize);
         BMSfreeBlockMemorySize(blkmem, domchg, sizeof(SCIP_DOMCHGDYN));
         break;
      default:
         SCIPerrorMessage("invalid domain change type\n");
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** converts a static domain change data into a dynamic one */
static
SCIP_RETCODE domchgMakeDynamic(
   SCIP_DOMCHG**         domchg,             /**< pointer to domain change data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(domchg != NULL);
   assert(blkmem != NULL);

   SCIPdebugMessage("making domain change data %p pointing to %p dynamic\n", (void*)domchg, (void*)*domchg);

   if( *domchg == NULL )
   {
      SCIP_CALL( domchgCreate(domchg, blkmem) );
   }
   else
   {
      switch( (*domchg)->domchgdyn.domchgtype )
      {
      case SCIP_DOMCHGTYPE_BOUND:
         SCIP_ALLOC( BMSreallocBlockMemorySize(blkmem, domchg, sizeof(SCIP_DOMCHGBOUND), sizeof(SCIP_DOMCHGDYN)) );
         (*domchg)->domchgdyn.nholechgs = 0;
         (*domchg)->domchgdyn.holechgs = NULL;
         (*domchg)->domchgdyn.boundchgssize = (*domchg)->domchgdyn.nboundchgs;
         (*domchg)->domchgdyn.holechgssize = 0;
         (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_DYNAMIC; /*lint !e641*/
         break;
      case SCIP_DOMCHGTYPE_BOTH:
         SCIP_ALLOC( BMSreallocBlockMemorySize(blkmem, domchg, sizeof(SCIP_DOMCHGBOTH), sizeof(SCIP_DOMCHGDYN)) );
         (*domchg)->domchgdyn.boundchgssize = (*domchg)->domchgdyn.nboundchgs;
         (*domchg)->domchgdyn.holechgssize = (*domchg)->domchgdyn.nholechgs;
         (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_DYNAMIC; /*lint !e641*/
         break;
      case SCIP_DOMCHGTYPE_DYNAMIC:
         break;
      default:
         SCIPerrorMessage("invalid domain change type\n");
         return SCIP_INVALIDDATA;
      }
   }
#ifndef NDEBUG
   {
      int i;
      for( i = 0; i < (int)(*domchg)->domchgbound.nboundchgs; ++i )
         assert(SCIPvarGetType((*domchg)->domchgbound.boundchgs[i].var) == SCIP_VARTYPE_CONTINUOUS
            || EPSISINT((*domchg)->domchgbound.boundchgs[i].newbound, 1e-06));
   }
#endif

   return SCIP_OKAY;
}

/** converts a dynamic domain change data into a static one, using less memory than for a dynamic one */
SCIP_RETCODE SCIPdomchgMakeStatic(
   SCIP_DOMCHG**         domchg,             /**< pointer to domain change data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(domchg != NULL);
   assert(blkmem != NULL);

   SCIPsetDebugMsg(set, "making domain change data %p pointing to %p static\n", (void*)domchg, (void*)*domchg);

   if( *domchg != NULL )
   {
      switch( (*domchg)->domchgdyn.domchgtype )
      {
      case SCIP_DOMCHGTYPE_BOUND:
         if( (*domchg)->domchgbound.nboundchgs == 0 )
         {
            SCIP_CALL( SCIPdomchgFree(domchg, blkmem, set, eventqueue, lp) );
         }
         break;
      case SCIP_DOMCHGTYPE_BOTH:
         if( (*domchg)->domchgboth.nholechgs == 0 )
         {
            if( (*domchg)->domchgbound.nboundchgs == 0 )
            {
               SCIP_CALL( SCIPdomchgFree(domchg, blkmem, set, eventqueue, lp) );
            }
            else
            {
               SCIP_ALLOC( BMSreallocBlockMemorySize(blkmem, domchg, sizeof(SCIP_DOMCHGBOTH), sizeof(SCIP_DOMCHGBOUND)) );
               (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_BOUND; /*lint !e641*/
            }
         }
         break;
      case SCIP_DOMCHGTYPE_DYNAMIC:
         if( (*domchg)->domchgboth.nholechgs == 0 )
         {
            if( (*domchg)->domchgbound.nboundchgs == 0 )
            {
               SCIP_CALL( SCIPdomchgFree(domchg, blkmem, set, eventqueue, lp) );
            }
            else
            {
               /* shrink dynamic size arrays to their minimal sizes */
               SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*domchg)->domchgdyn.boundchgs, \
                     (*domchg)->domchgdyn.boundchgssize, (*domchg)->domchgdyn.nboundchgs) ); /*lint !e571*/
               BMSfreeBlockMemoryArrayNull(blkmem, &(*domchg)->domchgdyn.holechgs, (*domchg)->domchgdyn.holechgssize);

               /* convert into static domain change */
               SCIP_ALLOC( BMSreallocBlockMemorySize(blkmem, domchg, sizeof(SCIP_DOMCHGDYN), sizeof(SCIP_DOMCHGBOUND)) );
               (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_BOUND; /*lint !e641*/
            }
         }
         else
         {
            /* shrink dynamic size arrays to their minimal sizes */
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*domchg)->domchgdyn.boundchgs, \
                  (*domchg)->domchgdyn.boundchgssize, (*domchg)->domchgdyn.nboundchgs) ); /*lint !e571*/
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*domchg)->domchgdyn.holechgs, \
                  (*domchg)->domchgdyn.holechgssize, (*domchg)->domchgdyn.nholechgs) );

            /* convert into static domain change */
            SCIP_ALLOC( BMSreallocBlockMemorySize(blkmem, domchg, sizeof(SCIP_DOMCHGDYN), sizeof(SCIP_DOMCHGBOTH)) );
            (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_BOTH; /*lint !e641*/
         }
         break;
      default:
         SCIPerrorMessage("invalid domain change type\n");
         return SCIP_INVALIDDATA;
      }
#ifndef NDEBUG
      if( *domchg != NULL )
      {
         int i;
         for( i = 0; i < (int)(*domchg)->domchgbound.nboundchgs; ++i )
            assert(SCIPvarGetType((*domchg)->domchgbound.boundchgs[i].var) == SCIP_VARTYPE_CONTINUOUS
               || SCIPsetIsFeasIntegral(set, (*domchg)->domchgbound.boundchgs[i].newbound));
      }
#endif
   }

   return SCIP_OKAY;
}

/** ensures, that boundchgs array can store at least num entries */
static
SCIP_RETCODE domchgEnsureBoundchgsSize(
   SCIP_DOMCHG*          domchg,             /**< domain change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(domchg != NULL);
   assert(domchg->domchgdyn.domchgtype == SCIP_DOMCHGTYPE_DYNAMIC); /*lint !e641*/

   if( num > domchg->domchgdyn.boundchgssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &domchg->domchgdyn.boundchgs, domchg->domchgdyn.boundchgssize, newsize) );
      domchg->domchgdyn.boundchgssize = newsize;
   }
   assert(num <= domchg->domchgdyn.boundchgssize);

   return SCIP_OKAY;
}

/** ensures, that holechgs array can store at least num additional entries */
static
SCIP_RETCODE domchgEnsureHolechgsSize(
   SCIP_DOMCHG*          domchg,             /**< domain change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of additional entries to store */
   )
{
   assert(domchg != NULL);
   assert(domchg->domchgdyn.domchgtype == SCIP_DOMCHGTYPE_DYNAMIC); /*lint !e641*/

   if( num > domchg->domchgdyn.holechgssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &domchg->domchgdyn.holechgs, domchg->domchgdyn.holechgssize, newsize) );
      domchg->domchgdyn.holechgssize = newsize;
   }
   assert(num <= domchg->domchgdyn.holechgssize);

   return SCIP_OKAY;
}

/** applies domain change */
SCIP_RETCODE SCIPdomchgApply(
   SCIP_DOMCHG*          domchg,             /**< domain change to apply */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   depth,              /**< depth in the tree, where the domain change takes place */
   SCIP_Bool*            cutoff              /**< pointer to store whether an infeasible domain change was detected */
   )
{
   int i;

   assert(cutoff != NULL);

   *cutoff = FALSE;

   SCIPsetDebugMsg(set, "applying domain changes at %p in depth %d\n", (void*)domchg, depth);

   if( domchg == NULL )
      return SCIP_OKAY;

   /* apply bound changes */
   for( i = 0; i < (int)domchg->domchgbound.nboundchgs; ++i )
   {
      SCIP_CALL( SCIPboundchgApply(&domchg->domchgbound.boundchgs[i], blkmem, set, stat, lp,
            branchcand, eventqueue, depth, i, cutoff) );
      if( *cutoff )
         break;
   }
   SCIPsetDebugMsg(set, " -> %u bound changes (cutoff %u)\n", domchg->domchgbound.nboundchgs, *cutoff);

   /* mark all bound changes after a cutoff redundant */
   for( ; i < (int)domchg->domchgbound.nboundchgs; ++i )
      domchg->domchgbound.boundchgs[i].redundant = TRUE;

   /* apply holelist changes */
   if( domchg->domchgdyn.domchgtype != SCIP_DOMCHGTYPE_BOUND ) /*lint !e641*/
   {
      for( i = 0; i < domchg->domchgboth.nholechgs; ++i )
         *(domchg->domchgboth.holechgs[i].ptr) = domchg->domchgboth.holechgs[i].newlist;
      SCIPsetDebugMsg(set, " -> %d hole changes\n", domchg->domchgboth.nholechgs);
   }

   return SCIP_OKAY;
}

/** undoes domain change */
SCIP_RETCODE SCIPdomchgUndo(
   SCIP_DOMCHG*          domchg,             /**< domain change to remove */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   int i;

   SCIPsetDebugMsg(set, "undoing domain changes at %p\n", (void*)domchg);
   if( domchg == NULL )
      return SCIP_OKAY;

   /* undo holelist changes */
   if( domchg->domchgdyn.domchgtype != SCIP_DOMCHGTYPE_BOUND ) /*lint !e641*/
   {
      for( i = domchg->domchgboth.nholechgs-1; i >= 0; --i )
         *(domchg->domchgboth.holechgs[i].ptr) = domchg->domchgboth.holechgs[i].oldlist;
      SCIPsetDebugMsg(set, " -> %d hole changes\n", domchg->domchgboth.nholechgs);
   }

   /* undo bound changes */
   for( i = domchg->domchgbound.nboundchgs-1; i >= 0; --i )
   {
      SCIP_CALL( SCIPboundchgUndo(&domchg->domchgbound.boundchgs[i], blkmem, set, stat, lp, branchcand, eventqueue) );
   }
   SCIPsetDebugMsg(set, " -> %u bound changes\n", domchg->domchgbound.nboundchgs);

   return SCIP_OKAY;
}

/** applies domain change to the global problem */
SCIP_RETCODE SCIPdomchgApplyGlobal(
   SCIP_DOMCHG*          domchg,             /**< domain change to apply */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether an infeasible domain change was detected */
   )
{
   int i;

   assert(cutoff != NULL);

   *cutoff = FALSE;

   if( domchg == NULL )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "applying domain changes at %p to the global problem\n", (void*)domchg);

   /* apply bound changes */
   for( i = 0; i < (int)domchg->domchgbound.nboundchgs; ++i )
   {
      SCIP_CALL( boundchgApplyGlobal(&domchg->domchgbound.boundchgs[i], blkmem, set, stat, lp,
            branchcand, eventqueue, cliquetable, cutoff) );
      if( *cutoff )
         break;
   }
   SCIPsetDebugMsg(set, " -> %u global bound changes\n", domchg->domchgbound.nboundchgs);

   /**@todo globally apply holelist changes - how can this be done without confusing pointer updates? */

   return SCIP_OKAY;
}

/** adds bound change to domain changes */
SCIP_RETCODE SCIPdomchgAddBoundchg(
   SCIP_DOMCHG**         domchg,             /**< pointer to domain change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound for var: lower or upper bound */
   SCIP_BOUNDCHGTYPE     boundchgtype,       /**< type of bound change: branching decision or inference */
   SCIP_Real             lpsolval,           /**< solval of variable in last LP on path to node, or SCIP_INVALID if unknown */
   SCIP_VAR*             infervar,           /**< variable that was changed (parent of var, or var itself), or NULL */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change, or NULL */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_BOUNDTYPE        inferboundtype      /**< type of bound for inference var: lower or upper bound */
   )
{
   SCIP_BOUNDCHG* boundchg;

   assert(domchg != NULL);
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));
   assert(!SCIPvarIsBinary(var) || SCIPsetIsEQ(set, newbound, boundtype == SCIP_BOUNDTYPE_LOWER ? 1.0 : 0.0));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING || infervar != NULL);
   assert((boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER) == (infercons != NULL));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER || inferprop == NULL);

   SCIPsetDebugMsg(set, "adding %s bound change <%s: %g> of variable <%s> to domain change at %p pointing to %p\n",
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING ? "branching" : "inference",
      newbound, var->name, (void*)domchg, (void*)*domchg);

   /* if domain change data doesn't exist, create it;
    * if domain change is static, convert it into dynamic change
    */
   if( *domchg == NULL )
   {
      SCIP_CALL( domchgCreate(domchg, blkmem) );
   }
   else if( (*domchg)->domchgdyn.domchgtype != SCIP_DOMCHGTYPE_DYNAMIC ) /*lint !e641*/
   {
      SCIP_CALL( domchgMakeDynamic(domchg, blkmem) );
   }
   assert(*domchg != NULL && (*domchg)->domchgdyn.domchgtype == SCIP_DOMCHGTYPE_DYNAMIC); /*lint !e641*/

   /* get memory for additional bound change */
   SCIP_CALL( domchgEnsureBoundchgsSize(*domchg, blkmem, set, (*domchg)->domchgdyn.nboundchgs+1) );

   /* fill in the bound change data */
   boundchg = &(*domchg)->domchgdyn.boundchgs[(*domchg)->domchgdyn.nboundchgs];
   boundchg->var = var;
   switch( boundchgtype )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
      boundchg->data.branchingdata.lpsolval = lpsolval;
      break;
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      assert(infercons != NULL);
      boundchg->data.inferencedata.var = infervar;
      boundchg->data.inferencedata.reason.cons = infercons;
      boundchg->data.inferencedata.info = inferinfo; 
      break;
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      boundchg->data.inferencedata.var = infervar;
      boundchg->data.inferencedata.reason.prop = inferprop;
      boundchg->data.inferencedata.info = inferinfo; 
      break;
   default:
      SCIPerrorMessage("invalid bound change type %d\n", boundchgtype);
      return SCIP_INVALIDDATA;
   }

   boundchg->newbound = newbound;
   boundchg->boundchgtype = boundchgtype; /*lint !e641*/
   boundchg->boundtype = boundtype; /*lint !e641*/
   boundchg->inferboundtype = inferboundtype; /*lint !e641*/
   boundchg->applied = FALSE;
   boundchg->redundant = FALSE;
   (*domchg)->domchgdyn.nboundchgs++;

   /* capture branching and inference data associated with the bound changes */
   SCIP_CALL( boundchgCaptureData(boundchg) );

#ifdef SCIP_DISABLED_CODE /* expensive debug check */
#ifdef SCIP_MORE_DEBUG
   {
      int i;
      for( i = 0; i < (int)(*domchg)->domchgbound.nboundchgs; ++i )
         assert(SCIPvarGetType((*domchg)->domchgbound.boundchgs[i].var) == SCIP_VARTYPE_CONTINUOUS
            || SCIPsetIsFeasIntegral(set, (*domchg)->domchgbound.boundchgs[i].newbound));
   }
#endif
#endif

   return SCIP_OKAY;
}

/** adds hole change to domain changes */
SCIP_RETCODE SCIPdomchgAddHolechg(
   SCIP_DOMCHG**         domchg,             /**< pointer to domain change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_HOLELIST**       ptr,                /**< changed list pointer */
   SCIP_HOLELIST*        newlist,            /**< new value of list pointer */
   SCIP_HOLELIST*        oldlist             /**< old value of list pointer */
   )
{
   SCIP_HOLECHG* holechg;

   assert(domchg != NULL);
   assert(ptr != NULL);

   /* if domain change data doesn't exist, create it;
    * if domain change is static, convert it into dynamic change
    */
   if( *domchg == NULL )
   {
      SCIP_CALL( domchgCreate(domchg, blkmem) );
   }
   else if( (*domchg)->domchgdyn.domchgtype != SCIP_DOMCHGTYPE_DYNAMIC ) /*lint !e641*/
   {
      SCIP_CALL( domchgMakeDynamic(domchg, blkmem) );
   }
   assert(*domchg != NULL && (*domchg)->domchgdyn.domchgtype == SCIP_DOMCHGTYPE_DYNAMIC); /*lint !e641*/

   /* get memory for additional hole change */
   SCIP_CALL( domchgEnsureHolechgsSize(*domchg, blkmem, set, (*domchg)->domchgdyn.nholechgs+1) );

   /* fill in the hole change data */
   holechg = &(*domchg)->domchgdyn.holechgs[(*domchg)->domchgdyn.nholechgs];
   holechg->ptr = ptr;
   holechg->newlist = newlist;
   holechg->oldlist = oldlist;
   (*domchg)->domchgdyn.nholechgs++;

   return SCIP_OKAY;
}




/*
 * methods for variables 
 */

/** returns adjusted lower bound value, which is rounded for integral variable types */
static
SCIP_Real adjustedLb(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Real             lb                  /**< lower bound to adjust */
   )
{
   if( lb < 0 && SCIPsetIsInfinity(set, -lb) )
      return -SCIPsetInfinity(set);
   else if( lb > 0 && SCIPsetIsInfinity(set, lb) )
      return SCIPsetInfinity(set);
   else if( vartype != SCIP_VARTYPE_CONTINUOUS )
      return SCIPsetFeasCeil(set, lb);
   else if( SCIPsetIsZero(set, lb) )
      return 0.0;
   else
      return lb;
}

/** returns adjusted upper bound value, which is rounded for integral variable types */
static
SCIP_Real adjustedUb(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Real             ub                  /**< upper bound to adjust */
   )
{
   if( ub > 0 && SCIPsetIsInfinity(set, ub) )
      return SCIPsetInfinity(set);
   else if( ub < 0 && SCIPsetIsInfinity(set, -ub) )
      return -SCIPsetInfinity(set);
   else if( vartype != SCIP_VARTYPE_CONTINUOUS )
      return SCIPsetFeasFloor(set, ub);
   else if( SCIPsetIsZero(set, ub) )
      return 0.0;
   else
      return ub;
}

/** removes (redundant) cliques, implications and variable bounds of variable from all other variables' implications and variable
 *  bounds arrays, and optionally removes them also from the variable itself
 */
SCIP_RETCODE SCIPvarRemoveCliquesImplicsVbs(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             irrelevantvar,      /**< has the variable become irrelevant? */
   SCIP_Bool             onlyredundant,      /**< should only the redundant implications and variable bounds be removed? */
   SCIP_Bool             removefromvar       /**< should the implications and variable bounds be removed from the var itself? */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarIsActive(var) || SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);

   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);

   SCIPsetDebugMsg(set, "removing %s implications and vbounds of %s<%s>[%g,%g]\n",
      onlyredundant ? "redundant" : "all", irrelevantvar ? "irrelevant " : "", SCIPvarGetName(var), lb, ub);

   /* remove implications of (fixed) binary variable */
   if( var->implics != NULL && (!onlyredundant || lb > 0.5 || ub < 0.5) )
   {
      SCIP_Bool varfixing;

      assert(SCIPvarIsBinary(var));

      varfixing = FALSE;
      do
      {
         SCIP_VAR** implvars;
         SCIP_BOUNDTYPE* impltypes;
         int nimpls;
         int i;

         nimpls = SCIPimplicsGetNImpls(var->implics, varfixing);
         implvars = SCIPimplicsGetVars(var->implics, varfixing);
         impltypes = SCIPimplicsGetTypes(var->implics, varfixing);

         for( i = 0; i < nimpls; i++ )
         {
            SCIP_VAR* implvar;
            SCIP_BOUNDTYPE impltype;

            implvar = implvars[i];
            impltype = impltypes[i];
            assert(implvar != var);

            /* remove for all implications z == 0 / 1  ==>  x <= p / x >= p (x not binary)
             * the following variable bound from x's variable bounds 
             *   x <= b*z+d (z in vubs of x)            , for z == 0 / 1  ==>  x <= p
             *   x >= b*z+d (z in vlbs of x)            , for z == 0 / 1  ==>  x >= p
             */
            if( impltype == SCIP_BOUNDTYPE_UPPER )
            {
               if( implvar->vubs != NULL ) /* implvar may have been aggregated in the mean time */
               {
                  SCIPsetDebugMsg(set, "deleting variable bound: <%s> == %u  ==>  <%s> <= %g\n",
                     SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar),
                     SCIPimplicsGetBounds(var->implics, varfixing)[i]);
                  SCIP_CALL( SCIPvboundsDel(&implvar->vubs, blkmem, var, varfixing) );
                  implvar->closestvblpcount = -1;
                  var->closestvblpcount = -1;
               }
            }
            else
            {
               if( implvar->vlbs != NULL ) /* implvar may have been aggregated in the mean time */
               {
                  SCIPsetDebugMsg(set, "deleting variable bound: <%s> == %u  ==>  <%s> >= %g\n",
                     SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar),
                     SCIPimplicsGetBounds(var->implics, varfixing)[i]);
                  SCIP_CALL( SCIPvboundsDel(&implvar->vlbs, blkmem, var, !varfixing) );
                  implvar->closestvblpcount = -1;
                  var->closestvblpcount = -1;
               }
            }
         }
         varfixing = !varfixing;
      }
      while( varfixing == TRUE );

      if( removefromvar )
      {
         /* free the implications data structures */
         SCIPimplicsFree(&var->implics, blkmem);
      }
   }

   /* remove the (redundant) variable lower bounds */
   if( var->vlbs != NULL )
   {
      SCIP_VAR** vars;
      SCIP_Real* coefs;
      SCIP_Real* constants;
      int nvbds;
      int newnvbds;
      int i;

      nvbds = SCIPvboundsGetNVbds(var->vlbs);
      vars = SCIPvboundsGetVars(var->vlbs);
      coefs = SCIPvboundsGetCoefs(var->vlbs);
      constants = SCIPvboundsGetConstants(var->vlbs);

      /* remove for all variable bounds x >= b*z+d the following implication from z's implications 
       *   z == ub  ==>  x >= b*ub + d           , if b > 0
       *   z == lb  ==>  x >= b*lb + d           , if b < 0
       */
      newnvbds = 0;
      for( i = 0; i < nvbds; i++ )
      {
         SCIP_VAR* implvar;
         SCIP_Real coef;

         assert(newnvbds <= i);

         implvar = vars[i];
         assert(implvar != NULL);

         coef = coefs[i];
         assert(!SCIPsetIsZero(set, coef));

         /* check, if we want to remove the variable bound */
         if( onlyredundant )
         {
            SCIP_Real vbound;

            vbound = MAX(coef * SCIPvarGetUbGlobal(implvar), coef * SCIPvarGetLbGlobal(implvar)) + constants[i];  /*lint !e666*/
            if( SCIPsetIsFeasGT(set, vbound, lb) )
            {
               /* the variable bound is not redundant: keep it */
               if( removefromvar )
               {
                  if( newnvbds < i )
                  {
                     vars[newnvbds] = implvar;
                     coefs[newnvbds] = coef;
                     constants[newnvbds] = constants[i];
                  }
                  newnvbds++;
               }
               continue;
            }
         }

         /* remove the corresponding implication */
         if( implvar->implics != NULL ) /* variable may have been aggregated in the mean time */
         {
            SCIPsetDebugMsg(set, "deleting implication: <%s> == %d  ==>  <%s> >= %g\n",
               SCIPvarGetName(implvar), (coef > 0.0), SCIPvarGetName(var), MAX(coef, 0.0) + constants[i]);
            SCIP_CALL( SCIPimplicsDel(&implvar->implics, blkmem, set, (coef > 0.0), var, SCIP_BOUNDTYPE_LOWER) );
         }
         if( coef > 0.0 && implvar->vubs != NULL ) /* implvar may have been aggregated in the mean time */
         {
            SCIPsetDebugMsg(set, "deleting variable upper bound from <%s> involving variable %s\n",
               SCIPvarGetName(implvar), SCIPvarGetName(var));
            SCIP_CALL( SCIPvboundsDel(&implvar->vubs, blkmem, var, FALSE) );
            implvar->closestvblpcount = -1;
            var->closestvblpcount = -1;
         }
         else if( coef < 0.0 && implvar->vlbs != NULL ) /* implvar may have been aggregated in the mean time */
         {
            SCIPsetDebugMsg(set, "deleting variable lower bound from <%s> involving variable %s\n",
               SCIPvarGetName(implvar), SCIPvarGetName(var));
            SCIP_CALL( SCIPvboundsDel(&implvar->vlbs, blkmem, var, TRUE) );
            implvar->closestvblpcount = -1;
            var->closestvblpcount = -1;
         }
      }

      if( removefromvar )
      {
         /* update the number of variable bounds */
         SCIPvboundsShrink(&var->vlbs, blkmem, newnvbds);
         var->closestvblpcount = -1;
      }
   }

   /**@todo in general, variable bounds like x >= b*z + d corresponding to an implication like z = ub ==> x >= b*ub + d 
    *       might be missing because we only add variable bounds with reasonably small value of b. thus, we currently 
    *       cannot remove such variables x from z's implications.
    */

   /* remove the (redundant) variable upper bounds */
   if( var->vubs != NULL )
   {
      SCIP_VAR** vars;
      SCIP_Real* coefs;
      SCIP_Real* constants;
      int nvbds;
      int newnvbds;
      int i;

      nvbds = SCIPvboundsGetNVbds(var->vubs);
      vars = SCIPvboundsGetVars(var->vubs);
      coefs = SCIPvboundsGetCoefs(var->vubs);
      constants = SCIPvboundsGetConstants(var->vubs);

      /* remove for all variable bounds x <= b*z+d the following implication from z's implications 
       *   z == lb  ==>  x <= b*lb + d           , if b > 0
       *   z == ub  ==>  x <= b*ub + d           , if b < 0
       */
      newnvbds = 0;
      for( i = 0; i < nvbds; i++ )
      {
         SCIP_VAR* implvar;
         SCIP_Real coef;

         assert(newnvbds <= i);

         implvar = vars[i];
         assert(implvar != NULL);

         coef = coefs[i];
         assert(!SCIPsetIsZero(set, coef));

         /* check, if we want to remove the variable bound */
         if( onlyredundant )
         {
            SCIP_Real vbound;

            vbound = MIN(coef * SCIPvarGetUbGlobal(implvar), coef * SCIPvarGetLbGlobal(implvar)) + constants[i];  /*lint !e666*/
            if( SCIPsetIsFeasLT(set, vbound, ub) )
            {
               /* the variable bound is not redundant: keep it */
               if( removefromvar )
               {
                  if( newnvbds < i )
                  {
                     vars[newnvbds] = implvar;
                     coefs[newnvbds] = coefs[i];
                     constants[newnvbds] = constants[i];
                  }
                  newnvbds++;
               }
               continue;
            }
         }

         /* remove the corresponding implication */
         if( implvar->implics != NULL ) /* variable may have been aggregated in the mean time */
         {
            SCIPsetDebugMsg(set, "deleting implication: <%s> == %d  ==>  <%s> <= %g\n",
               SCIPvarGetName(implvar), (coef < 0.0), SCIPvarGetName(var), MIN(coef, 0.0) + constants[i]);
            SCIP_CALL( SCIPimplicsDel(&implvar->implics, blkmem, set, (coef < 0.0), var, SCIP_BOUNDTYPE_UPPER) );
         }
         if( coef < 0.0 && implvar->vubs != NULL ) /* implvar may have been aggregated in the mean time */
         {
            SCIPsetDebugMsg(set, "deleting variable upper bound from <%s> involving variable %s\n",
               SCIPvarGetName(implvar), SCIPvarGetName(var));
            SCIP_CALL( SCIPvboundsDel(&implvar->vubs, blkmem, var, TRUE) );
            implvar->closestvblpcount = -1;
            var->closestvblpcount = -1;
         }
         else if( coef > 0.0 && implvar->vlbs != NULL ) /* implvar may have been aggregated in the mean time */
         {
            SCIPsetDebugMsg(set, "deleting variable lower bound from <%s> involving variable %s\n",
               SCIPvarGetName(implvar), SCIPvarGetName(var));
            SCIP_CALL( SCIPvboundsDel(&implvar->vlbs, blkmem, var, FALSE) );
            implvar->closestvblpcount = -1;
            var->closestvblpcount = -1;
         }
      }

      if( removefromvar )
      {
         /* update the number of variable bounds */
         SCIPvboundsShrink(&var->vubs, blkmem, newnvbds);
         var->closestvblpcount = -1;
      }
   }

   /* remove the variable from all cliques */
   if( SCIPvarIsBinary(var) )
   {
      SCIPcliquelistRemoveFromCliques(var->cliquelist, cliquetable, var, irrelevantvar);
      SCIPcliquelistFree(&var->cliquelist, blkmem);
   }

   /**@todo variable bounds like x <= b*z + d with z general integer are not removed from x's vbd arrays, because
    *       z has no link (like in the binary case) to x
    */

   return SCIP_OKAY;
}

/** sets the variable name */
static
SCIP_RETCODE varSetName(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< problem statistics, or NULL */
   const char*           name                /**< name of variable, or NULL for automatic name creation */
   )
{
   assert(blkmem != NULL);
   assert(var != NULL);

   if( name == NULL )
   {
      char s[SCIP_MAXSTRLEN];

      assert(stat != NULL);

      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "_var%d_", stat->nvaridx);
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &var->name, s, strlen(s)+1) );
   }
   else
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &var->name, name, strlen(name)+1) );
   }

   return SCIP_OKAY;
}


/** creates variable; if variable is of integral type, fractional bounds are automatically rounded; an integer variable
 *  with bounds zero and one is automatically converted into a binary variable
 */
static
SCIP_RETCODE varCreate(
   SCIP_VAR**            var,                /**< pointer to variable data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable, or NULL */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data, or NULL */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable, or NULL */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   )
{
   assert(var != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);

   /* adjust bounds of variable */
   lb = adjustedLb(set, vartype, lb);
   ub = adjustedUb(set, vartype, ub);

   /* convert [0,1]-integers into binary variables and check that binary variables have correct bounds */
   if( (SCIPsetIsEQ(set, lb, 0.0) || SCIPsetIsEQ(set, lb, 1.0))
      && (SCIPsetIsEQ(set, ub, 0.0) || SCIPsetIsEQ(set, ub, 1.0)) )
   {
      if( vartype == SCIP_VARTYPE_INTEGER )
         vartype = SCIP_VARTYPE_BINARY;
   }
   else
   {
      if( vartype == SCIP_VARTYPE_BINARY )
      {
         SCIPerrorMessage("invalid bounds [%.2g,%.2g] for binary variable <%s>\n", lb, ub, name);
         return SCIP_INVALIDDATA;
      }
   }

   assert(vartype != SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, lb, 0.0) || SCIPsetIsEQ(set, lb, 1.0));
   assert(vartype != SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, ub, 0.0) || SCIPsetIsEQ(set, ub, 1.0));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, var) );

   /* set variable's name */
   SCIP_CALL( varSetName(*var, blkmem, stat, name) );

#ifndef NDEBUG
   (*var)->scip = set->scip;
#endif
   (*var)->obj = obj;
   (*var)->unchangedobj = obj;
   (*var)->branchfactor = 1.0;
   (*var)->rootsol = 0.0;
   (*var)->bestrootsol = 0.0;
   (*var)->bestrootredcost = 0.0;
   (*var)->bestrootlpobjval = SCIP_INVALID;
   (*var)->relaxsol = 0.0;
   (*var)->nlpsol = 0.0;
   (*var)->primsolavg = 0.5 * (lb + ub);
   (*var)->conflictlb = SCIP_REAL_MIN;
   (*var)->conflictub = SCIP_REAL_MAX;
   (*var)->conflictrelaxedlb = (*var)->conflictlb;
   (*var)->conflictrelaxedub = (*var)->conflictub;
   (*var)->lazylb = -SCIPsetInfinity(set);
   (*var)->lazyub = SCIPsetInfinity(set);
   (*var)->glbdom.holelist = NULL;
   (*var)->glbdom.lb = lb;
   (*var)->glbdom.ub = ub;
   (*var)->locdom.holelist = NULL;
   (*var)->locdom.lb = lb;
   (*var)->locdom.ub = ub;
   (*var)->varcopy = varcopy;
   (*var)->vardelorig = vardelorig;
   (*var)->vartrans = vartrans;
   (*var)->vardeltrans = vardeltrans;
   (*var)->vardata = vardata;
   (*var)->parentvars = NULL;
   (*var)->negatedvar = NULL;
   (*var)->vlbs = NULL;
   (*var)->vubs = NULL;
   (*var)->implics = NULL;
   (*var)->cliquelist = NULL;
   (*var)->eventfilter = NULL;
   (*var)->lbchginfos = NULL;
   (*var)->ubchginfos = NULL;
   (*var)->index = stat->nvaridx;
   (*var)->probindex = -1;
   (*var)->pseudocandindex = -1;
   (*var)->eventqueueindexobj = -1;
   (*var)->eventqueueindexlb = -1;
   (*var)->eventqueueindexub = -1;
   (*var)->parentvarssize = 0;
   (*var)->nparentvars = 0;
   (*var)->nuses = 0;
   (*var)->nlocksdown = 0;
   (*var)->nlocksup = 0;
   (*var)->branchpriority = 0;
   (*var)->branchdirection = SCIP_BRANCHDIR_AUTO; /*lint !e641*/
   (*var)->lbchginfossize = 0;
   (*var)->nlbchginfos = 0;
   (*var)->ubchginfossize = 0;
   (*var)->nubchginfos = 0;
   (*var)->conflictlbcount = 0;
   (*var)->conflictubcount = 0;
   (*var)->closestvlbidx = -1;
   (*var)->closestvubidx = -1;
   (*var)->closestvblpcount = -1;
   (*var)->initial = initial;
   (*var)->removable = removable;
   (*var)->deleted = FALSE;
   (*var)->donotmultaggr = FALSE;
   (*var)->vartype = vartype; /*lint !e641*/
   (*var)->pseudocostflag = FALSE;
   (*var)->eventqueueimpl = FALSE;
   (*var)->deletable = FALSE;
   (*var)->delglobalstructs = FALSE;

   stat->nvaridx++;

   /* create branching and inference history entries */
   SCIP_CALL( SCIPhistoryCreate(&(*var)->history, blkmem) );
   SCIP_CALL( SCIPhistoryCreate(&(*var)->historycrun, blkmem) );

   /* the value based history is only created on demand */
   (*var)->valuehistory = NULL;

   return SCIP_OKAY;
}

/** creates and captures an original problem variable; an integer variable with bounds
 *  zero and one is automatically converted into a binary variable
 */
SCIP_RETCODE SCIPvarCreateOriginal(
   SCIP_VAR**            var,                /**< pointer to variable data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable, or NULL */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data, or NULL */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable, or NULL */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   )
{
   assert(var != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);

   /* create variable */
   SCIP_CALL( varCreate(var, blkmem, set, stat, name, lb, ub, obj, vartype, initial, removable,
         varcopy, vardelorig, vartrans, vardeltrans, vardata) );

   /* set variable status and data */
   (*var)->varstatus = SCIP_VARSTATUS_ORIGINAL; /*lint !e641*/
   (*var)->data.original.origdom.holelist = NULL;
   (*var)->data.original.origdom.lb = lb;
   (*var)->data.original.origdom.ub = ub;
   (*var)->data.original.transvar = NULL;

   /* capture variable */
   SCIPvarCapture(*var);

   return SCIP_OKAY;
}

/** creates and captures a loose variable belonging to the transformed problem; an integer variable with bounds
 *  zero and one is automatically converted into a binary variable
 */
SCIP_RETCODE SCIPvarCreateTransformed(
   SCIP_VAR**            var,                /**< pointer to variable data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable, or NULL */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data, or NULL */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable, or NULL */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   )
{
   assert(var != NULL);
   assert(blkmem != NULL);

   /* create variable */
   SCIP_CALL( varCreate(var, blkmem, set, stat, name, lb, ub, obj, vartype, initial, removable,
         varcopy, vardelorig, vartrans, vardeltrans, vardata) );

   /* create event filter for transformed variable */
   SCIP_CALL( SCIPeventfilterCreate(&(*var)->eventfilter, blkmem) );

   /* set variable status and data */
   (*var)->varstatus = SCIP_VARSTATUS_LOOSE; /*lint !e641*/

   /* capture variable */
   SCIPvarCapture(*var);

   return SCIP_OKAY;   
}

/** copies and captures a variable from source to target SCIP; an integer variable with bounds zero and one is
 *  automatically converted into a binary variable; in case the variable data cannot be copied the variable is not
 *  copied at all
 */
SCIP_RETCODE SCIPvarCopy(
   SCIP_VAR**            var,                /**< pointer to store the target variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP_VAR*             sourcevar,          /**< source variable */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints */
   SCIP_Bool             global              /**< should global or local bounds be used? */
   )
{
   SCIP_VARDATA* targetdata;
   SCIP_RESULT result;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(sourcescip != NULL);
   assert(sourcevar != NULL);
   assert(var != NULL);
   assert(set->stage == SCIP_STAGE_PROBLEM);
   assert(varmap != NULL);
   assert(consmap != NULL);

   /** @todo copy hole lists */
   assert(global || SCIPvarGetHolelistLocal(sourcevar) == NULL);
   assert(!global || SCIPvarGetHolelistGlobal(sourcevar) == NULL);

   result = SCIP_DIDNOTRUN;
   targetdata = NULL;

   if( SCIPvarGetStatus(sourcevar) == SCIP_VARSTATUS_ORIGINAL )
   {
      lb = SCIPvarGetLbOriginal(sourcevar);
      ub = SCIPvarGetUbOriginal(sourcevar);
   }
   else
   {
      lb = global ? SCIPvarGetLbGlobal(sourcevar) : SCIPvarGetLbLocal(sourcevar);
      ub = global ? SCIPvarGetUbGlobal(sourcevar) : SCIPvarGetUbLocal(sourcevar);
   }

   /* creates and captures the variable in the target SCIP and initialize callback methods and variable data to NULL */
   SCIP_CALL( SCIPvarCreateOriginal(var, blkmem, set, stat, SCIPvarGetName(sourcevar), 
         lb, ub, SCIPvarGetObj(sourcevar), SCIPvarGetType(sourcevar),
         SCIPvarIsInitial(sourcevar), SCIPvarIsRemovable(sourcevar), 
         NULL, NULL, NULL, NULL, NULL) );
   assert(*var != NULL);

   /* directly copy donotmultaggr flag */
   (*var)->donotmultaggr = sourcevar->donotmultaggr;

   /* insert variable into mapping between source SCIP and the target SCIP */
   assert(!SCIPhashmapExists(varmap, sourcevar));
   SCIP_CALL( SCIPhashmapInsert(varmap, sourcevar, *var) );

   /* in case there exists variable data and the variable data copy callback, try to copy variable data */
   if( sourcevar->vardata != NULL && sourcevar->varcopy != NULL )
   {
      SCIP_CALL( sourcevar->varcopy(set->scip, sourcescip, sourcevar, sourcevar->vardata, 
            varmap, consmap, (*var), &targetdata, &result) );

      /* evaluate result */
      if( result != SCIP_DIDNOTRUN && result != SCIP_SUCCESS )
      {
         SCIPerrorMessage("variable data copying method returned invalid result <%d>\n", result);
         return SCIP_INVALIDRESULT;
      }

      assert(targetdata == NULL || result == SCIP_SUCCESS);

      /* if copying was successful, add the created variable data to the variable as well as all callback methods */
      if( result == SCIP_SUCCESS )
      {
         (*var)->varcopy = sourcevar->varcopy;
         (*var)->vardelorig = sourcevar->vardelorig;
         (*var)->vartrans = sourcevar->vartrans;
         (*var)->vardeltrans = sourcevar->vardeltrans;
         (*var)->vardata = targetdata;
      }
   }

   /* we initialize histories of the variables by copying the source variable-information */
   if( set->history_allowtransfer )
   {
      SCIPvarMergeHistories((*var), sourcevar, stat);
   }

   /* in case the copying was successfully, add the created variable data to the variable as well as all callback
    * methods 
    */
   if( result == SCIP_SUCCESS )
   {
      (*var)->varcopy = sourcevar->varcopy;
      (*var)->vardelorig = sourcevar->vardelorig;
      (*var)->vartrans = sourcevar->vartrans;
      (*var)->vardeltrans = sourcevar->vardeltrans;
      (*var)->vardata = targetdata;
   }

   SCIPsetDebugMsg(set, "created copy <%s> of variable <%s>\n", SCIPvarGetName(*var), SCIPvarGetName(sourcevar));

   return SCIP_OKAY;
}

/** parse given string for a SCIP_Real bound */
static
SCIP_RETCODE parseValue(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           str,                /**< string to parse */
   SCIP_Real*            value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed */
   )
{
   /* first check for infinity value */
   if( strncmp(str, "+inf", 4) == 0 )
   {
      *value = SCIPsetInfinity(set);
      (*endptr) = (char*)str + 4;
   }
   else if( strncmp(str, "-inf", 4) == 0 )
   {
      *value = -SCIPsetInfinity(set);
      (*endptr) = (char*)str + 4;
   }
   else
   {
      if( !SCIPstrToRealValue(str, value, endptr) )
         return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** parse the characters as bounds */
static
SCIP_RETCODE parseBounds(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           str,                /**< string to parse */
   char*                 type,               /**< bound type (global, local, or lazy) */
   SCIP_Real*            lb,                 /**< pointer to store the lower bound */
   SCIP_Real*            ub,                 /**< pointer to store the upper bound */
   char**                endptr              /**< pointer to store the final string position if successfully parsed (or NULL if an error occured) */
   )
{
   char token[SCIP_MAXSTRLEN];
   char* tmpend;

   SCIPsetDebugMsg(set, "parsing bounds: '%s'\n", str);

   /* get bound type */
   SCIPstrCopySection(str, ' ', ' ', type, SCIP_MAXSTRLEN, endptr);
   if ( strncmp(type, "original", 8) != 0 && strncmp(type, "global", 6) != 0 && strncmp(type, "local", 5) != 0 && strncmp(type, "lazy", 4) != 0 )
   {
      SCIPsetDebugMsg(set, "unkown bound type <%s>\n", type);
      *endptr = NULL;
      return SCIP_OKAY;
   }

   SCIPsetDebugMsg(set, "parsed bound type <%s>\n", type);

   /* get lower bound */
   SCIPstrCopySection(str, '[', ',', token, SCIP_MAXSTRLEN, endptr);
   str = *endptr;
   SCIP_CALL( parseValue(set, token, lb, &tmpend) );

   /* get upper bound */
   SCIP_CALL( parseValue(set, str, ub, endptr) );

   SCIPsetDebugMsg(set, "parsed bounds: [%g,%g]\n", *lb, *ub);

   /* skip end of bounds */
   while ( **endptr != '\0' && (**endptr == ']' || **endptr == ',') )
      ++(*endptr);

   return SCIP_OKAY;
}

/** parses a given string for a variable informations */
static
SCIP_RETCODE varParse(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           str,                /**< string to parse */
   char*                 name,               /**< pointer to store the variable name */
   SCIP_Real*            lb,                 /**< pointer to store the lower bound */
   SCIP_Real*            ub,                 /**< pointer to store the upper bound */
   SCIP_Real*            obj,                /**< pointer to store the objective coefficient */
   SCIP_VARTYPE*         vartype,            /**< pointer to store the variable type */
   SCIP_Real*            lazylb,             /**< pointer to store if the lower bound is lazy */
   SCIP_Real*            lazyub,             /**< pointer to store if the upper bound is lazy */
   SCIP_Bool             local,              /**< should the local bound be applied */
   char**                endptr,             /**< pointer to store the final string position if successfully */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   )
{
   SCIP_Real parsedlb;
   SCIP_Real parsedub;
   char token[SCIP_MAXSTRLEN];
   char* strptr;
   int i;

   assert(lb != NULL);
   assert(ub != NULL);
   assert(obj != NULL);
   assert(vartype != NULL);
   assert(lazylb != NULL);
   assert(lazyub != NULL);
   assert(success != NULL);

   (*success) = TRUE;

   /* copy variable type */
   SCIPstrCopySection(str, '[', ']', token, SCIP_MAXSTRLEN, endptr);
   assert(str != *endptr);
   SCIPsetDebugMsg(set, "parsed variable type <%s>\n", token);

   /* get variable type */
   if( strncmp(token, "binary", 3) == 0 )
      (*vartype) = SCIP_VARTYPE_BINARY;
   else if( strncmp(token, "integer", 3) == 0 )
      (*vartype) = SCIP_VARTYPE_INTEGER;
   else if( strncmp(token, "implicit", 3) == 0 )
      (*vartype) = SCIP_VARTYPE_IMPLINT;
   else if( strncmp(token, "continuous", 3) == 0 )
      (*vartype) = SCIP_VARTYPE_CONTINUOUS;
   else
   {
      SCIPmessagePrintWarning(messagehdlr, "unknown variable type\n");
      (*success) = FALSE;
      return SCIP_OKAY;
   }

   /* move string pointer behind variable type */
   str = *endptr;

   /* get variable name */
   SCIPstrCopySection(str, '<', '>', name, SCIP_MAXSTRLEN, endptr);
   assert(endptr != NULL);
   SCIPsetDebugMsg(set, "parsed variable name <%s>\n", name);

   /* move string pointer behind variable name */
   str = *endptr;

   /* cut out objective coefficient */
   SCIPstrCopySection(str, '=', ',', token, SCIP_MAXSTRLEN, endptr);

   /* move string pointer behind objective coefficient */
   str = *endptr;

   /* get objective coefficient */
   if( !SCIPstrToRealValue(token, obj, endptr) )
   {
      *endptr = NULL;
      return SCIP_READERROR;
   }

   SCIPsetDebugMsg(set, "parsed objective coefficient <%g>\n", *obj);

   /* parse global/original bounds */
   SCIP_CALL( parseBounds(set, str, token, lb, ub, endptr) );
   assert(strncmp(token, "global", 6) == 0 || strncmp(token, "original", 8) == 0);

   /* initialize the lazy bound */
   *lazylb = -SCIPsetInfinity(set);
   *lazyub =  SCIPsetInfinity(set);

   /* store pointer */
   strptr = *endptr;

   /* possibly parse optional local and lazy bounds */
   for( i = 0; i < 2 && *endptr != NULL && **endptr != '\0'; ++i )
   {
      /* start after previous bounds */
      strptr = *endptr;

      /* parse global bounds */
      SCIP_CALL( parseBounds(set, strptr, token, &parsedlb, &parsedub, endptr) );

      /* stop if parsing of bounds failed */
      if( *endptr == NULL )
         break;

      if( strncmp(token, "local", 5) == 0 && local )
      {
         *lb = parsedlb;
         *ub = parsedub;
      }
      else if( strncmp(token, "lazy", 4) == 0 )
      {
         *lazylb = parsedlb;
         *lazyub = parsedub;
      }
   }

   /* restore pointer */
   if ( *endptr == NULL )
      *endptr = strptr;

   /* check bounds for binary variables */
   if ( (*vartype) == SCIP_VARTYPE_BINARY )
   {
      if ( SCIPsetIsLT(set, *lb, 0.0) || SCIPsetIsGT(set, *ub, 1.0) )
      {
         SCIPerrorMessage("Parsed invalid bounds for binary variable <%s>: [%f, %f].\n", name, *lb, *ub);
         return SCIP_READERROR;
      }
      if ( !SCIPsetIsInfinity(set, -(*lazylb)) && !SCIPsetIsInfinity(set, *lazyub) && 
           ( SCIPsetIsLT(set, *lazylb, 0.0) || SCIPsetIsGT(set, *lazyub, 1.0) ) )
      {
         SCIPerrorMessage("Parsed invalid lazy bounds for binary variable <%s>: [%f, %f].\n", name, *lazylb, *lazyub);
         return SCIP_READERROR;
      }
   }

   return SCIP_OKAY;
}

/** parses variable information (in cip format) out of a string; if the parsing process was successful an original
 *  variable is created and captured; if variable is of integral type, fractional bounds are automatically rounded; an
 *  integer variable with bounds zero and one is automatically converted into a binary variable
 */
SCIP_RETCODE SCIPvarParseOriginal(
   SCIP_VAR**            var,                /**< pointer to variable data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           str,                /**< string to parse */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   SCIP_VARDATA*         vardata,            /**< user data for this specific variable */
   char**                endptr,             /**< pointer to store the final string position if successfully */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   SCIP_VARTYPE vartype;
   SCIP_Real lazylb;
   SCIP_Real lazyub;

   assert(var != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(endptr != NULL);
   assert(success != NULL);

   /* parse string in cip format for variable information */
   SCIP_CALL( varParse(set, messagehdlr, str, name, &lb, &ub, &obj, &vartype, &lazylb, &lazyub, FALSE, endptr, success) );

   if( *success )
   {
      /* create variable */
      SCIP_CALL( varCreate(var, blkmem, set, stat, name, lb, ub, obj, vartype, initial, removable,
            varcopy, vardelorig, vartrans, vardeltrans, vardata) );

      /* set variable status and data */
      (*var)->varstatus = SCIP_VARSTATUS_ORIGINAL; /*lint !e641*/
      (*var)->data.original.origdom.holelist = NULL;
      (*var)->data.original.origdom.lb = lb;
      (*var)->data.original.origdom.ub = ub;
      (*var)->data.original.transvar = NULL;

      /* set lazy status of variable bounds */
      (*var)->lazylb = lazylb;
      (*var)->lazyub = lazyub;

      /* capture variable */
      SCIPvarCapture(*var);
   }

   return SCIP_OKAY;
}

/** parses variable information (in cip format) out of a string; if the parsing process was successful a loose variable
 *  belonging to the transformed problem is created and captured; if variable is of integral type, fractional bounds are
 *  automatically rounded; an integer variable with bounds zero and one is automatically converted into a binary
 *  variable
 */
SCIP_RETCODE SCIPvarParseTransformed(
   SCIP_VAR**            var,                /**< pointer to variable data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           str,                /**< string to parse */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   SCIP_VARDATA*         vardata,            /**< user data for this specific variable */
   char**                endptr,             /**< pointer to store the final string position if successfully */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   SCIP_VARTYPE vartype;
   SCIP_Real lazylb;
   SCIP_Real lazyub;

   assert(var != NULL);
   assert(blkmem != NULL);
   assert(endptr != NULL);
   assert(success != NULL);

   /* parse string in cip format for variable information */
   SCIP_CALL( varParse(set, messagehdlr, str, name, &lb, &ub, &obj, &vartype, &lazylb, &lazyub, TRUE, endptr, success) );

   if( *success )
   {
      /* create variable */
      SCIP_CALL( varCreate(var, blkmem, set, stat, name, lb, ub, obj, vartype, initial, removable,
            varcopy, vardelorig, vartrans, vardeltrans, vardata) );

      /* create event filter for transformed variable */
      SCIP_CALL( SCIPeventfilterCreate(&(*var)->eventfilter, blkmem) );

      /* set variable status and data */
      (*var)->varstatus = SCIP_VARSTATUS_LOOSE; /*lint !e641*/

      /* set lazy status of variable bounds */
      (*var)->lazylb = lazylb;
      (*var)->lazyub = lazyub;

      /* capture variable */
      SCIPvarCapture(*var);
   }

   return SCIP_OKAY;   
}

/** ensures, that parentvars array of var can store at least num entries */
static
SCIP_RETCODE varEnsureParentvarsSize(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(var->nparentvars <= var->parentvarssize);

   if( num > var->parentvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &var->parentvars, var->parentvarssize, newsize) );
      var->parentvarssize = newsize;
   }
   assert(num <= var->parentvarssize);

   return SCIP_OKAY;
}

/** adds variable to parent list of a variable and captures parent variable */
static
SCIP_RETCODE varAddParent(
   SCIP_VAR*             var,                /**< variable to add parent to */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             parentvar           /**< parent variable to add */
   )
{
   assert(var != NULL);
   assert(parentvar != NULL);

   /* the direct original counterpart must be stored as first parent */
   assert(var->nparentvars == 0 || SCIPvarGetStatus(parentvar) != SCIP_VARSTATUS_ORIGINAL);

   SCIPsetDebugMsg(set, "adding parent <%s>[%p] to variable <%s>[%p] in slot %d\n",
      parentvar->name, (void*)parentvar, var->name, (void*)var, var->nparentvars);

   SCIP_CALL( varEnsureParentvarsSize(var, blkmem, set, var->nparentvars+1) );

   var->parentvars[var->nparentvars] = parentvar;
   var->nparentvars++;

   SCIPvarCapture(parentvar);

   return SCIP_OKAY;
}

/** deletes and releases all variables from the parent list of a variable, frees the memory of parents array */
static
SCIP_RETCODE varFreeParents(
   SCIP_VAR**            var,                /**< pointer to variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue (or NULL, if it's an original variable) */
   SCIP_LP*              lp                  /**< current LP data (or NULL, if it's an original variable) */
   )
{
   SCIP_VAR* parentvar;
   int i;

   SCIPsetDebugMsg(set, "free parents of <%s>\n", (*var)->name);

   /* release the parent variables and remove the link from the parent variable to the child */
   for( i = 0; i < (*var)->nparentvars; ++i )
   {
      assert((*var)->parentvars != NULL);
      parentvar = (*var)->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         assert(parentvar->data.original.transvar == *var);
         assert(&parentvar->data.original.transvar != var);
         parentvar->data.original.transvar = NULL;
         break;

      case SCIP_VARSTATUS_AGGREGATED:
         assert(parentvar->data.aggregate.var == *var);
         assert(&parentvar->data.aggregate.var != var);
         parentvar->data.aggregate.var = NULL;
         break;

#if 0
      /* The following code is unclear: should the current variable be removed from its parents? */
      case SCIP_VARSTATUS_MULTAGGR:
         assert(parentvar->data.multaggr.vars != NULL);
         for( v = 0; v < parentvar->data.multaggr.nvars && parentvar->data.multaggr.vars[v] != *var; ++v )
         {}
         assert(v < parentvar->data.multaggr.nvars && parentvar->data.multaggr.vars[v] == *var);
         if( v < parentvar->data.multaggr.nvars-1 )
         {
            parentvar->data.multaggr.vars[v] = parentvar->data.multaggr.vars[parentvar->data.multaggr.nvars-1];
            parentvar->data.multaggr.scalars[v] = parentvar->data.multaggr.scalars[parentvar->data.multaggr.nvars-1];
         }
         parentvar->data.multaggr.nvars--;
         break;
#endif

      case SCIP_VARSTATUS_NEGATED:
         assert(parentvar->negatedvar == *var);
         assert((*var)->negatedvar == parentvar);
         parentvar->negatedvar = NULL;
         (*var)->negatedvar = NULL;
         break;

      default:
         SCIPerrorMessage("parent variable is neither ORIGINAL, AGGREGATED nor NEGATED\n");
         return SCIP_INVALIDDATA;
      }  /*lint !e788*/

      SCIP_CALL( SCIPvarRelease(&(*var)->parentvars[i], blkmem, set, eventqueue, lp) );
   }

   /* free parentvars array */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*var)->parentvars, (*var)->parentvarssize);

   return SCIP_OKAY;
}

/** frees a variable */
static
SCIP_RETCODE varFree(
   SCIP_VAR**            var,                /**< pointer to variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue (may be NULL, if it's not a column variable) */
   SCIP_LP*              lp                  /**< current LP data (may be NULL, if it's not a column variable) */
   )
{
   assert(var != NULL);
   assert(*var != NULL);
   assert(SCIPvarGetStatus(*var) != SCIP_VARSTATUS_COLUMN || &(*var)->data.col->var != var);
   assert((*var)->nuses == 0);
   assert((*var)->probindex == -1);

   SCIPsetDebugMsg(set, "free variable <%s> with status=%d\n", (*var)->name, SCIPvarGetStatus(*var));

   switch( SCIPvarGetStatus(*var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert((*var)->data.original.transvar == NULL); /* cannot free variable, if transformed variable is still existing */
      holelistFree(&(*var)->data.original.origdom.holelist, blkmem);
      assert((*var)->data.original.origdom.holelist == NULL);
      break;
   case SCIP_VARSTATUS_LOOSE:
      break;
   case SCIP_VARSTATUS_COLUMN:
      SCIP_CALL( SCIPcolFree(&(*var)->data.col, blkmem, set, eventqueue, lp) );  /* free corresponding LP column */
      break;
   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
      break;
   case SCIP_VARSTATUS_MULTAGGR:
      BMSfreeBlockMemoryArray(blkmem, &(*var)->data.multaggr.vars, (*var)->data.multaggr.varssize);
      BMSfreeBlockMemoryArray(blkmem, &(*var)->data.multaggr.scalars, (*var)->data.multaggr.varssize);
      break;
   case SCIP_VARSTATUS_NEGATED:
      break;
   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   /* release all parent variables and free the parentvars array */
   SCIP_CALL( varFreeParents(var, blkmem, set, eventqueue, lp) );

   /* free user data */
   if( SCIPvarGetStatus(*var) == SCIP_VARSTATUS_ORIGINAL )
   {
      if( (*var)->vardelorig != NULL )
      {
         SCIP_CALL( (*var)->vardelorig(set->scip, *var, &(*var)->vardata) );
      }
   }
   else
   {
      if( (*var)->vardeltrans != NULL )
      {
         SCIP_CALL( (*var)->vardeltrans(set->scip, *var, &(*var)->vardata) );
      }
   }

   /* free event filter */
   if( (*var)->eventfilter != NULL )
   {
      SCIP_CALL( SCIPeventfilterFree(&(*var)->eventfilter, blkmem, set) );
   }
   assert((*var)->eventfilter == NULL);

   /* free hole lists */
   holelistFree(&(*var)->glbdom.holelist, blkmem);
   holelistFree(&(*var)->locdom.holelist, blkmem);
   assert((*var)->glbdom.holelist == NULL);
   assert((*var)->locdom.holelist == NULL);

   /* free variable bounds data structures */
   SCIPvboundsFree(&(*var)->vlbs, blkmem);
   SCIPvboundsFree(&(*var)->vubs, blkmem);

   /* free implications data structures */
   SCIPimplicsFree(&(*var)->implics, blkmem);

   /* free clique list data structures */
   SCIPcliquelistFree(&(*var)->cliquelist, blkmem);

   /* free bound change information arrays */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*var)->lbchginfos, (*var)->lbchginfossize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*var)->ubchginfos, (*var)->ubchginfossize);

   /* free branching and inference history entries */
   SCIPhistoryFree(&(*var)->history, blkmem);
   SCIPhistoryFree(&(*var)->historycrun, blkmem);
   SCIPvaluehistoryFree(&(*var)->valuehistory, blkmem);

   /* free variable data structure */
   BMSfreeBlockMemoryArray(blkmem, &(*var)->name, strlen((*var)->name)+1);
   BMSfreeBlockMemory(blkmem, var);

   return SCIP_OKAY;
}

/** increases usage counter of variable */
void SCIPvarCapture(
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert(var != NULL);
   assert(var->nuses >= 0);

   SCIPdebugMessage("capture variable <%s> with nuses=%d\n", var->name, var->nuses);
   var->nuses++;
}

/** decreases usage counter of variable, and frees memory if necessary */
SCIP_RETCODE SCIPvarRelease(
   SCIP_VAR**            var,                /**< pointer to variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data (or NULL, if it's an original variable) */
   )
{
   assert(var != NULL);
   assert(*var != NULL);
   assert((*var)->nuses >= 1);
   assert(blkmem != NULL);
   assert((*var)->scip == set->scip);

   SCIPsetDebugMsg(set, "release variable <%s> with nuses=%d\n", (*var)->name, (*var)->nuses);
   (*var)->nuses--;
   if( (*var)->nuses == 0 )
   {
      SCIP_CALL( varFree(var, blkmem, set, eventqueue, lp) );
   }

   *var = NULL;

   return SCIP_OKAY;
}

/** change variable name */
SCIP_RETCODE SCIPvarChgName(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name                /**< name of variable */
   )
{
   assert(name != NULL);

   /* remove old variable name */
   BMSfreeBlockMemoryArray(blkmem, &var->name, strlen(var->name)+1);

   /* set new variable name */
   SCIP_CALL( varSetName(var, blkmem, NULL, name) );

   return SCIP_OKAY;
}

/** initializes variable data structure for solving */
void SCIPvarInitSolve(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   SCIPhistoryReset(var->historycrun);
   var->conflictlbcount = 0;
   var->conflictubcount = 0;
}

/** outputs the given bounds into the file stream */
static
void printBounds(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub,                 /**< upper bound */
   const char*           name                /**< bound type name */
   )
{
   assert(set != NULL);

   SCIPmessageFPrintInfo(messagehdlr, file, ", %s=", name);
   if( SCIPsetIsInfinity(set, lb) )
      SCIPmessageFPrintInfo(messagehdlr, file, "[+inf,");
   else if( SCIPsetIsInfinity(set, -lb) )
      SCIPmessageFPrintInfo(messagehdlr, file, "[-inf,");
   else
      SCIPmessageFPrintInfo(messagehdlr, file, "[%.15g,", lb);
   if( SCIPsetIsInfinity(set, ub) )
      SCIPmessageFPrintInfo(messagehdlr, file, "+inf]");
   else if( SCIPsetIsInfinity(set, -ub) )
      SCIPmessageFPrintInfo(messagehdlr, file, "-inf]");
   else
      SCIPmessageFPrintInfo(messagehdlr, file, "%.15g]", ub);
}

/** prints hole list to file stream */
static
void printHolelist(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_HOLELIST*        holelist,           /**< hole list pointer to hole of interest */
   const char*           name                /**< hole type name */
   )
{  /*lint --e{715}*/
   SCIP_Real left;
   SCIP_Real right;

   if( holelist == NULL )
      return;

   left = SCIPholelistGetLeft(holelist);
   right = SCIPholelistGetRight(holelist);

   /* display first hole */
   SCIPmessageFPrintInfo(messagehdlr, file, ", %s=(%g,%g)", name, left, right);
   holelist = SCIPholelistGetNext(holelist);

   while(holelist != NULL  )
   {
      left = SCIPholelistGetLeft(holelist);
      right = SCIPholelistGetRight(holelist);

      /* display hole */
      SCIPmessageFPrintInfo(messagehdlr, file, "(%g,%g)", left, right);

      /* get next hole */
      holelist = SCIPholelistGetNext(holelist);
   }
}

/** outputs variable information into file stream */
SCIP_RETCODE SCIPvarPrint(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_HOLELIST* holelist;
   SCIP_Real lb;
   SCIP_Real ub;
   int i;

   assert(var != NULL);
   assert(var->scip == set->scip);

   /* type of variable */
   switch( SCIPvarGetType(var) )
   {
   case SCIP_VARTYPE_BINARY:
      SCIPmessageFPrintInfo(messagehdlr, file, "  [binary]");
      break;
   case SCIP_VARTYPE_INTEGER:
      SCIPmessageFPrintInfo(messagehdlr, file, "  [integer]");
      break;
   case SCIP_VARTYPE_IMPLINT:
      SCIPmessageFPrintInfo(messagehdlr, file, "  [implicit]");
      break;
   case SCIP_VARTYPE_CONTINUOUS:
      SCIPmessageFPrintInfo(messagehdlr, file, "  [continuous]");
      break;
   default:
      SCIPerrorMessage("unknown variable type\n");
      SCIPABORT();
      return SCIP_ERROR; /*lint !e527*/
   }

   /* name */
   SCIPmessageFPrintInfo(messagehdlr, file, " <%s>:", var->name);

   /* objective value */
   SCIPmessageFPrintInfo(messagehdlr, file, " obj=%.15g", var->obj);

   /* bounds (global bounds for transformed variables, original bounds for original variables) */
   if( !SCIPvarIsTransformed(var) )
   {
      /* output original bound */
      lb = SCIPvarGetLbOriginal(var);
      ub = SCIPvarGetUbOriginal(var);
      printBounds(set, messagehdlr, file, lb, ub, "original bounds");

      /* output lazy bound */
      lb = SCIPvarGetLbLazy(var);
      ub = SCIPvarGetUbLazy(var);

      /* only display the lazy bounds if they are different from [-infinity,infinity] */
      if( !SCIPsetIsInfinity(set, -lb) || !SCIPsetIsInfinity(set, ub) )
         printBounds(set, messagehdlr, file, lb, ub, "lazy bounds");

      holelist = SCIPvarGetHolelistOriginal(var);
      printHolelist(set, messagehdlr, file, holelist, "original holes");
   }
   else
   {
      /* output global bound */
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);
      printBounds(set, messagehdlr, file, lb, ub, "global bounds");

      /* output local bound */
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      printBounds(set, messagehdlr, file, lb, ub, "local bounds");

      /* output lazy bound */
      lb = SCIPvarGetLbLazy(var);
      ub = SCIPvarGetUbLazy(var);

      /* only display the lazy bounds if they are different from [-infinity,infinity] */
      if( !SCIPsetIsInfinity(set, -lb) || !SCIPsetIsInfinity(set, ub) )
         printBounds(set, messagehdlr, file, lb, ub, "lazy bounds");

      /* global hole list */
      holelist = SCIPvarGetHolelistGlobal(var);
      printHolelist(set, messagehdlr, file, holelist, "global holes");

      /* local hole list */
      holelist = SCIPvarGetHolelistLocal(var);
      printHolelist(set, messagehdlr, file, holelist, "local holes");
   }

   /* fixings and aggregations */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      break;

   case SCIP_VARSTATUS_FIXED:
      SCIPmessageFPrintInfo(messagehdlr, file, ", fixed:");
      if( SCIPsetIsInfinity(set, var->glbdom.lb) )
         SCIPmessageFPrintInfo(messagehdlr, file, "+inf");
      else if( SCIPsetIsInfinity(set, -var->glbdom.lb) )
         SCIPmessageFPrintInfo(messagehdlr, file, "-inf");
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "%.15g", var->glbdom.lb);
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      SCIPmessageFPrintInfo(messagehdlr, file, ", aggregated:");
      if( !SCIPsetIsZero(set, var->data.aggregate.constant) )
         SCIPmessageFPrintInfo(messagehdlr, file, " %.15g", var->data.aggregate.constant);
      SCIPmessageFPrintInfo(messagehdlr, file, " %+.15g<%s>", var->data.aggregate.scalar, SCIPvarGetName(var->data.aggregate.var));
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPmessageFPrintInfo(messagehdlr, file, ", aggregated:");
      if( var->data.multaggr.nvars == 0 || !SCIPsetIsZero(set, var->data.multaggr.constant) )
         SCIPmessageFPrintInfo(messagehdlr, file, " %.15g", var->data.multaggr.constant);
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         SCIPmessageFPrintInfo(messagehdlr, file, " %+.15g<%s>", var->data.multaggr.scalars[i], SCIPvarGetName(var->data.multaggr.vars[i]));
      break;

   case SCIP_VARSTATUS_NEGATED:
      SCIPmessageFPrintInfo(messagehdlr, file, ", negated: %.15g - <%s>", var->data.negate.constant, SCIPvarGetName(var->negatedvar));
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_ERROR; /*lint !e527*/
   }

   SCIPmessageFPrintInfo(messagehdlr, file, "\n");

   return SCIP_OKAY;
}

/** issues a VARUNLOCKED event on the given variable */
static
SCIP_RETCODE varEventVarUnlocked(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_EVENT* event;

   assert(var != NULL);
   assert(var->nlocksdown <= 1 && var->nlocksup <= 1);
   assert(var->scip == set->scip);

   /* issue VARUNLOCKED event on variable */
   SCIP_CALL( SCIPeventCreateVarUnlocked(&event, blkmem, var) );
   SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, NULL, &event) );

   return SCIP_OKAY;
}

/** modifies lock numbers for rounding */
SCIP_RETCODE SCIPvarAddLocks(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   addnlocksdown,      /**< increase in number of rounding down locks */
   int                   addnlocksup         /**< increase in number of rounding up locks */
   )
{
   SCIP_VAR* lockvar;

   assert(var != NULL);
   assert(var->nlocksup >= 0);
   assert(var->nlocksdown >= 0);
   assert(var->scip == set->scip);

   if( addnlocksdown == 0 && addnlocksup == 0 )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "add rounding locks %d/%d to variable <%s> (locks=%d/%d)\n",
      addnlocksdown, addnlocksup, var->name, var->nlocksdown, var->nlocksup);

   lockvar = var;

   while( TRUE ) /*lint !e716 */
   {
      assert(lockvar != NULL);

      switch( SCIPvarGetStatus(lockvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( lockvar->data.original.transvar != NULL )
         {
            lockvar = lockvar->data.original.transvar;
            break;
         }
         else
         {
            lockvar->nlocksdown += addnlocksdown;
            lockvar->nlocksup += addnlocksup;

            assert(lockvar->nlocksdown >= 0);
            assert(lockvar->nlocksup >= 0);

            return SCIP_OKAY;
         }
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_FIXED:
         lockvar->nlocksdown += addnlocksdown;
         lockvar->nlocksup += addnlocksup;

         assert(lockvar->nlocksdown >= 0);
         assert(lockvar->nlocksup >= 0);

         if( lockvar->nlocksdown <= 1 && lockvar->nlocksup <= 1 )
         {
            SCIP_CALL( varEventVarUnlocked(lockvar, blkmem, set, eventqueue) );
         }

         return SCIP_OKAY;
      case SCIP_VARSTATUS_AGGREGATED:
         if( lockvar->data.aggregate.scalar < 0.0 )
         {
            int tmp = addnlocksup;

            addnlocksup = addnlocksdown;
            addnlocksdown = tmp;
         }

         lockvar = lockvar->data.aggregate.var;
         break;
      case SCIP_VARSTATUS_MULTAGGR:
      {
         int v;

         assert(!lockvar->donotmultaggr);

         for( v = lockvar->data.multaggr.nvars - 1; v >= 0; --v )
         {
            if( lockvar->data.multaggr.scalars[v] > 0.0 )
            {
               SCIP_CALL( SCIPvarAddLocks(lockvar->data.multaggr.vars[v], blkmem, set, eventqueue, addnlocksdown,
                     addnlocksup) );
            }
            else
            {
               SCIP_CALL( SCIPvarAddLocks(lockvar->data.multaggr.vars[v], blkmem, set, eventqueue, addnlocksup,
                     addnlocksdown) );
            }
         }
         return SCIP_OKAY;
      }
      case SCIP_VARSTATUS_NEGATED:
      {
         int tmp = addnlocksup;

         assert(lockvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(lockvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(lockvar->negatedvar->negatedvar == lockvar);

         addnlocksup = addnlocksdown;
         addnlocksdown = tmp;

         lockvar = lockvar->negatedvar;
         break;
      }
      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }
   }
}

/** gets number of locks for rounding down */
int SCIPvarGetNLocksDown(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   int nlocks;
   int i;

   assert(var != NULL);
   assert(var->nlocksdown >= 0);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar != NULL )
         return SCIPvarGetNLocksDown(var->data.original.transvar);
      else
         return var->nlocksdown;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      return var->nlocksdown;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNLocksDown(var->data.aggregate.var);
      else
         return SCIPvarGetNLocksUp(var->data.aggregate.var);

   case SCIP_VARSTATUS_MULTAGGR:
      assert(!var->donotmultaggr);
      nlocks = 0;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         if( var->data.multaggr.scalars[i] > 0.0 )
            nlocks += SCIPvarGetNLocksDown(var->data.multaggr.vars[i]);
         else
            nlocks += SCIPvarGetNLocksUp(var->data.multaggr.vars[i]);
      }
      return nlocks;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return SCIPvarGetNLocksUp(var->negatedvar);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return INT_MAX; /*lint !e527*/
   }
}

/** gets number of locks for rounding up */
int SCIPvarGetNLocksUp(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   int nlocks;
   int i;

   assert(var != NULL);
   assert(var->nlocksup >= 0);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar != NULL )
         return SCIPvarGetNLocksUp(var->data.original.transvar);
      else
         return var->nlocksup;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      return var->nlocksup;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNLocksUp(var->data.aggregate.var);
      else
         return SCIPvarGetNLocksDown(var->data.aggregate.var);

   case SCIP_VARSTATUS_MULTAGGR:
      assert(!var->donotmultaggr);
      nlocks = 0;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         if( var->data.multaggr.scalars[i] > 0.0 )
            nlocks += SCIPvarGetNLocksUp(var->data.multaggr.vars[i]);
         else
            nlocks += SCIPvarGetNLocksDown(var->data.multaggr.vars[i]);
      }
      return nlocks;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return SCIPvarGetNLocksDown(var->negatedvar);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return INT_MAX; /*lint !e527*/
   }
}

/** is it possible, to round variable down and stay feasible? */
SCIP_Bool SCIPvarMayRoundDown(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   return (SCIPvarGetNLocksDown(var) == 0);
}

/** is it possible, to round variable up and stay feasible? */
SCIP_Bool SCIPvarMayRoundUp(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   return (SCIPvarGetNLocksUp(var) == 0);
}

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 */
SCIP_RETCODE SCIPvarTransform(
   SCIP_VAR*             origvar,            /**< original problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_OBJSENSE         objsense,           /**< objective sense of original problem; transformed is always MINIMIZE */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable */
   )
{
   char name[SCIP_MAXSTRLEN];

   assert(origvar != NULL);
   assert(origvar->scip == set->scip);
   assert(SCIPvarGetStatus(origvar) == SCIP_VARSTATUS_ORIGINAL);
   assert(SCIPsetIsEQ(set, origvar->glbdom.lb, origvar->locdom.lb));
   assert(SCIPsetIsEQ(set, origvar->glbdom.ub, origvar->locdom.ub));
   assert(origvar->vlbs == NULL);
   assert(origvar->vubs == NULL);
   assert(transvar != NULL);

   /* check if variable is already transformed */
   if( origvar->data.original.transvar != NULL )
   {
      *transvar = origvar->data.original.transvar;
      SCIPvarCapture(*transvar);
   }
   else
   {
      /* create transformed variable */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "t_%s", origvar->name);
      SCIP_CALL( SCIPvarCreateTransformed(transvar, blkmem, set, stat, name,
            origvar->glbdom.lb, origvar->glbdom.ub, (SCIP_Real)objsense * origvar->obj,
            SCIPvarGetType(origvar), origvar->initial, origvar->removable,
            origvar->vardelorig, origvar->vartrans, origvar->vardeltrans, origvar->varcopy, NULL) );

      /* copy the branch factor and priority */
      (*transvar)->branchfactor = origvar->branchfactor;
      (*transvar)->branchpriority = origvar->branchpriority;
      (*transvar)->branchdirection = origvar->branchdirection; /*lint !e732*/

      /* duplicate hole lists */
      SCIP_CALL( holelistDuplicate(&(*transvar)->glbdom.holelist, blkmem, set, origvar->glbdom.holelist) );
      SCIP_CALL( holelistDuplicate(&(*transvar)->locdom.holelist, blkmem, set, origvar->locdom.holelist) );

      /* link original and transformed variable */
      origvar->data.original.transvar = *transvar;
      SCIP_CALL( varAddParent(*transvar, blkmem, set, origvar) );

      /* copy rounding locks */
      (*transvar)->nlocksdown = origvar->nlocksdown;
      (*transvar)->nlocksup = origvar->nlocksup;
      assert((*transvar)->nlocksdown >= 0);
      assert((*transvar)->nlocksup >= 0);

      /* copy doNotMultiaggr status */
      (*transvar)->donotmultaggr = origvar->donotmultaggr;

      /* copy lazy bounds */
      (*transvar)->lazylb = origvar->lazylb;
      (*transvar)->lazyub = origvar->lazyub;

      /* transfer eventual variable statistics; do not update global statistics, because this has been done
       * when original variable was created
       */
      SCIPhistoryUnite((*transvar)->history, origvar->history, FALSE);

      /* transform user data */
      if( origvar->vartrans != NULL )
      {
         SCIP_CALL( origvar->vartrans(set->scip, origvar, origvar->vardata, *transvar, &(*transvar)->vardata) );
      }
      else
         (*transvar)->vardata = origvar->vardata;
   }

   SCIPsetDebugMsg(set, "transformed variable: <%s>[%p] -> <%s>[%p]\n", origvar->name, (void*)origvar, (*transvar)->name, (void*)*transvar);

   return SCIP_OKAY;
}

/** gets corresponding transformed variable of an original or negated original variable */
SCIP_RETCODE SCIPvarGetTransformed(
   SCIP_VAR*             origvar,            /**< original problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable, or NULL if not existing yet */
   )
{
   assert(origvar != NULL);
   assert(SCIPvarGetStatus(origvar) == SCIP_VARSTATUS_ORIGINAL || SCIPvarGetStatus(origvar) == SCIP_VARSTATUS_NEGATED);
   assert(origvar->scip == set->scip);

   if( SCIPvarGetStatus(origvar) == SCIP_VARSTATUS_NEGATED )
   {
      assert(origvar->negatedvar != NULL);
      assert(SCIPvarGetStatus(origvar->negatedvar) == SCIP_VARSTATUS_ORIGINAL);

      if( origvar->negatedvar->data.original.transvar == NULL )
         *transvar = NULL;
      else
      {
         SCIP_CALL( SCIPvarNegate(origvar->negatedvar->data.original.transvar, blkmem, set, stat, transvar) );
      }
   }
   else 
      *transvar = origvar->data.original.transvar;

   return SCIP_OKAY;
}

/** converts loose transformed variable into column variable, creates LP column */
SCIP_RETCODE SCIPvarColumn(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(var->scip == set->scip);

   SCIPsetDebugMsg(set, "creating column for variable <%s>\n", var->name);

   /* switch variable status */
   var->varstatus = SCIP_VARSTATUS_COLUMN; /*lint !e641*/

   /* create column of variable */
   SCIP_CALL( SCIPcolCreate(&var->data.col, blkmem, set, stat, var, 0, NULL, NULL, var->removable) );

   if( var->probindex != -1 )
   {
      /* inform problem about the variable's status change */
      SCIP_CALL( SCIPprobVarChangedStatus(prob, blkmem, set, NULL, NULL, var) );

      /* inform LP, that problem variable is now a column variable and no longer loose */
      SCIP_CALL( SCIPlpUpdateVarColumn(lp, set, var) );
   }

   return SCIP_OKAY;
}

/** converts column transformed variable back into loose variable, frees LP column */
SCIP_RETCODE SCIPvarLoose(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(var->scip == set->scip);
   assert(var->data.col != NULL);
   assert(var->data.col->lppos == -1);
   assert(var->data.col->lpipos == -1);

   SCIPsetDebugMsg(set, "deleting column for variable <%s>\n", var->name);

   /* free column of variable */
   SCIP_CALL( SCIPcolFree(&var->data.col, blkmem, set, eventqueue, lp) );

   /* switch variable status */
   var->varstatus = SCIP_VARSTATUS_LOOSE; /*lint !e641*/

   if( var->probindex != -1 )
   {
      /* inform problem about the variable's status change */
      SCIP_CALL( SCIPprobVarChangedStatus(prob, blkmem, set, NULL, NULL, var) );

      /* inform LP, that problem variable is now a loose variable and no longer a column */
      SCIP_CALL( SCIPlpUpdateVarLoose(lp, set, var) );
   }

   return SCIP_OKAY;
}

/** issues a VARFIXED event on the given variable and all its parents (except ORIGINAL parents);
 *  the event issuing on the parents is necessary, because unlike with bound changes, the parent variables
 *  are not informed about a fixing of an active variable they are pointing to
 */
static
SCIP_RETCODE varEventVarFixed(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   fixeventtype        /**< is this event a fixation(0), an aggregation(1), or a
					      *   multi-aggregation(2)
					      */
   )
{
   SCIP_EVENT* event;
   SCIP_VARSTATUS varstatus;
   int i;

   assert(var != NULL);
   assert(var->scip == set->scip);
   assert(0 <= fixeventtype && fixeventtype <= 2);

   /* issue VARFIXED event on variable */
   SCIP_CALL( SCIPeventCreateVarFixed(&event, blkmem, var) );
   SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, NULL, &event) );

#ifndef NDEBUG
   for( i = var->nparentvars -1; i >= 0; --i )
   {
      assert(SCIPvarGetStatus(var->parentvars[i]) != SCIP_VARSTATUS_MULTAGGR);
   }
#endif

   switch( fixeventtype )
   {
   case 0:
      /* process all parents of a fixed variable */
      for( i = var->nparentvars - 1; i >= 0; --i )
      {
	 varstatus = SCIPvarGetStatus(var->parentvars[i]);

	 assert(varstatus != SCIP_VARSTATUS_FIXED);

	 /* issue event on all not yet fixed parent variables, (that should already issued this event) except the original
	  * one
	  */
	 if( varstatus != SCIP_VARSTATUS_ORIGINAL )
	 {
	    SCIP_CALL( varEventVarFixed(var->parentvars[i], blkmem, set, eventqueue, fixeventtype) );
	 }
      }
      break;
   case 1:
      /* process all parents of a aggregated variable */
      for( i = var->nparentvars - 1; i >= 0; --i )
      {
	 varstatus = SCIPvarGetStatus(var->parentvars[i]);

	 assert(varstatus != SCIP_VARSTATUS_FIXED);

	 /* issue event for not aggregated parent variable, because for these and its parents the var event was already
          * issued(, except the original one)
          *
          * @note that even before an aggregated parent variable, there might be variables, for which the vent was not
          *       yet issued
	  */
         if( varstatus == SCIP_VARSTATUS_AGGREGATED )
            continue;

	 if( varstatus != SCIP_VARSTATUS_ORIGINAL )
	 {
	    SCIP_CALL( varEventVarFixed(var->parentvars[i], blkmem, set, eventqueue, fixeventtype) );
	 }
      }
      break;
   case 2:
      /* process all parents of a aggregated variable */
      for( i = var->nparentvars - 1; i >= 0; --i )
      {
	 varstatus = SCIPvarGetStatus(var->parentvars[i]);

	 assert(varstatus != SCIP_VARSTATUS_FIXED);

	 /* issue event on all parent variables except the original one */
	 if( varstatus != SCIP_VARSTATUS_ORIGINAL )
	 {
	    SCIP_CALL( varEventVarFixed(var->parentvars[i], blkmem, set, eventqueue, fixeventtype) );
	 }
      }
      break;
   default:
      SCIPerrorMessage("unknown variable fixation event origin\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** converts variable into fixed variable */
SCIP_RETCODE SCIPvarFix(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real             fixedval,           /**< value to fix variable at */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            fixed               /**< pointer to store whether the fixing was performed (variable was unfixed) */
   )
{
   SCIP_Real obj;
   SCIP_Real childfixedval;

   assert(var != NULL);
   assert(var->scip == set->scip);
   assert(SCIPsetIsEQ(set, var->glbdom.lb, var->locdom.lb));
   assert(SCIPsetIsEQ(set, var->glbdom.ub, var->locdom.ub));
   assert(infeasible != NULL);
   assert(fixed != NULL);

   SCIPsetDebugMsg(set, "fix variable <%s>[%g,%g] to %g\n", var->name, var->glbdom.lb, var->glbdom.ub, fixedval);

   *infeasible = FALSE;
   *fixed = FALSE;

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
   {
      *infeasible = !SCIPsetIsFeasEQ(set, fixedval, var->locdom.lb);
      SCIPsetDebugMsg(set, " -> variable already fixed to %g (fixedval=%g): infeasible=%u\n", var->locdom.lb, fixedval, *infeasible);
      return SCIP_OKAY;
   }
   else if( (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && !SCIPsetIsFeasIntegral(set, fixedval))
      || SCIPsetIsFeasLT(set, fixedval, var->locdom.lb)
      || SCIPsetIsFeasGT(set, fixedval, var->locdom.ub) )
   {
      SCIPsetDebugMsg(set, " -> fixing infeasible: locdom=[%g,%g], fixedval=%g\n", var->locdom.lb, var->locdom.ub, fixedval);
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
      {
         SCIPerrorMessage("cannot fix an untransformed original variable\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarFix(var->data.original.transvar, blkmem, set, stat, transprob, origprob, primal, tree, reopt,
            lp, branchcand, eventqueue, cliquetable, fixedval, infeasible, fixed) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      assert(!SCIPeventqueueIsDelayed(eventqueue)); /* otherwise, the pseudo objective value update gets confused */

      /* set the fixed variable's objective value to 0.0 */
      obj = var->obj;
      SCIP_CALL( SCIPvarChgObj(var, blkmem, set, transprob, primal, lp, eventqueue, 0.0) );

      /* since we change the variable type form loose to fixed, we have to adjust the number of loose
       * variables in the LP data structure; the loose objective value (looseobjval) in the LP data structure, however,
       * gets adjusted automatically, due to the event SCIP_EVENTTYPE_OBJCHANGED which dropped in the moment where the
       * objective of this variable is set to zero
       */
      SCIPlpDecNLoosevars(lp);

      /* change variable's bounds to fixed value (thereby removing redundant implications and variable bounds) */
      holelistFree(&var->glbdom.holelist, blkmem);
      holelistFree(&var->locdom.holelist, blkmem);
      SCIP_CALL( SCIPvarChgLbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, fixedval) );
      SCIP_CALL( SCIPvarChgUbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, fixedval) );

      /* explicitly set variable's bounds, even if the fixed value is in epsilon range of the old bound */
      var->glbdom.lb = fixedval;
      var->glbdom.ub = fixedval;
      var->locdom.lb = fixedval;
      var->locdom.ub = fixedval;

      /* delete implications and variable bounds information */
      SCIP_CALL( SCIPvarRemoveCliquesImplicsVbs(var, blkmem, cliquetable, set, FALSE, FALSE, TRUE) );
      assert(var->vlbs == NULL);
      assert(var->vubs == NULL);
      assert(var->implics == NULL);
      assert(var->cliquelist == NULL);

      /* clear the history of the variable */
      SCIPhistoryReset(var->history);
      SCIPhistoryReset(var->historycrun);

      /* convert variable into fixed variable */
      var->varstatus = SCIP_VARSTATUS_FIXED; /*lint !e641*/

      /* inform problem about the variable's status change */
      if( var->probindex != -1 )
      {
         SCIP_CALL( SCIPprobVarChangedStatus(transprob, blkmem, set, branchcand, cliquetable, var) );
      }

      /* reset the objective value of the fixed variable, thus adjusting the problem's objective offset */
      SCIP_CALL( SCIPvarAddObj(var, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, eventqueue, obj) );

      /* issue VARFIXED event */
      SCIP_CALL( varEventVarFixed(var, blkmem, set, eventqueue, 0) );

      *fixed = TRUE;
      break;

   case SCIP_VARSTATUS_COLUMN:
      SCIPerrorMessage("cannot fix a column variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot fix a fixed variable again\n");  /*lint !e527*/
      SCIPABORT(); /* case is already handled in earlier if condition */
      return SCIP_INVALIDDATA;  /*lint !e527*/

   case SCIP_VARSTATUS_AGGREGATED:
      /* fix aggregation variable y in x = a*y + c, instead of fixing x directly */
      assert(SCIPsetIsZero(set, var->obj));
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      if( SCIPsetIsInfinity(set, fixedval) || SCIPsetIsInfinity(set, -fixedval) )
         childfixedval = (var->data.aggregate.scalar < 0.0 ? -fixedval : fixedval);
      else
         childfixedval = (fixedval - var->data.aggregate.constant)/var->data.aggregate.scalar;
      SCIP_CALL( SCIPvarFix(var->data.aggregate.var, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp,
            branchcand, eventqueue, cliquetable, childfixedval, infeasible, fixed) );
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot fix a multiple aggregated variable\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;  /*lint !e527*/

   case SCIP_VARSTATUS_NEGATED:
      /* fix negation variable x in x' = offset - x, instead of fixing x' directly */
      assert(SCIPsetIsZero(set, var->obj));
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarFix(var->negatedvar, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp,
            branchcand, eventqueue, cliquetable, var->data.negate.constant - fixedval, infeasible, fixed) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** transforms given variables, scalars and constant to the corresponding active variables, scalars and constant
 *
 * If the number of needed active variables is greater than the available slots in the variable array, nothing happens except
 * that the required size is stored in the corresponding variable; hence, if afterwards the required size is greater than the
 * available slots (varssize), nothing happens; otherwise, the active variable representation is stored in the arrays.
 *
 * The reason for this approach is that we cannot reallocate memory, since we do not know how the
 * memory has been allocated (e.g., by a C++ 'new' or SCIP functions).
 */
SCIP_RETCODE SCIPvarGetActiveRepresentatives(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< variable array to get active variables */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and scalars array */
   int                   varssize,           /**< available slots in vars and scalars array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   int*                  requiredsize,       /**< pointer to store the required array size for the active variables */
   SCIP_Bool             mergemultiples      /**< should multiple occurrences of a var be replaced by a single coeff? */
   )
{
   SCIP_VAR** activevars;
   SCIP_Real* activescalars;
   int nactivevars;
   SCIP_Real activeconstant;
   SCIP_Bool activeconstantinf;
   int activevarssize;

   SCIP_VAR* var;
   SCIP_Real scalar;
   int v;
   int k;

   SCIP_VAR** tmpvars;
   SCIP_VAR** multvars;
   SCIP_Real* tmpscalars;
   SCIP_Real* multscalars;
   int tmpvarssize;
   int ntmpvars;
   int nmultvars;

   SCIP_VAR* multvar;
   SCIP_Real multscalar;
   SCIP_Real multconstant;
   int pos;

   int noldtmpvars;

   SCIP_VAR** tmpvars2;
   SCIP_Real* tmpscalars2;
   int tmpvarssize2;
   int ntmpvars2;

   SCIP_Bool sortagain = FALSE;

   assert(set != NULL);
   assert(nvars != NULL);
   assert(scalars != NULL || *nvars == 0);
   assert(constant != NULL);
   assert(requiredsize != NULL);
   assert(*nvars <= varssize);

   *requiredsize = 0;

   if( *nvars == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);

   /* handle the "easy" case of just one variable and avoid memory allocation if the variable is already active */
   if( *nvars == 1 && (vars[0]->varstatus == ((int) SCIP_VARSTATUS_COLUMN) || vars[0]->varstatus == ((int) SCIP_VARSTATUS_LOOSE)) )
   {
      *requiredsize = 1;

      return SCIP_OKAY;
   }

   nactivevars = 0;
   activeconstant = 0.0;
   activeconstantinf = FALSE;
   activevarssize = (*nvars) * 2;
   ntmpvars = *nvars;
   tmpvarssize = *nvars;

   tmpvarssize2 = 1;

   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &tmpvars2, tmpvarssize2) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &tmpscalars2, tmpvarssize2) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activevars, activevarssize) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activescalars, activevarssize) );
   SCIP_CALL( SCIPsetDuplicateBufferArray(set, &tmpvars, vars, ntmpvars) );
   SCIP_CALL( SCIPsetDuplicateBufferArray(set, &tmpscalars, scalars, ntmpvars) );

   /* to avoid unnecessary expanding of variable arrays while disaggregating several variables multiple times combine same variables
    * first, first get all corresponding variables with status loose, column, multaggr or fixed
    */
   for( v = ntmpvars - 1; v >= 0; --v )
   {
      var = tmpvars[v];
      scalar = tmpscalars[v];

      assert(var != NULL);
      /* transforms given variable, scalar and constant to the corresponding active, fixed, or
       * multi-aggregated variable, scalar and constant; if the variable resolves to a fixed
       * variable, "scalar" will be 0.0 and the value of the sum will be stored in "constant".
       */
      SCIP_CALL( SCIPvarGetProbvarSum(&var, set, &scalar, &activeconstant) );
      assert(var != NULL);

      assert(SCIPsetIsInfinity(set, activeconstant) == (activeconstant == SCIPsetInfinity(set))); /*lint !e777*/
      assert(SCIPsetIsInfinity(set, -activeconstant) == (activeconstant == -SCIPsetInfinity(set))); /*lint !e777*/

      activeconstantinf = SCIPsetIsInfinity(set, activeconstant) || SCIPsetIsInfinity(set, -activeconstant);

      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED);

      tmpvars[v] = var;
      tmpscalars[v] = scalar;
   }
   noldtmpvars = ntmpvars;

   /* sort all variables to combine equal variables easily */
   SCIPsortPtrReal((void**)tmpvars, tmpscalars, SCIPvarComp, noldtmpvars);
   ntmpvars = 0;
   for( v = 1; v < noldtmpvars; ++v )
   {
      /* combine same variables */
      if( SCIPvarCompare(tmpvars[v], tmpvars[ntmpvars]) == 0 )
      {
         tmpscalars[ntmpvars] += tmpscalars[v];
      }
      else
      {
         ++ntmpvars;
         if( v > ntmpvars )
         {
            tmpscalars[ntmpvars] = tmpscalars[v];
            tmpvars[ntmpvars] = tmpvars[v];
         }
      }
   }
   ++ntmpvars;

#ifdef SCIP_MORE_DEBUG
   for( v = 1; v < ntmpvars; ++v )
      assert(SCIPvarCompare(tmpvars[v], tmpvars[v-1]) > 0);
#endif


   /* collect for each variable the representation in active variables */
   while( ntmpvars >= 1 )
   {
      --ntmpvars;
      ntmpvars2 = 0;
      var = tmpvars[ntmpvars];
      scalar = tmpscalars[ntmpvars];

      assert(var != NULL);

      /* TODO: maybe we should test here on SCIPsetIsZero() instead of 0.0 */
      if( scalar == 0.0 )
         continue;

      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED);

      switch( SCIPvarGetStatus(var) )
      {
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
         /* x = a*y + c */
         if( nactivevars >= activevarssize )
         {
            activevarssize *= 2;
            SCIP_CALL( SCIPsetReallocBufferArray(set, &activevars, activevarssize) );
            SCIP_CALL( SCIPsetReallocBufferArray(set, &activescalars, activevarssize) );
            assert(nactivevars < activevarssize);
         }
         activevars[nactivevars] = var;
         activescalars[nactivevars] = scalar;
         nactivevars++;
         break;

      case SCIP_VARSTATUS_MULTAGGR:
         /* x = a_1*y_1 + ... + a_n*y_n + c */
         nmultvars = var->data.multaggr.nvars;
         multvars = var->data.multaggr.vars;
         multscalars = var->data.multaggr.scalars;
         sortagain = TRUE;

         if( nmultvars + ntmpvars > tmpvarssize )
         {
            while( nmultvars + ntmpvars > tmpvarssize )
               tmpvarssize *= 2;
            SCIP_CALL( SCIPsetReallocBufferArray(set, &tmpvars, tmpvarssize) );
            SCIP_CALL( SCIPsetReallocBufferArray(set, &tmpscalars, tmpvarssize) );
            assert(nmultvars + ntmpvars <= tmpvarssize);
         }

         if( nmultvars > tmpvarssize2 )
         {
            while( nmultvars > tmpvarssize2 )
               tmpvarssize2 *= 2;
            SCIP_CALL( SCIPsetReallocBufferArray(set, &tmpvars2, tmpvarssize2) );
            SCIP_CALL( SCIPsetReallocBufferArray(set, &tmpscalars2, tmpvarssize2) );
            assert(nmultvars <= tmpvarssize2);
         }

         --nmultvars;

         for( ; nmultvars >= 0; --nmultvars )
         {
            multvar = multvars[nmultvars];
            multscalar = multscalars[nmultvars];
            multconstant = 0;

            assert(multvar != NULL);
            SCIP_CALL( SCIPvarGetProbvarSum(&multvar, set, &multscalar, &multconstant) );
            assert(multvar != NULL);

            assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
               || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
               || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR
               || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED);

            if( !activeconstantinf )
            {
               assert(!SCIPsetIsInfinity(set, scalar) && !SCIPsetIsInfinity(set, -scalar));

               if( SCIPsetIsInfinity(set, multconstant) || SCIPsetIsInfinity(set, -multconstant) )
               {
                  assert(scalar != 0.0);
                  if( scalar * multconstant > 0.0 )
                  {
                     activeconstant = SCIPsetInfinity(set);
                     activeconstantinf = TRUE;
                  }
                  else
                  {
                     activeconstant = -SCIPsetInfinity(set);
                     activeconstantinf = TRUE;
                  }
               }
               else
                  activeconstant += scalar * multconstant;
            }
#ifndef NDEBUG
            else
            {
               assert(!SCIPsetIsInfinity(set, activeconstant) || !(scalar * multconstant < 0.0 &&
                     (SCIPsetIsInfinity(set, multconstant) || SCIPsetIsInfinity(set, -multconstant))));
               assert(!SCIPsetIsInfinity(set, -activeconstant) || !(scalar * multconstant > 0.0 &&
                     (SCIPsetIsInfinity(set, multconstant) || SCIPsetIsInfinity(set, -multconstant))));
            }
#endif

            if( SCIPsortedvecFindPtr((void**)tmpvars, SCIPvarComp, multvar, ntmpvars, &pos) )
            {
               assert(SCIPvarCompare(tmpvars[pos], multvar) == 0);
               tmpscalars[pos] += scalar * multscalar;
            }
            else
            {
               tmpvars2[ntmpvars2] = multvar;
               tmpscalars2[ntmpvars2] = scalar * multscalar;
               ++(ntmpvars2);
               assert(ntmpvars2 <= tmpvarssize2);
            }
         }

         if( ntmpvars2 > 0 )
         {
            /* sort all variables to combine equal variables easily */
            SCIPsortPtrReal((void**)tmpvars2, tmpscalars2, SCIPvarComp, ntmpvars2);
            pos = 0;
            for( v = 1; v < ntmpvars2; ++v )
            {
               /* combine same variables */
               if( SCIPvarCompare(tmpvars2[v], tmpvars2[pos]) == 0 )
               {
                  tmpscalars2[pos] += tmpscalars2[v];
               }
               else
               {
                  ++pos;
                  if( v > pos )
                  {
                     tmpscalars2[pos] = tmpscalars2[v];
                     tmpvars2[pos] = tmpvars2[v];
                  }
               }
            }
            ntmpvars2 = pos + 1;
#ifdef SCIP_MORE_DEBUG
            for( v = 1; v < ntmpvars2; ++v )
            {
               assert(SCIPvarCompare(tmpvars2[v], tmpvars2[v-1]) > 0);
            }
            for( v = 1; v < ntmpvars; ++v )
            {
               assert(SCIPvarCompare(tmpvars[v], tmpvars[v-1]) > 0);
            }
#endif
            v = ntmpvars - 1;
            k = ntmpvars2 - 1;
            pos = ntmpvars + ntmpvars2 - 1;
            ntmpvars += ntmpvars2;

            while( v >= 0 && k >= 0 )
            {
               assert(pos >= 0);
               assert(SCIPvarCompare(tmpvars[v], tmpvars2[k]) != 0);
               if( SCIPvarCompare(tmpvars[v], tmpvars2[k]) >= 0 )
               {
                  tmpvars[pos] = tmpvars[v];
                  tmpscalars[pos] = tmpscalars[v];
                  --v;
               }
               else
               {
                  tmpvars[pos] = tmpvars2[k];
                  tmpscalars[pos] = tmpscalars2[k];
                  --k;
               }
               --pos;
               assert(pos >= 0);
            }
            while( v >= 0 )
            {
               assert(pos >= 0);
               tmpvars[pos] = tmpvars[v];
               tmpscalars[pos] = tmpscalars[v];
               --v;
               --pos;
            }
            while( k >= 0 )
            {
               assert(pos >= 0);
               tmpvars[pos] = tmpvars2[k];
               tmpscalars[pos] = tmpscalars2[k];
               --k;
               --pos;
            }
         }
#ifdef SCIP_MORE_DEBUG
         for( v = 1; v < ntmpvars; ++v )
         {
            assert(SCIPvarCompare(tmpvars[v], tmpvars[v-1]) > 0);
         }
#endif

         if( !activeconstantinf )
         {
            assert(!SCIPsetIsInfinity(set, scalar) && !SCIPsetIsInfinity(set, -scalar));

            multconstant = SCIPvarGetMultaggrConstant(var);

            if( SCIPsetIsInfinity(set, multconstant) || SCIPsetIsInfinity(set, -multconstant) )
            {
               assert(scalar != 0.0);
               if( scalar * multconstant > 0.0 )
               {
                  activeconstant = SCIPsetInfinity(set);
                  activeconstantinf = TRUE;
               }
               else
               {
                  activeconstant = -SCIPsetInfinity(set);
                  activeconstantinf = TRUE;
               }
            }
            else
               activeconstant += scalar * multconstant;
         }
#ifndef NDEBUG
         else
         {
            multconstant = SCIPvarGetMultaggrConstant(var);
            assert(!SCIPsetIsInfinity(set, activeconstant) || !(scalar * multconstant < 0.0 &&
                  (SCIPsetIsInfinity(set, multconstant) || SCIPsetIsInfinity(set, -multconstant))));
            assert(!SCIPsetIsInfinity(set, -activeconstant) || !(scalar * multconstant > 0.0 &&
                  (SCIPsetIsInfinity(set, multconstant) || SCIPsetIsInfinity(set, -multconstant))));
         }
#endif
         break;

      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_ORIGINAL:
      case SCIP_VARSTATUS_AGGREGATED:
      case SCIP_VARSTATUS_NEGATED:
      default:
         /* case x = c, but actually we should not be here, since SCIPvarGetProbvarSum() returns a scalar of 0.0 for
          * fixed variables and is handled already
          */
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED);
         assert(SCIPsetIsZero(set, var->glbdom.lb) && SCIPsetIsEQ(set, var->glbdom.lb, var->glbdom.ub));
      }
   }

   if( mergemultiples )
   {
      if( sortagain )
      {
         /* sort variable and scalar array by variable index */
         SCIPsortPtrReal((void**)activevars, activescalars, SCIPvarComp, nactivevars);

         /* eliminate duplicates and count required size */
         v = nactivevars - 1;
         while( v > 0 )
         {
            /* combine both variable since they are the same */
            if( SCIPvarCompare(activevars[v - 1], activevars[v]) == 0 )
            {
               if( activescalars[v - 1] + activescalars[v] != 0.0 )
               {
                  activescalars[v - 1] += activescalars[v];
                  --nactivevars;
                  activevars[v] = activevars[nactivevars];
                  activescalars[v] = activescalars[nactivevars];
               }
               else
               {
                  --nactivevars;
                  activevars[v] = activevars[nactivevars];
                  activescalars[v] = activescalars[nactivevars];
                  --nactivevars;
                  --v;
                  activevars[v] = activevars[nactivevars];
                  activescalars[v] = activescalars[nactivevars];
               }
            }
            --v;
         }
      }
      /* the variables were added in reverse order, we revert the order now;
       * this should not be necessary, but not doing this changes the behavior sometimes
       */
      else
      {
         SCIP_VAR* tmpvar;
         SCIP_Real tmpscalar;

         for( v = 0; v < nactivevars / 2; ++v )
         {
            tmpvar = activevars[v];
            tmpscalar = activescalars[v];
            activevars[v] = activevars[nactivevars - 1 - v];
            activescalars[v] = activescalars[nactivevars - 1 - v];
            activevars[nactivevars - 1 - v] = tmpvar;
            activescalars[nactivevars - 1 - v] = tmpscalar;
         }
      }
   }
   *requiredsize = nactivevars;

   if( varssize >= *requiredsize )
   {
      assert(vars != NULL);

      *nvars = *requiredsize;

      if( !SCIPsetIsInfinity(set, *constant) && !SCIPsetIsInfinity(set, -(*constant)) )
      {
         /* if the activeconstant is infinite, the constant pointer gets the same value, otherwise add the value */
         if( activeconstantinf )
            (*constant) = activeconstant;
         else
            (*constant) += activeconstant;
      }
#ifndef NDEBUG
      else
      {
         assert(!SCIPsetIsInfinity(set, (*constant)) || !SCIPsetIsInfinity(set, -activeconstant));
         assert(!SCIPsetIsInfinity(set, -(*constant)) || !SCIPsetIsInfinity(set, activeconstant));
      }
#endif

      /* copy active variable and scalar array to the given arrays */
      for( v = 0; v < *nvars; ++v )
      {
         vars[v] = activevars[v];
         scalars[v] = activescalars[v]; /*lint !e613*/
      }
   }

   assert(SCIPsetIsInfinity(set, *constant) == ((*constant) == SCIPsetInfinity(set))); /*lint !e777*/
   assert(SCIPsetIsInfinity(set, -(*constant)) == ((*constant) == -SCIPsetInfinity(set))); /*lint !e777*/

   SCIPsetFreeBufferArray(set, &tmpscalars);
   SCIPsetFreeBufferArray(set, &tmpvars);
   SCIPsetFreeBufferArray(set, &activescalars);
   SCIPsetFreeBufferArray(set, &activevars);
   SCIPsetFreeBufferArray(set, &tmpscalars2);
   SCIPsetFreeBufferArray(set, &tmpvars2);

   return SCIP_OKAY;
}


/** flattens aggregation graph of multi-aggregated variable in order to avoid exponential recursion later on */
SCIP_RETCODE SCIPvarFlattenAggregationGraph(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_Real multconstant;
   int multvarssize;
   int nmultvars;
   int multrequiredsize;

   assert( var != NULL );
   assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR );
   assert(var->scip == set->scip);

   multconstant = var->data.multaggr.constant;
   nmultvars = var->data.multaggr.nvars;
   multvarssize = var->data.multaggr.varssize;

   SCIP_CALL( SCIPvarGetActiveRepresentatives(set, var->data.multaggr.vars, var->data.multaggr.scalars, &nmultvars, multvarssize, &multconstant, &multrequiredsize, TRUE) );

   if( multrequiredsize > multvarssize )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(var->data.multaggr.vars), multvarssize, multrequiredsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(var->data.multaggr.scalars), multvarssize, multrequiredsize) );
      multvarssize = multrequiredsize;
      SCIP_CALL( SCIPvarGetActiveRepresentatives(set, var->data.multaggr.vars, var->data.multaggr.scalars, &nmultvars, multvarssize, &multconstant, &multrequiredsize, TRUE) );
      assert( multrequiredsize <= multvarssize );
   }
   /**@note After the flattening the multi aggregation might resolve to be in fact an aggregation (or even a fixing?).
    * This issue is not resolved right now, since var->data.multaggr.nvars < 2 should not cause troubles. However, one
    * may loose performance hereby, since aggregated variables are easier to handle.
    * 
    * Note, that there are two cases where SCIPvarFlattenAggregationGraph() is called: The easier one is that it is
    * called while installing the multi-aggregation. in principle, the described issue could be handled straightforward
    * in this case by aggregating or fixing the variable instead.  The more complicated case is the one, when the
    * multi-aggregation is used, e.g., in linear presolving (and the variable is already declared to be multi-aggregated).
    *
    * By now, it is not allowed to fix or aggregate multi-aggregated variables which would be necessary in this case.
    *
    * The same issue appears in the SCIPvarGetProbvar...() methods.
    */

   var->data.multaggr.constant = multconstant;
   var->data.multaggr.nvars = nmultvars;
   var->data.multaggr.varssize = multvarssize;

   return SCIP_OKAY;
}

/** merge two variable histories together; a typical use case is that \p othervar is an image of the target variable
 *  in a SCIP copy. Method should be applied with care, especially because no internal checks are performed whether
 *  the history merge is reasonable
 *
 *  @note Do not use this method if the two variables originate from two SCIP's with different objective functions, since
 *        this corrupts the variable pseudo costs
 *  @note Apply with care; no internal checks are performed if the two variables should be merged
 */
void SCIPvarMergeHistories(
   SCIP_VAR*             targetvar,          /**< the variable that should contain both histories afterwards */
   SCIP_VAR*             othervar,           /**< the variable whose history is to be merged with that of the target variable */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   /* merge only the history of the current run into the target history */
   SCIPhistoryUnite(targetvar->history, othervar->historycrun, FALSE);

   /* apply the changes also to the global history */
   SCIPhistoryUnite(stat->glbhistory, othervar->historycrun, FALSE);
}

/** sets the history of a variable; this method is typically used within reoptimization to keep and update the variable
 *  history over several iterations
 */
void SCIPvarSetHistory(
   SCIP_VAR*             var,                /**< variable */
   SCIP_HISTORY*         history,            /**< the history which is to set */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   /* merge only the history of the current run into the target history */
   SCIPhistoryUnite(var->history, history, FALSE);

   /* apply the changes also to the global history */
   SCIPhistoryUnite(stat->glbhistory, history, FALSE);
}

/** tightens the bounds of both variables in aggregation x = a*y + c */
static
SCIP_RETCODE varUpdateAggregationBounds(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             aggvar,             /**< variable y in aggregation x = a*y + c */
   SCIP_Real             scalar,             /**< multiplier a in aggregation x = a*y + c */
   SCIP_Real             constant,           /**< constant shift c in aggregation x = a*y + c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            fixed               /**< pointer to store whether the variables were fixed */
   )
{
   SCIP_Real varlb;
   SCIP_Real varub;
   SCIP_Real aggvarlb;
   SCIP_Real aggvarub;
   SCIP_Bool aggvarbdschanged;

   assert(var != NULL);
   assert(var->scip == set->scip);
   assert(aggvar != NULL);
   assert(!SCIPsetIsZero(set, scalar));
   assert(infeasible != NULL);
   assert(fixed != NULL);

   *infeasible = FALSE;
   *fixed = FALSE;

   SCIPsetDebugMsg(set, "updating bounds of variables in aggregation <%s> == %g*<%s> %+g\n", var->name, scalar, aggvar->name, constant);
   SCIPsetDebugMsg(set, "  old bounds: <%s> [%g,%g]   <%s> [%g,%g]\n",
      var->name, var->glbdom.lb, var->glbdom.ub, aggvar->name, aggvar->glbdom.lb, aggvar->glbdom.ub);

   /* loop as long additional changes may be found */
   do
   {
      aggvarbdschanged = FALSE;

      /* update the bounds of the aggregated variable x in x = a*y + c */
      if( scalar > 0.0 )
      {
         if( SCIPsetIsInfinity(set, -aggvar->glbdom.lb) )
            varlb = -SCIPsetInfinity(set);
         else
            varlb = aggvar->glbdom.lb * scalar + constant;
         if( SCIPsetIsInfinity(set, aggvar->glbdom.ub) )
            varub = SCIPsetInfinity(set);
         else
            varub = aggvar->glbdom.ub * scalar + constant;
      }
      else
      {
         if( SCIPsetIsInfinity(set, -aggvar->glbdom.lb) )
            varub = SCIPsetInfinity(set);
         else
            varub = aggvar->glbdom.lb * scalar + constant;
         if( SCIPsetIsInfinity(set, aggvar->glbdom.ub) )
            varlb = -SCIPsetInfinity(set);
         else
            varlb = aggvar->glbdom.ub * scalar + constant;
      }
      varlb = MAX(varlb, var->glbdom.lb);
      varub = MIN(varub, var->glbdom.ub);
      SCIPvarAdjustLb(var, set, &varlb);
      SCIPvarAdjustUb(var, set, &varub);

      /* check the new bounds */
      if( SCIPsetIsGT(set, varlb, varub) )
      {
         /* the aggregation is infeasible */
         *infeasible = TRUE;
         return SCIP_OKAY;
      }
      else if( SCIPsetIsEQ(set, varlb, varub) )
      {
         /* the aggregated variable is fixed -> fix both variables */
         SCIP_CALL( SCIPvarFix(var, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand,
               eventqueue, cliquetable, varlb, infeasible, fixed) );
         if( !(*infeasible) )
         {
            SCIP_Bool aggfixed;

            SCIP_CALL( SCIPvarFix(aggvar, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand,
                  eventqueue, cliquetable, (varlb-constant)/scalar, infeasible, &aggfixed) );
            assert(*fixed == aggfixed);
         }
         return SCIP_OKAY;
      }
      else
      {
         if( SCIPsetIsGT(set, varlb, var->glbdom.lb) )
         {
            SCIP_CALL( SCIPvarChgLbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, varlb) );
         }
         if( SCIPsetIsLT(set, varub, var->glbdom.ub) )
         {
            SCIP_CALL( SCIPvarChgUbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, varub) );
         }

         /* update the hole list of the aggregation variable */
         /**@todo update hole list of aggregation variable */
      }

      /* update the bounds of the aggregation variable y in x = a*y + c  ->  y = (x-c)/a */
      if( scalar > 0.0 )
      {
         if( SCIPsetIsInfinity(set, -var->glbdom.lb) )
            aggvarlb = -SCIPsetInfinity(set);
         else
            aggvarlb = (var->glbdom.lb - constant) / scalar;
         if( SCIPsetIsInfinity(set, var->glbdom.ub) )
            aggvarub = SCIPsetInfinity(set);
         else
            aggvarub = (var->glbdom.ub - constant) / scalar;
      }
      else
      {
         if( SCIPsetIsInfinity(set, -var->glbdom.lb) )
            aggvarub = SCIPsetInfinity(set);
         else
            aggvarub = (var->glbdom.lb - constant) / scalar;
         if( SCIPsetIsInfinity(set, var->glbdom.ub) )
            aggvarlb = -SCIPsetInfinity(set);
         else
            aggvarlb = (var->glbdom.ub - constant) / scalar;
      }
      aggvarlb = MAX(aggvarlb, aggvar->glbdom.lb);
      aggvarub = MIN(aggvarub, aggvar->glbdom.ub);
      SCIPvarAdjustLb(aggvar, set, &aggvarlb);
      SCIPvarAdjustUb(aggvar, set, &aggvarub);

      /* check the new bounds */
      if( SCIPsetIsGT(set, aggvarlb, aggvarub) )
      {
         /* the aggregation is infeasible */
         *infeasible = TRUE;
         return SCIP_OKAY;
      }
      else if( SCIPsetIsEQ(set, aggvarlb, aggvarub) )
      {
         /* the aggregation variable is fixed -> fix both variables */
         SCIP_CALL( SCIPvarFix(aggvar, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand,
               eventqueue, cliquetable, aggvarlb, infeasible, fixed) );
         if( !(*infeasible) )
         {
            SCIP_Bool varfixed;

            SCIP_CALL( SCIPvarFix(var, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand,
                  eventqueue, cliquetable, aggvarlb * scalar + constant, infeasible, &varfixed) );
            assert(*fixed == varfixed);
         }
         return SCIP_OKAY;
      }
      else
      {
         SCIP_Real oldbd;
         if( SCIPsetIsGT(set, aggvarlb, aggvar->glbdom.lb) )
         {
            oldbd = aggvar->glbdom.lb;
            SCIP_CALL( SCIPvarChgLbGlobal(aggvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, aggvarlb) );
            aggvarbdschanged = !SCIPsetIsEQ(set, oldbd, aggvar->glbdom.lb);
         }
         if( SCIPsetIsLT(set, aggvarub, aggvar->glbdom.ub) )
         {
            oldbd = aggvar->glbdom.ub;
            SCIP_CALL( SCIPvarChgUbGlobal(aggvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, aggvarub) );
            aggvarbdschanged = aggvarbdschanged || !SCIPsetIsEQ(set, oldbd, aggvar->glbdom.ub);
         }

         /* update the hole list of the aggregation variable */
         /**@todo update hole list of aggregation variable */
      }
   }
   while( aggvarbdschanged );

   SCIPsetDebugMsg(set, "  new bounds: <%s> [%g,%g]   <%s> [%g,%g]\n",
      var->name, var->glbdom.lb, var->glbdom.ub, aggvar->name, aggvar->glbdom.lb, aggvar->glbdom.ub);

   return SCIP_OKAY;
}

/** converts loose variable into aggregated variable */
SCIP_RETCODE SCIPvarAggregate(
   SCIP_VAR*             var,                /**< loose problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             aggvar,             /**< loose variable y in aggregation x = a*y + c */
   SCIP_Real             scalar,             /**< multiplier a in aggregation x = a*y + c */
   SCIP_Real             constant,           /**< constant shift c in aggregation x = a*y + c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Real* constants;
   SCIP_Real obj;
   SCIP_Real branchfactor;
   SCIP_Bool fixed;
   int branchpriority;
   int nlocksdown;
   int nlocksup;
   int nvbds;
   int i;
   int j;

   assert(var != NULL);
   assert(aggvar != NULL);
   assert(var->scip == set->scip);
   assert(var->glbdom.lb == var->locdom.lb); /*lint !e777*/
   assert(var->glbdom.ub == var->locdom.ub); /*lint !e777*/
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(!SCIPeventqueueIsDelayed(eventqueue)); /* otherwise, the pseudo objective value update gets confused */
   assert(infeasible != NULL);
   assert(aggregated != NULL);

   *infeasible = FALSE;
   *aggregated = FALSE;

   /* get active problem variable of aggregation variable */
   SCIP_CALL( SCIPvarGetProbvarSum(&aggvar, set, &scalar, &constant) );

   /* aggregation is a fixing, if the scalar is zero */
   if( SCIPsetIsZero(set, scalar) )
   {
      SCIP_CALL( SCIPvarFix(var, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, eventqueue,
            cliquetable, constant, infeasible, aggregated) );
      return SCIP_OKAY;
   }

   /* don't perform the aggregation if the aggregation variable is multi-aggregated itself */
   if( SCIPvarGetStatus(aggvar) == SCIP_VARSTATUS_MULTAGGR )
      return SCIP_OKAY;

   /**@todo currently we don't perform the aggregation if the aggregation variable has a non-empty hole list; this
    *  should be changed in the future
    */
   if( SCIPvarGetHolelistGlobal(var) != NULL )
      return SCIP_OKAY;

   assert(aggvar->glbdom.lb == aggvar->locdom.lb); /*lint !e777*/
   assert(aggvar->glbdom.ub == aggvar->locdom.ub); /*lint !e777*/
   assert(SCIPvarGetStatus(aggvar) == SCIP_VARSTATUS_LOOSE);

   SCIPsetDebugMsg(set, "aggregate variable <%s>[%g,%g] == %g*<%s>[%g,%g] %+g\n", var->name, var->glbdom.lb, var->glbdom.ub,
      scalar, aggvar->name, aggvar->glbdom.lb, aggvar->glbdom.ub, constant);

   /* if variable and aggregation variable are equal, the variable can be fixed: x == a*x + c  =>  x == c/(1-a) */
   if( var == aggvar )
   {
      if( SCIPsetIsEQ(set, scalar, 1.0) )
         *infeasible = !SCIPsetIsZero(set, constant);
      else
      {
         SCIP_CALL( SCIPvarFix(var, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand,
               eventqueue, cliquetable, constant/(1.0-scalar), infeasible, aggregated) );
      }
      return SCIP_OKAY;
   }

   /* tighten the bounds of aggregated and aggregation variable */
   SCIP_CALL( varUpdateAggregationBounds(var, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp,
         branchcand, eventqueue, cliquetable, aggvar, scalar, constant, infeasible, &fixed) );
   if( *infeasible || fixed )
   {
      *aggregated = fixed;
      return SCIP_OKAY;
   }

   /* delete implications and variable bounds of the aggregated variable from other variables, but keep them in the
    * aggregated variable
    */
   SCIP_CALL( SCIPvarRemoveCliquesImplicsVbs(var, blkmem, cliquetable, set, FALSE, FALSE, FALSE) );
   assert(var->cliquelist == NULL);

   /* set the aggregated variable's objective value to 0.0 */
   obj = var->obj;
   SCIP_CALL( SCIPvarChgObj(var, blkmem, set, transprob, primal, lp, eventqueue, 0.0) );

   /* unlock all rounding locks */
   nlocksdown = var->nlocksdown;
   nlocksup = var->nlocksup;
   var->nlocksdown = 0;
   var->nlocksup = 0;

   /* check, if variable should be used as NEGATED variable of the aggregation variable */
   if( SCIPvarIsBinary(var) && SCIPvarIsBinary(aggvar)
      && var->negatedvar == NULL && aggvar->negatedvar == NULL
      && SCIPsetIsEQ(set, scalar, -1.0) && SCIPsetIsEQ(set, constant, 1.0) )
   {
      /* link both variables as negation pair */
      var->varstatus = SCIP_VARSTATUS_NEGATED; /*lint !e641*/
      var->data.negate.constant = 1.0;
      var->negatedvar = aggvar;
      aggvar->negatedvar = var;

      /* copy doNotMultiaggr status */
      aggvar->donotmultaggr |= var->donotmultaggr;

      /* mark both variables to be non-deletable */
      SCIPvarMarkNotDeletable(var);
      SCIPvarMarkNotDeletable(aggvar);
   }
   else
   {
      /* convert variable into aggregated variable */
      var->varstatus = SCIP_VARSTATUS_AGGREGATED; /*lint !e641*/
      var->data.aggregate.var = aggvar;
      var->data.aggregate.scalar = scalar;
      var->data.aggregate.constant = constant;

      /* copy doNotMultiaggr status */
      aggvar->donotmultaggr |= var->donotmultaggr;

      /* mark both variables to be non-deletable */
      SCIPvarMarkNotDeletable(var);
      SCIPvarMarkNotDeletable(aggvar);
   }

   /* make aggregated variable a parent of the aggregation variable */
   SCIP_CALL( varAddParent(aggvar, blkmem, set, var) );

   /* relock the rounding locks of the variable, thus increasing the locks of the aggregation variable */
   SCIP_CALL( SCIPvarAddLocks(var, blkmem, set, eventqueue, nlocksdown, nlocksup) );

   /* move the variable bounds to the aggregation variable:
    *  - add all variable bounds again to the variable, thus adding it to the aggregation variable
    *  - free the variable bounds data structures
    */
   if( var->vlbs != NULL )
   {
      nvbds = SCIPvboundsGetNVbds(var->vlbs);
      vars = SCIPvboundsGetVars(var->vlbs);
      coefs = SCIPvboundsGetCoefs(var->vlbs);
      constants = SCIPvboundsGetConstants(var->vlbs);
      for( i = 0; i < nvbds && !(*infeasible); ++i )
      {
         SCIP_CALL( SCIPvarAddVlb(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable, branchcand,
               eventqueue, vars[i], coefs[i], constants[i], FALSE, infeasible, NULL) );
      }
   }
   if( var->vubs != NULL )
   {
      nvbds = SCIPvboundsGetNVbds(var->vubs);
      vars = SCIPvboundsGetVars(var->vubs);
      coefs = SCIPvboundsGetCoefs(var->vubs);
      constants = SCIPvboundsGetConstants(var->vubs);
      for( i = 0; i < nvbds && !(*infeasible); ++i )
      {
         SCIP_CALL( SCIPvarAddVub(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable, branchcand,
               eventqueue, vars[i], coefs[i], constants[i], FALSE, infeasible, NULL) );
      }
   }
   SCIPvboundsFree(&var->vlbs, blkmem);
   SCIPvboundsFree(&var->vubs, blkmem);

   /* move the implications to the aggregation variable:
    *  - add all implications again to the variable, thus adding it to the aggregation variable
    *  - free the implications data structures
    */
   if( var->implics != NULL && SCIPvarGetType(aggvar) == SCIP_VARTYPE_BINARY )
   {
      assert(SCIPvarIsBinary(var));
      for( i = 0; i < 2; ++i )
      {
         SCIP_VAR** implvars;
         SCIP_BOUNDTYPE* impltypes;
         SCIP_Real* implbounds;
         int nimpls;

         nimpls = SCIPimplicsGetNImpls(var->implics, (SCIP_Bool)i);
         implvars = SCIPimplicsGetVars(var->implics, (SCIP_Bool)i);
         impltypes = SCIPimplicsGetTypes(var->implics, (SCIP_Bool)i);
         implbounds = SCIPimplicsGetBounds(var->implics, (SCIP_Bool)i);

         for( j = 0; j < nimpls && !(*infeasible); ++j )
         {
            /* @todo can't we omit transitive closure, because it should already have been done when adding the
             *       implication to the aggregated variable?
             */
            SCIP_CALL( SCIPvarAddImplic(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable,
                  branchcand, eventqueue, (SCIP_Bool)i, implvars[j], impltypes[j], implbounds[j], FALSE, infeasible,
                  NULL) );
            assert(nimpls == SCIPimplicsGetNImpls(var->implics, (SCIP_Bool)i));
         }
      }
   }
   SCIPimplicsFree(&var->implics, blkmem);

   /* add the history entries to the aggregation variable and clear the history of the aggregated variable */
   SCIPhistoryUnite(aggvar->history, var->history, scalar < 0.0);
   SCIPhistoryUnite(aggvar->historycrun, var->historycrun, scalar < 0.0);
   SCIPhistoryReset(var->history);
   SCIPhistoryReset(var->historycrun);

   /* update flags of aggregation variable */
   aggvar->removable &= var->removable;

   /* update branching factors and priorities of both variables to be the maximum of both variables */
   branchfactor = MAX(aggvar->branchfactor, var->branchfactor);
   branchpriority = MAX(aggvar->branchpriority, var->branchpriority);
   SCIP_CALL( SCIPvarChgBranchFactor(aggvar, set, branchfactor) );
   SCIP_CALL( SCIPvarChgBranchPriority(aggvar, branchpriority) );
   SCIP_CALL( SCIPvarChgBranchFactor(var, set, branchfactor) );
   SCIP_CALL( SCIPvarChgBranchPriority(var, branchpriority) );

   /* update branching direction of both variables to agree to a single direction */
   if( scalar >= 0.0 )
   {
      if( (SCIP_BRANCHDIR)var->branchdirection == SCIP_BRANCHDIR_AUTO )
      {
         SCIP_CALL( SCIPvarChgBranchDirection(var, (SCIP_BRANCHDIR)aggvar->branchdirection) );
      }
      else if( (SCIP_BRANCHDIR)aggvar->branchdirection == SCIP_BRANCHDIR_AUTO )
      {
         SCIP_CALL( SCIPvarChgBranchDirection(aggvar, (SCIP_BRANCHDIR)var->branchdirection) );
      }
      else if( var->branchdirection != aggvar->branchdirection )
      {
         SCIP_CALL( SCIPvarChgBranchDirection(var, SCIP_BRANCHDIR_AUTO) );
      }
   }
   else
   {
      if( (SCIP_BRANCHDIR)var->branchdirection == SCIP_BRANCHDIR_AUTO )
      {
         SCIP_CALL( SCIPvarChgBranchDirection(var, SCIPbranchdirOpposite((SCIP_BRANCHDIR)aggvar->branchdirection)) );
      }
      else if( (SCIP_BRANCHDIR)aggvar->branchdirection == SCIP_BRANCHDIR_AUTO )
      {
         SCIP_CALL( SCIPvarChgBranchDirection(aggvar, SCIPbranchdirOpposite((SCIP_BRANCHDIR)var->branchdirection)) );
      }
      else if( var->branchdirection != aggvar->branchdirection )
      {
         SCIP_CALL( SCIPvarChgBranchDirection(var, SCIP_BRANCHDIR_AUTO) );
      }
   }

   if( var->probindex != -1 )
   {
      /* inform problem about the variable's status change */
      SCIP_CALL( SCIPprobVarChangedStatus(transprob, blkmem, set, branchcand, cliquetable, var) );
   }

   /* reset the objective value of the aggregated variable, thus adjusting the objective value of the aggregation
    * variable and the problem's objective offset
    */
   SCIP_CALL( SCIPvarAddObj(var, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, eventqueue, obj) );

   /* issue VARFIXED event */
   SCIP_CALL( varEventVarFixed(var, blkmem, set, eventqueue, 1) );

   *aggregated = TRUE;

   return SCIP_OKAY;
}

/** Tries to aggregate an equality a*x + b*y == c consisting of two (implicit) integral active problem variables x and
 *  y. An integer aggregation (i.e. integral coefficients a' and b', such that a'*x + b'*y == c') is searched.
 *
 *  This can lead to the detection of infeasibility (e.g. if c' is fractional), or to a rejection of the aggregation
 *  (denoted by aggregated == FALSE), if the resulting integer coefficients are too large and thus numerically instable.
 */
static
SCIP_RETCODE tryAggregateIntVars(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             varx,               /**< integral variable x in equality a*x + b*y == c */
   SCIP_VAR*             vary,               /**< integral variable y in equality a*x + b*y == c */
   SCIP_Real             scalarx,            /**< multiplier a in equality a*x + b*y == c */
   SCIP_Real             scalary,            /**< multiplier b in equality a*x + b*y == c */
   SCIP_Real             rhs,                /**< right hand side c in equality a*x + b*y == c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   )
{
   SCIP_VAR* aggvar;
   char aggvarname[SCIP_MAXSTRLEN];
   SCIP_Longint scalarxn = 0;
   SCIP_Longint scalarxd = 0;
   SCIP_Longint scalaryn = 0;
   SCIP_Longint scalaryd = 0;
   SCIP_Longint a;
   SCIP_Longint b;
   SCIP_Longint c;
   SCIP_Longint scm;
   SCIP_Longint gcd;
   SCIP_Longint currentclass;
   SCIP_Longint classstep;
   SCIP_Longint xsol;
   SCIP_Longint ysol;
   SCIP_Bool success;
   SCIP_VARTYPE vartype;

#define MAXDNOM 1000000LL

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(cliquetable != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(varx != NULL);
   assert(vary != NULL);
   assert(varx != vary);
   assert(infeasible != NULL);
   assert(aggregated != NULL);
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PRESOLVING);
   assert(SCIPvarGetStatus(varx) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetType(varx) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(varx) == SCIP_VARTYPE_IMPLINT);
   assert(SCIPvarGetStatus(vary) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetType(vary) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(vary) == SCIP_VARTYPE_IMPLINT);
   assert(!SCIPsetIsZero(set, scalarx));
   assert(!SCIPsetIsZero(set, scalary));

   *infeasible = FALSE;
   *aggregated = FALSE;

   /* get rational representation of coefficients */
   success = SCIPrealToRational(scalarx, -SCIPsetEpsilon(set), SCIPsetEpsilon(set), MAXDNOM, &scalarxn, &scalarxd);
   if( success )
      success = SCIPrealToRational(scalary, -SCIPsetEpsilon(set), SCIPsetEpsilon(set), MAXDNOM, &scalaryn, &scalaryd);
   if( !success )
      return SCIP_OKAY;
   assert(scalarxd >= 1);
   assert(scalaryd >= 1);

   /* multiply equality with smallest common denominator */
   scm = SCIPcalcSmaComMul(scalarxd, scalaryd);
   a = (scm/scalarxd)*scalarxn;
   b = (scm/scalaryd)*scalaryn;
   rhs *= scm;

   /* divide equality by the greatest common divisor of a and b */
   gcd = SCIPcalcGreComDiv(ABS(a), ABS(b));
   a /= gcd;
   b /= gcd;
   rhs /= gcd;
   assert(a != 0);
   assert(b != 0);

   /* check, if right hand side is integral */
   if( !SCIPsetIsFeasIntegral(set, rhs) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   c = (SCIP_Longint)(SCIPsetFeasFloor(set, rhs));

   if( REALABS((SCIP_Real)(c/a)) > SCIPsetGetHugeValue(set) * SCIPsetFeastol(set) ) /*lint !e653*/
      return SCIP_OKAY;

   /* check, if we are in an easy case with either |a| = 1 or |b| = 1 */
   if( (a == 1 || a == -1) && SCIPvarGetType(vary) == SCIP_VARTYPE_INTEGER )
   {
      /* aggregate x = - b/a*y + c/a */
      /*lint --e{653}*/
      SCIP_CALL( SCIPvarAggregate(varx, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, cliquetable,
            branchcand, eventqueue, vary, (SCIP_Real)(-b/a), (SCIP_Real)(c/a), infeasible, aggregated) );
      assert(*aggregated);
      return SCIP_OKAY;
   }
   if( (b == 1 || b == -1) && SCIPvarGetType(varx) == SCIP_VARTYPE_INTEGER )
   {
      /* aggregate y = - a/b*x + c/b */
      /*lint --e{653}*/
      SCIP_CALL( SCIPvarAggregate(vary, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, cliquetable,
            branchcand, eventqueue, varx, (SCIP_Real)(-a/b), (SCIP_Real)(c/b), infeasible, aggregated) );
      assert(*aggregated);
      return SCIP_OKAY;
   }

   /* Both variables are integers, their coefficients are not multiples of each other, and they don't have any
    * common divisor. Let (x',y') be a solution of the equality
    *   a*x + b*y == c    ->   a*x == c - b*y
    * Then x = -b*z + x', y = a*z + y' with z integral gives all solutions to the equality.
    */

   /* find initial solution (x',y'):
    *  - find y' such that c - b*y' is a multiple of a
    *    - start in equivalence class c%a
    *    - step through classes, where each step increases class number by (-b)%a, until class 0 is visited
    *    - if equivalence class 0 is visited, we are done: y' equals the number of steps taken
    *    - because a and b don't have a common divisor, each class is visited at most once, and at most a-1 steps are needed
    *  - calculate x' with x' = (c - b*y')/a (which must be integral)
    *
    * Algorithm works for a > 0 only.
    */
   if( a < 0 )
   {
      a = -a;
      b = -b;
      c = -c;
   }
   assert(0 <= a);

   /* search upwards from ysol = 0 */
   ysol = 0;
   currentclass = c%a;
   if( currentclass < 0 )
      currentclass += a;
   assert(0 <= currentclass && currentclass < a);

   classstep = (-b)%a;

   if( classstep < 0 )
      classstep += a;
   assert(0 <= classstep && classstep < a);

   while( currentclass != 0 )
   {
      assert(0 <= currentclass && currentclass < a);
      currentclass += classstep;
      if( currentclass >= a )
         currentclass -= a;
      ysol++;
   }
   assert(ysol < a);
   assert(((c - b*ysol)%a) == 0);

   xsol = (c - b*ysol)/a;

   /* determine variable type for new artificial variable:
    *
    * if both variables are implicit integer the new variable can be implicit too, because the integer implication on
    * these both variables should be enforced by some other variables, otherwise the new variable needs to be of
    * integral type
    */
   vartype = ((SCIPvarGetType(varx) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(vary) == SCIP_VARTYPE_INTEGER)
      ? SCIP_VARTYPE_INTEGER : SCIP_VARTYPE_IMPLINT);

   /* feasible solutions are (x,y) = (x',y') + z*(-b,a)
    * - create new integer variable z with infinite bounds
    * - aggregate variable x = -b*z + x'
    * - aggregate variable y =  a*z + y'
    * - the bounds of z are calculated automatically during aggregation
    */
   (void) SCIPsnprintf(aggvarname, SCIP_MAXSTRLEN, "agg%d", stat->nvaridx);
   SCIP_CALL( SCIPvarCreateTransformed(&aggvar, blkmem, set, stat,
         aggvarname, -SCIPsetInfinity(set), SCIPsetInfinity(set), 0.0, vartype,
         SCIPvarIsInitial(varx) || SCIPvarIsInitial(vary), SCIPvarIsRemovable(varx) && SCIPvarIsRemovable(vary),
         NULL, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPprobAddVar(transprob, blkmem, set, lp, branchcand, eventfilter, eventqueue, aggvar) );

   SCIP_CALL( SCIPvarAggregate(varx, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, cliquetable,
         branchcand, eventqueue, aggvar, (SCIP_Real)(-b), (SCIP_Real)xsol, infeasible, aggregated) );
   assert(*aggregated || *infeasible);

   if( !(*infeasible) )
   {
      SCIP_CALL( SCIPvarAggregate(vary, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, cliquetable,
            branchcand, eventqueue, aggvar, (SCIP_Real)a, (SCIP_Real)ysol, infeasible, aggregated) );
      assert(*aggregated || *infeasible);
   }

   /* release z */
   SCIP_CALL( SCIPvarRelease(&aggvar, blkmem, set, eventqueue, lp) );

   return SCIP_OKAY;
}

/** performs second step of SCIPaggregateVars():
 *  the variable to be aggregated is chosen among active problem variables x' and y', preferring a less strict variable
 *  type as aggregation variable (i.e. continuous variables are preferred over implicit integers, implicit integers
 *  or integers over binaries). If none of the variables is continuous, it is tried to find an integer
 *  aggregation (i.e. integral coefficients a'' and b'', such that a''*x' + b''*y' == c''). This can lead to
 *  the detection of infeasibility (e.g. if c'' is fractional), or to a rejection of the aggregation (denoted by
 *  aggregated == FALSE), if the resulting integer coefficients are too large and thus numerically instable.
 *
 *  @todo check for fixings, infeasibility, bound changes, or domain holes:
 *     a) if there is no easy aggregation and we have one binary variable and another integer/implicit/binary variable
 *     b) for implicit integer variables with fractional aggregation scalar (we cannot (for technical reasons) and do
 *        not want to aggregate implicit integer variables, since we loose the corresponding divisibility property)
 */
SCIP_RETCODE SCIPvarTryAggregateVars(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             varx,               /**< variable x in equality a*x + b*y == c */
   SCIP_VAR*             vary,               /**< variable y in equality a*x + b*y == c */
   SCIP_Real             scalarx,            /**< multiplier a in equality a*x + b*y == c */
   SCIP_Real             scalary,            /**< multiplier b in equality a*x + b*y == c */
   SCIP_Real             rhs,                /**< right hand side c in equality a*x + b*y == c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   )
{
   SCIP_Bool easyaggr;
   SCIP_Real maxscalar;
   SCIP_Real absquot;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(cliquetable != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(varx != NULL);
   assert(vary != NULL);
   assert(varx != vary);
   assert(infeasible != NULL);
   assert(aggregated != NULL);
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PRESOLVING);
   assert(SCIPvarGetStatus(varx) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetStatus(vary) == SCIP_VARSTATUS_LOOSE);
   assert(!SCIPsetIsZero(set, scalarx));
   assert(!SCIPsetIsZero(set, scalary));

   *infeasible = FALSE;
   *aggregated = FALSE;

   absquot = REALABS(scalarx / scalary);
   maxscalar = SCIPsetFeastol(set) / SCIPsetEpsilon(set);
   maxscalar = MAX(maxscalar, 1.0);

   if( absquot > maxscalar || absquot < 1 / maxscalar )
      return SCIP_OKAY;

   /* prefer aggregating the variable of more general type (preferred aggregation variable is varx) */
   if( SCIPvarGetType(vary) > SCIPvarGetType(varx) ||
         (SCIPvarGetType(vary) == SCIPvarGetType(varx) && !SCIPvarIsBinary(vary) && SCIPvarIsBinary(varx))  )
   {
      SCIP_VAR* var;
      SCIP_Real scalar;

      /* switch the variables, such that varx is the variable of more general type (cont > implint > int > bin) */
      var = vary;
      vary = varx;
      varx = var;
      scalar = scalary;
      scalary = scalarx;
      scalarx = scalar;
   }

   /* don't aggregate if the aggregation would lead to a binary variable aggregated to a non-binary variable */
   if( SCIPvarIsBinary(varx) && !SCIPvarIsBinary(vary) )
      return SCIP_OKAY;

   assert(SCIPvarGetType(varx) >= SCIPvarGetType(vary));

   /* figure out, which variable should be aggregated */
   easyaggr = FALSE;

   /* check if it is an easy aggregation that means:
    *
    *   a*x + b*y == c -> x == -b/a * y + c/a iff |b/a| > feastol and |a/b| > feastol
    */
   if( !SCIPsetIsFeasZero(set, scalary/scalarx) && !SCIPsetIsFeasZero(set, scalarx/scalary) )
   {
      if( SCIPvarGetType(varx) == SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(vary) < SCIP_VARTYPE_CONTINUOUS )
      {
         easyaggr = TRUE;
      }
      else if( SCIPsetIsFeasIntegral(set, scalary/scalarx) )
      {
         easyaggr = TRUE;
      }
      else if( SCIPsetIsFeasIntegral(set, scalarx/scalary) && SCIPvarGetType(vary) == SCIPvarGetType(varx) )
      {
         /* we have an easy aggregation if we flip the variables x and y */
         SCIP_VAR* var;
         SCIP_Real scalar;

         /* switch the variables, such that varx is the aggregated variable */
         var = vary;
         vary = varx;
         varx = var;
         scalar = scalary;
         scalary = scalarx;
         scalarx = scalar;
         easyaggr = TRUE;
      }
      else if( SCIPvarGetType(varx) == SCIP_VARTYPE_CONTINUOUS )
      {
         /* the aggregation is still easy if both variables are continuous */
         assert(SCIPvarGetType(vary) == SCIP_VARTYPE_CONTINUOUS); /* otherwise we are in the first case */
         easyaggr = TRUE;
      }
   }

   /* did we find an "easy" aggregation? */
   if( easyaggr )
   {
      SCIP_Real scalar;
      SCIP_Real constant;

      assert(SCIPvarGetType(varx) >= SCIPvarGetType(vary));

      /* calculate aggregation scalar and constant: a*x + b*y == c  =>  x == -b/a * y + c/a */
      scalar = -scalary/scalarx;
      constant = rhs/scalarx;

      if( REALABS(constant) > SCIPsetGetHugeValue(set) * SCIPsetFeastol(set) ) /*lint !e653*/
         return SCIP_OKAY;

      /* check aggregation for integer feasibility */
      if( SCIPvarGetType(varx) != SCIP_VARTYPE_CONTINUOUS
         && SCIPvarGetType(vary) != SCIP_VARTYPE_CONTINUOUS
         && SCIPsetIsFeasIntegral(set, scalar) && !SCIPsetIsFeasIntegral(set, constant) )
      {
         *infeasible = TRUE;
         return SCIP_OKAY;
      }

      /* if the aggregation scalar is fractional, we cannot (for technical reasons) and do not want to aggregate implicit integer variables,
       * since then we would loose the corresponding divisibility property
       */
      assert(SCIPvarGetType(varx) != SCIP_VARTYPE_IMPLINT || SCIPsetIsFeasIntegral(set, scalar));

      /* aggregate the variable */
      SCIP_CALL( SCIPvarAggregate(varx, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, cliquetable,
            branchcand, eventqueue, vary, scalar, constant, infeasible, aggregated) );
      assert(*aggregated || *infeasible);
   }
   else if( (SCIPvarGetType(varx) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(varx) == SCIP_VARTYPE_IMPLINT)
      && (SCIPvarGetType(vary) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(vary) == SCIP_VARTYPE_IMPLINT) )
   {
      /* the variables are both integral: we have to try to find an integer aggregation */
      SCIP_CALL( tryAggregateIntVars(set, blkmem, stat, transprob, origprob, primal, tree, reopt, lp, cliquetable,
            branchcand, eventfilter, eventqueue, varx, vary, scalarx, scalary, rhs, infeasible, aggregated) );
   }

   return SCIP_OKAY;
}

/** converts variable into multi-aggregated variable */
SCIP_RETCODE SCIPvarMultiaggregate(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   naggvars,           /**< number n of variables in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_VAR**            aggvars,            /**< variables y_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real*            scalars,            /**< multipliers a_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real             constant,           /**< constant shift c in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   )
{
   SCIP_VAR** tmpvars;
   SCIP_Real* tmpscalars;
   SCIP_Real obj;
   SCIP_Real branchfactor;
   int branchpriority;
   SCIP_BRANCHDIR branchdirection;
   int nlocksdown;
   int nlocksup;
   int v;
   SCIP_Real tmpconstant;
   SCIP_Real tmpscalar;
   int ntmpvars;
   int tmpvarssize;
   int tmprequiredsize;

   assert(var != NULL);
   assert(var->scip == set->scip);
   assert(var->glbdom.lb == var->locdom.lb); /*lint !e777*/
   assert(var->glbdom.ub == var->locdom.ub); /*lint !e777*/
   assert(naggvars == 0 || aggvars != NULL);
   assert(naggvars == 0 || scalars != NULL);
   assert(infeasible != NULL);
   assert(aggregated != NULL);

   SCIPsetDebugMsg(set, "trying multi-aggregating variable <%s> == ...%d vars... %+g\n", var->name, naggvars, constant);

   *infeasible = FALSE;
   *aggregated = FALSE;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
      {
         SCIPerrorMessage("cannot multi-aggregate an untransformed original variable\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarMultiaggregate(var->data.original.transvar, blkmem, set, stat, transprob, origprob, primal, tree,
            reopt, lp, cliquetable, branchcand, eventfilter, eventqueue, naggvars, aggvars, scalars, constant, infeasible, aggregated) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      assert(!SCIPeventqueueIsDelayed(eventqueue)); /* otherwise, the pseudo objective value update gets confused */

      /* check if we would create a self-reference */
      ntmpvars = naggvars;
      tmpvarssize = naggvars;
      tmpconstant = constant;
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &tmpvars, aggvars, ntmpvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &tmpscalars, scalars, ntmpvars) );

      /* get all active variables for multi-aggregation */
      SCIP_CALL( SCIPvarGetActiveRepresentatives(set, tmpvars, tmpscalars, &ntmpvars, tmpvarssize, &tmpconstant, &tmprequiredsize, FALSE) );
      if( tmprequiredsize > tmpvarssize )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &tmpvars, tmpvarssize, tmprequiredsize) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &tmpscalars, tmpvarssize, tmprequiredsize) );
         tmpvarssize = tmprequiredsize;
         SCIP_CALL( SCIPvarGetActiveRepresentatives(set, tmpvars, tmpscalars, &ntmpvars, tmpvarssize, &tmpconstant, &tmprequiredsize, FALSE) );
         assert( tmprequiredsize <= tmpvarssize );
      }

      tmpscalar = 0.0;

      /* iterate over all active variables of the multi-aggregation and filter all variables which are equal to the
       * possible multi-aggregated variable
       */
      for( v = ntmpvars - 1; v >= 0; --v )
      {
         assert(tmpvars[v] != NULL);
         assert(SCIPvarGetStatus(tmpvars[v]) == SCIP_VARSTATUS_LOOSE);

         if( tmpvars[v]->index == var->index )
         {
            tmpscalar += tmpscalars[v];
            tmpvars[v] = tmpvars[ntmpvars - 1];
            tmpscalars[v] = tmpscalars[ntmpvars - 1];
            --ntmpvars;
         }
      }

      /* this means that x = x + a_1*y_1 + ... + a_n*y_n + c */
      if( SCIPsetIsEQ(set, tmpscalar, 1.0) )
      {
         if( ntmpvars == 0 )
         {
            if( SCIPsetIsZero(set, tmpconstant) ) /* x = x */
            {
               SCIPsetDebugMsg(set, "Possible multi-aggregation was completely resolved and detected to be redundant.\n");
               goto TERMINATE;
            }
            else /* 0 = c and c != 0 */
            {
               SCIPsetDebugMsg(set, "Multi-aggregation was completely resolved and led to infeasibility.\n");
               *infeasible = TRUE;
               goto TERMINATE;
            }
         }
         else if( ntmpvars == 1 ) /* 0 = a*y + c => y = -c/a */
         {
            assert(tmpscalars[0] != 0.0);
            assert(tmpvars[0] != NULL);

            SCIPsetDebugMsg(set, "Possible multi-aggregation led to fixing of variable <%s> to %g.\n", SCIPvarGetName(tmpvars[0]), -constant/tmpscalars[0]);
            SCIP_CALL( SCIPvarFix(tmpvars[0], blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp,
                  branchcand, eventqueue, cliquetable, -constant/tmpscalars[0], infeasible, aggregated) );
            goto TERMINATE;
         }
         else if( ntmpvars == 2 ) /* 0 = a_1*y_1 + a_2*y_2 + c => y_1 = -a_2/a_1 * y_2 - c/a_1 */
         {
	    /* both variables are different active problem variables, and both scalars are non-zero: try to aggregate them */

            SCIPsetDebugMsg(set, "Possible multi-aggregation led to aggregation of variables <%s> and <%s> with scalars %g and %g and constant %g.\n",
                  SCIPvarGetName(tmpvars[0]), SCIPvarGetName(tmpvars[1]), tmpscalars[0], tmpscalars[1], -tmpconstant);

            SCIP_CALL( SCIPvarTryAggregateVars(set, blkmem, stat, transprob, origprob, primal, tree, reopt, lp,
                  cliquetable, branchcand, eventfilter, eventqueue, tmpvars[0], tmpvars[1], tmpscalars[0],
                  tmpscalars[1], -tmpconstant, infeasible, aggregated) );

            goto TERMINATE;
         }
         else
            /* @todo: it is possible to multi-aggregate another variable, does it make sense?,
             *        rest looks like 0 = a_1*y_1 + ... + a_n*y_n + c and has at least three variables
             */
            goto TERMINATE;
      }
      /* this means that x = b*x + a_1*y_1 + ... + a_n*y_n + c */
      else if( !SCIPsetIsZero(set, tmpscalar) )
      {
         tmpscalar = 1 - tmpscalar;
         tmpconstant /= tmpscalar;
         for( v = ntmpvars - 1; v >= 0; --v )
            tmpscalars[v] /= tmpscalar;
      }

      /* check, if we are in one of the simple cases */
      if( ntmpvars == 0 )
      {
         SCIPsetDebugMsg(set, "Possible multi-aggregation led to fixing of variable <%s> to %g.\n", SCIPvarGetName(var), tmpconstant);
         SCIP_CALL( SCIPvarFix(var, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand,
               eventqueue, cliquetable, tmpconstant, infeasible, aggregated) );
         goto TERMINATE;
      }

      /* if only one aggregation variable is left, we perform a normal aggregation instead of a multi-aggregation */
      if( ntmpvars == 1 )
      {
            SCIPsetDebugMsg(set, "Possible multi-aggregation led to aggregation of variables <%s> and <%s> with scalars %g and %g and constant %g.\n",
                  SCIPvarGetName(var), SCIPvarGetName(tmpvars[0]), 1.0, -tmpscalars[0], tmpconstant);

            SCIP_CALL( SCIPvarTryAggregateVars(set, blkmem, stat, transprob, origprob, primal, tree, reopt, lp,
                  cliquetable, branchcand, eventfilter, eventqueue, var, tmpvars[0], 1.0, -tmpscalars[0], tmpconstant,
                  infeasible, aggregated) );

         goto TERMINATE;
      }

      /**@todo currently we don't perform the multi aggregation if the multi aggregation variable has a non
       *  empty hole list; this should be changed in the future  */
      if( SCIPvarGetHolelistGlobal(var) != NULL )
         goto TERMINATE;

      /* if the variable is not allowed to be multi-aggregated */
      if( SCIPvarDoNotMultaggr(var) )
      {
         SCIPsetDebugMsg(set, "variable is not allowed to be multi-aggregated.\n");
         goto TERMINATE;
      }

      /* if the variable to be multi-aggregated has implications or variable bounds (i.e. is the implied variable or
       * variable bound variable of another variable), we have to remove it from the other variables implications or
       * variable bounds
       */
      SCIP_CALL( SCIPvarRemoveCliquesImplicsVbs(var, blkmem, cliquetable, set, FALSE, FALSE, TRUE) );
      assert(var->vlbs == NULL);
      assert(var->vubs == NULL);
      assert(var->implics == NULL);
      assert(var->cliquelist == NULL);

      /* set the aggregated variable's objective value to 0.0 */
      obj = var->obj;
      SCIP_CALL( SCIPvarChgObj(var, blkmem, set, transprob, primal, lp, eventqueue, 0.0) );

      /* since we change the variable type form loose to multi aggregated, we have to adjust the number of loose
       * variables in the LP data structure; the loose objective value (looseobjval) in the LP data structure, however,
       * gets adjusted automatically, due to the event SCIP_EVENTTYPE_OBJCHANGED which dropped in the moment where the
       * objective of this variable is set to zero
       */
      SCIPlpDecNLoosevars(lp);

      /* unlock all rounding locks */
      nlocksdown = var->nlocksdown;
      nlocksup = var->nlocksup;
      var->nlocksdown = 0;
      var->nlocksup = 0;

      /* convert variable into multi-aggregated variable */
      var->varstatus = SCIP_VARSTATUS_MULTAGGR; /*lint !e641*/
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &var->data.multaggr.vars, tmpvars, ntmpvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &var->data.multaggr.scalars, tmpscalars, ntmpvars) );
      var->data.multaggr.constant = tmpconstant;
      var->data.multaggr.nvars = ntmpvars;
      var->data.multaggr.varssize = ntmpvars;

      /* mark variable to be non-deletable */
      SCIPvarMarkNotDeletable(var);

      /* relock the rounding locks of the variable, thus increasing the locks of the aggregation variables */
      SCIP_CALL( SCIPvarAddLocks(var, blkmem, set, eventqueue, nlocksdown, nlocksup) );

      /* update flags and branching factors and priorities of aggregation variables;
       * update preferred branching direction of all aggregation variables that don't have a preferred direction yet
       */
      branchfactor = var->branchfactor;
      branchpriority = var->branchpriority;
      branchdirection = (SCIP_BRANCHDIR)var->branchdirection;

      for( v = 0; v < ntmpvars; ++v )
      {
         assert(tmpvars[v] != NULL);
         tmpvars[v]->removable &= var->removable;
         branchfactor = MAX(tmpvars[v]->branchfactor, branchfactor);
         branchpriority = MAX(tmpvars[v]->branchpriority, branchpriority);

         /* mark variable to be non-deletable */
         SCIPvarMarkNotDeletable(tmpvars[v]);
      }
      for( v = 0; v < ntmpvars; ++v )
      {
         SCIP_CALL( SCIPvarChgBranchFactor(tmpvars[v], set, branchfactor) );
         SCIP_CALL( SCIPvarChgBranchPriority(tmpvars[v], branchpriority) );
         if( (SCIP_BRANCHDIR)tmpvars[v]->branchdirection == SCIP_BRANCHDIR_AUTO )
         {
            if( tmpscalars[v] >= 0.0 )
            {
               SCIP_CALL( SCIPvarChgBranchDirection(tmpvars[v], branchdirection) );
            }
            else
            {
               SCIP_CALL( SCIPvarChgBranchDirection(tmpvars[v], SCIPbranchdirOpposite(branchdirection)) );
            }
         }
      }
      SCIP_CALL( SCIPvarChgBranchFactor(var, set, branchfactor) );
      SCIP_CALL( SCIPvarChgBranchPriority(var, branchpriority) );

      if( var->probindex != -1 )
      {
         /* inform problem about the variable's status change */
         SCIP_CALL( SCIPprobVarChangedStatus(transprob, blkmem, set, branchcand, cliquetable, var) );
      }

      /* issue VARFIXED event */
      SCIP_CALL( varEventVarFixed(var, blkmem, set, eventqueue, 2) );

      /* reset the objective value of the aggregated variable, thus adjusting the objective value of the aggregation
       * variables and the problem's objective offset
       */
      SCIP_CALL( SCIPvarAddObj(var, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, eventqueue, obj) );

      *aggregated = TRUE;

   TERMINATE:      
      BMSfreeBlockMemoryArray(blkmem, &tmpscalars, tmpvarssize);
      BMSfreeBlockMemoryArray(blkmem, &tmpvars, tmpvarssize);

      break;

   case SCIP_VARSTATUS_COLUMN:
      SCIPerrorMessage("cannot multi-aggregate a column variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot multi-aggregate a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      SCIPerrorMessage("cannot multi-aggregate an aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot multi-aggregate a multiple aggregated variable again\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      /* aggregate negation variable x in x' = offset - x, instead of aggregating x' directly:
       *   x' = a_1*y_1 + ... + a_n*y_n + c  ->  x = offset - x' = offset - a_1*y_1 - ... - a_n*y_n - c
       */
      assert(SCIPsetIsZero(set, var->obj));
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);

      /* switch the signs of the aggregation scalars */
      for( v = 0; v < naggvars; ++v )
         scalars[v] *= -1.0;

      /* perform the multi aggregation on the negation variable */
      SCIP_CALL( SCIPvarMultiaggregate(var->negatedvar, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp,
            cliquetable, branchcand, eventfilter, eventqueue, naggvars, aggvars, scalars,
	    var->data.negate.constant - constant, infeasible, aggregated) );

      /* switch the signs of the aggregation scalars again, to reset them to their original values */
      for( v = 0; v < naggvars; ++v )
         scalars[v] *= -1.0;
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** transformed variables are resolved to their active, fixed, or multi-aggregated problem variable of a variable,
 * or for original variables the same variable is returned
 */
static
SCIP_VAR* varGetActiveVar(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_VAR* retvar;

   assert(var != NULL);

   retvar = var;

   SCIPdebugMessage("get active variable of <%s>\n", var->name);

   while( TRUE ) /*lint !e716 */
   {
      assert(retvar != NULL);

      switch( SCIPvarGetStatus(retvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_FIXED:
	 return retvar;

      case SCIP_VARSTATUS_MULTAGGR:
	 /* handle multi-aggregated variables depending on one variable only (possibly caused by SCIPvarFlattenAggregationGraph()) */
	 if ( retvar->data.multaggr.nvars == 1 )
	    retvar = retvar->data.multaggr.vars[0];
	 else
	    return retvar;
	 break;

      case SCIP_VARSTATUS_AGGREGATED:
	 retvar = retvar->data.aggregate.var;
	 break;

      case SCIP_VARSTATUS_NEGATED:
	 retvar = retvar->negatedvar;
	 break;

      default:
	 SCIPerrorMessage("unknown variable status\n");
	 SCIPABORT();
	 return NULL; /*lint !e527*/
      }
   }
}

/** returns whether variable is not allowed to be multi-aggregated */
SCIP_Bool SCIPvarDoNotMultaggr(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_VAR* retvar;

   assert(var != NULL);

   retvar = varGetActiveVar(var);
   assert(retvar != NULL);

   switch( SCIPvarGetStatus(retvar) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      return retvar->donotmultaggr;

   case SCIP_VARSTATUS_MULTAGGR:
      return FALSE;

   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_NEGATED:
   default:
      SCIPerrorMessage("wrong variable status\n");
      SCIPABORT();
      return FALSE; /*lint !e527 */
   }
}

/** gets negated variable x' = offset - x of problem variable x; the negated variable is created if not yet existing;
 *  the negation offset of binary variables is always 1, the offset of other variables is fixed to lb + ub when the
 *  negated variable is created
 */
SCIP_RETCODE SCIPvarNegate(
   SCIP_VAR*             var,                /**< problem variable to negate */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR**            negvar              /**< pointer to store the negated variable */
   )
{
   assert(var != NULL);
   assert(var->scip == set->scip);
   assert(negvar != NULL);

   /* check, if we already created the negated variable */
   if( var->negatedvar == NULL )
   {
      char negvarname[SCIP_MAXSTRLEN];

      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_NEGATED);

      SCIPsetDebugMsg(set, "creating negated variable of <%s>\n", var->name);

      /* negation is only possible for bounded variables */
      if( SCIPsetIsInfinity(set, -var->glbdom.lb) || SCIPsetIsInfinity(set, var->glbdom.ub) )
      {
         SCIPerrorMessage("cannot negate unbounded variable\n");
         return SCIP_INVALIDDATA;
      }

      (void) SCIPsnprintf(negvarname, SCIP_MAXSTRLEN, "%s_neg", var->name);

      /* create negated variable */
      SCIP_CALL( varCreate(negvar, blkmem, set, stat, negvarname, var->glbdom.lb, var->glbdom.ub, 0.0,
            SCIPvarGetType(var), var->initial, var->removable, NULL, NULL, NULL, NULL, NULL) );
      (*negvar)->varstatus = SCIP_VARSTATUS_NEGATED; /*lint !e641*/
      if( SCIPvarIsBinary(var) )
         (*negvar)->data.negate.constant = 1.0;
      else
         (*negvar)->data.negate.constant = var->glbdom.lb + var->glbdom.ub;

      /* create event filter for transformed variable */
      if( SCIPvarIsTransformed(var) )
      {
         SCIP_CALL( SCIPeventfilterCreate(&(*negvar)->eventfilter, blkmem) );
      }

      /* set the bounds corresponding to the negation variable */
      (*negvar)->glbdom.lb = (*negvar)->data.negate.constant - var->glbdom.ub;
      (*negvar)->glbdom.ub = (*negvar)->data.negate.constant - var->glbdom.lb;
      (*negvar)->locdom.lb = (*negvar)->data.negate.constant - var->locdom.ub;
      (*negvar)->locdom.ub = (*negvar)->data.negate.constant - var->locdom.lb;
      /**@todo create holes in the negated variable corresponding to the holes of the negation variable */

      /* link the variables together */
      var->negatedvar = *negvar;
      (*negvar)->negatedvar = var;

      /* mark both variables to be non-deletable */
      SCIPvarMarkNotDeletable(var);
      SCIPvarMarkNotDeletable(*negvar);

      /* copy the branch factor and priority, and use the negative preferred branching direction */
      (*negvar)->branchfactor = var->branchfactor;
      (*negvar)->branchpriority = var->branchpriority;
      (*negvar)->branchdirection = SCIPbranchdirOpposite((SCIP_BRANCHDIR)var->branchdirection); /*lint !e641*/

      /* copy doNotMultiaggr status */
      (*negvar)->donotmultaggr = var->donotmultaggr;

      /* copy lazy bounds (they have to be flipped) */
      (*negvar)->lazylb = (*negvar)->data.negate.constant - var->lazyub;
      (*negvar)->lazyub = (*negvar)->data.negate.constant - var->lazylb;

      /* make negated variable a parent of the negation variable (negated variable is captured as a parent) */
      SCIP_CALL( varAddParent(var, blkmem, set, *negvar) );
      assert((*negvar)->nuses == 1);
   }
   assert(var->negatedvar != NULL);

   /* return the negated variable */
   *negvar = var->negatedvar;

   /* exactly one variable of the negation pair has to be marked as negated variable */
   assert((SCIPvarGetStatus(*negvar) == SCIP_VARSTATUS_NEGATED) != (SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED));

   return SCIP_OKAY;
}

/** informs variable that its position in problem's vars array changed */
static
void varSetProbindex(
   SCIP_VAR*             var,                /**< problem variable */
   int                   probindex           /**< new problem index of variable (-1 for removal) */
   )
{
   assert(var != NULL);
   assert(probindex >= 0 || var->vlbs == NULL);
   assert(probindex >= 0 || var->vubs == NULL);
   assert(probindex >= 0 || var->implics == NULL);

   var->probindex = probindex;
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
   {
      assert(var->data.col != NULL);
      var->data.col->var_probindex = probindex;
   }
}

/** informs variable that its position in problem's vars array changed */
void SCIPvarSetProbindex(
   SCIP_VAR*             var,                /**< problem variable */
   int                   probindex           /**< new problem index of variable */
   )
{
   assert(var != NULL);
   assert(probindex >= 0);

   varSetProbindex(var, probindex);
}

/** gives the variable a new name
 *
 *  @note the old pointer is overwritten, which might result in a memory leakage
 */
void SCIPvarSetNamePointer(
   SCIP_VAR*             var,                /**< problem variable */
   const char*           name                /**< new name of variable */
   )
{
   assert(var != NULL);
   assert(name != NULL);

   var->name = (char*)name;
}

/** informs variable that it will be removed from the problem; adjusts probindex and removes variable from the
 *  implication graph;
 *  If 'final' is TRUE, the thorough implication graph removal is not performed. Instead, only the
 *  variable bounds and implication data structures of the variable are freed. Since in the final removal
 *  of all variables from the transformed problem, this deletes the implication graph completely and is faster
 *  than removing the variables one by one, each time updating all lists of the other variables.
 */
SCIP_RETCODE SCIPvarRemove(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             final               /**< is this the final removal of all problem variables? */
   )
{
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(var->scip == set->scip);

   /* if the variable is active in the transformed problem, remove it from the implication graph */
   if( SCIPvarIsTransformed(var)
      && (SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN) )
   {
      if( final )
      {
         /* just destroy the data structures */
         SCIPvboundsFree(&var->vlbs, blkmem);
         SCIPvboundsFree(&var->vubs, blkmem);
         SCIPimplicsFree(&var->implics, blkmem);
      }
      else
      {
         /* unlink the variable from all other variables' lists and free the data structures */
         SCIP_CALL( SCIPvarRemoveCliquesImplicsVbs(var, blkmem, cliquetable, set, FALSE, FALSE, TRUE) );
      }
   }

   /* mark the variable to be no longer a member of the problem */
   varSetProbindex(var, -1);

   return SCIP_OKAY;
}

/** marks the variable to be deleted from the problem */
void SCIPvarMarkDeleted(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->probindex != -1);

   var->deleted = TRUE;
}

/** marks the variable to not to be multi-aggregated */
SCIP_RETCODE SCIPvarMarkDoNotMultaggr(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_VAR* retvar;

   assert(var != NULL);

   retvar = varGetActiveVar(var);
   assert(retvar != NULL);

   switch( SCIPvarGetStatus(retvar) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      retvar->donotmultaggr = TRUE;
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot mark a multi-aggregated variable to not be multi-aggregated.\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_NEGATED:
   default:
      SCIPerrorMessage("wrong variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** changes type of variable; cannot be called, if var belongs to a problem */
SCIP_RETCODE SCIPvarChgType(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_VARTYPE          vartype             /**< new type of variable */
   )
{
   assert(var != NULL);

   SCIPdebugMessage("change type of <%s> from %d to %d\n", var->name, SCIPvarGetType(var), vartype);

   if( var->probindex >= 0 )
   {
      SCIPerrorMessage("cannot change type of variable already in the problem\n");
      return SCIP_INVALIDDATA;
   }

   var->vartype = vartype; /*lint !e641*/
   if( var->negatedvar != NULL )
      var->negatedvar->vartype = vartype; /*lint !e641*/

   return SCIP_OKAY;
}

/** appends OBJCHANGED event to the event queue */
static
SCIP_RETCODE varEventObjChanged(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             oldobj,             /**< old objective value for variable */
   SCIP_Real             newobj              /**< new objective value for variable */
   )
{
   SCIP_EVENT* event;

   assert(var != NULL);
   assert(var->scip == set->scip);
   assert(var->eventfilter != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarIsTransformed(var));
   assert(!SCIPsetIsEQ(set, oldobj, newobj));

   SCIP_CALL( SCIPeventCreateObjChanged(&event, blkmem, var, oldobj, newobj) );
   SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, primal, lp, NULL, NULL, &event) );

   return SCIP_OKAY;
}

/** changes objective value of variable */
SCIP_RETCODE SCIPvarChgObj(
   SCIP_VAR*             var,                /**< variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             newobj              /**< new objective value for variable */
   )
{
   SCIP_Real oldobj;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   SCIPsetDebugMsg(set, "changing objective value of <%s> from %g to %g\n", var->name, var->obj, newobj);

   if( !SCIPsetIsEQ(set, var->obj, newobj) )
   {
      switch( SCIPvarGetStatus(var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.original.transvar != NULL )
         {
            assert(SCIPprobIsTransformed(prob));

            SCIP_CALL( SCIPvarChgObj(var->data.original.transvar, blkmem, set, prob, primal, lp, eventqueue,
                  (SCIP_Real) prob->objsense * newobj/prob->objscale) );
         }
         else
            assert(set->stage == SCIP_STAGE_PROBLEM);

         var->obj = newobj;
         var->unchangedobj = newobj;

         break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
         oldobj = var->obj;
         var->obj = newobj;

         /* update unchanged objective value of variable */
         if( !lp->divingobjchg )
            var->unchangedobj = newobj;

         /* update the number of variables with non-zero objective coefficient;
          * we only want to do the update, if the variable is added to the problem;
          * since the objective of inactive variables cannot be changed, this corresponds to probindex != -1
          */
         if( SCIPvarIsActive(var) )
            SCIPprobUpdateNObjVars(prob, set, oldobj, var->obj);

         SCIP_CALL( varEventObjChanged(var, blkmem, set, primal, lp, eventqueue, oldobj, var->obj) );
         break;

      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_AGGREGATED:
      case SCIP_VARSTATUS_MULTAGGR:
      case SCIP_VARSTATUS_NEGATED:
         SCIPerrorMessage("cannot change objective value of a fixed, aggregated, multi-aggregated, or negated variable\n");
         return SCIP_INVALIDDATA;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** adds value to objective value of variable */
SCIP_RETCODE SCIPvarAddObj(
   SCIP_VAR*             var,                /**< variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             addobj              /**< additional objective value for variable */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(set->stage < SCIP_STAGE_INITSOLVE);

   SCIPsetDebugMsg(set, "adding %g to objective value %g of <%s>\n", addobj, var->obj, var->name);

   if( !SCIPsetIsZero(set, addobj) )
   {
      SCIP_Real oldobj;
      int i;

      switch( SCIPvarGetStatus(var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.original.transvar != NULL )
         {
            SCIP_CALL( SCIPvarAddObj(var->data.original.transvar, blkmem, set, stat, transprob, origprob, primal, tree,
                  reopt, lp, eventqueue, (SCIP_Real) transprob->objsense * addobj/transprob->objscale) );
         }
         else
            assert(set->stage == SCIP_STAGE_PROBLEM);

         var->obj += addobj;
         var->unchangedobj += addobj;
         assert(SCIPsetIsEQ(set, var->obj, var->unchangedobj));

         break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
         oldobj = var->obj;
         var->obj += addobj;

         /* update unchanged objective value of variable */
         if( !lp->divingobjchg )
         {
            var->unchangedobj += addobj;
            assert(SCIPsetIsEQ(set, var->obj, var->unchangedobj));
         }

         /* update the number of variables with non-zero objective coefficient;
          * we only want to do the update, if the variable is added to the problem;
          * since the objective of inactive variables cannot be changed, this corresponds to probindex != -1
          */
         if( SCIPvarIsActive(var) )
	    SCIPprobUpdateNObjVars(transprob, set, oldobj, var->obj);

         SCIP_CALL( varEventObjChanged(var, blkmem, set, primal, lp, eventqueue, oldobj, var->obj) );
         break;

      case SCIP_VARSTATUS_FIXED:
         assert(SCIPsetIsEQ(set, var->locdom.lb, var->locdom.ub));
         SCIPprobAddObjoffset(transprob, var->locdom.lb * addobj);
         SCIP_CALL( SCIPprimalUpdateObjoffset(primal, blkmem, set, stat, eventqueue, transprob, origprob, tree, reopt, lp) );
         break;

      case SCIP_VARSTATUS_AGGREGATED:
         /* x = a*y + c  ->  add a*addobj to obj. val. of y, and c*addobj to obj. offset of problem */
         SCIPprobAddObjoffset(transprob, var->data.aggregate.constant * addobj);
         SCIP_CALL( SCIPprimalUpdateObjoffset(primal, blkmem, set, stat, eventqueue, transprob, origprob, tree, reopt, lp) );
         SCIP_CALL( SCIPvarAddObj(var->data.aggregate.var, blkmem, set, stat, transprob, origprob, primal, tree, reopt,
               lp, eventqueue, var->data.aggregate.scalar * addobj) );
         break;

      case SCIP_VARSTATUS_MULTAGGR:
         assert(!var->donotmultaggr);
         /* x = a_1*y_1 + ... + a_n*y_n  + c  ->  add a_i*addobj to obj. val. of y_i, and c*addobj to obj. offset */
         SCIPprobAddObjoffset(transprob, var->data.multaggr.constant * addobj);
         SCIP_CALL( SCIPprimalUpdateObjoffset(primal, blkmem, set, stat, eventqueue, transprob, origprob, tree, reopt, lp) );
         for( i = 0; i < var->data.multaggr.nvars; ++i )
         {
            SCIP_CALL( SCIPvarAddObj(var->data.multaggr.vars[i], blkmem, set, stat, transprob, origprob, primal, tree,
                  reopt, lp, eventqueue, var->data.multaggr.scalars[i] * addobj) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED:
         /* x' = offset - x  ->  add -addobj to obj. val. of x and offset*addobj to obj. offset of problem */
         assert(var->negatedvar != NULL);
         assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(var->negatedvar->negatedvar == var);
         SCIPprobAddObjoffset(transprob, var->data.negate.constant * addobj);
         SCIP_CALL( SCIPprimalUpdateObjoffset(primal, blkmem, set, stat, eventqueue, transprob, origprob, tree, reopt, lp) );
         SCIP_CALL( SCIPvarAddObj(var->negatedvar, blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp,
               eventqueue, -addobj) );
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** changes objective value of variable in current dive */
SCIP_RETCODE SCIPvarChgObjDive(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newobj              /**< new objective value for variable */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(lp != NULL);

   SCIPsetDebugMsg(set, "changing objective of <%s> to %g in current dive\n", var->name, newobj);

   if( SCIPsetIsZero(set, newobj) )
      newobj = 0.0;

   /* change objective value of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      SCIP_CALL( SCIPvarChgObjDive(var->data.original.transvar, set, lp, newobj) );
      break;

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      SCIP_CALL( SCIPcolChgObj(var->data.col, set, lp, newobj) );
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      /* nothing to do here: only the constant shift in objective function would change */
      break;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      SCIP_CALL( SCIPvarChgObjDive(var->data.aggregate.var, set, lp, newobj / var->data.aggregate.scalar) );
      /* the constant can be ignored, because it would only affect the objective shift */
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot change diving objective value of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarChgObjDive(var->negatedvar, set, lp, -newobj) );
      /* the offset can be ignored, because it would only affect the objective shift */
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** adjust lower bound to integral value, if variable is integral */
void SCIPvarAdjustLb(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            lb                  /**< pointer to lower bound to adjust */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(lb != NULL);

   SCIPsetDebugMsg(set, "adjust lower bound %g of <%s>\n", *lb, var->name);

   *lb = adjustedLb(set, SCIPvarGetType(var), *lb);
}

/** adjust upper bound to integral value, if variable is integral */
void SCIPvarAdjustUb(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            ub                  /**< pointer to upper bound to adjust */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(ub != NULL);

   SCIPsetDebugMsg(set, "adjust upper bound %g of <%s>\n", *ub, var->name);

   *ub = adjustedUb(set, SCIPvarGetType(var), *ub);
}

/** adjust lower or upper bound to integral value, if variable is integral */
void SCIPvarAdjustBd(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound to adjust */
   SCIP_Real*            bd                  /**< pointer to bound to adjust */
   )
{
   assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

   if( boundtype == SCIP_BOUNDTYPE_LOWER )
      SCIPvarAdjustLb(var, set, bd);
   else
      SCIPvarAdjustUb(var, set, bd);
}

/** changes lower bound of original variable in original problem */
SCIP_RETCODE SCIPvarChgLbOriginal(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   int i;

   assert(var != NULL);
   assert(!SCIPvarIsTransformed(var));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(set->stage == SCIP_STAGE_PROBLEM);

   /* check that the bound is feasible */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || SCIPsetIsLE(set, newbound, SCIPvarGetUbOriginal(var)));
   /* adjust bound to integral value if variable is of integral type */
   newbound = adjustedLb(set, SCIPvarGetType(var), newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   /* original domains are only stored for ORIGINAL variables, not for NEGATED */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
   {
      SCIPsetDebugMsg(set, "changing original lower bound of <%s> from %g to %g\n",
         var->name, var->data.original.origdom.lb, newbound);

      if( SCIPsetIsEQ(set, var->data.original.origdom.lb, newbound) )
         return SCIP_OKAY;

      /* change the bound */
      var->data.original.origdom.lb = newbound;
   }
   else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
   {
      assert( var->negatedvar != NULL );
      SCIP_CALL( SCIPvarChgUbOriginal(var->negatedvar, set, var->data.negate.constant - newbound) );
   }

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      SCIP_VAR* parentvar;

      parentvar = var->parentvars[i];
      assert(parentvar != NULL);
      assert(SCIPvarGetStatus(parentvar) == SCIP_VARSTATUS_NEGATED);
      assert(parentvar->negatedvar == var);
      assert(var->negatedvar == parentvar);

      SCIP_CALL( SCIPvarChgUbOriginal(parentvar, set, parentvar->data.negate.constant - newbound) );
   }

   return SCIP_OKAY;
}

/** changes upper bound of original variable in original problem */
SCIP_RETCODE SCIPvarChgUbOriginal(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   int i;

   assert(var != NULL);
   assert(!SCIPvarIsTransformed(var));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(set->stage == SCIP_STAGE_PROBLEM);

   /* check that the bound is feasible */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || SCIPsetIsGE(set, newbound, SCIPvarGetLbOriginal(var)));
   /* adjust bound to integral value if variable is of integral type */
   newbound = adjustedUb(set, SCIPvarGetType(var), newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   /* original domains are only stored for ORIGINAL variables, not for NEGATED */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
   {
      SCIPsetDebugMsg(set, "changing original upper bound of <%s> from %g to %g\n",
         var->name, var->data.original.origdom.ub, newbound);

      if( SCIPsetIsEQ(set, var->data.original.origdom.ub, newbound) )
         return SCIP_OKAY;

      /* change the bound */
      var->data.original.origdom.ub = newbound;
   }
   else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
   {
      assert( var->negatedvar != NULL );
      SCIP_CALL( SCIPvarChgLbOriginal(var->negatedvar, set, var->data.negate.constant - newbound) );
   }

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      SCIP_VAR* parentvar;

      parentvar = var->parentvars[i];
      assert(parentvar != NULL);
      assert(SCIPvarGetStatus(parentvar) == SCIP_VARSTATUS_NEGATED);
      assert(parentvar->negatedvar == var);
      assert(var->negatedvar == parentvar);

      SCIP_CALL( SCIPvarChgLbOriginal(parentvar, set, parentvar->data.negate.constant - newbound) );
   }

   return SCIP_OKAY;
}

/** appends GLBCHANGED event to the event queue */
static
SCIP_RETCODE varEventGlbChanged(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             oldbound,           /**< old lower bound for variable */
   SCIP_Real             newbound            /**< new lower bound for variable */
   )
{
   assert(var != NULL);
   assert(var->eventfilter != NULL);
   assert(SCIPvarIsTransformed(var));
   assert(!SCIPsetIsEQ(set, oldbound, newbound));
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* check, if the variable is being tracked for bound changes
    * COLUMN and LOOSE variables are tracked always, because global/root pseudo objective value has to be updated
    */
   if( (var->eventfilter->len > 0 && (var->eventfilter->eventmask & SCIP_EVENTTYPE_GLBCHANGED) != 0)
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      SCIP_EVENT* event;

      SCIPsetDebugMsg(set, "issue GLBCHANGED event for variable <%s>: %g -> %g\n", var->name, oldbound, newbound);

      SCIP_CALL( SCIPeventCreateGlbChanged(&event, blkmem, var, oldbound, newbound) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, lp, branchcand, NULL, &event) );
   }

   return SCIP_OKAY;
}

/** appends GUBCHANGED event to the event queue */
static
SCIP_RETCODE varEventGubChanged(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             oldbound,           /**< old lower bound for variable */
   SCIP_Real             newbound            /**< new lower bound for variable */
   )
{
   assert(var != NULL);
   assert(var->eventfilter != NULL);
   assert(SCIPvarIsTransformed(var));
   assert(!SCIPsetIsEQ(set, oldbound, newbound));
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* check, if the variable is being tracked for bound changes
    * COLUMN and LOOSE variables are tracked always, because global/root pseudo objective value has to be updated
    */
   if( (var->eventfilter->len > 0 && (var->eventfilter->eventmask & SCIP_EVENTTYPE_GUBCHANGED) != 0)
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      SCIP_EVENT* event;

      SCIPsetDebugMsg(set, "issue GUBCHANGED event for variable <%s>: %g -> %g\n", var->name, oldbound, newbound);

      SCIP_CALL( SCIPeventCreateGubChanged(&event, blkmem, var, oldbound, newbound) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, lp, branchcand, NULL, &event) );
   }

   return SCIP_OKAY;
}

/** appends GHOLEADDED event to the event queue */
static
SCIP_RETCODE varEventGholeAdded(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right               /**< right bound of open interval in new hole */
   )
{
   assert(var != NULL);
   assert(var->eventfilter != NULL);
   assert(SCIPvarIsTransformed(var));
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(SCIPsetIsLT(set, left, right));

   /* check, if the variable is being tracked for bound changes */
   if( (var->eventfilter->len > 0 && (var->eventfilter->eventmask & SCIP_EVENTTYPE_GHOLEADDED) != 0) )
   {
      SCIP_EVENT* event;

      SCIPsetDebugMsg(set, "issue GHOLEADDED event for variable <%s>: (%.15g,%.15g)\n", var->name, left, right);

      SCIP_CALL( SCIPeventCreateGholeAdded(&event, blkmem, var, left, right) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, NULL, &event) );
   }

   return SCIP_OKAY;
}

/** increases root bound change statistics after a global bound change */
static
void varIncRootboundchgs(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(stat != NULL);

   if( SCIPvarIsActive(var) && SCIPvarIsTransformed(var) && set->stage == SCIP_STAGE_SOLVING )
   {
      stat->nrootboundchgs++;
      stat->nrootboundchgsrun++;
      if( SCIPvarIsIntegral(var) && SCIPvarGetLbGlobal(var) + 0.5 > SCIPvarGetUbGlobal(var) )
      {
         stat->nrootintfixings++;
         stat->nrootintfixingsrun++;
      }
   }
}

/* forward declaration, because both methods call each other recursively */

/* performs the current change in upper bound, changes all parents accordingly */
static
SCIP_RETCODE varProcessChgUbGlobal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real             newbound            /**< new bound for variable */
   );

/** performs the current change in lower bound, changes all parents accordingly */
static
SCIP_RETCODE varProcessChgLbGlobal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   SCIP_VAR* parentvar;
   SCIP_Real oldbound;
   int i;

   assert(var != NULL);
   /* local domains can violate global bounds but not more than feasibility epsilon */
   assert(SCIPsetIsFeasLE(set, var->glbdom.lb, var->locdom.lb));
   assert(SCIPsetIsFeasLE(set, var->locdom.ub, var->glbdom.ub));
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(stat != NULL);

   /* adjust bound to integral value if variable is of integral type */
   newbound = adjustedLb(set, SCIPvarGetType(var), newbound);

   /* check that the bound is feasible */
   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM && newbound > var->glbdom.ub )
   {
      /* due to numerics we only want to be feasible in feasibility tolerance */
      assert(SCIPsetIsFeasLE(set, newbound, var->glbdom.ub));
      newbound = var->glbdom.ub;
   }
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));

   assert(var->vartype != SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, newbound, 0.0) || SCIPsetIsEQ(set, newbound, 1.0));  /*lint !e641*/

   SCIPsetDebugMsg(set, "process changing global lower bound of <%s> from %f to %f\n", var->name, var->glbdom.lb, newbound);

   if( SCIPsetIsEQ(set, newbound, var->glbdom.lb) )
      return SCIP_OKAY;

   /* check bound on debugging solution */
   SCIP_CALL( SCIPdebugCheckLbGlobal(set->scip, var, newbound) ); /*lint !e506 !e774*/

   /* change the bound */
   oldbound = var->glbdom.lb;
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || SCIPsetIsFeasLE(set, newbound, var->glbdom.ub));
   var->glbdom.lb = newbound;
   assert( SCIPsetIsFeasLE(set, var->glbdom.lb, var->locdom.lb) );
   assert( SCIPsetIsFeasLE(set, var->locdom.ub, var->glbdom.ub) );

   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM )
   {
      /* merges overlapping holes into single holes, moves bounds respectively */
      domMerge(&var->glbdom, blkmem, set, &newbound, NULL);
   }

   /* update the root bound changes counters */
   varIncRootboundchgs(var, set, stat);

   /* update the lbchginfos array by replacing worse local bounds with the new global bound and changing the
    * redundant bound changes to be branching decisions
    */
   for( i = 0; i < var->nlbchginfos; ++i )
   {
      assert(var->lbchginfos[i].var == var);

      if( var->lbchginfos[i].oldbound < var->glbdom.lb )
      {
         SCIPsetDebugMsg(set, " -> adjust lower bound change <%s>: %g -> %g due to new global lower bound %g\n",
            SCIPvarGetName(var), var->lbchginfos[i].oldbound, var->lbchginfos[i].newbound, var->glbdom.lb);
         var->lbchginfos[i].oldbound = var->glbdom.lb;
         if( SCIPsetIsLE(set, var->lbchginfos[i].newbound, var->glbdom.lb) )
         {
            /* this bound change is redundant due to the new global bound */
            var->lbchginfos[i].newbound = var->glbdom.lb;
            var->lbchginfos[i].boundchgtype = SCIP_BOUNDCHGTYPE_BRANCHING; /*lint !e641*/
            var->lbchginfos[i].redundant = TRUE;
         }
         else
            break; /* from now on, the remaining local bound changes are not redundant */
      }
      else
         break; /* from now on, the remaining local bound changes are not redundant */
   }

   /* remove redundant implications and variable bounds */
   if( (SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE)
      && (!set->reopt_enable || set->stage == SCIP_STAGE_PRESOLVING) )
   {
      SCIP_CALL( SCIPvarRemoveCliquesImplicsVbs(var, blkmem, cliquetable, set, FALSE, TRUE, TRUE) );
   }

   /* issue bound change event */
   assert(SCIPvarIsTransformed(var) == (var->eventfilter != NULL));
   if( var->eventfilter != NULL )
   {
      SCIP_CALL( varEventGlbChanged(var, blkmem, set, lp, branchcand, eventqueue, oldbound, newbound) );
   }

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         SCIP_CALL( varProcessChgLbGlobal(parentvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newbound) );
         break;

      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         SCIPerrorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         return SCIP_INVALIDDATA;

      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            SCIP_Real parentnewbound;

            /* a > 0 -> change lower bound of y */
            assert((SCIPsetIsInfinity(set, -parentvar->glbdom.lb) && SCIPsetIsInfinity(set, -oldbound))
               || SCIPsetIsFeasEQ(set, parentvar->glbdom.lb, oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant)
               || (SCIPsetIsZero(set, parentvar->glbdom.lb / parentvar->data.aggregate.scalar) && SCIPsetIsZero(set, oldbound)));

            if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
               parentnewbound = parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant;
            else
               parentnewbound = newbound;
            SCIP_CALL( varProcessChgLbGlobal(parentvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, parentnewbound) );
         }
         else
         {
            SCIP_Real parentnewbound;

            /* a < 0 -> change upper bound of y */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert((SCIPsetIsInfinity(set, parentvar->glbdom.ub) && SCIPsetIsInfinity(set, -oldbound))
               || SCIPsetIsFeasEQ(set, parentvar->glbdom.ub, oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant)
               || (SCIPsetIsZero(set, parentvar->glbdom.ub / parentvar->data.aggregate.scalar) && SCIPsetIsZero(set, oldbound)));

            if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
               parentnewbound = parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant;
            else
               parentnewbound = -newbound;
            SCIP_CALL( varProcessChgUbGlobal(parentvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, parentnewbound) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         SCIP_CALL( varProcessChgUbGlobal(parentvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable,
               parentvar->data.negate.constant - newbound) );
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** performs the current change in upper bound, changes all parents accordingly */
static
SCIP_RETCODE varProcessChgUbGlobal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   SCIP_VAR* parentvar;
   SCIP_Real oldbound;
   int i;

   assert(var != NULL);
   /* local domains can violate global bounds but not more than feasibility epsilon */
   assert(SCIPsetIsFeasLE(set, var->glbdom.lb , var->locdom.lb));
   assert(SCIPsetIsFeasLE(set, var->locdom.ub, var->glbdom.ub));
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(stat != NULL);

   /* adjust bound to integral value if variable is of integral type */
   newbound = adjustedUb(set, SCIPvarGetType(var), newbound);

   /* check that the bound is feasible */
   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM && newbound < var->glbdom.lb )
   {
      /* due to numerics we only want to be feasible in feasibility tolerance */
      assert(SCIPsetIsFeasGE(set, newbound, var->glbdom.lb));
      newbound = var->glbdom.lb;
   }
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));

   assert(var->vartype != SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, newbound, 0.0) || SCIPsetIsEQ(set, newbound, 1.0));  /*lint !e641*/

   SCIPsetDebugMsg(set, "process changing global upper bound of <%s> from %f to %f\n", var->name, var->glbdom.ub, newbound);

   if( SCIPsetIsEQ(set, newbound, var->glbdom.ub) )
      return SCIP_OKAY;

   /* check bound on debugging solution */
   SCIP_CALL( SCIPdebugCheckUbGlobal(set->scip, var, newbound) ); /*lint !e506 !e774*/

   /* change the bound */
   oldbound = var->glbdom.ub;
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || SCIPsetIsFeasGE(set, newbound, var->glbdom.lb));
   var->glbdom.ub = newbound;
   assert( SCIPsetIsFeasLE(set, var->glbdom.lb, var->locdom.lb) );
   assert( SCIPsetIsFeasLE(set, var->locdom.ub, var->glbdom.ub) );

   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM )
   {
      /* merges overlapping holes into single holes, moves bounds respectively */
      domMerge(&var->glbdom, blkmem, set, NULL, &newbound);
   }

   /* update the root bound changes counters */
   varIncRootboundchgs(var, set, stat);

   /* update the ubchginfos array by replacing worse local bounds with the new global bound and changing the
    * redundant bound changes to be branching decisions
    */
   for( i = 0; i < var->nubchginfos; ++i )
   {
      assert(var->ubchginfos[i].var == var);
      if( var->ubchginfos[i].oldbound > var->glbdom.ub )
      {
         SCIPsetDebugMsg(set, " -> adjust upper bound change <%s>: %g -> %g due to new global upper bound %g\n",
            SCIPvarGetName(var), var->ubchginfos[i].oldbound, var->ubchginfos[i].newbound, var->glbdom.ub);
         var->ubchginfos[i].oldbound = var->glbdom.ub;
         if( SCIPsetIsGE(set, var->ubchginfos[i].newbound, var->glbdom.ub) )
         {
            /* this bound change is redundant due to the new global bound */
            var->ubchginfos[i].newbound = var->glbdom.ub;
            var->ubchginfos[i].boundchgtype = SCIP_BOUNDCHGTYPE_BRANCHING; /*lint !e641*/
            var->ubchginfos[i].redundant = TRUE;
         }
         else
            break; /* from now on, the remaining local bound changes are not redundant */
      }
      else
         break; /* from now on, the remaining local bound changes are not redundant */
   }

   /* remove redundant implications and variable bounds */
   if( (SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE)
      && (!set->reopt_enable || set->stage == SCIP_STAGE_PRESOLVING) )
   {
      SCIP_CALL( SCIPvarRemoveCliquesImplicsVbs(var, blkmem, cliquetable, set, FALSE, TRUE, TRUE) );
   }

   /* issue bound change event */
   assert(SCIPvarIsTransformed(var) == (var->eventfilter != NULL));
   if( var->eventfilter != NULL )
   {
      SCIP_CALL( varEventGubChanged(var, blkmem, set, lp, branchcand, eventqueue, oldbound, newbound) );
   }

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         SCIP_CALL( varProcessChgUbGlobal(parentvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newbound) );
         break;

      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         SCIPerrorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         return SCIP_INVALIDDATA;

      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            SCIP_Real parentnewbound;

            /* a > 0 -> change upper bound of y */
            assert((SCIPsetIsInfinity(set, parentvar->glbdom.ub) && SCIPsetIsInfinity(set, oldbound))
               || SCIPsetIsFeasEQ(set, parentvar->glbdom.ub,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
               parentnewbound = parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant;
            else
               parentnewbound = newbound;
            SCIP_CALL( varProcessChgUbGlobal(parentvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, parentnewbound) );
         }
         else
         {
            SCIP_Real parentnewbound;

            /* a < 0 -> change lower bound of y */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert((SCIPsetIsInfinity(set, -parentvar->glbdom.lb) && SCIPsetIsInfinity(set, oldbound))
               || SCIPsetIsFeasEQ(set, parentvar->glbdom.lb,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
               parentnewbound = parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant;
            else
               parentnewbound = -newbound;
            SCIP_CALL( varProcessChgLbGlobal(parentvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, parentnewbound) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         SCIP_CALL( varProcessChgLbGlobal(parentvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable,
               parentvar->data.negate.constant - newbound) );
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** changes global lower bound of variable; if possible, adjusts bound to integral value;
 *  updates local lower bound if the global bound is tighter
 */
SCIP_RETCODE SCIPvarChgLbGlobal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* check that the bound is feasible; this must be w.r.t. feastol because SCIPvarFix() allows fixings that are outside
    * of the domain within feastol
    */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || !SCIPsetIsFeasGT(set, newbound, var->glbdom.ub));

   /* adjust bound to integral value if variable is of integral type */
   newbound = adjustedLb(set, SCIPvarGetType(var), newbound);

   /* check that the adjusted bound is feasible
    * @todo this does not have to be the case if the original problem was infeasible due to bounds and we are called
    *       here because we reset bounds to their original value!
    */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || !SCIPsetIsFeasGT(set, newbound, var->glbdom.ub));

   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM )
   {
      /* we do not want to exceed the upperbound, which could have happened due to numerics */
      newbound = MIN(newbound, var->glbdom.ub);
   }
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));

   /* the new global bound has to be tighter except we are in the original problem; this must be w.r.t. feastol because
    * SCIPvarFix() allows fixings that are outside of the domain within feastol
    */
   assert(lp == NULL || SCIPsetIsFeasLE(set, var->glbdom.lb, newbound) || (set->reopt_enable && set->stage == SCIP_STAGE_PRESOLVED));

   SCIPsetDebugMsg(set, "changing global lower bound of <%s> from %g to %g\n", var->name, var->glbdom.lb, newbound);

   if( SCIPsetIsEQ(set, var->glbdom.lb, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar != NULL )
      {
         SCIP_CALL( SCIPvarChgLbGlobal(var->data.original.transvar, blkmem, set, stat, lp, branchcand, eventqueue,
               cliquetable, newbound) );
      }
      else
      {
         assert(set->stage == SCIP_STAGE_PROBLEM);
         if( newbound > SCIPvarGetLbLocal(var) )
         {
            SCIP_CALL( SCIPvarChgLbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, newbound) );
         }
         SCIP_CALL( varProcessChgLbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newbound) );
      }
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      if( newbound > SCIPvarGetLbLocal(var) )
      {
         SCIP_CALL( SCIPvarChgLbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, newbound) );
      }
      SCIP_CALL( varProcessChgLbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newbound) );
      break;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a > 0 -> change lower bound of y */
         assert((SCIPsetIsInfinity(set, -var->glbdom.lb) && SCIPsetIsInfinity(set, -var->data.aggregate.var->glbdom.lb))
            || SCIPsetIsFeasEQ(set, var->glbdom.lb,
               var->data.aggregate.var->glbdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = newbound;
         SCIP_CALL( SCIPvarChgLbGlobal(var->data.aggregate.var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable,
               childnewbound) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a < 0 -> change upper bound of y */
         assert((SCIPsetIsInfinity(set, -var->glbdom.lb) && SCIPsetIsInfinity(set, var->data.aggregate.var->glbdom.ub))
            || SCIPsetIsFeasEQ(set, var->glbdom.lb,
               var->data.aggregate.var->glbdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = -newbound;
         SCIP_CALL( SCIPvarChgUbGlobal(var->data.aggregate.var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable,
               childnewbound) );
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot change the bounds of a multi-aggregated variable.\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarChgUbGlobal(var->negatedvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable,
            var->data.negate.constant - newbound) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** changes global upper bound of variable; if possible, adjusts bound to integral value;
 *  updates local upper bound if the global bound is tighter
 */
SCIP_RETCODE SCIPvarChgUbGlobal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* check that the bound is feasible; this must be w.r.t. feastol because SCIPvarFix() allows fixings that are outside
    * of the domain within feastol
    */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || !SCIPsetIsFeasLT(set, newbound, var->glbdom.lb));

   /* adjust bound to integral value if variable is of integral type */
   newbound = adjustedUb(set, SCIPvarGetType(var), newbound);

   /* check that the adjusted bound is feasible
    * @todo this does not have to be the case if the original problem was infeasible due to bounds and we are called
    *       here because we reset bounds to their original value!
    */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || !SCIPsetIsFeasLT(set, newbound, var->glbdom.lb));

   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM )
   {
      /* we do not want to undercut the lowerbound, which could have happened due to numerics */
      newbound = MAX(newbound, var->glbdom.lb);
   }
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));

   /* the new global bound has to be tighter except we are in the original problem; this must be w.r.t. feastol because
    * SCIPvarFix() allows fixings that are outside of the domain within feastol
    */
   assert(lp == NULL || SCIPsetIsFeasGE(set, var->glbdom.ub, newbound) || (set->reopt_enable && set->stage == SCIP_STAGE_PRESOLVED));

   SCIPsetDebugMsg(set, "changing global upper bound of <%s> from %g to %g\n", var->name, var->glbdom.ub, newbound);

   if( SCIPsetIsEQ(set, var->glbdom.ub, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar != NULL )
      {
         SCIP_CALL( SCIPvarChgUbGlobal(var->data.original.transvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable,
               newbound) );
      }
      else
      {
         assert(set->stage == SCIP_STAGE_PROBLEM);
         if( newbound < SCIPvarGetUbLocal(var) )
         {
            SCIP_CALL( SCIPvarChgUbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, newbound) );
         }
         SCIP_CALL( varProcessChgUbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newbound) );
      }
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      if( newbound < SCIPvarGetUbLocal(var) )
      {
         SCIP_CALL( SCIPvarChgUbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, newbound) );
      }
      SCIP_CALL( varProcessChgUbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newbound) );
      break;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a > 0 -> change lower bound of y */
         assert((SCIPsetIsInfinity(set, var->glbdom.ub) && SCIPsetIsInfinity(set, var->data.aggregate.var->glbdom.ub))
            || SCIPsetIsFeasEQ(set, var->glbdom.ub,
               var->data.aggregate.var->glbdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = newbound;
         SCIP_CALL( SCIPvarChgUbGlobal(var->data.aggregate.var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable,
               childnewbound) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a < 0 -> change upper bound of y */
         assert((SCIPsetIsInfinity(set, var->glbdom.ub) && SCIPsetIsInfinity(set, -var->data.aggregate.var->glbdom.lb))
            || SCIPsetIsFeasEQ(set, var->glbdom.ub,
               var->data.aggregate.var->glbdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = -newbound;
         SCIP_CALL( SCIPvarChgLbGlobal(var->data.aggregate.var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable,
               childnewbound) );
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot change the bounds of a multi-aggregated variable.\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarChgLbGlobal(var->negatedvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable,
            var->data.negate.constant - newbound) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** changes lazy lower bound of the variable, this is only possible if the variable is not in the LP yet */
SCIP_RETCODE SCIPvarChgLbLazy(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             lazylb              /**< the lazy lower bound to be set */
   )
{
   assert(var != NULL);
   assert(var->probindex != -1);
   assert(SCIPsetIsFeasGE(set, var->glbdom.ub, lazylb));
   assert(SCIPsetIsFeasGE(set, var->lazyub, lazylb));
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* variable should not be in the LP */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      return SCIP_INVALIDCALL;

   var->lazylb = lazylb;

   return SCIP_OKAY;
}

/** changes lazy upper bound of the variable, this is only possible if the variable is not in the LP yet */
SCIP_RETCODE SCIPvarChgUbLazy(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             lazyub              /**< the lazy lower bound to be set */
   )
{
   assert(var != NULL);
   assert(var->probindex != -1);
   assert(SCIPsetIsFeasGE(set, lazyub, var->glbdom.lb));
   assert(SCIPsetIsFeasGE(set, lazyub, var->lazylb));
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* variable should not be in the LP */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      return SCIP_INVALIDCALL;

   var->lazyub = lazyub;

   return SCIP_OKAY;
}


/** changes global bound of variable; if possible, adjusts bound to integral value;
 *  updates local bound if the global bound is tighter
 */
SCIP_RETCODE SCIPvarChgBdGlobal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real             newbound,           /**< new bound for variable */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   /* apply bound change to the LP data */
   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      return SCIPvarChgLbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newbound);
   case SCIP_BOUNDTYPE_UPPER:
      return SCIPvarChgUbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newbound);
   default:
      SCIPerrorMessage("unknown bound type\n");
      return SCIP_INVALIDDATA;
   }
}

/** appends LBTIGHTENED or LBRELAXED event to the event queue */
static
SCIP_RETCODE varEventLbChanged(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             oldbound,           /**< old lower bound for variable */
   SCIP_Real             newbound            /**< new lower bound for variable */
   )
{
   assert(var != NULL);
   assert(var->eventfilter != NULL);
   assert(SCIPvarIsTransformed(var));
   assert(!SCIPsetIsEQ(set, oldbound, newbound) || newbound == var->glbdom.lb); /*lint !e777*/
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* check, if the variable is being tracked for bound changes
    * COLUMN and LOOSE variables are tracked always, because row activities and LP changes have to be updated
    */
   if( (var->eventfilter->len > 0 && (var->eventfilter->eventmask & SCIP_EVENTTYPE_LBCHANGED) != 0)
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      SCIP_EVENT* event;

      SCIPsetDebugMsg(set, "issue LBCHANGED event for variable <%s>: %g -> %g\n", var->name, oldbound, newbound);

      SCIP_CALL( SCIPeventCreateLbChanged(&event, blkmem, var, oldbound, newbound) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, lp, branchcand, NULL, &event) );
   }

   return SCIP_OKAY;
}

/** appends UBTIGHTENED or UBRELAXED event to the event queue */
static
SCIP_RETCODE varEventUbChanged(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             oldbound,           /**< old upper bound for variable */
   SCIP_Real             newbound            /**< new upper bound for variable */
   )
{
   assert(var != NULL);
   assert(var->eventfilter != NULL);
   assert(SCIPvarIsTransformed(var));
   assert(!SCIPsetIsEQ(set, oldbound, newbound) || newbound == var->glbdom.ub); /*lint !e777*/
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* check, if the variable is being tracked for bound changes
    * COLUMN and LOOSE variables are tracked always, because row activities and LP changes have to be updated
    */
   if( (var->eventfilter->len > 0 && (var->eventfilter->eventmask & SCIP_EVENTTYPE_UBCHANGED) != 0)
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      SCIP_EVENT* event;

      SCIPsetDebugMsg(set, "issue UBCHANGED event for variable <%s>: %g -> %g\n", var->name, oldbound, newbound);

      SCIP_CALL( SCIPeventCreateUbChanged(&event, blkmem, var, oldbound, newbound) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, lp, branchcand, NULL, &event) );
   }

   return SCIP_OKAY;
}

/* forward declaration, because both methods call each other recursively */

/* performs the current change in upper bound, changes all parents accordingly */
static
SCIP_RETCODE varProcessChgUbLocal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics, or NULL if the bound change belongs to updating the parent variables */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             newbound            /**< new bound for variable */
   );

/** performs the current change in lower bound, changes all parents accordingly */
static
SCIP_RETCODE varProcessChgLbLocal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics, or NULL if the bound change belongs to updating the parent variables */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   SCIP_VAR* parentvar;
   SCIP_Real oldbound;
   int i;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert((SCIPvarGetType(var) == SCIP_VARTYPE_BINARY && (SCIPsetIsZero(set, newbound) || SCIPsetIsEQ(set, newbound, 1.0)
            || SCIPsetIsEQ(set, newbound, var->locdom.ub)))
      || (SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS && (SCIPsetIsIntegral(set, newbound)
            || SCIPsetIsEQ(set, newbound, var->locdom.ub)))
      || SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

   /* check that the bound is feasible */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || SCIPsetIsLE(set, newbound, var->glbdom.ub));
   /* adjust bound to integral value if variable is of integral type */
   newbound = adjustedLb(set, SCIPvarGetType(var), newbound);

   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM )
   {
      /* we do not want to exceed the upper bound, which could have happened due to numerics */
      newbound = MIN(newbound, var->locdom.ub);

      /* we do not want to undercut the global lower bound, which could have happened due to numerics */
      newbound = MAX(newbound, var->glbdom.lb);
   }
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));

   SCIPsetDebugMsg(set, "process changing lower bound of <%s> from %g to %g\n", var->name, var->locdom.lb, newbound);

   if( SCIPsetIsEQ(set, newbound, var->glbdom.lb) && var->glbdom.lb != var->locdom.lb ) /*lint !e777*/
      newbound = var->glbdom.lb;
   else if( SCIPsetIsEQ(set, newbound, var->locdom.lb) )
      return SCIP_OKAY;

   /* change the bound */
   oldbound = var->locdom.lb;
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || SCIPsetIsFeasLE(set, newbound, var->locdom.ub));
   var->locdom.lb = newbound;

   /* update statistic; during the update steps of the parent variable we pass a NULL pointer to ensure that we only
    * once update the statistic
    */
   if( stat != NULL )
      SCIPstatIncrement(stat, set, domchgcount);

   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM )
   {
      /* merges overlapping holes into single holes, moves bounds respectively */
      domMerge(&var->locdom, blkmem, set, &newbound, NULL);
   }

   /* issue bound change event */
   assert(SCIPvarIsTransformed(var) == (var->eventfilter != NULL));
   if( var->eventfilter != NULL )
   {
      SCIP_CALL( varEventLbChanged(var, blkmem, set, lp, branchcand, eventqueue, oldbound, newbound) );
   }

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         SCIP_CALL( varProcessChgLbLocal(parentvar, blkmem, set, NULL, lp, branchcand, eventqueue, newbound) );
         break;

      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         SCIPerrorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         return SCIP_INVALIDDATA;

      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            SCIP_Real parentnewbound;

            /* a > 0 -> change lower bound of y */
            assert((SCIPsetIsInfinity(set, -parentvar->locdom.lb) && SCIPsetIsInfinity(set, -oldbound))
               || SCIPsetIsFeasEQ(set, parentvar->locdom.lb, oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant)
               || (SCIPsetIsZero(set, parentvar->locdom.lb / parentvar->data.aggregate.scalar) && SCIPsetIsZero(set, oldbound)));

            if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            {
               parentnewbound = parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant;
               /* if parent's new lower bound exceeds its upper bound, then this could be due to numerical difficulties, e.g., if numbers are large
                * thus, at least a relative comparision of the new lower bound and the current upper bound should proof consistency
                * as a result, the parent's lower bound is set to it's upper bound, and not above
                */
               if( parentnewbound > parentvar->glbdom.ub )
               {
                  /* due to numerics we only need to be feasible w.r.t. feasibility tolerance */
                  assert(SCIPsetIsFeasLE(set, parentnewbound, parentvar->glbdom.ub));
                  parentnewbound = parentvar->glbdom.ub;
               }
            }
            else
               parentnewbound = newbound;
            SCIP_CALL( varProcessChgLbLocal(parentvar, blkmem, set, NULL, lp, branchcand, eventqueue, parentnewbound) );
         }
         else
         {
            SCIP_Real parentnewbound;

            /* a < 0 -> change upper bound of y */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert((SCIPsetIsInfinity(set, parentvar->locdom.ub) && SCIPsetIsInfinity(set, -oldbound))
               || SCIPsetIsFeasEQ(set, parentvar->locdom.ub, oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant)
               || (SCIPsetIsZero(set, parentvar->locdom.ub / parentvar->data.aggregate.scalar) && SCIPsetIsZero(set, oldbound)));

            if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            {
               parentnewbound = parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant;
               /* if parent's new upper bound is below its lower bound, then this could be due to numerical difficulties, e.g., if numbers are large
                * thus, at least a relative comparision of the new upper bound and the current lower bound should proof consistency
                * as a result, the parent's upper bound is set to it's lower bound, and not below
                */
               if( parentnewbound < parentvar->glbdom.lb )
               {
                  /* due to numerics we only need to be feasible w.r.t. feasibility tolerance */
                  assert(SCIPsetIsFeasGE(set, parentnewbound, parentvar->glbdom.lb));
                  parentnewbound = parentvar->glbdom.lb;
               }
            }
            else
               parentnewbound = -newbound;
            SCIP_CALL( varProcessChgUbLocal(parentvar, blkmem, set, NULL, lp, branchcand, eventqueue, parentnewbound) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x = offset - x'  ->  x' = offset - x */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         SCIP_CALL( varProcessChgUbLocal(parentvar, blkmem, set, NULL, lp, branchcand, eventqueue,
               parentvar->data.negate.constant - newbound) );
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** performs the current change in upper bound, changes all parents accordingly */
static
SCIP_RETCODE varProcessChgUbLocal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics, or NULL if the bound change belongs to updating the parent variables */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   SCIP_VAR* parentvar;
   SCIP_Real oldbound;
   int i;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert((SCIPvarGetType(var) == SCIP_VARTYPE_BINARY && (SCIPsetIsZero(set, newbound) || SCIPsetIsEQ(set, newbound, 1.0)
            || SCIPsetIsEQ(set, newbound, var->locdom.lb)))
      || (SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS && (SCIPsetIsIntegral(set, newbound)
            || SCIPsetIsEQ(set, newbound, var->locdom.lb)))
      || SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

   /* check that the bound is feasible */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || SCIPsetIsGE(set, newbound, var->glbdom.lb));
   /* adjust bound to integral value if variable is of integral type */
   newbound = adjustedUb(set, SCIPvarGetType(var), newbound);

   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM )
   {
      /* we do not want to undercut the lower bound, which could have happened due to numerics */
      newbound = MAX(newbound, var->locdom.lb);

      /* we do not want to exceed the global upper bound, which could have happened due to numerics */
      newbound = MIN(newbound, var->glbdom.ub);
   }
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));

   SCIPsetDebugMsg(set, "process changing upper bound of <%s> from %g to %g\n", var->name, var->locdom.ub, newbound);

   if( SCIPsetIsEQ(set, newbound, var->glbdom.ub) && var->glbdom.ub != var->locdom.ub  ) /*lint !e777*/
      newbound = var->glbdom.ub;
   else if( SCIPsetIsEQ(set, newbound, var->locdom.ub) )
      return SCIP_OKAY;

   /* change the bound */
   oldbound = var->locdom.ub;
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || SCIPsetIsFeasGE(set, newbound, var->locdom.lb));
   var->locdom.ub = newbound;

   /* update statistic; during the update steps of the parent variable we pass a NULL pointer to ensure that we only
    * once update the statistic
    */
   if( stat != NULL )
      SCIPstatIncrement(stat, set, domchgcount);

   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM )
   {
      /* merges overlapping holes into single holes, moves bounds respectively */
      domMerge(&var->locdom, blkmem, set, NULL, &newbound);
   }

   /* issue bound change event */
   assert(SCIPvarIsTransformed(var) == (var->eventfilter != NULL));
   if( var->eventfilter != NULL )
   {
      SCIP_CALL( varEventUbChanged(var, blkmem, set, lp, branchcand, eventqueue, oldbound, newbound) );
   }

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         SCIP_CALL( varProcessChgUbLocal(parentvar, blkmem, set, NULL, lp, branchcand, eventqueue, newbound) );
         break;

      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         SCIPerrorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         return SCIP_INVALIDDATA;

      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            SCIP_Real parentnewbound;

            /* a > 0 -> change upper bound of x */
            assert((SCIPsetIsInfinity(set, parentvar->locdom.ub) && SCIPsetIsInfinity(set, oldbound))
               || SCIPsetIsFeasEQ(set, parentvar->locdom.ub,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            {
               parentnewbound = parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant;
               /* if parent's new upper bound is below its lower bound, then this could be due to numerical difficulties, e.g., if numbers are large
                * thus, at least a relative comparision of the new upper bound and the current lower bound should proof consistency
                * as a result, the parent's upper bound is set to it's lower bound, and not below
                */
               if( parentnewbound < parentvar->glbdom.lb )
               {
                  /* due to numerics we only need to be feasible w.r.t. feasibility tolerance */
                  assert(SCIPsetIsFeasGE(set, parentnewbound, parentvar->glbdom.lb));
                  parentnewbound = parentvar->glbdom.lb;
               }
            }
            else
               parentnewbound = newbound;
            SCIP_CALL( varProcessChgUbLocal(parentvar, blkmem, set, NULL, lp, branchcand, eventqueue, parentnewbound) );
         }
         else
         {
            SCIP_Real parentnewbound;

            /* a < 0 -> change lower bound of x */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert((SCIPsetIsInfinity(set, -parentvar->locdom.lb) && SCIPsetIsInfinity(set, oldbound))
               || SCIPsetIsFeasEQ(set, parentvar->locdom.lb,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            {
               parentnewbound = parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant;
               /* if parent's new lower bound exceeds its upper bound, then this could be due to numerical difficulties, e.g., if numbers are large
                * thus, at least a relative comparision of the new lower bound and the current upper bound should proof consistency
                * as a result, the parent's lower bound is set to it's upper bound, and not above
                */
               if( parentnewbound > parentvar->glbdom.ub )
               {
                  /* due to numerics we only need to be feasible w.r.t. feasibility tolerance */
                  assert(SCIPsetIsFeasLE(set, parentnewbound, parentvar->glbdom.ub));
                  parentnewbound = parentvar->glbdom.ub;
               }
            }
            else
               parentnewbound = -newbound;
            SCIP_CALL( varProcessChgLbLocal(parentvar, blkmem, set, NULL, lp, branchcand, eventqueue, parentnewbound) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x = offset - x'  ->  x' = offset - x */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         SCIP_CALL( varProcessChgLbLocal(parentvar, blkmem, set, NULL, lp, branchcand, eventqueue,
               parentvar->data.negate.constant - newbound) );
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** changes current local lower bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
SCIP_RETCODE SCIPvarChgLbLocal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* check that the bound is feasible; this must be w.r.t. feastol because SCIPvarFix() allows fixings that are outside
    * of the domain within feastol
    */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || !SCIPsetIsFeasGT(set, newbound, var->locdom.ub));

   /* adjust bound to integral value if variable is of integral type */
   newbound = adjustedLb(set, SCIPvarGetType(var), newbound);

   /* check that the adjusted bound is feasible */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || !SCIPsetIsFeasGT(set, newbound, var->locdom.ub));

   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM )
   {
      /* we do not want to exceed the upperbound, which could have happened due to numerics */
      newbound = MIN(newbound, var->locdom.ub);
   }
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));

   SCIPsetDebugMsg(set, "changing lower bound of <%s>[%g,%g] to %g\n", var->name, var->locdom.lb, var->locdom.ub, newbound);

   if( SCIPsetIsEQ(set, var->locdom.lb, newbound) && (!SCIPsetIsEQ(set, var->glbdom.lb, newbound) || var->locdom.lb == newbound) ) /*lint !e777*/
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar != NULL )
      {
         SCIP_CALL( SCIPvarChgLbLocal(var->data.original.transvar, blkmem, set, stat, lp, branchcand, eventqueue,
               newbound) );
      }
      else
      {
         assert(set->stage == SCIP_STAGE_PROBLEM);
         SCIP_CALL( varProcessChgLbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, newbound) );
      }
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      SCIP_CALL( varProcessChgLbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, newbound) );
      break;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a > 0 -> change lower bound of y */
         assert((SCIPsetIsInfinity(set, -var->locdom.lb) && SCIPsetIsInfinity(set, -var->data.aggregate.var->locdom.lb))
            || SCIPsetIsFeasEQ(set, var->locdom.lb,
               var->data.aggregate.var->locdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = newbound;
         SCIP_CALL( SCIPvarChgLbLocal(var->data.aggregate.var, blkmem, set, stat, lp, branchcand, eventqueue,
               childnewbound) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a < 0 -> change upper bound of y */
         assert((SCIPsetIsInfinity(set, -var->locdom.lb) && SCIPsetIsInfinity(set, var->data.aggregate.var->locdom.ub))
            || SCIPsetIsFeasEQ(set, var->locdom.lb,
               var->data.aggregate.var->locdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = -newbound;
         SCIP_CALL( SCIPvarChgUbLocal(var->data.aggregate.var, blkmem, set, stat, lp, branchcand, eventqueue,
               childnewbound) );
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot change the bounds of a multi-aggregated variable.\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarChgUbLocal(var->negatedvar, blkmem, set, stat, lp, branchcand, eventqueue,
            var->data.negate.constant - newbound) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** changes current local upper bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
SCIP_RETCODE SCIPvarChgUbLocal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* check that the bound is feasible; this must be w.r.t. feastol because SCIPvarFix() allows fixings that are outside
    * of the domain within feastol
    */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || !SCIPsetIsFeasLT(set, newbound, var->locdom.lb));

   /* adjust bound to integral value if variable is of integral type */
   newbound = adjustedUb(set, SCIPvarGetType(var), newbound);

   /* check that the adjusted bound is feasible */
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM || !SCIPsetIsFeasLT(set, newbound, var->locdom.lb));

   if( SCIPsetGetStage(set) != SCIP_STAGE_PROBLEM )
   {
      /* we do not want to undercut the lowerbound, which could have happened due to numerics */
      newbound = MAX(newbound, var->locdom.lb);
   }
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));

   SCIPsetDebugMsg(set, "changing upper bound of <%s>[%g,%g] to %g\n", var->name, var->locdom.lb, var->locdom.ub, newbound);

   if( SCIPsetIsEQ(set, var->locdom.ub, newbound) && (!SCIPsetIsEQ(set, var->glbdom.ub, newbound) || var->locdom.ub == newbound) ) /*lint !e777*/
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar != NULL )
      {
         SCIP_CALL( SCIPvarChgUbLocal(var->data.original.transvar, blkmem, set, stat, lp, branchcand, eventqueue, newbound) );
      }
      else
      {
         assert(set->stage == SCIP_STAGE_PROBLEM);
         SCIP_CALL( varProcessChgUbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, newbound) );
      }
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      SCIP_CALL( varProcessChgUbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, newbound) );
      break;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a > 0 -> change upper bound of y */
         assert((SCIPsetIsInfinity(set, var->locdom.ub) && SCIPsetIsInfinity(set, var->data.aggregate.var->locdom.ub))
            || SCIPsetIsFeasEQ(set, var->locdom.ub,
               var->data.aggregate.var->locdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = newbound;
         SCIP_CALL( SCIPvarChgUbLocal(var->data.aggregate.var, blkmem, set, stat, lp, branchcand, eventqueue,
               childnewbound) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a < 0 -> change lower bound of y */
         assert((SCIPsetIsInfinity(set, var->locdom.ub) && SCIPsetIsInfinity(set, -var->data.aggregate.var->locdom.lb))
            || SCIPsetIsFeasEQ(set, var->locdom.ub,
               var->data.aggregate.var->locdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = -newbound;
         SCIP_CALL( SCIPvarChgLbLocal(var->data.aggregate.var, blkmem, set, stat, lp, branchcand, eventqueue,
               childnewbound) );
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot change the bounds of a multi-aggregated variable.\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarChgLbLocal(var->negatedvar, blkmem, set, stat, lp, branchcand, eventqueue,
            var->data.negate.constant - newbound) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** changes current local bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
SCIP_RETCODE SCIPvarChgBdLocal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             newbound,           /**< new bound for variable */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   /* apply bound change to the LP data */
   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      return SCIPvarChgLbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, newbound);
   case SCIP_BOUNDTYPE_UPPER:
      return SCIPvarChgUbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, newbound);
   default:
      SCIPerrorMessage("unknown bound type\n");
      return SCIP_INVALIDDATA;
   }
}

/** changes lower bound of variable in current dive; if possible, adjusts bound to integral value */
SCIP_RETCODE SCIPvarChgLbDive(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(lp != NULL);
   assert(SCIPlpDiving(lp));

   /* adjust bound for integral variables */
   SCIPvarAdjustLb(var, set, &newbound);

   SCIPsetDebugMsg(set, "changing lower bound of <%s> to %g in current dive\n", var->name, newbound);

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      SCIP_CALL( SCIPvarChgLbDive(var->data.original.transvar, set, lp, newbound) );
      break;

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      SCIP_CALL( SCIPcolChgLb(var->data.col, set, lp, newbound) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      SCIPerrorMessage("cannot change variable's bounds in dive for LOOSE variables\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a > 0 -> change lower bound of y */
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = newbound;
         SCIP_CALL( SCIPvarChgLbDive(var->data.aggregate.var, set, lp, childnewbound) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a < 0 -> change upper bound of y */
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = -newbound;
         SCIP_CALL( SCIPvarChgUbDive(var->data.aggregate.var, set, lp, childnewbound) );
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot change the bounds of a multi-aggregated variable.\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarChgUbDive(var->negatedvar, set, lp, var->data.negate.constant - newbound) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** changes upper bound of variable in current dive; if possible, adjusts bound to integral value */
SCIP_RETCODE SCIPvarChgUbDive(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(lp != NULL);
   assert(SCIPlpDiving(lp));

   /* adjust bound for integral variables */
   SCIPvarAdjustUb(var, set, &newbound);

   SCIPsetDebugMsg(set, "changing upper bound of <%s> to %g in current dive\n", var->name, newbound);

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      SCIP_CALL( SCIPvarChgUbDive(var->data.original.transvar, set, lp, newbound) );
      break;

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      SCIP_CALL( SCIPcolChgUb(var->data.col, set, lp, newbound) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      SCIPerrorMessage("cannot change variable's bounds in dive for LOOSE variables\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a > 0 -> change upper bound of y */
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = newbound;
         SCIP_CALL( SCIPvarChgUbDive(var->data.aggregate.var, set, lp, childnewbound) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         SCIP_Real childnewbound;

         /* a < 0 -> change lower bound of y */
         if( !SCIPsetIsInfinity(set, -newbound) && !SCIPsetIsInfinity(set, newbound) )
            childnewbound = (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar;
         else
            childnewbound = -newbound;
         SCIP_CALL( SCIPvarChgLbDive(var->data.aggregate.var, set, lp, childnewbound) );
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot change the bounds of a multi-aggregated variable.\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarChgLbDive(var->negatedvar, set, lp, var->data.negate.constant - newbound) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** for a multi-aggregated variable, gives the local lower bound computed by adding the local bounds from all
 *  aggregation variables, this lower bound may be tighter than the one given by SCIPvarGetLbLocal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
SCIP_Real SCIPvarGetMultaggrLbLocal(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;
   SCIP_Real lb;
   SCIP_Real bnd;
   SCIP_VAR* aggrvar;
   SCIP_Bool posinf;
   SCIP_Bool neginf;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert((SCIP_VARSTATUS) var->varstatus == SCIP_VARSTATUS_MULTAGGR);

   posinf = FALSE;
   neginf = FALSE;
   lb = var->data.multaggr.constant;
   for( i = var->data.multaggr.nvars-1 ; i >= 0 ; --i )
   {
      aggrvar = var->data.multaggr.vars[i];
      if( var->data.multaggr.scalars[i] > 0.0 )
      {
         bnd = SCIPvarGetStatus(aggrvar) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrLbLocal(aggrvar, set) : SCIPvarGetLbLocal(aggrvar);

         if( SCIPsetIsInfinity(set, bnd) )
            posinf = TRUE;
         else if( SCIPsetIsInfinity(set, -bnd) )
            neginf = TRUE;
         else
            lb += var->data.multaggr.scalars[i] * bnd;
      }
      else
      {
         bnd = SCIPvarGetStatus(aggrvar) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrUbLocal(aggrvar, set) : SCIPvarGetUbLocal(aggrvar);

         if( SCIPsetIsInfinity(set, -bnd) )
            posinf = TRUE;
         else if( SCIPsetIsInfinity(set, bnd) )
            neginf = TRUE;
         else
            lb += var->data.multaggr.scalars[i] * bnd;
      }

      /* stop if two diffrent infinities (or a -infinity) were found and return local lower bound of multi aggregated
       * variable
       */
      if( neginf )
         return SCIPvarGetLbLocal(var);
   }

   /* if positive infinity flag was set to true return infinity */
   if( posinf )
      return SCIPsetInfinity(set);

   return (MAX(lb, SCIPvarGetLbLocal(var))); /*lint !e666*/
}

/** for a multi-aggregated variable, gives the local upper bound computed by adding the local bounds from all
 *  aggregation variables, this upper bound may be tighter than the one given by SCIPvarGetUbLocal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
SCIP_Real SCIPvarGetMultaggrUbLocal(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;
   SCIP_Real ub;
   SCIP_Real bnd;
   SCIP_VAR* aggrvar;
   SCIP_Bool posinf;
   SCIP_Bool neginf;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert((SCIP_VARSTATUS) var->varstatus == SCIP_VARSTATUS_MULTAGGR);

   posinf = FALSE;
   neginf = FALSE;
   ub = var->data.multaggr.constant;
   for( i = var->data.multaggr.nvars-1 ; i >= 0 ; --i )
   {
      aggrvar = var->data.multaggr.vars[i];
      if( var->data.multaggr.scalars[i] > 0.0 )
      {
         bnd = SCIPvarGetStatus(aggrvar) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrUbLocal(aggrvar, set) : SCIPvarGetUbLocal(aggrvar);

         if( SCIPsetIsInfinity(set, bnd) )
            posinf = TRUE;
         else if( SCIPsetIsInfinity(set, -bnd) )
            neginf = TRUE;
         else
            ub += var->data.multaggr.scalars[i] * bnd;
      }
      else
      {
         bnd = SCIPvarGetStatus(aggrvar) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrLbLocal(aggrvar, set) : SCIPvarGetLbLocal(aggrvar);

         if( SCIPsetIsInfinity(set, -bnd) )
            posinf = TRUE;
         else if( SCIPsetIsInfinity(set, bnd) )
            neginf = TRUE;
         else
            ub += var->data.multaggr.scalars[i] * bnd;
      }

      /* stop if two diffrent infinities (or a -infinity) were found and return local upper bound of multi aggregated
       * variable
       */
      if( posinf )
         return SCIPvarGetUbLocal(var);
   }

   /* if negative infinity flag was set to true return -infinity */
   if( neginf )
      return -SCIPsetInfinity(set);

   return (MIN(ub, SCIPvarGetUbLocal(var))); /*lint !e666*/
}

/** for a multi-aggregated variable, gives the global lower bound computed by adding the global bounds from all
 *  aggregation variables, this global bound may be tighter than the one given by SCIPvarGetLbGlobal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
SCIP_Real SCIPvarGetMultaggrLbGlobal(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;
   SCIP_Real lb;
   SCIP_Real bnd;
   SCIP_VAR* aggrvar;
   SCIP_Bool posinf;
   SCIP_Bool neginf;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert((SCIP_VARSTATUS) var->varstatus == SCIP_VARSTATUS_MULTAGGR);

   posinf = FALSE;
   neginf = FALSE;
   lb = var->data.multaggr.constant;
   for( i = var->data.multaggr.nvars-1 ; i >= 0 ; --i )
   {
      aggrvar = var->data.multaggr.vars[i];
      if( var->data.multaggr.scalars[i] > 0.0 )
      {
         bnd = SCIPvarGetStatus(aggrvar) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrLbGlobal(aggrvar, set) : SCIPvarGetLbGlobal(aggrvar);

         if( SCIPsetIsInfinity(set, bnd) )
            posinf = TRUE;
         else if( SCIPsetIsInfinity(set, -bnd) )
            neginf = TRUE;
         else
            lb += var->data.multaggr.scalars[i] * bnd;
      }
      else
      {
         bnd = SCIPvarGetStatus(aggrvar) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrUbGlobal(aggrvar, set) : SCIPvarGetUbGlobal(aggrvar);

         if( SCIPsetIsInfinity(set, -bnd) )
            posinf = TRUE;
         else if( SCIPsetIsInfinity(set, bnd) )
            neginf = TRUE;
         else
            lb += var->data.multaggr.scalars[i] * bnd;
      }

      /* stop if two diffrent infinities (or a -infinity) were found and return global lower bound of multi aggregated
       * variable
       */
      if( neginf )
         return SCIPvarGetLbGlobal(var);
   }

   /* if positive infinity flag was set to true return infinity */
   if( posinf )
      return SCIPsetInfinity(set);

   return (MAX(lb, SCIPvarGetLbGlobal(var))); /*lint !e666*/
}

/** for a multi-aggregated variable, gives the global upper bound computed by adding the global bounds from all
 *  aggregation variables, this upper bound may be tighter than the one given by SCIPvarGetUbGlobal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
SCIP_Real SCIPvarGetMultaggrUbGlobal(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;
   SCIP_Real ub;
   SCIP_Real bnd;
   SCIP_VAR* aggrvar;
   SCIP_Bool posinf;
   SCIP_Bool neginf;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert((SCIP_VARSTATUS) var->varstatus == SCIP_VARSTATUS_MULTAGGR);

   posinf = FALSE;
   neginf = FALSE;
   ub = var->data.multaggr.constant;
   for( i = var->data.multaggr.nvars-1 ; i >= 0 ; --i )
   {
      aggrvar = var->data.multaggr.vars[i];
      if( var->data.multaggr.scalars[i] > 0.0 )
      {
         bnd = SCIPvarGetStatus(aggrvar) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrUbGlobal(aggrvar, set) : SCIPvarGetUbGlobal(aggrvar);

         if( SCIPsetIsInfinity(set, bnd) )
            posinf = TRUE;
         else if( SCIPsetIsInfinity(set, -bnd) )
            neginf = TRUE;
         else
            ub += var->data.multaggr.scalars[i] * bnd;
      }
      else
      {
         bnd = SCIPvarGetStatus(aggrvar) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrLbGlobal(aggrvar, set) : SCIPvarGetLbGlobal(aggrvar);

         if( SCIPsetIsInfinity(set, -bnd) )
            posinf = TRUE;
         else if( SCIPsetIsInfinity(set, bnd) )
            neginf = TRUE;
         else
            ub += var->data.multaggr.scalars[i] * bnd;
      }

      /* stop if two diffrent infinities (or a -infinity) were found and return local upper bound of multi aggregated
       * variable
       */
      if( posinf )
         return SCIPvarGetUbGlobal(var);
   }

   /* if negative infinity flag was set to true return -infinity */
   if( neginf )
      return -SCIPsetInfinity(set);

   return (MIN(ub, SCIPvarGetUbGlobal(var))); /*lint !e666*/
}

/** adds a hole to the original domain of the variable */
SCIP_RETCODE SCIPvarAddHoleOriginal(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right               /**< right bound of open interval in new hole */
   )
{
   SCIP_Bool added;

   assert(var != NULL);
   assert(!SCIPvarIsTransformed(var));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(set->stage == SCIP_STAGE_PROBLEM);

   SCIPsetDebugMsg(set, "adding original hole (%g,%g) to <%s>\n", left, right, var->name);

   if( SCIPsetIsEQ(set, left, right) )
      return SCIP_OKAY;

   /* the interval should not be empty */
   assert(SCIPsetIsLT(set, left, right));

   /* the the interval bound should already be adjusted */
   assert(SCIPsetIsEQ(set, left, adjustedUb(set, SCIPvarGetType(var), left)));
   assert(SCIPsetIsEQ(set, right, adjustedLb(set, SCIPvarGetType(var), right)));

   /* the the interval should lay between the lower and upper bound */
   assert(SCIPsetIsGE(set, left, SCIPvarGetLbOriginal(var)));
   assert(SCIPsetIsLE(set, right, SCIPvarGetUbOriginal(var)));

   /* add domain hole */
   SCIP_CALL( domAddHole(&var->data.original.origdom, blkmem, set, left, right, &added) );

   /* merges overlapping holes into single holes, moves bounds respectively if hole was added */
   if( added )
   {
      domMerge(&var->data.original.origdom, blkmem, set, NULL, NULL);
   }

   /**@todo add hole in parent and child variables (just like with bound changes);
    *       warning! original vars' holes are in original blkmem, transformed vars' holes in transformed blkmem
    */

   return SCIP_OKAY;
}

/** performs the current add of domain, changes all parents accordingly */
static
SCIP_RETCODE varProcessAddHoleGlobal(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right,              /**< right bound of open interval in new hole */
   SCIP_Bool*            added               /**< pointer to store whether the hole was added */
   )
{
   SCIP_VAR* parentvar;
   SCIP_Real newlb;
   SCIP_Real newub;
   int i;

   assert(var != NULL);
   assert(added != NULL);
   assert(blkmem != NULL);

   /* the interval should not be empty */
   assert(SCIPsetIsLT(set, left, right));

   /* the interval bound should already be adjusted */
   assert(SCIPsetIsEQ(set, left, adjustedUb(set, SCIPvarGetType(var), left)));
   assert(SCIPsetIsEQ(set, right, adjustedLb(set, SCIPvarGetType(var), right)));

   /* the interval should lay between the lower and upper bound */
   assert(SCIPsetIsGE(set, left, SCIPvarGetLbGlobal(var)));
   assert(SCIPsetIsLE(set, right, SCIPvarGetUbGlobal(var)));

   /* @todo add debugging mechanism for holes when using a debugging solution */

   /* add hole to hole list */
   SCIP_CALL( domAddHole(&var->glbdom, blkmem, set, left, right, added) );

   /* check if the hole is redundant */
   if( !(*added) )
      return SCIP_OKAY;

   /* current bounds */
   newlb = var->glbdom.lb;
   newub = var->glbdom.ub;

   /* merge domain holes */
   domMerge(&var->glbdom, blkmem, set, &newlb, &newub);

   /* the bound should not be changed */
   assert(SCIPsetIsEQ(set, newlb, var->glbdom.lb));
   assert(SCIPsetIsEQ(set, newub, var->glbdom.ub));

   /* issue bound change event */
   assert(SCIPvarIsTransformed(var) == (var->eventfilter != NULL));
   if( var->eventfilter != NULL )
   {
      SCIP_CALL( varEventGholeAdded(var, blkmem, set, eventqueue, left, right) );
   }

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      SCIP_Real parentnewleft;
      SCIP_Real parentnewright;
      SCIP_Bool localadded;

      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         parentnewleft = left;
         parentnewright = right;
         break;

      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         SCIPerrorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         return SCIP_INVALIDDATA;

      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);

         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change upper bound of x */
            parentnewleft = parentvar->data.aggregate.scalar * left + parentvar->data.aggregate.constant;
            parentnewright = parentvar->data.aggregate.scalar * right + parentvar->data.aggregate.constant;
         }
         else
         {
            /* a < 0 -> change lower bound of x */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));

            parentnewright = parentvar->data.aggregate.scalar * left + parentvar->data.aggregate.constant;
            parentnewleft = parentvar->data.aggregate.scalar * right + parentvar->data.aggregate.constant;
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x = offset - x'  ->  x' = offset - x */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);

         parentnewright = -left + parentvar->data.negate.constant;
         parentnewleft = -right + parentvar->data.negate.constant;
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }

      SCIPsetDebugMsg(set, "add global hole (%g,%g) to parent variable <%s>\n", parentnewleft, parentnewright, SCIPvarGetName(parentvar));

      /* perform hole added for parent variable */
      assert(blkmem != NULL);
      assert(SCIPsetIsLT(set, parentnewleft, parentnewright));
      SCIP_CALL( varProcessAddHoleGlobal(parentvar, blkmem, set, stat, eventqueue, 
            parentnewleft, parentnewright, &localadded) );
      assert(localadded);
   }

   return SCIP_OKAY;
}

/** adds a hole to the variable's global and local domain */
SCIP_RETCODE SCIPvarAddHoleGlobal(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right,              /**< right bound of open interval in new hole */
   SCIP_Bool*            added               /**< pointer to store whether the hole was added */
   )
{
   SCIP_Real childnewleft;
   SCIP_Real childnewright;

   assert(var != NULL);
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
   assert(blkmem != NULL);
   assert(added != NULL);

   SCIPsetDebugMsg(set, "adding global hole (%g,%g) to <%s>\n", left, right, var->name);

   /* the interval should not be empty */
   assert(SCIPsetIsLT(set, left, right));

   /* the the interval bound should already be adjusted */
   assert(SCIPsetIsEQ(set, left, adjustedUb(set, SCIPvarGetType(var), left)));
   assert(SCIPsetIsEQ(set, right, adjustedLb(set, SCIPvarGetType(var), right)));

   /* the the interval should lay between the lower and upper bound */
   assert(SCIPsetIsGE(set, left, SCIPvarGetLbGlobal(var)));
   assert(SCIPsetIsLE(set, right, SCIPvarGetUbGlobal(var)));

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar != NULL )
      {
         SCIP_CALL( SCIPvarAddHoleGlobal(var->data.original.transvar, blkmem, set, stat, eventqueue,
               left, right, added) );
      }
      else
      {
         assert(set->stage == SCIP_STAGE_PROBLEM);

         SCIP_CALL( varProcessAddHoleGlobal(var, blkmem, set, stat, eventqueue, left, right, added) );
         if( *added )
         {
            SCIP_Bool localadded;

            SCIP_CALL( SCIPvarAddHoleLocal(var, blkmem, set, stat, eventqueue, left, right, &localadded) );
         }
      }
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      SCIP_CALL( varProcessAddHoleGlobal(var, blkmem, set, stat, eventqueue, left, right, added) );
      if( *added )
      {
         SCIP_Bool localadded;

         SCIP_CALL( SCIPvarAddHoleLocal(var, blkmem, set, stat, eventqueue, left, right, &localadded) );
      }
      break;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot add hole of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);

      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change lower bound of y */
         childnewleft = (left - var->data.aggregate.constant)/var->data.aggregate.scalar;
         childnewright = (right - var->data.aggregate.constant)/var->data.aggregate.scalar;
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         childnewright = (left - var->data.aggregate.constant)/var->data.aggregate.scalar;
         childnewleft = (right - var->data.aggregate.constant)/var->data.aggregate.scalar;
      }
      else
      { 
         SCIPerrorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarAddHoleGlobal(var->data.aggregate.var, blkmem, set, stat, eventqueue, 
            childnewleft, childnewright, added) );
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot add a hole of a multi-aggregated variable.\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);

      childnewright = -left + var->data.negate.constant;
      childnewleft = -right + var->data.negate.constant;

      SCIP_CALL( SCIPvarAddHoleGlobal(var->negatedvar, blkmem, set, stat, eventqueue, 
            childnewleft, childnewright, added) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** performs the current add of domain, changes all parents accordingly */
static
SCIP_RETCODE varProcessAddHoleLocal(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right,              /**< right bound of open interval in new hole */
   SCIP_Bool*            added               /**< pointer to store whether the hole was added, or NULL */
   )
{
   SCIP_VAR* parentvar;
   SCIP_Real newlb;
   SCIP_Real newub;
   int i;

   assert(var != NULL);
   assert(added != NULL);
   assert(blkmem != NULL);

   /* the interval should not be empty */
   assert(SCIPsetIsLT(set, left, right));

   /* the the interval bound should already be adjusted */
   assert(SCIPsetIsEQ(set, left, adjustedUb(set, SCIPvarGetType(var), left)));
   assert(SCIPsetIsEQ(set, right, adjustedLb(set, SCIPvarGetType(var), right)));

   /* the the interval should lay between the lower and upper bound */
   assert(SCIPsetIsGE(set, left, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsLE(set, right, SCIPvarGetUbLocal(var)));

   /* add hole to hole list */
   SCIP_CALL( domAddHole(&var->locdom, blkmem, set, left, right, added) );

   /* check if the hole is redundant */
   if( !(*added) )
      return SCIP_OKAY;

   /* current bounds */
   newlb = var->locdom.lb;
   newub = var->locdom.ub;

   /* merge domain holes */
   domMerge(&var->locdom, blkmem, set, &newlb, &newub);

   /* the bound should not be changed */
   assert(SCIPsetIsEQ(set, newlb, var->locdom.lb));
   assert(SCIPsetIsEQ(set, newub, var->locdom.ub));

#if 0
   /* issue bound change event */
   assert(SCIPvarIsTransformed(var) == (var->eventfilter != NULL));
   if( var->eventfilter != NULL )
   {
      SCIP_CALL( varEventLholeAdded(var, blkmem, set, lp, branchcand, eventqueue, left, right) );
   }
#endif

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      SCIP_Real parentnewleft;
      SCIP_Real parentnewright;
      SCIP_Bool localadded;

      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         parentnewleft = left;
         parentnewright = right;
         break;

      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         SCIPerrorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         return SCIP_INVALIDDATA;

      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);

         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change upper bound of x */
            parentnewleft = parentvar->data.aggregate.scalar * left + parentvar->data.aggregate.constant;
            parentnewright = parentvar->data.aggregate.scalar * right + parentvar->data.aggregate.constant;
         }
         else
         {
            /* a < 0 -> change lower bound of x */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));

            parentnewright = parentvar->data.aggregate.scalar * left + parentvar->data.aggregate.constant;
            parentnewleft = parentvar->data.aggregate.scalar * right + parentvar->data.aggregate.constant;
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x = offset - x'  ->  x' = offset - x */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);

         parentnewright = -left + parentvar->data.negate.constant;
         parentnewleft = -right + parentvar->data.negate.constant;
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }

      SCIPsetDebugMsg(set, "add local hole (%g,%g) to parent variable <%s>\n", parentnewleft, parentnewright, SCIPvarGetName(parentvar));

      /* perform hole added for parent variable */
      assert(blkmem != NULL);
      assert(SCIPsetIsLT(set, parentnewleft, parentnewright));
      SCIP_CALL( varProcessAddHoleLocal(parentvar, blkmem, set, stat, eventqueue, 
            parentnewleft, parentnewright, &localadded) );
      assert(localadded);
   }

   return SCIP_OKAY;
}

/** adds a hole to the variable's current local domain */
SCIP_RETCODE SCIPvarAddHoleLocal(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right,              /**< right bound of open interval in new hole */
   SCIP_Bool*            added               /**< pointer to store whether the hole was added */
   )
{
   SCIP_Real childnewleft;
   SCIP_Real childnewright;

   assert(var != NULL);

   SCIPsetDebugMsg(set, "adding local hole (%g,%g) to <%s>\n", left, right, var->name);

   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
   assert(blkmem != NULL);
   assert(added != NULL);

   /* the interval should not be empty */
   assert(SCIPsetIsLT(set, left, right));

   /* the the interval bound should already be adjusted */
   assert(SCIPsetIsEQ(set, left, adjustedUb(set, SCIPvarGetType(var), left)));
   assert(SCIPsetIsEQ(set, right, adjustedLb(set, SCIPvarGetType(var), right)));

   /* the the interval should lay between the lower and upper bound */
   assert(SCIPsetIsGE(set, left, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsLE(set, right, SCIPvarGetUbLocal(var)));

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar != NULL )
      {
         SCIP_CALL( SCIPvarAddHoleLocal(var->data.original.transvar, blkmem, set, stat, eventqueue, 
               left, right, added) );
      }
      else
      {
         assert(set->stage == SCIP_STAGE_PROBLEM);
         SCIPstatIncrement(stat, set, domchgcount);
         SCIP_CALL( varProcessAddHoleLocal(var, blkmem, set, stat, eventqueue, left, right, added) );
      }
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      SCIPstatIncrement(stat, set, domchgcount);
      SCIP_CALL( varProcessAddHoleLocal(var, blkmem, set, stat, eventqueue, left, right, added) );
      break;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot add domain hole to a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);

      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change lower bound of y */
         childnewleft = (left - var->data.aggregate.constant)/var->data.aggregate.scalar;
         childnewright = (right - var->data.aggregate.constant)/var->data.aggregate.scalar;
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         childnewright = (left - var->data.aggregate.constant)/var->data.aggregate.scalar;
         childnewleft = (right - var->data.aggregate.constant)/var->data.aggregate.scalar;
      }
      else
      {         
         SCIPerrorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarAddHoleLocal(var->data.aggregate.var, blkmem, set, stat, eventqueue, 
            childnewleft, childnewright, added) );
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot add domain hole to a multi-aggregated variable.\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);

      childnewright = -left + var->data.negate.constant;
      childnewleft = -right + var->data.negate.constant;

      SCIP_CALL( SCIPvarAddHoleLocal(var->negatedvar, blkmem, set, stat, eventqueue, childnewleft, childnewright, added) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** resets the global and local bounds of original variable to their original values */
SCIP_RETCODE SCIPvarResetBounds(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(SCIPvarIsOriginal(var));
   /* resetting of bounds on original variables which have a transformed counterpart easily fails if, e.g.,
    * the transformed variable has been fixed */ 
   assert(SCIPvarGetTransVar(var) == NULL);

   /* copy the original bounds back to the global and local bounds */
   SCIP_CALL( SCIPvarChgLbGlobal(var, blkmem, set, stat, NULL, NULL, NULL, NULL, var->data.original.origdom.lb) );
   SCIP_CALL( SCIPvarChgUbGlobal(var, blkmem, set, stat, NULL, NULL, NULL, NULL, var->data.original.origdom.ub) );
   SCIP_CALL( SCIPvarChgLbLocal(var, blkmem, set, stat, NULL, NULL, NULL, var->data.original.origdom.lb) );
   SCIP_CALL( SCIPvarChgUbLocal(var, blkmem, set, stat, NULL, NULL, NULL, var->data.original.origdom.ub) );

   /* free the global and local holelists and duplicate the original ones */
   /**@todo this has also to be called recursively with methods similar to SCIPvarChgLbGlobal() */
   holelistFree(&var->glbdom.holelist, blkmem);
   holelistFree(&var->locdom.holelist, blkmem);
   SCIP_CALL( holelistDuplicate(&var->glbdom.holelist, blkmem, set, var->data.original.origdom.holelist) );
   SCIP_CALL( holelistDuplicate(&var->locdom.holelist, blkmem, set, var->data.original.origdom.holelist) );

   return SCIP_OKAY;
}

/** issues a IMPLADDED event on the given variable */
static
SCIP_RETCODE varEventImplAdded(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_EVENT* event;

   assert(var != NULL);

   /* issue IMPLADDED event on variable */
   SCIP_CALL( SCIPeventCreateImplAdded(&event, blkmem, var) );
   SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, NULL, &event) );

   return SCIP_OKAY;
}

/** actually performs the addition of a variable bound to the variable's vbound arrays */
static
SCIP_RETCODE varAddVbound(
   SCIP_VAR*             var,                /**< problem variable x in x <= b*z + d  or  x >= b*z + d */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_BOUNDTYPE        vbtype,             /**< type of variable bound (LOWER or UPPER) */
   SCIP_VAR*             vbvar,              /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             vbcoef,             /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             vbconstant          /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   )
{
   SCIP_Bool added;

   /* It can happen that the variable "var" and the variable "vbvar" are the same variable. For example if a variable
    * gets aggregated, the variable bounds (vbound) of that variable are copied to the other variable. A variable bound
    * variable of the aggregated variable might be the same as the one its gets aggregated too.
    * 
    * If the variable "var" and the variable "vbvar" are the same, the variable bound which should be added here has to
    * be redundant. This is the case since an infeasibility should have be detected in the previous methods. As well as
    * the bounds of the variable which should be also already be tightened in the previous methods. Therefore, the
    * variable bound can be ignored.
    *
    * From the way the the variable bound system is implemented (detecting infeasibility, tighten bounds), the
    * equivalence of the variables should be checked here.
    */
   if( var == vbvar )
   {
      /* in this case the variable bound has to be redundant, this means for possible assignments to this variable; this
       * can be checked via the global bounds of the variable */
#ifndef NDEBUG
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      if(vbtype == SCIP_BOUNDTYPE_LOWER)
      {
         if( vbcoef > 0.0 )
         {
            assert(SCIPsetIsGE(set, lb,  lb * vbcoef + vbconstant) );
            assert(SCIPsetIsGE(set, ub,  ub * vbcoef + vbconstant) );
         }
         else
         {
            assert(SCIPsetIsGE(set, lb,  ub * vbcoef + vbconstant) );
            assert(SCIPsetIsGE(set, ub,  lb * vbcoef + vbconstant) );
         }
      }
      else
      {
         assert(vbtype == SCIP_BOUNDTYPE_UPPER);
         if( vbcoef > 0.0 )
         {
            assert(SCIPsetIsLE(set, lb,  lb * vbcoef + vbconstant) );
            assert(SCIPsetIsLE(set, ub,  ub * vbcoef + vbconstant) );
         }
         else
         {
            assert(SCIPsetIsLE(set, lb,  ub * vbcoef + vbconstant) );
            assert(SCIPsetIsLE(set, ub,  lb * vbcoef + vbconstant) );
         }
      }
#endif
      SCIPsetDebugMsg(set, "redundant variable bound: <%s> %s %g<%s> %+g\n",
         SCIPvarGetName(var), vbtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", vbcoef, SCIPvarGetName(vbvar), vbconstant);

      return SCIP_OKAY;
   }

   SCIPsetDebugMsg(set, "adding variable bound: <%s> %s %g<%s> %+g\n",
      SCIPvarGetName(var), vbtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", vbcoef, SCIPvarGetName(vbvar), vbconstant);

   /* check variable bound on debugging solution */
   SCIP_CALL( SCIPdebugCheckVbound(set, var, vbtype, vbvar, vbcoef, vbconstant) ); /*lint !e506 !e774*/

   /* perform the addition */
   if( vbtype == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_CALL( SCIPvboundsAdd(&var->vlbs, blkmem, set, vbtype, vbvar, vbcoef, vbconstant, &added) );
   }
   else
   {
      SCIP_CALL( SCIPvboundsAdd(&var->vubs, blkmem, set, vbtype, vbvar, vbcoef, vbconstant, &added) );
   }
   var->closestvblpcount = -1;

   if( added )
   {
      /* issue IMPLADDED event */
      SCIP_CALL( varEventImplAdded(var, blkmem, set, eventqueue) );
   }

   return SCIP_OKAY;
}

/** checks whether the given implication is redundant or infeasible w.r.t. the implied variables global bounds */
static
void checkImplic(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool*            redundant,          /**< pointer to store whether the implication is redundant */
   SCIP_Bool*            infeasible          /**< pointer to store whether the implication is infeasible */
   )
{
   SCIP_Real impllb;
   SCIP_Real implub;

   assert(redundant != NULL);
   assert(infeasible != NULL);

   impllb = SCIPvarGetLbGlobal(implvar);
   implub = SCIPvarGetUbGlobal(implvar);
   if( impltype == SCIP_BOUNDTYPE_LOWER )
   {
      *infeasible = SCIPsetIsFeasGT(set, implbound, implub);
      *redundant = SCIPsetIsFeasLE(set, implbound, impllb);
   }
   else
   {
      *infeasible = SCIPsetIsFeasLT(set, implbound, impllb);
      *redundant = SCIPsetIsFeasGE(set, implbound, implub);
   }
}

/** applies the given implication, if it is not redundant */
static
SCIP_RETCODE applyImplic(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   )
{
   SCIP_Real implub;
   SCIP_Real impllb;

   assert(infeasible != NULL);

   *infeasible = FALSE;

   implub = SCIPvarGetUbGlobal(implvar);
   impllb = SCIPvarGetLbGlobal(implvar);
   if( impltype == SCIP_BOUNDTYPE_LOWER )
   {
      if( SCIPsetIsFeasGT(set, implbound, implub) )
      {
         /* the implication produces a conflict: the problem is infeasible */
         *infeasible = TRUE;
      }
      else if( SCIPsetIsFeasGT(set, implbound, impllb) )
      {
         /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
          * with the local bound, in this case we need to store the bound change as pending bound change
          */
         if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
         {
            assert(tree != NULL);
            assert(transprob != NULL);
            assert(SCIPprobIsTransformed(transprob));

            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, implvar, implbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
         }
         else
         {
            SCIP_CALL( SCIPvarChgLbGlobal(implvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, implbound) );
         }

         if( nbdchgs != NULL )
            (*nbdchgs)++;
      }
   }
   else
   {
      if( SCIPsetIsFeasLT(set, implbound, impllb) )
      {
         /* the implication produces a conflict: the problem is infeasible */
         *infeasible = TRUE;
      }
      else if( SCIPsetIsFeasLT(set, implbound, implub) )
      {
         /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
          * with the local bound, in this case we need to store the bound change as pending bound change
          */
         if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
         {
            assert(tree != NULL);
            assert(transprob != NULL);
            assert(SCIPprobIsTransformed(transprob));

            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, implvar, implbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
         else
         {
            SCIP_CALL( SCIPvarChgUbGlobal(implvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, implbound) );
         }

         if( nbdchgs != NULL )
            (*nbdchgs)++;
      }
   }

   return SCIP_OKAY;
}

/** actually performs the addition of an implication to the variable's implication arrays,
 *  and adds the corresponding implication or variable bound to the implied variable;
 *  if the implication is conflicting, the variable is fixed to the opposite value;
 *  if the variable is already fixed to the given value, the implication is performed immediately;
 *  if the implication is redundant with respect to the variables' global bounds, it is ignored
 */
static
SCIP_RETCODE varAddImplic(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool             isshortcut,         /**< is the implication a shortcut, i.e., added as part of the transitive closure of another implication? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs,            /**< pointer to count the number of performed bound changes, or NULL */
   SCIP_Bool*            added               /**< pointer to store whether an implication was added */
   )
{
   SCIP_Bool redundant;
   SCIP_Bool conflict;

   assert(var != NULL);
   assert(SCIPvarIsActive(var));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(SCIPvarIsActive(implvar) || SCIPvarGetStatus(implvar) == SCIP_VARSTATUS_FIXED);
   assert(infeasible != NULL);
   assert(added != NULL);

   /* check implication on debugging solution */
   SCIP_CALL( SCIPdebugCheckImplic(set, var, varfixing, implvar, impltype, implbound) ); /*lint !e506 !e774*/

   *infeasible = FALSE;
   *added = FALSE;

   /* check, if the implication is redundant or infeasible */
   checkImplic(set, implvar, impltype, implbound, &redundant, &conflict);
   assert(!redundant || !conflict);
   if( redundant )
      return SCIP_OKAY;

   if( var == implvar )
   {
      /* special cases appear were a bound to a variable implies itself to be outside the bounds:
       * x == varfixing  =>  x < 0 or x > 1
       */
      if( SCIPsetIsLT(set, implbound, 0.0) || SCIPsetIsGT(set, implbound, 1.0) )
         conflict = TRUE;
      else
      {
         /* variable implies itself: x == varfixing  =>  x == (impltype == SCIP_BOUNDTYPE_LOWER) */
         assert(SCIPsetIsZero(set, implbound) || SCIPsetIsEQ(set, implbound, 1.0));
         assert(SCIPsetIsZero(set, implbound) == (impltype == SCIP_BOUNDTYPE_UPPER));
         assert(SCIPsetIsEQ(set, implbound, 1.0) == (impltype == SCIP_BOUNDTYPE_LOWER));
         conflict = conflict || ((varfixing == TRUE) == (impltype == SCIP_BOUNDTYPE_UPPER));
         if( !conflict )
            return SCIP_OKAY;
      }
   }

   /* check, if the variable is already fixed */
   if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5 )
   {
      /* if the variable is fixed to the given value, perform the implication; otherwise, ignore the implication */
      if( varfixing == (SCIPvarGetLbGlobal(var) > 0.5) )
      {
         SCIP_CALL( applyImplic(blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               cliquetable, implvar, impltype, implbound, infeasible, nbdchgs) );
      }
      return SCIP_OKAY;
   }

   assert((impltype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsGT(set, implbound, SCIPvarGetLbGlobal(implvar)))
      || (impltype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsLT(set, implbound, SCIPvarGetUbGlobal(implvar))));

   if( !conflict )
   {
      assert(SCIPvarIsActive(implvar)); /* a fixed implvar would either cause a redundancy or infeasibility */

      if( SCIPvarIsBinary(implvar) )
      {
         SCIP_VAR* vars[2];
         SCIP_Bool vals[2];

         assert(SCIPsetIsFeasEQ(set, implbound, 1.0) || SCIPsetIsFeasZero(set, implbound));
         assert((impltype == SCIP_BOUNDTYPE_UPPER) == SCIPsetIsFeasZero(set, implbound));

         vars[0] = var;
         vars[1] = implvar;
         vals[0] = varfixing;
         vals[1] = (impltype == SCIP_BOUNDTYPE_UPPER);

         /* add the clique to the clique table */
         SCIP_CALL( SCIPcliquetableAdd(cliquetable, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
               eventqueue, vars, vals, 2, FALSE, &conflict, nbdchgs) );

         if( !conflict )
            return SCIP_OKAY;
      }
      else
      {
         /* add implication x == 0/1 -> y <= b / y >= b to the implications list of x */
         SCIPsetDebugMsg(set, "adding implication: <%s> == %u  ==>  <%s> %s %g\n",
            SCIPvarGetName(var), varfixing,
            SCIPvarGetName(implvar), impltype == SCIP_BOUNDTYPE_UPPER ? "<=" : ">=", implbound);
         SCIP_CALL( SCIPimplicsAdd(&var->implics, blkmem, set, stat, varfixing, implvar, impltype, implbound,
               isshortcut, &conflict, added) );
      }
   }
   assert(!conflict || !(*added));

   /* on conflict, fix the variable to the opposite value */
   if( conflict )
   {
      SCIPsetDebugMsg(set, " -> implication yields a conflict: fix <%s> == %d\n", SCIPvarGetName(var), !varfixing);

      /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
       * with the local bound, in this case we need to store the bound change as pending bound change
       */
      if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
      {
         assert(tree != NULL);
         assert(transprob != NULL);
         assert(SCIPprobIsTransformed(transprob));

         if( varfixing )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, var, 0.0, SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
         else
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, var, 1.0, SCIP_BOUNDTYPE_LOWER, FALSE) );
         }
      }
      else
      {
         if( varfixing )
         {
            SCIP_CALL( SCIPvarChgUbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, 0.0) );
         }
         else
         {
            SCIP_CALL( SCIPvarChgLbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, 1.0) );
         }
      }
      if( nbdchgs != NULL )
         (*nbdchgs)++;

      return SCIP_OKAY;
   }
   else if( *added )
   {
      /* issue IMPLADDED event */
      SCIP_CALL( varEventImplAdded(var, blkmem, set, eventqueue) );
   }
   else
   {
      /* the implication was redundant: the inverse is also redundant */
      return SCIP_OKAY;
   }

   assert(SCIPvarIsActive(implvar)); /* a fixed implvar would either cause a redundancy or infeasibility */

   /* check, whether implied variable is binary */
   if( !SCIPvarIsBinary(implvar) )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      /* add inverse variable bound to the variable bounds of y with global bounds y \in [lb,ub]:
       *   x == 0 -> y <= b  <->  y <= (ub - b)*x + b
       *   x == 1 -> y <= b  <->  y <= (b - ub)*x + ub
       *   x == 0 -> y >= b  <->  y >= (lb - b)*x + b
       *   x == 1 -> y >= b  <->  y >= (b - lb)*x + lb
       * for numerical reasons, ignore variable bounds with large absolute coefficient
       */
      lb = SCIPvarGetLbGlobal(implvar);
      ub = SCIPvarGetUbGlobal(implvar);
      if( impltype == SCIP_BOUNDTYPE_UPPER )
      {
         if( REALABS(implbound - ub) <= MAXABSVBCOEF ) 
         {
            SCIP_CALL( varAddVbound(implvar, blkmem, set, eventqueue, SCIP_BOUNDTYPE_UPPER, var,
                  varfixing ? implbound - ub : ub - implbound, varfixing ? ub : implbound) );
         }
      }
      else
      {
         if( REALABS(implbound - lb) <= MAXABSVBCOEF ) 
         {
            SCIP_CALL( varAddVbound(implvar, blkmem, set, eventqueue, SCIP_BOUNDTYPE_LOWER, var,
                  varfixing ? implbound - lb : lb - implbound, varfixing ? lb : implbound) );
         }
      }
   }

   return SCIP_OKAY;
}

/** adds transitive closure for binary implication x = a -> y = b */
static
SCIP_RETCODE varAddTransitiveBinaryClosureImplic(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_Bool             implvarfixing,      /**< fixing b in implication */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   )
{
   SCIP_VAR** implvars;
   SCIP_BOUNDTYPE* impltypes;
   SCIP_Real* implbounds;
   int nimpls;
   int i;

   *infeasible = FALSE;

   /* binary variable: implications of implvar */
   nimpls = SCIPimplicsGetNImpls(implvar->implics, implvarfixing);
   implvars = SCIPimplicsGetVars(implvar->implics, implvarfixing);
   impltypes = SCIPimplicsGetTypes(implvar->implics, implvarfixing);
   implbounds = SCIPimplicsGetBounds(implvar->implics, implvarfixing);

   /* if variable has too many implications, the implication graph may become too dense */
   i = MIN(nimpls, MAXIMPLSCLOSURE) - 1;

   /* we have to iterate from back to front, because in varAddImplic() it may happen that a conflict is detected and
    * implvars[i] is fixed, s.t. the implication y == varfixing -> z <= b / z >= b is deleted; this affects the
    * array over which we currently iterate; the only thing that can happen, is that elements of the array are
    * deleted; in this case, the subsequent elements are moved to the front; if we iterate from back to front, the
    * only thing that can happen is that we add the same implication twice - this does no harm
    */
   while ( i >= 0 && !(*infeasible) )
   {
      SCIP_Bool added;

      assert(implvars[i] != implvar);

      /* we have x == varfixing -> y == implvarfixing -> z <= b / z >= b:
       * add implication x == varfixing -> z <= b / z >= b to the implications list of x 
       */
      if( SCIPvarIsActive(implvars[i]) )
      {
         SCIP_CALL( varAddImplic(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable, branchcand,
               eventqueue, varfixing, implvars[i], impltypes[i], implbounds[i], TRUE, infeasible, nbdchgs, &added) );
         assert(SCIPimplicsGetNImpls(implvar->implics, implvarfixing) <= nimpls);
         nimpls = SCIPimplicsGetNImpls(implvar->implics, implvarfixing);
         i = MIN(i, nimpls); /* some elements from the array could have been removed */
      }
      --i;
   }

   return SCIP_OKAY;
}

/** adds given implication to the variable's implication list, and adds all implications directly implied by this
 *  implication to the variable's implication list;
 *  if the implication is conflicting, the variable is fixed to the opposite value;
 *  if the variable is already fixed to the given value, the implication is performed immediately;
 *  if the implication is redundant with respect to the variables' global bounds, it is ignored
 */
static
SCIP_RETCODE varAddTransitiveImplic(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool             transitive,         /**< should transitive closure of implication also be added? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   )
{
   SCIP_Bool added;

   assert(var != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(SCIPvarIsActive(var));
   assert(implvar != NULL);
   assert(SCIPvarIsActive(implvar) || SCIPvarGetStatus(implvar) == SCIP_VARSTATUS_FIXED);
   assert(infeasible != NULL);

   /* add implication x == varfixing -> y <= b / y >= b to the implications list of x */
   SCIP_CALL( varAddImplic(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable, branchcand,
         eventqueue, varfixing, implvar, impltype, implbound, FALSE, infeasible, nbdchgs, &added) );

   if( *infeasible || var == implvar || !transitive || !added )
      return SCIP_OKAY;

   assert(SCIPvarIsActive(implvar)); /* a fixed implvar would either cause a redundancy or infeasibility */

   /* add transitive closure */
   if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
   {
      SCIP_Bool implvarfixing;

      implvarfixing = (impltype == SCIP_BOUNDTYPE_LOWER);

      /* binary variable: implications of implvar */
      SCIP_CALL( varAddTransitiveBinaryClosureImplic(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
            cliquetable, branchcand, eventqueue, varfixing, implvar, implvarfixing, infeasible, nbdchgs) );

      /* inverse implication */
      if( !(*infeasible) )
      {
         SCIP_CALL( varAddTransitiveBinaryClosureImplic(implvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
               cliquetable, branchcand, eventqueue, !implvarfixing, var, !varfixing, infeasible, nbdchgs) );
      }
   }
   else
   {
      /* non-binary variable: variable lower bounds of implvar */
      if( impltype == SCIP_BOUNDTYPE_UPPER && implvar->vlbs != NULL )
      {
         SCIP_VAR** vlbvars;
         SCIP_Real* vlbcoefs;
         SCIP_Real* vlbconstants;
         int nvlbvars;
         int i;

         nvlbvars = SCIPvboundsGetNVbds(implvar->vlbs);
         vlbvars = SCIPvboundsGetVars(implvar->vlbs);
         vlbcoefs = SCIPvboundsGetCoefs(implvar->vlbs);
         vlbconstants = SCIPvboundsGetConstants(implvar->vlbs);

         /* we have to iterate from back to front, because in varAddImplic() it may happen that a conflict is detected and
          * vlbvars[i] is fixed, s.t. the variable bound is deleted; this affects the array over which we currently
          * iterate; the only thing that can happen, is that elements of the array are deleted; in this case, the
          * subsequent elements are moved to the front; if we iterate from back to front, the only thing that can happen
          * is that we add the same implication twice - this does no harm
          */
         i = nvlbvars-1;
         while ( i >= 0 && !(*infeasible) )
         {
            assert(vlbvars[i] != implvar);
            assert(!SCIPsetIsZero(set, vlbcoefs[i]));

            /* we have x == varfixing -> y <= b and y >= c*z + d:
             *   c > 0: add implication x == varfixing -> z <= (b-d)/c to the implications list of x
             *   c < 0: add implication x == varfixing -> z >= (b-d)/c to the implications list of x
             *
             * @note during an aggregation the aggregated variable "aggrvar" (the one which will have the status
             *       SCIP_VARSTATUS_AGGREGATED afterwards) copies its variable lower and uppers bounds to the
             *       aggregation variable (the one which will stay active); 
             *
             *       W.l.o.g. we consider the variable upper bounds for now. Let "vubvar" be a variable upper bound of
             *       the aggregated variable "aggvar"; During that copying of that variable upper bound variable
             *       "vubvar" the variable lower and upper bounds of this variable "vubvar" are also considered; note
             *       that the "aggvar" can be a variable lower bound variable of the variable "vubvar"; Due to that
             *       situation it can happen that we reach that code place where "vlbvars[i] == aggvar". In particular
             *       the "aggvar" has already the variable status SCIP_VARSTATUS_AGGREGATED or SCIP_VARSTATUS_NEGATED
             *       but is still active since the aggregation is not finished yet (in SCIPvarAggregate()); therefore we
             *       have to explicitly check that the active variable has not a variable status
             *       SCIP_VARSTATUS_AGGREGATED or SCIP_VARSTATUS_NEGATED;
             */
            if( SCIPvarIsActive(vlbvars[i]) && SCIPvarGetStatus(vlbvars[i]) != SCIP_VARSTATUS_AGGREGATED && SCIPvarGetStatus(vlbvars[i]) != SCIP_VARSTATUS_NEGATED )
            {
               SCIP_Real vbimplbound;

               vbimplbound = (implbound - vlbconstants[i])/vlbcoefs[i];
               if( vlbcoefs[i] >= 0.0 )
               {
                  vbimplbound = adjustedUb(set, SCIPvarGetType(vlbvars[i]), vbimplbound);
                  SCIP_CALL( varAddImplic(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable,
                        branchcand, eventqueue, varfixing, vlbvars[i], SCIP_BOUNDTYPE_UPPER, vbimplbound, TRUE,
                        infeasible, nbdchgs, &added) );
               }
               else
               {
                  vbimplbound = adjustedLb(set, SCIPvarGetType(vlbvars[i]), vbimplbound);
                  SCIP_CALL( varAddImplic(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable,
                        branchcand, eventqueue, varfixing, vlbvars[i], SCIP_BOUNDTYPE_LOWER, vbimplbound, TRUE,
                        infeasible, nbdchgs, &added) );
               }
               nvlbvars = SCIPvboundsGetNVbds(implvar->vlbs);
               i = MIN(i, nvlbvars); /* some elements from the array could have been removed */
            }
            --i;
         }
      }

      /* non-binary variable: variable upper bounds of implvar */
      if( impltype == SCIP_BOUNDTYPE_LOWER && implvar->vubs != NULL )
      {
         SCIP_VAR** vubvars;
         SCIP_Real* vubcoefs;
         SCIP_Real* vubconstants;
         int nvubvars;
         int i;

         nvubvars = SCIPvboundsGetNVbds(implvar->vubs);
         vubvars = SCIPvboundsGetVars(implvar->vubs);
         vubcoefs = SCIPvboundsGetCoefs(implvar->vubs);
         vubconstants = SCIPvboundsGetConstants(implvar->vubs);

         /* we have to iterate from back to front, because in varAddImplic() it may happen that a conflict is detected and
          * vubvars[i] is fixed, s.t. the variable bound is deleted; this affects the array over which we currently
          * iterate; the only thing that can happen, is that elements of the array are deleted; in this case, the
          * subsequent elements are moved to the front; if we iterate from back to front, the only thing that can happen
          * is that we add the same implication twice - this does no harm
          */
         i = nvubvars-1;
         while ( i >= 0 && !(*infeasible) )
         {
            assert(vubvars[i] != implvar);
            assert(!SCIPsetIsZero(set, vubcoefs[i]));

            /* we have x == varfixing -> y >= b and y <= c*z + d:
             *   c > 0: add implication x == varfixing -> z >= (b-d)/c to the implications list of x
             *   c < 0: add implication x == varfixing -> z <= (b-d)/c to the implications list of x
             *
             * @note during an aggregation the aggregated variable "aggrvar" (the one which will have the status
             *       SCIP_VARSTATUS_AGGREGATED afterwards) copies its variable lower and uppers bounds to the
             *       aggregation variable (the one which will stay active); 
             *
             *       W.l.o.g. we consider the variable lower bounds for now. Let "vlbvar" be a variable lower bound of
             *       the aggregated variable "aggvar"; During that copying of that variable lower bound variable
             *       "vlbvar" the variable lower and upper bounds of this variable "vlbvar" are also considered; note
             *       that the "aggvar" can be a variable upper bound variable of the variable "vlbvar"; Due to that
             *       situation it can happen that we reach that code place where "vubvars[i] == aggvar". In particular
             *       the "aggvar" has already the variable status SCIP_VARSTATUS_AGGREGATED or SCIP_VARSTATUS_NEGATED
             *       but is still active since the aggregation is not finished yet (in SCIPvarAggregate()); therefore we
             *       have to explicitly check that the active variable has not a variable status
             *       SCIP_VARSTATUS_AGGREGATED or SCIP_VARSTATUS_NEGATED;
             */
            if( SCIPvarIsActive(vubvars[i]) && SCIPvarGetStatus(vubvars[i]) != SCIP_VARSTATUS_AGGREGATED && SCIPvarGetStatus(vubvars[i]) != SCIP_VARSTATUS_NEGATED )
            {
               SCIP_Real vbimplbound;

               vbimplbound = (implbound - vubconstants[i])/vubcoefs[i];
               if( vubcoefs[i] >= 0.0 )
               {
                  vbimplbound = adjustedLb(set, SCIPvarGetType(vubvars[i]), vbimplbound);
                  SCIP_CALL( varAddImplic(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable,
                        branchcand, eventqueue, varfixing, vubvars[i], SCIP_BOUNDTYPE_LOWER, vbimplbound, TRUE,
                        infeasible, nbdchgs, &added) );
               }
               else
               {
                  vbimplbound = adjustedUb(set, SCIPvarGetType(vubvars[i]), vbimplbound);
                  SCIP_CALL( varAddImplic(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable,
                        branchcand, eventqueue, varfixing, vubvars[i], SCIP_BOUNDTYPE_UPPER, vbimplbound, TRUE,
                        infeasible, nbdchgs, &added) );
               }
               nvubvars = SCIPvboundsGetNVbds(implvar->vubs);
               i = MIN(i, nvubvars); /* some elements from the array could have been removed */
            }
            --i;
         }
      }
   }

   return SCIP_OKAY;
}

/** informs variable x about a globally valid variable lower bound x >= b*z + d with integer variable z;
 *  if z is binary, the corresponding valid implication for z is also added;
 *  improves the global bounds of the variable and the vlb variable if possible
 */
SCIP_RETCODE SCIPvarAddVlb(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             vlbvar,             /**< variable z    in x >= b*z + d */
   SCIP_Real             vlbcoef,            /**< coefficient b in x >= b*z + d */
   SCIP_Real             vlbconstant,        /**< constant d    in x >= b*z + d */
   SCIP_Bool             transitive,         /**< should transitive closure of implication also be added? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(SCIPvarGetType(vlbvar) != SCIP_VARTYPE_CONTINUOUS);
   assert(infeasible != NULL);

   SCIPsetDebugMsg(set, "adding variable lower bound <%s> >= %g<%s> + %g\n", SCIPvarGetName(var), vlbcoef, SCIPvarGetName(vlbvar), vlbconstant);

   *infeasible = FALSE;
   if( nbdchgs != NULL )
      *nbdchgs = 0;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      SCIP_CALL( SCIPvarAddVlb(var->data.original.transvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
            cliquetable, branchcand, eventqueue, vlbvar, vlbcoef, vlbconstant, transitive, infeasible, nbdchgs) );
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      /* transform b*z + d into the corresponding sum after transforming z to an active problem variable */
      SCIP_CALL( SCIPvarGetProbvarSum(&vlbvar, set, &vlbcoef, &vlbconstant) );
      SCIPsetDebugMsg(set, " -> transformed to variable lower bound <%s> >= %g<%s> + %g\n", SCIPvarGetName(var), vlbcoef, SCIPvarGetName(vlbvar), vlbconstant);

      /* if the vlb coefficient is zero, just update the lower bound of the variable */
      if( SCIPsetIsZero(set, vlbcoef) )
      {
         if( SCIPsetIsFeasGT(set, vlbconstant, SCIPvarGetUbGlobal(var)) )
            *infeasible = TRUE;
         else if( SCIPsetIsFeasGT(set, vlbconstant, SCIPvarGetLbGlobal(var)) )
         {
            /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
             * with the local bound, in this case we need to store the bound change as pending bound change
             */
            if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
            {
               assert(tree != NULL);
               assert(transprob != NULL);
               assert(SCIPprobIsTransformed(transprob));

               SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                     tree, reopt, lp, branchcand, eventqueue, cliquetable, var, vlbconstant, SCIP_BOUNDTYPE_LOWER, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPvarChgLbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, vlbconstant) );
            }

            if( nbdchgs != NULL )
               (*nbdchgs)++;
         }
      }
      else if( SCIPvarIsActive(vlbvar) )
      {
         SCIP_Real xlb;
         SCIP_Real xub;
         SCIP_Real zlb;
         SCIP_Real zub;
         SCIP_Real minvlb;
         SCIP_Real maxvlb;

         assert(SCIPvarGetStatus(vlbvar) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(vlbvar) == SCIP_VARSTATUS_COLUMN);
         assert(vlbcoef != 0.0);

         minvlb = -SCIPsetInfinity(set);
         maxvlb = SCIPsetInfinity(set);

         xlb = SCIPvarGetLbGlobal(var);
         xub = SCIPvarGetUbGlobal(var);
         zlb = SCIPvarGetLbGlobal(vlbvar);
         zub = SCIPvarGetUbGlobal(vlbvar);

         /* improve global bounds of vlb variable, and calculate minimal and maximal value of variable bound */
         if( vlbcoef >= 0.0 )
         {
            SCIP_Real newzub;

            if( !SCIPsetIsInfinity(set, xub) )
            {
               /* x >= b*z + d  ->  z <= (x-d)/b */
               newzub = (xub - vlbconstant)/vlbcoef;

               /* return if the new bound is less than -infinity */
               if( SCIPsetIsInfinity(set, REALABS(newzub)) )
                  return SCIP_OKAY;

               if( SCIPsetIsFeasLT(set, newzub, zlb) )
               {
                  *infeasible = TRUE;
                  return SCIP_OKAY;
               }
               if( SCIPsetIsFeasLT(set, newzub, zub) )
               {
                  /* bound might be adjusted due to integrality condition */
                  newzub = adjustedUb(set, SCIPvarGetType(vlbvar), newzub);

                  /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
                   * with the local bound, in this case we need to store the bound change as pending bound change
                   */
                  if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
                  {
                     assert(tree != NULL);
                     assert(transprob != NULL);
                     assert(SCIPprobIsTransformed(transprob));

                     SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                           tree, reopt, lp, branchcand, eventqueue, cliquetable, vlbvar, newzub, SCIP_BOUNDTYPE_UPPER, FALSE) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPvarChgUbGlobal(vlbvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newzub) );
                  }
                  zub = newzub;

                  if( nbdchgs != NULL )
                     (*nbdchgs)++;
               }
               maxvlb = vlbcoef * zub + vlbconstant;
               if( !SCIPsetIsInfinity(set, -zlb) )
                  minvlb = vlbcoef * zlb + vlbconstant;
            }
            else
            {
               if( !SCIPsetIsInfinity(set, zub) )
                  maxvlb = vlbcoef * zub + vlbconstant;
               if( !SCIPsetIsInfinity(set, -zlb) )
                  minvlb = vlbcoef * zlb + vlbconstant;
            }
         }
         else
         {
            SCIP_Real newzlb;

            if( !SCIPsetIsInfinity(set, xub) )
            {
               /* x >= b*z + d  ->  z >= (x-d)/b */
               newzlb = (xub - vlbconstant)/vlbcoef;

               /* return if the new bound is larger than infinity */
               if( SCIPsetIsInfinity(set, REALABS(newzlb)) )
                  return SCIP_OKAY;

               if( SCIPsetIsFeasGT(set, newzlb, zub) )
               {
                  *infeasible = TRUE;
                  return SCIP_OKAY;
               }
               if( SCIPsetIsFeasGT(set, newzlb, zlb) )
               {
                  /* bound might be adjusted due to integrality condition */
                  newzlb = adjustedLb(set, SCIPvarGetType(vlbvar), newzlb);

                  /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
                   * with the local bound, in this case we need to store the bound change as pending bound change
                   */
                  if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
                  {
                     assert(tree != NULL);
                     assert(transprob != NULL);
                     assert(SCIPprobIsTransformed(transprob));

                     SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                           tree, reopt, lp, branchcand, eventqueue, cliquetable, vlbvar, newzlb, SCIP_BOUNDTYPE_LOWER, FALSE) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPvarChgLbGlobal(vlbvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newzlb) );
                  }
                  zlb = newzlb;

                  if( nbdchgs != NULL )
                     (*nbdchgs)++;
               }
               maxvlb = vlbcoef * zlb + vlbconstant;
               if( !SCIPsetIsInfinity(set, zub) )
                  minvlb = vlbcoef * zub + vlbconstant;
            }
            else
            {
               if( !SCIPsetIsInfinity(set, -zlb) )
                  maxvlb = vlbcoef * zlb + vlbconstant;
               if( !SCIPsetIsInfinity(set, zub) )
                  minvlb = vlbcoef * zub + vlbconstant;
            }
         }
         if( maxvlb < minvlb )
            maxvlb = minvlb;

         /* adjust bounds due to integrality of variable */
         minvlb = adjustedLb(set, SCIPvarGetType(var), minvlb);
         maxvlb = adjustedLb(set, SCIPvarGetType(var), maxvlb);

         /* check bounds for feasibility */
         if( SCIPsetIsFeasGT(set, minvlb, xub) || (var == vlbvar && SCIPsetIsEQ(set, vlbcoef, 1.0) && SCIPsetIsFeasPositive(set, vlbconstant))  )
         {
            *infeasible = TRUE;
            return SCIP_OKAY;
         }
         /* improve global lower bound of variable */
         if( SCIPsetIsFeasGT(set, minvlb, xlb) )
         {
            /* bound might be adjusted due to integrality condition */
            minvlb = adjustedLb(set, SCIPvarGetType(var), minvlb);

            /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
             * with the local bound, in this case we need to store the bound change as pending bound change
             */
            if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
            {
               assert(tree != NULL);
               assert(transprob != NULL);
               assert(SCIPprobIsTransformed(transprob));

               SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                     tree, reopt, lp, branchcand, eventqueue, cliquetable, var, minvlb, SCIP_BOUNDTYPE_LOWER, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPvarChgLbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, minvlb) );
            }
            xlb = minvlb;

            if( nbdchgs != NULL )
               (*nbdchgs)++;
         }
         minvlb = xlb;

         /* improve variable bound for binary z by moving the variable's global bound to the vlb constant */
         if( SCIPvarGetType(vlbvar) == SCIP_VARTYPE_BINARY )
         {
            /* b > 0: x >= (maxvlb - minvlb) * z + minvlb
             * b < 0: x >= (minvlb - maxvlb) * z + maxvlb
             */

            assert(!SCIPsetIsInfinity(set, maxvlb) && !SCIPsetIsInfinity(set, -minvlb));

            if( vlbcoef >= 0.0 )
            {
               vlbcoef = maxvlb - minvlb;
               vlbconstant = minvlb;
            }
            else
            {
               vlbcoef = minvlb - maxvlb;
               vlbconstant = maxvlb;
            }
         }

         /* add variable bound to the variable bounds list */
         if( SCIPsetIsFeasGT(set, maxvlb, xlb) )
         {
            assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_FIXED);
            assert(!SCIPsetIsZero(set, vlbcoef));

            /* if one of the variables is binary, add the corresponding implication to the variable's implication
             * list, thereby also adding the variable bound (or implication) to the other variable
             */
            if( SCIPvarGetType(vlbvar) == SCIP_VARTYPE_BINARY )
            {
               /* add corresponding implication:
                *   b > 0, x >= b*z + d  <->  z == 1 -> x >= b+d
                *   b < 0, x >= b*z + d  <->  z == 0 -> x >= d
                */
               SCIP_CALL( varAddTransitiveImplic(vlbvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                     cliquetable, branchcand, eventqueue, (vlbcoef >= 0.0), var, SCIP_BOUNDTYPE_LOWER, maxvlb, transitive,
                     infeasible, nbdchgs) );
            }
            else if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
            {
               /* add corresponding implication:
                *   b > 0, x >= b*z + d  <->  x == 0 -> z <= -d/b
                *   b < 0, x >= b*z + d  <->  x == 0 -> z >= -d/b
                */
               SCIP_Real implbound;
               implbound = -vlbconstant/vlbcoef;

               /* tighten the implication bound if the variable is integer */
               if( SCIPvarIsIntegral(vlbvar) )
               {
                  if( vlbcoef >= 0 )
                     implbound = SCIPsetFloor(set, implbound);
                  else
                     implbound = SCIPsetCeil(set, implbound);
               }
               SCIP_CALL( varAddTransitiveImplic(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                     cliquetable, branchcand, eventqueue, FALSE, vlbvar, (vlbcoef >= 0.0 ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER),
                     implbound, transitive, infeasible, nbdchgs) );
            }
            else
            {
               SCIP_CALL( varAddVbound(var, blkmem, set, eventqueue, SCIP_BOUNDTYPE_LOWER, vlbvar, vlbcoef, vlbconstant) );
            }
         }
      }
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      /* x = a*y + c:  x >= b*z + d  <=>  a*y + c >= b*z + d  <=>  y >= b/a * z + (d-c)/a, if a > 0
       *                                                           y <= b/a * z + (d-c)/a, if a < 0
       */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> add variable lower bound */
         SCIP_CALL( SCIPvarAddVlb(var->data.aggregate.var, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
               cliquetable, branchcand, eventqueue, vlbvar, vlbcoef/var->data.aggregate.scalar,
               (vlbconstant - var->data.aggregate.constant)/var->data.aggregate.scalar, transitive, infeasible, nbdchgs) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> add variable upper bound */
         SCIP_CALL( SCIPvarAddVub(var->data.aggregate.var, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
               cliquetable, branchcand, eventqueue, vlbvar, vlbcoef/var->data.aggregate.scalar,
               (vlbconstant - var->data.aggregate.constant)/var->data.aggregate.scalar, transitive, infeasible, nbdchgs) );
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      /* nothing to do here */
      break;

   case SCIP_VARSTATUS_NEGATED:
      /* x = offset - x':  x >= b*z + d  <=>  offset - x' >= b*z + d  <=>  x' <= -b*z + (offset-d) */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarAddVub(var->negatedvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable,
            branchcand, eventqueue, vlbvar, -vlbcoef, var->data.negate.constant - vlbconstant, transitive, infeasible,
            nbdchgs) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** informs variable x about a globally valid variable upper bound x <= b*z + d with integer variable z;
 *  if z is binary, the corresponding valid implication for z is also added;
 *  updates the global bounds of the variable and the vub variable correspondingly
 */
SCIP_RETCODE SCIPvarAddVub(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             vubvar,             /**< variable z    in x <= b*z + d */
   SCIP_Real             vubcoef,            /**< coefficient b in x <= b*z + d */
   SCIP_Real             vubconstant,        /**< constant d    in x <= b*z + d */
   SCIP_Bool             transitive,         /**< should transitive closure of implication also be added? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(SCIPvarGetType(vubvar) != SCIP_VARTYPE_CONTINUOUS);
   assert(infeasible != NULL);

   SCIPsetDebugMsg(set, "adding variable upper bound <%s> <= %g<%s> + %g\n", SCIPvarGetName(var), vubcoef, SCIPvarGetName(vubvar), vubconstant);

   *infeasible = FALSE;
   if( nbdchgs != NULL )
      *nbdchgs = 0;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      SCIP_CALL( SCIPvarAddVub(var->data.original.transvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
            cliquetable, branchcand, eventqueue, vubvar, vubcoef, vubconstant, transitive, infeasible, nbdchgs) );
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      /* transform b*z + d into the corresponding sum after transforming z to an active problem variable */
      SCIP_CALL( SCIPvarGetProbvarSum(&vubvar, set, &vubcoef, &vubconstant) );
      SCIPsetDebugMsg(set, " -> transformed to variable upper bound <%s> <= %g<%s> + %g\n",
         SCIPvarGetName(var), vubcoef, SCIPvarGetName(vubvar), vubconstant);

      /* if the vub coefficient is zero, just update the upper bound of the variable */
      if( SCIPsetIsZero(set, vubcoef) )
      {
         if( SCIPsetIsFeasLT(set, vubconstant, SCIPvarGetLbGlobal(var)) )
            *infeasible = TRUE;
         else if( SCIPsetIsFeasLT(set, vubconstant, SCIPvarGetUbGlobal(var)) )
         {
            /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
             * with the local bound, in this case we need to store the bound change as pending bound change
             */
            if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
            {
               assert(tree != NULL);
               assert(transprob != NULL);
               assert(SCIPprobIsTransformed(transprob));

               SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                     tree, reopt, lp, branchcand, eventqueue, cliquetable, var, vubconstant, SCIP_BOUNDTYPE_UPPER, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPvarChgUbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, vubconstant) );
            }

            if( nbdchgs != NULL )
               (*nbdchgs)++;
         }
      }
      else if( SCIPvarIsActive(vubvar) )
      {
         SCIP_Real xlb;
         SCIP_Real xub;
         SCIP_Real zlb;
         SCIP_Real zub;
         SCIP_Real minvub;
         SCIP_Real maxvub;

         assert(SCIPvarGetStatus(vubvar) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(vubvar) == SCIP_VARSTATUS_COLUMN);
         assert(vubcoef != 0.0);

         minvub = -SCIPsetInfinity(set);
         maxvub = SCIPsetInfinity(set);

         xlb = SCIPvarGetLbGlobal(var);
         xub = SCIPvarGetUbGlobal(var);
         zlb = SCIPvarGetLbGlobal(vubvar);
         zub = SCIPvarGetUbGlobal(vubvar);

         /* improve global bounds of vub variable, and calculate minimal and maximal value of variable bound */
         if( vubcoef >= 0.0 )
         {
            SCIP_Real newzlb;

            if( !SCIPsetIsInfinity(set, -xlb) )
            {
               /* x <= b*z + d  ->  z >= (x-d)/b */
               newzlb = (xlb - vubconstant)/vubcoef;
               if( SCIPsetIsFeasGT(set, newzlb, zub) )
               {
                  *infeasible = TRUE;
                  return SCIP_OKAY;
               }
               if( SCIPsetIsFeasGT(set, newzlb, zlb) )
               {
                  /* bound might be adjusted due to integrality condition */
                  newzlb = adjustedLb(set, SCIPvarGetType(vubvar), newzlb);

                  /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
                   * with the local bound, in this case we need to store the bound change as pending bound change
                   */
                  if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
                  {
                     assert(tree != NULL);
                     assert(transprob != NULL);
                     assert(SCIPprobIsTransformed(transprob));

                     SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                           tree, reopt, lp, branchcand, eventqueue, cliquetable, vubvar, newzlb, SCIP_BOUNDTYPE_LOWER, FALSE) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPvarChgLbGlobal(vubvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newzlb) );
                  }
                  zlb = newzlb;

                  if( nbdchgs != NULL )
                     (*nbdchgs)++;
               }
               minvub = vubcoef * zlb + vubconstant;
               if( !SCIPsetIsInfinity(set, zub) )
                  maxvub = vubcoef * zub + vubconstant;
            }
            else
            {
               if( !SCIPsetIsInfinity(set, zub) )
                  maxvub = vubcoef * zub + vubconstant;
               if( !SCIPsetIsInfinity(set, -zlb) )
                  minvub = vubcoef * zlb + vubconstant;
            }
         }
         else
         {
            SCIP_Real newzub;

            if( !SCIPsetIsInfinity(set, -xlb) )
            {
               /* x <= b*z + d  ->  z <= (x-d)/b */
               newzub = (xlb - vubconstant)/vubcoef;
               if( SCIPsetIsFeasLT(set, newzub, zlb) )
               {
                  *infeasible = TRUE;
                  return SCIP_OKAY;
               }
               if( SCIPsetIsFeasLT(set, newzub, zub) )
               {
                  /* bound might be adjusted due to integrality condition */
                  newzub = adjustedUb(set, SCIPvarGetType(vubvar), newzub);

                  /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
                   * with the local bound, in this case we need to store the bound change as pending bound change
                   */
                  if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
                  {
                     assert(tree != NULL);
                     assert(transprob != NULL);
                     assert(SCIPprobIsTransformed(transprob));

                     SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                           tree, reopt, lp, branchcand, eventqueue, cliquetable, vubvar, newzub, SCIP_BOUNDTYPE_UPPER, FALSE) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPvarChgUbGlobal(vubvar, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, newzub) );
                  }
                  zub = newzub;

                  if( nbdchgs != NULL )
                     (*nbdchgs)++;
               }
               minvub = vubcoef * zub + vubconstant;
               if( !SCIPsetIsInfinity(set, -zlb) )
                  maxvub = vubcoef * zlb + vubconstant;
            }
            else
            {
               if( !SCIPsetIsInfinity(set, zub) )
                  minvub = vubcoef * zub + vubconstant;
               if( !SCIPsetIsInfinity(set, -zlb) )
                  maxvub = vubcoef * zlb + vubconstant;
            }

         }
         if( minvub > maxvub )
            minvub = maxvub;

         /* adjust bounds due to integrality of vub variable */
         minvub = adjustedUb(set, SCIPvarGetType(var), minvub);
         maxvub = adjustedUb(set, SCIPvarGetType(var), maxvub);

         /* check bounds for feasibility */
         if( SCIPsetIsFeasLT(set, maxvub, xlb) || (var == vubvar && SCIPsetIsEQ(set, vubcoef, 1.0) && SCIPsetIsFeasNegative(set, vubconstant))  )
         {
            *infeasible = TRUE;
            return SCIP_OKAY;
         }
         /* improve global upper bound of variable */
         if( SCIPsetIsFeasLT(set, maxvub, xub) )
         {
            /* bound might be adjusted due to integrality condition */
            maxvub = adjustedUb(set, SCIPvarGetType(var), maxvub);

            /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
             * with the local bound, in this case we need to store the bound change as pending bound change
             */
            if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
            {
               assert(tree != NULL);
               assert(transprob != NULL);
               assert(SCIPprobIsTransformed(transprob));

               SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                     tree, reopt, lp, branchcand, eventqueue, cliquetable, var, maxvub, SCIP_BOUNDTYPE_UPPER, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPvarChgUbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, maxvub) );
            }
            xub = maxvub;

            if( nbdchgs != NULL )
               (*nbdchgs)++;
         }
         maxvub = xub;

         /* improve variable bound for binary z by moving the variable's global bound to the vub constant */
         if( SCIPvarIsBinary(vubvar) )
         {
            /* b > 0: x <= (maxvub - minvub) * z + minvub
             * b < 0: x <= (minvub - maxvub) * z + maxvub
             */

            assert(!SCIPsetIsInfinity(set, maxvub) && !SCIPsetIsInfinity(set, -minvub));

            if( vubcoef >= 0.0 )
            {
               vubcoef = maxvub - minvub;
               vubconstant = minvub;
            }
            else
            {
               vubcoef = minvub - maxvub;
               vubconstant = maxvub;
            }
         }

         /* add variable bound to the variable bounds list */
         if( SCIPsetIsFeasLT(set, minvub, xub) )
         {
            assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_FIXED);
            assert(!SCIPsetIsZero(set, vubcoef));

            /* if one of the variables is binary, add the corresponding implication to the variable's implication
             * list, thereby also adding the variable bound (or implication) to the other variable
             */
            if( SCIPvarGetType(vubvar) == SCIP_VARTYPE_BINARY )
            {
               /* add corresponding implication:
                *   b > 0, x <= b*z + d  <->  z == 0 -> x <= d
                *   b < 0, x <= b*z + d  <->  z == 1 -> x <= b+d
                */
               SCIP_CALL( varAddTransitiveImplic(vubvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                     cliquetable, branchcand, eventqueue, (vubcoef < 0.0), var, SCIP_BOUNDTYPE_UPPER, minvub, transitive,
                     infeasible, nbdchgs) );
            }
            else if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
            {
               /* add corresponding implication:
                *   b > 0, x <= b*z + d  <->  x == 1 -> z >= (1-d)/b
                *   b < 0, x <= b*z + d  <->  x == 1 -> z <= (1-d)/b
                */
               SCIP_CALL( varAddTransitiveImplic(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                     cliquetable, branchcand, eventqueue, TRUE, vubvar, (vubcoef >= 0.0 ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER),
                     (1.0-vubconstant)/vubcoef, transitive, infeasible, nbdchgs) );
            }
            else
            {
               SCIP_CALL( varAddVbound(var, blkmem, set, eventqueue, SCIP_BOUNDTYPE_UPPER, vubvar, vubcoef, vubconstant) );
            }
         }
      }
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      /* x = a*y + c:  x <= b*z + d  <=>  a*y + c <= b*z + d  <=>  y <= b/a * z + (d-c)/a, if a > 0
       *                                                           y >= b/a * z + (d-c)/a, if a < 0
       */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> add variable upper bound */
         SCIP_CALL( SCIPvarAddVub(var->data.aggregate.var, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
               cliquetable, branchcand, eventqueue, vubvar, vubcoef/var->data.aggregate.scalar,
               (vubconstant - var->data.aggregate.constant)/var->data.aggregate.scalar, transitive, infeasible, nbdchgs) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> add variable lower bound */
         SCIP_CALL( SCIPvarAddVlb(var->data.aggregate.var, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
               cliquetable, branchcand, eventqueue, vubvar, vubcoef/var->data.aggregate.scalar,
               (vubconstant - var->data.aggregate.constant)/var->data.aggregate.scalar, transitive, infeasible, nbdchgs) );
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      /* nothing to do here */
      break;

   case SCIP_VARSTATUS_NEGATED:
      /* x = offset - x':  x <= b*z + d  <=>  offset - x' <= b*z + d  <=>  x' >= -b*z + (offset-d) */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarAddVlb(var->negatedvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable,
            branchcand, eventqueue, vubvar, -vubcoef, var->data.negate.constant - vubconstant, transitive, infeasible,
            nbdchgs) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** informs binary variable x about a globally valid implication:  x == 0 or x == 1  ==>  y <= b  or  y >= b;
 *  also adds the corresponding implication or variable bound to the implied variable;
 *  if the implication is conflicting, the variable is fixed to the opposite value;
 *  if the variable is already fixed to the given value, the implication is performed immediately;
 *  if the implication is redundant with respect to the variables' global bounds, it is ignored
 */
SCIP_RETCODE SCIPvarAddImplic(
   SCIP_VAR*             var,                /**< problem variable  */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool             transitive,         /**< should transitive closure of implication also be added? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(infeasible != NULL);

   *infeasible = FALSE;
   if( nbdchgs != NULL )
      *nbdchgs = 0;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      SCIP_CALL( SCIPvarAddImplic(var->data.original.transvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
            cliquetable, branchcand, eventqueue, varfixing, implvar, impltype, implbound, transitive, infeasible,
            nbdchgs) );
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      /* if the variable is fixed (although it has no FIXED status), and varfixing corresponds to the fixed value of
       * the variable, the implication can be applied directly;
       * otherwise, add implication to the implications list (and add inverse of implication to the implied variable)
       */
      if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5 )
      {
         if( varfixing == (SCIPvarGetLbGlobal(var) > 0.5) )
         {
            SCIP_CALL( applyImplic(blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
                  cliquetable, implvar, impltype, implbound, infeasible, nbdchgs) );
         }
      }
      else
      {
         SCIP_CALL( SCIPvarGetProbvarBound(&implvar, &implbound, &impltype) );
         SCIPvarAdjustBd(implvar, set, impltype, &implbound);
         if( SCIPvarIsActive(implvar) || SCIPvarGetStatus(implvar) == SCIP_VARSTATUS_FIXED )
         {
            SCIP_CALL( varAddTransitiveImplic(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable,
                  branchcand, eventqueue, varfixing, implvar, impltype, implbound, transitive, infeasible, nbdchgs) );
         }
      }
      break;

   case SCIP_VARSTATUS_FIXED:
      /* if varfixing corresponds to the fixed value of the variable, the implication can be applied directly */
      if( varfixing == (SCIPvarGetLbGlobal(var) > 0.5) )
      {
         SCIP_CALL( applyImplic(blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               cliquetable, implvar, impltype, implbound, infeasible, nbdchgs) );
      }
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      /* implication added for x == 1:
       *   x == 1 && x =  1*z + 0  ==>  y <= b or y >= b    <==>    z >= 1  ==>  y <= b or y >= b 
       *   x == 1 && x = -1*z + 1  ==>  y <= b or y >= b    <==>    z <= 0  ==>  y <= b or y >= b
       * implication added for x == 0:
       *   x == 0 && x =  1*z + 0  ==>  y <= b or y >= b    <==>    z <= 0  ==>  y <= b or y >= b
       *   x == 0 && x = -1*z + 1  ==>  y <= b or y >= b    <==>    z >= 1  ==>  y <= b or y >= b
       *
       * use only binary variables z 
       */
      assert(var->data.aggregate.var != NULL);
      if( SCIPvarIsBinary(var->data.aggregate.var) )
      {
         assert( (SCIPsetIsEQ(set, var->data.aggregate.scalar, 1.0) && SCIPsetIsZero(set, var->data.aggregate.constant))
            || (SCIPsetIsEQ(set, var->data.aggregate.scalar, -1.0) && SCIPsetIsEQ(set, var->data.aggregate.constant, 1.0)) );

         if( var->data.aggregate.scalar > 0 )
         {
            SCIP_CALL( SCIPvarAddImplic(var->data.aggregate.var, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                  cliquetable, branchcand, eventqueue, varfixing, implvar, impltype, implbound, transitive, infeasible,
                  nbdchgs) );
         }
         else
         {
            SCIP_CALL( SCIPvarAddImplic(var->data.aggregate.var, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                  cliquetable, branchcand, eventqueue, !varfixing, implvar, impltype, implbound, transitive, infeasible,
                  nbdchgs) );
         }
      }
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      /* nothing to do here */
      break;

   case SCIP_VARSTATUS_NEGATED:
      /* implication added for x == 1:
       *   x == 1 && x = -1*z + 1  ==>  y <= b or y >= b    <==>    z <= 0  ==>  y <= b or y >= b
       * implication added for x == 0:
       *   x == 0 && x = -1*z + 1  ==>  y <= b or y >= b    <==>    z >= 1  ==>  y <= b or y >= b
       */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      assert(SCIPvarIsBinary(var->negatedvar));

      if( SCIPvarGetType(var->negatedvar) == SCIP_VARTYPE_BINARY )
      {
         SCIP_CALL( SCIPvarAddImplic(var->negatedvar, blkmem, set, stat,  transprob, origprob, tree, reopt, lp,
               cliquetable, branchcand, eventqueue, !varfixing, implvar, impltype, implbound, transitive, infeasible, nbdchgs) );
      }
      /* in case one both variables are not of binary type we have to add the implication as variable bounds */
      else
      {
         /* if the implied variable is of binary type exchange the variables */
         if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
         {
            SCIP_CALL( SCIPvarAddImplic(implvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable,
                  branchcand, eventqueue, (impltype == SCIP_BOUNDTYPE_UPPER) ? TRUE : FALSE, var->negatedvar,
                  varfixing ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER, varfixing ? 1.0 : 0.0, transitive,
                  infeasible, nbdchgs) );
         }
         else
         {
            /* both variables are not of binary type but are implicit binary; in that case we can only add this
             * implication as variable bounds
             */

            /* add variable lower bound on the negation of var */
            if( varfixing )
            {
               /* (x = 1 => i) z = 0 ii) z = 1) <=> ( i) z = 1 ii) z = 0 => ~x = 1), this is done by adding ~x >= b*z + d
                * as variable lower bound
                */
               SCIP_CALL( SCIPvarAddVlb(var->negatedvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                     cliquetable, branchcand, eventqueue, implvar, (impltype == SCIP_BOUNDTYPE_UPPER) ? 1.0 : -1.0,
                     (impltype == SCIP_BOUNDTYPE_UPPER) ? 0.0 : 1.0, transitive, infeasible, nbdchgs) );
            }
            else
            {
               /* (x = 0 => i) z = 0 ii) z = 1) <=> ( i) z = 1 ii) z = 0 => ~x = 0), this is done by adding ~x <= b*z + d
                * as variable upper bound
                */
               SCIP_CALL( SCIPvarAddVub(var->negatedvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                     cliquetable, branchcand, eventqueue, implvar, (impltype == SCIP_BOUNDTYPE_UPPER) ? -1.0 : 1.0,
                     (impltype == SCIP_BOUNDTYPE_UPPER) ? 1.0 : 0.0, transitive, infeasible, nbdchgs) );
            }

            /* add variable bound on implvar */
            if( impltype == SCIP_BOUNDTYPE_UPPER )
            {
               /* (z = 1 => i) x = 0 ii) x = 1) <=> ( i) ~x = 0 ii) ~x = 1 => z = 0), this is done by adding z <= b*~x + d
                * as variable upper bound
                */
               SCIP_CALL( SCIPvarAddVub(implvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable,
                     branchcand, eventqueue, var->negatedvar, (varfixing) ? 1.0 : -1.0,
                     (varfixing) ? 0.0 : 1.0, transitive, infeasible, nbdchgs) );
            }
            else
            {
               /* (z = 0 => i) x = 0 ii) x = 1) <=> ( i) ~x = 0 ii) ~x = 1 => z = 1), this is done by adding z >= b*~x + d
                * as variable upper bound
                */
               SCIP_CALL( SCIPvarAddVlb(implvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp, cliquetable,
                     branchcand, eventqueue, var->negatedvar, (varfixing) ? -1.0 : 1.0, (varfixing) ? 1.0 : 0.0,
                     transitive, infeasible, nbdchgs) );
            }
         }
      }
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** returns whether there is an implication x == varfixing -> y <= b or y >= b in the implication graph;
 *  implications that are represented as cliques in the clique table are not regarded (use SCIPvarsHaveCommonClique());
 *  both variables must be active, variable x must be binary
 */
SCIP_Bool SCIPvarHasImplic(
   SCIP_VAR*             var,                /**< problem variable x */
   SCIP_Bool             varfixing,          /**< FALSE if y should be searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   SCIP_BOUNDTYPE        impltype            /**< type of implication y <=/>= b to search for */
   )
{
   assert(var != NULL);
   assert(implvar != NULL);
   assert(SCIPvarIsActive(var));
   assert(SCIPvarIsActive(implvar));
   assert(SCIPvarIsBinary(var));

   return var->implics != NULL && SCIPimplicsContainsImpl(var->implics, varfixing, implvar, impltype);
}

/** returns whether there is an implication x == varfixing -> y == implvarfixing in the implication graph;
 *  implications that are represented as cliques in the clique table are not regarded (use SCIPvarsHaveCommonClique());
 *  both variables must be active binary variables
 */
SCIP_Bool SCIPvarHasBinaryImplic(
   SCIP_VAR*             var,                /**< problem variable x */
   SCIP_Bool             varfixing,          /**< FALSE if y should be searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   SCIP_Bool             implvarfixing       /**< value of the implied variable to search for */
   )
{
   assert(SCIPvarIsBinary(implvar));

   return SCIPvarHasImplic(var, varfixing, implvar, implvarfixing ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER);
}

/** fixes the bounds of a binary variable to the given value, counting bound changes and detecting infeasibility */
SCIP_RETCODE SCIPvarFixBinary(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             value,              /**< value to fix variable to */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   if( value == FALSE )
   {
      if( var->glbdom.lb > 0.5 )
         *infeasible = TRUE;
      else if( var->glbdom.ub > 0.5 )
      {
         /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
          * with the local bound, in this case we need to store the bound change as pending bound change
          */
         if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
         {
            assert(tree != NULL);
            assert(transprob != NULL);
            assert(SCIPprobIsTransformed(transprob));

            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, var, 0.0, SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
         else
         {
            SCIP_CALL( SCIPvarChgUbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, 0.0) );
         }

         if( nbdchgs != NULL )
            (*nbdchgs)++;
      }
   }
   else
   {
      if( var->glbdom.ub < 0.5 )
         *infeasible = TRUE;
      else if( var->glbdom.lb < 0.5 )
      {
         /* during solving stage it can happen that the global bound change cannot be applied directly because it conflicts
          * with the local bound, in this case we need to store the bound change as pending bound change
          */
         if( SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING )
         {
            assert(tree != NULL);
            assert(transprob != NULL);
            assert(SCIPprobIsTransformed(transprob));

            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, var, 1.0, SCIP_BOUNDTYPE_LOWER, FALSE) );
         }
         else
         {
            SCIP_CALL( SCIPvarChgLbGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, cliquetable, 1.0) );
         }

         if( nbdchgs != NULL )
            (*nbdchgs)++;
      }
   }

   /* during presolving, the variable should have been removed immediately from all its cliques */
   assert(SCIPsetGetStage(set) >= SCIP_STAGE_SOLVING || var->cliquelist == NULL);

   return SCIP_OKAY;
}

/** adds the variable to the given clique and updates the list of cliques the binary variable is member of;
 *  if the variable now appears twice in the clique with the same value, it is fixed to the opposite value;
 *  if the variable now appears twice in the clique with opposite values, all other variables are fixed to
 *  the opposite of the value they take in the clique
 */
SCIP_RETCODE SCIPvarAddClique(
   SCIP_VAR*             var,                /**< problem variable  */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             value,              /**< value of the variable in the clique */
   SCIP_CLIQUE*          clique,             /**< clique the variable should be added to */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(SCIPvarIsBinary(var));
   assert(infeasible != NULL);

   *infeasible = FALSE;

   /* get corresponding active problem variable */
   SCIP_CALL( SCIPvarGetProbvarBinary(&var, &value) );
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
   assert(SCIPvarIsBinary(var));

   /* only column and loose variables may be member of a clique */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      SCIP_Bool doubleentry;
      SCIP_Bool oppositeentry;

      /* add variable to clique */
      SCIP_CALL( SCIPcliqueAddVar(clique, blkmem, set, var, value, &doubleentry, &oppositeentry) );

      /* add clique to variable's clique list */
      SCIP_CALL( SCIPcliquelistAdd(&var->cliquelist, blkmem, set, value, clique) );

      /* check consistency of cliquelist */
      SCIPcliquelistCheck(var->cliquelist, var);

      /* if the variable now appears twice with the same value in the clique, it can be fixed to the opposite value */
      if( doubleentry )
      {
         SCIP_CALL( SCIPvarFixBinary(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
               eventqueue, cliquetable, !value, infeasible, nbdchgs) );
      }

      /* if the variable appears with both values in the clique, all other variables of the clique can be fixed
       * to the opposite of the value they take in the clique
       */
      if( oppositeentry )
      {
         SCIP_VAR** vars;
         SCIP_Bool* values;
         int nvars;
         int i;

         nvars = SCIPcliqueGetNVars(clique);
         vars = SCIPcliqueGetVars(clique);
         values = SCIPcliqueGetValues(clique);
         for( i = 0; i < nvars && !(*infeasible); ++i )
         {
            if( vars[i] == var )
               continue;

            SCIP_CALL( SCIPvarFixBinary(vars[i], blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
                  eventqueue, cliquetable, !values[i], infeasible, nbdchgs) );
         }
      }
   }

   return SCIP_OKAY;
}

/** adds a filled clique to the cliquelists of all corresponding variables */
SCIP_RETCODE SCIPvarsAddClique(
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_Bool*            values,             /**< values of the variables in the clique */
   int                   nvars,              /**< number of problem variables */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CLIQUE*          clique              /**< clique that contains all given variables and values */
   )
{
   SCIP_VAR* var;
   int v;

   assert(vars != NULL);
   assert(values != NULL);
   assert(nvars > 0);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(clique != NULL);

   for( v = nvars - 1; v >= 0; --v )
   {
      var = vars[v];
      assert(SCIPvarIsBinary(var));
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

      /* add clique to variable's clique list */
      SCIP_CALL( SCIPcliquelistAdd(&var->cliquelist, blkmem, set, values[v], clique) );

      /* check consistency of cliquelist */
      SCIPcliquelistCheck(var->cliquelist, var);
   }

   return SCIP_OKAY;
}

/** adds a clique to the list of cliques of the given binary variable, but does not change the clique
 *  itself
 */
SCIP_RETCODE SCIPvarAddCliqueToList(
   SCIP_VAR*             var,                /**< problem variable  */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             value,              /**< value of the variable in the clique */
   SCIP_CLIQUE*          clique              /**< clique that should be removed from the variable's clique list */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsBinary(var));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

   /* add clique to variable's clique list */
   SCIP_CALL( SCIPcliquelistAdd(&var->cliquelist, blkmem, set, value, clique) );

   return SCIP_OKAY;
}


/** deletes a clique from the list of cliques the binary variable is member of, but does not change the clique
 *  itself
 */
SCIP_RETCODE SCIPvarDelCliqueFromList(
   SCIP_VAR*             var,                /**< problem variable  */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             value,              /**< value of the variable in the clique */
   SCIP_CLIQUE*          clique              /**< clique that should be removed from the variable's clique list */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsBinary(var));

   /* delete clique from variable's clique list */
   SCIP_CALL( SCIPcliquelistDel(&var->cliquelist, blkmem, value, clique) );

   return SCIP_OKAY;
}

/** deletes the variable from the given clique and updates the list of cliques the binary variable is member of */
SCIP_RETCODE SCIPvarDelClique(
   SCIP_VAR*             var,                /**< problem variable  */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             value,              /**< value of the variable in the clique */
   SCIP_CLIQUE*          clique              /**< clique the variable should be removed from */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsBinary(var));

   /* get corresponding active problem variable */
   SCIP_CALL( SCIPvarGetProbvarBinary(&var, &value) );
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
   assert(SCIPvarIsBinary(var));

   /* only column and loose variables may be member of a clique */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      /* delete clique from variable's clique list */
      SCIP_CALL( SCIPcliquelistDel(&var->cliquelist, blkmem, value, clique) );

      /* delete variable from clique */
      SCIPcliqueDelVar(clique, cliquetable, var, value);

      /* check consistency of cliquelist */
      SCIPcliquelistCheck(var->cliquelist, var);
   }

   return SCIP_OKAY;
}

/** returns whether there is a clique that contains both given variable/value pairs;
 *  the variables must be active binary variables;
 *  if regardimplics is FALSE, only the cliques in the clique table are looked at;
 *  if regardimplics is TRUE, both the cliques and the implications of the implication graph are regarded
 *
 *  @note a variable with it's negated variable are NOT! in a clique
 *  @note a variable with itself are in a clique
 */
SCIP_Bool SCIPvarsHaveCommonClique(
   SCIP_VAR*             var1,               /**< first variable */
   SCIP_Bool             value1,             /**< value of first variable */
   SCIP_VAR*             var2,               /**< second variable */
   SCIP_Bool             value2,             /**< value of second variable */
   SCIP_Bool             regardimplics       /**< should the implication graph also be searched for a clique? */
   )
{
   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(SCIPvarIsActive(var1));
   assert(SCIPvarIsActive(var2));
   assert(SCIPvarIsBinary(var1));
   assert(SCIPvarIsBinary(var2));

   return (SCIPcliquelistsHaveCommonClique(var1->cliquelist, value1, var2->cliquelist, value2)
      || (regardimplics && SCIPvarHasImplic(var1, value1, var2, value2 ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER)));
}

/** actually changes the branch factor of the variable and of all parent variables */
static
SCIP_RETCODE varProcessChgBranchFactor(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             branchfactor        /**< factor to weigh variable's branching score with */
   )
{
   SCIP_VAR* parentvar;
   SCIP_Real eps;
   int i;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* only use positive values */
   eps = SCIPsetEpsilon(set);
   branchfactor = MAX(branchfactor, eps);

   SCIPsetDebugMsg(set, "process changing branch factor of <%s> from %f to %f\n", var->name, var->branchfactor, branchfactor);

   if( SCIPsetIsEQ(set, branchfactor, var->branchfactor) )
      return SCIP_OKAY;

   /* change the branch factor */
   var->branchfactor = branchfactor;

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         /* do not change priorities across the border between transformed and original problem */
         break;

      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         SCIPerrorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/

      case SCIP_VARSTATUS_AGGREGATED:
      case SCIP_VARSTATUS_NEGATED:
         SCIP_CALL( varProcessChgBranchFactor(parentvar, set, branchfactor) );
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         SCIPABORT();
         return SCIP_ERROR; /*lint !e527*/
      }
   }

   return SCIP_OKAY;
}

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
SCIP_RETCODE SCIPvarChgBranchFactor(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             branchfactor        /**< factor to weigh variable's branching score with */
   )
{
   int v;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(branchfactor >= 0.0);

   SCIPdebugMessage("changing branch factor of <%s> from %g to %g\n", var->name, var->branchfactor, branchfactor);

   if( SCIPsetIsEQ(set, var->branchfactor, branchfactor) )
      return SCIP_OKAY;

   /* change priorities of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar != NULL )
      {
         SCIP_CALL( SCIPvarChgBranchFactor(var->data.original.transvar, set, branchfactor) );
      }
      else
      {
         assert(set->stage == SCIP_STAGE_PROBLEM);
         var->branchfactor = branchfactor;
      }
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      SCIP_CALL( varProcessChgBranchFactor(var, set, branchfactor) );
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      SCIP_CALL( SCIPvarChgBranchFactor(var->data.aggregate.var, set, branchfactor) );
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(!var->donotmultaggr);
      for( v = 0; v < var->data.multaggr.nvars; ++v )
      {
         SCIP_CALL( SCIPvarChgBranchFactor(var->data.multaggr.vars[v], set, branchfactor) );
      }
      break;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarChgBranchFactor(var->negatedvar, set, branchfactor) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_ERROR; /*lint !e527*/
   }

   return SCIP_OKAY;
}

/** actually changes the branch priority of the variable and of all parent variables */
static
SCIP_RETCODE varProcessChgBranchPriority(
   SCIP_VAR*             var,                /**< problem variable */
   int                   branchpriority      /**< branching priority of the variable */
   )
{
   SCIP_VAR* parentvar;
   int i;

   assert(var != NULL);

   SCIPdebugMessage("process changing branch priority of <%s> from %d to %d\n", 
      var->name, var->branchpriority, branchpriority);

   if( branchpriority == var->branchpriority )
      return SCIP_OKAY;

   /* change the branch priority */
   var->branchpriority = branchpriority;

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         /* do not change priorities across the border between transformed and original problem */
         break;

      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         SCIPerrorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/

      case SCIP_VARSTATUS_AGGREGATED:
      case SCIP_VARSTATUS_NEGATED:
         SCIP_CALL( varProcessChgBranchPriority(parentvar, branchpriority) );
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_ERROR;
      }
   }

   return SCIP_OKAY;
}

/** sets the branch priority of the variable; variables with higher branch priority are always preferred to variables
 *  with lower priority in selection of branching variable
 */
SCIP_RETCODE SCIPvarChgBranchPriority(
   SCIP_VAR*             var,                /**< problem variable */
   int                   branchpriority      /**< branching priority of the variable */
   )
{
   int v;

   assert(var != NULL);

   SCIPdebugMessage("changing branch priority of <%s> from %d to %d\n", var->name, var->branchpriority, branchpriority);

   if( var->branchpriority == branchpriority )
      return SCIP_OKAY;

   /* change priorities of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar != NULL )
      {
         SCIP_CALL( SCIPvarChgBranchPriority(var->data.original.transvar, branchpriority) );
      }
      else
         var->branchpriority = branchpriority;
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      SCIP_CALL( varProcessChgBranchPriority(var, branchpriority) );
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      SCIP_CALL( SCIPvarChgBranchPriority(var->data.aggregate.var, branchpriority) );
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(!var->donotmultaggr);      
      for( v = 0; v < var->data.multaggr.nvars; ++v )
      {
         SCIP_CALL( SCIPvarChgBranchPriority(var->data.multaggr.vars[v], branchpriority) );
      }
      break;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarChgBranchPriority(var->negatedvar, branchpriority) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_ERROR; /*lint !e527*/
   }

   return SCIP_OKAY;
}

/** actually changes the branch direction of the variable and of all parent variables */
static
SCIP_RETCODE varProcessChgBranchDirection(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        branchdirection     /**< preferred branch direction of the variable (downwards, upwards, auto) */
   )
{
   SCIP_VAR* parentvar;
   int i;

   assert(var != NULL);

   SCIPdebugMessage("process changing branch direction of <%s> from %u to %d\n", 
      var->name, var->branchdirection, branchdirection);

   if( branchdirection == (SCIP_BRANCHDIR)var->branchdirection )
      return SCIP_OKAY;

   /* change the branch direction */
   var->branchdirection = branchdirection; /*lint !e641*/

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         /* do not change directions across the border between transformed and original problem */
         break;

      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         SCIPerrorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/

      case SCIP_VARSTATUS_AGGREGATED:
         if( parentvar->data.aggregate.scalar > 0.0 )
         {
            SCIP_CALL( varProcessChgBranchDirection(parentvar, branchdirection) );
         }
         else
         {
            SCIP_CALL( varProcessChgBranchDirection(parentvar, SCIPbranchdirOpposite(branchdirection)) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED:
         SCIP_CALL( varProcessChgBranchDirection(parentvar, SCIPbranchdirOpposite(branchdirection)) );
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         SCIPABORT();
         return SCIP_ERROR; /*lint !e527*/
      }
   }

   return SCIP_OKAY;
}

/** sets the branch direction of the variable; variables with higher branch direction are always preferred to variables
 *  with lower direction in selection of branching variable
 */
SCIP_RETCODE SCIPvarChgBranchDirection(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        branchdirection     /**< preferred branch direction of the variable (downwards, upwards, auto) */
   )
{
   int v;

   assert(var != NULL);

   SCIPdebugMessage("changing branch direction of <%s> from %u to %d\n", var->name, var->branchdirection, branchdirection);

   if( (SCIP_BRANCHDIR)var->branchdirection == branchdirection )
      return SCIP_OKAY;

   /* change directions of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar != NULL )
      {
         SCIP_CALL( SCIPvarChgBranchDirection(var->data.original.transvar, branchdirection) );
      }
      else
         var->branchdirection = branchdirection; /*lint !e641*/
      break;

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      SCIP_CALL( varProcessChgBranchDirection(var, branchdirection) );
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      if( var->data.aggregate.scalar > 0.0 )
      {
         SCIP_CALL( SCIPvarChgBranchDirection(var->data.aggregate.var, branchdirection) );
      }
      else
      {
         SCIP_CALL( SCIPvarChgBranchDirection(var->data.aggregate.var, SCIPbranchdirOpposite(branchdirection)) );
      }
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(!var->donotmultaggr);
      for( v = 0; v < var->data.multaggr.nvars; ++v )
      {
         /* only update branching direction of aggregation variables, if they don't have a preferred direction yet */
         assert(var->data.multaggr.vars[v] != NULL);
         if( (SCIP_BRANCHDIR)var->data.multaggr.vars[v]->branchdirection == SCIP_BRANCHDIR_AUTO )
         {
            if( var->data.multaggr.scalars[v] > 0.0 )
            {
               SCIP_CALL( SCIPvarChgBranchDirection(var->data.multaggr.vars[v], branchdirection) );
            }
            else
            {
               SCIP_CALL( SCIPvarChgBranchDirection(var->data.multaggr.vars[v], SCIPbranchdirOpposite(branchdirection)) );
            }
         }
      }
      break;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarChgBranchDirection(var->negatedvar, SCIPbranchdirOpposite(branchdirection)) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_ERROR; /*lint !e527*/
   }

   return SCIP_OKAY;
}

/** compares the index of two variables, only active, fixed or negated variables are allowed, if a variable
 *  is negated then the index of the corresponding active variable is taken, returns -1 if first is
 *  smaller than, and +1 if first is greater than second variable index; returns 0 if both indices
 *  are equal, which means both variables are equal
 */
int SCIPvarCompareActiveAndNegated(
   SCIP_VAR*             var1,               /**< first problem variable */
   SCIP_VAR*             var2                /**< second problem variable */
   )
{
   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(SCIPvarIsActive(var1) || SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(var1) == SCIP_VARSTATUS_FIXED);
   assert(SCIPvarIsActive(var2) || SCIPvarGetStatus(var2) == SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(var2) == SCIP_VARSTATUS_FIXED);

   if( SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED )
      var1 = SCIPvarGetNegatedVar(var1);
   if( SCIPvarGetStatus(var2) == SCIP_VARSTATUS_NEGATED )
      var2 = SCIPvarGetNegatedVar(var2);

   assert(var1 != NULL);
   assert(var2 != NULL);

   if( SCIPvarGetIndex(var1) < SCIPvarGetIndex(var2) )
      return -1;
   else if( SCIPvarGetIndex(var1) > SCIPvarGetIndex(var2) )
      return +1;

   assert(var1 == var2);
   return 0;
}

/** comparison method for sorting active and negated variables by non-decreasing index, active and negated 
 *  variables are handled as the same variables
 */
SCIP_DECL_SORTPTRCOMP(SCIPvarCompActiveAndNegated)
{
   return SCIPvarCompareActiveAndNegated((SCIP_VAR*)elem1, (SCIP_VAR*)elem2);
}

/** compares the index of two variables, returns -1 if first is smaller than, and +1 if first is greater than second
 *  variable index; returns 0 if both indices are equal, which means both variables are equal
 */
int SCIPvarCompare(
   SCIP_VAR*             var1,               /**< first problem variable */
   SCIP_VAR*             var2                /**< second problem variable */
   )
{
   assert(var1 != NULL);
   assert(var2 != NULL);

   if( var1->index < var2->index )
      return -1;
   else if( var1->index > var2->index )
      return +1;
   else
   {
      assert(var1 == var2);
      return 0;
   }
}

/** comparison method for sorting variables by non-decreasing index */
SCIP_DECL_SORTPTRCOMP(SCIPvarComp)
{
   return SCIPvarCompare((SCIP_VAR*)elem1, (SCIP_VAR*)elem2);
}

/** comparison method for sorting variables by non-decreasing objective coefficient */
SCIP_DECL_SORTPTRCOMP(SCIPvarCompObj)
{
   SCIP_Real obj1;
   SCIP_Real obj2;

   obj1 = SCIPvarGetObj((SCIP_VAR*)elem1);
   obj2 = SCIPvarGetObj((SCIP_VAR*)elem2);

   if( obj1 < obj2 )
      return -1;
   else if( obj1 > obj2 )
      return +1;
   else
      return 0;
}

/** hash key retrieval function for variables */
SCIP_DECL_HASHGETKEY(SCIPvarGetHashkey)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff the indices of both variables are equal */
SCIP_DECL_HASHKEYEQ(SCIPvarIsHashkeyEq)
{  /*lint --e{715}*/
   if( key1 == key2 )
      return TRUE;
   return FALSE;
}

/** returns the hash value of the key */
SCIP_DECL_HASHKEYVAL(SCIPvarGetHashkeyVal)
{  /*lint --e{715}*/
   assert( SCIPvarGetIndex((SCIP_VAR*) key) >= 0 );
   return (unsigned int) SCIPvarGetIndex((SCIP_VAR*) key);
}

/** return for given variables all their active counterparts; all active variables will be pairwise different */
SCIP_RETCODE SCIPvarsGetActiveVars(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< variable array with given variables and as output all active
					      *   variables, if enough slots exist
					      */
   int*                  nvars,              /**< number of given variables, and as output number of active variables,
					      *   if enough slots exist
					      */
   int                   varssize,           /**< available slots in vars array */
   int*                  requiredsize        /**< pointer to store the required array size for the active variables */
   )
{
   SCIP_VAR** activevars;
   int nactivevars;
   int activevarssize;

   SCIP_VAR* var;
   int v;

   SCIP_VAR** tmpvars;
   SCIP_VAR** multvars;
   int tmpvarssize;
   int ntmpvars;
   int noldtmpvars;
   int nmultvars;

   assert(set != NULL);
   assert(nvars != NULL);
   assert(vars != NULL || *nvars == 0);
   assert(varssize >= *nvars);
   assert(requiredsize != NULL);

   *requiredsize = 0;

   if( *nvars == 0 )
      return SCIP_OKAY;

   nactivevars = 0;
   activevarssize = *nvars;
   ntmpvars = *nvars;
   tmpvarssize = *nvars;

   /* temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activevars, activevarssize) );
   SCIP_CALL( SCIPsetDuplicateBufferArray(set, &tmpvars, vars, ntmpvars) );

   noldtmpvars = ntmpvars;

   /* sort all variables to combine equal variables easily */
   SCIPsortPtr((void**)tmpvars, SCIPvarComp, ntmpvars);
   for( v = ntmpvars - 1; v > 0; --v )
   {
      /* combine same variables */
      if( SCIPvarCompare(tmpvars[v], tmpvars[v - 1]) == 0 )
      {
         --ntmpvars;
         tmpvars[v] = tmpvars[ntmpvars];
      }
   }
   /* sort all variables again to combine equal variables later on */
   if( noldtmpvars > ntmpvars )
      SCIPsortPtr((void**)tmpvars, SCIPvarComp, ntmpvars);

   /* collect for each variable the representation in active variables */
   while( ntmpvars >= 1 )
   {
      --ntmpvars;
      var = tmpvars[ntmpvars];
      assert( var != NULL );

      switch( SCIPvarGetStatus(var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
	 if( var->data.original.transvar == NULL )
	 {
	    SCIPerrorMessage("original variable has no transformed variable attached\n");
	    SCIPABORT();
	    return SCIP_INVALIDDATA; /*lint !e527*/
	 }
	 tmpvars[ntmpvars] = var->data.original.transvar;
	 ++ntmpvars;
	 break;

      case SCIP_VARSTATUS_AGGREGATED:
	 tmpvars[ntmpvars] = var->data.aggregate.var;
	 ++ntmpvars;
	 break;

      case SCIP_VARSTATUS_NEGATED:
	 tmpvars[ntmpvars] = var->negatedvar;
	 ++ntmpvars;
	 break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
	 /* check for space in temporary memory */
         if( nactivevars >= activevarssize )
         {
            activevarssize *= 2;
            SCIP_CALL( SCIPsetReallocBufferArray(set, &activevars, activevarssize) );
            assert(nactivevars < activevarssize);
         }
         activevars[nactivevars] = var;
         nactivevars++;
         break;

      case SCIP_VARSTATUS_MULTAGGR:
         /* x = a_1*y_1 + ... + a_n*y_n + c */
         nmultvars = var->data.multaggr.nvars;
         multvars = var->data.multaggr.vars;

	 /* check for space in temporary memory */
         if( nmultvars + ntmpvars > tmpvarssize )
         {
            while( nmultvars + ntmpvars > tmpvarssize )
               tmpvarssize *= 2;
            SCIP_CALL( SCIPsetReallocBufferArray(set, &tmpvars, tmpvarssize) );
            assert(nmultvars + ntmpvars <= tmpvarssize);
         }

	 /* copy all multi-aggregation variables into our working array */
	 BMScopyMemoryArray(&tmpvars[ntmpvars], multvars, nmultvars); /*lint !e866*/

	 /* get active, fixed or multi-aggregated corresponding variables for all new ones */
	 SCIPvarsGetProbvar(&tmpvars[ntmpvars], nmultvars);

	 ntmpvars += nmultvars;
	 noldtmpvars = ntmpvars;

	 /* sort all variables to combine equal variables easily */
	 SCIPsortPtr((void**)tmpvars, SCIPvarComp, ntmpvars);
	 for( v = ntmpvars - 1; v > 0; --v )
	 {
	    /* combine same variables */
	    if( SCIPvarCompare(tmpvars[v], tmpvars[v - 1]) == 0 )
	    {
	       --ntmpvars;
	       tmpvars[v] = tmpvars[ntmpvars];
	    }
	 }
	 /* sort all variables again to combine equal variables later on */
	 if( noldtmpvars > ntmpvars )
	    SCIPsortPtr((void**)tmpvars, SCIPvarComp, ntmpvars);

         break;

      case SCIP_VARSTATUS_FIXED:
	 /* no need for memorizing fixed variables */
         break;

      default:
	 SCIPerrorMessage("unknown variable status\n");
         SCIPABORT();
	 return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   /* sort variable array by variable index */
   SCIPsortPtr((void**)activevars, SCIPvarComp, nactivevars);

   /* eliminate duplicates and count required size */
   v = nactivevars - 1;
   while( v > 0 )
   {
      /* combine both variable since they are the same */
      if( SCIPvarCompare(activevars[v - 1], activevars[v]) == 0 )
      {
	 --nactivevars;
	 activevars[v] = activevars[nactivevars];
      }
      --v;
   }
   *requiredsize = nactivevars;

   if( varssize >= *requiredsize )
   {
      assert(vars != NULL);

      *nvars = *requiredsize;
      BMScopyMemoryArray(vars, activevars, nactivevars);
   }

   SCIPsetFreeBufferArray(set, &tmpvars);
   SCIPsetFreeBufferArray(set, &activevars);

   return SCIP_OKAY;
}

/** gets corresponding active, fixed, or multi-aggregated problem variables of given variables,
 *  @note the content of the given array will/might change
 */
void SCIPvarsGetProbvar(
   SCIP_VAR**            vars,               /**< array of problem variables */
   int                   nvars               /**< number of variables */
   )
{
   int v;

   assert(vars != NULL || nvars == 0);

   for( v = nvars - 1; v >= 0; --v )
   {
      assert(vars != NULL);
      assert(vars[v] != NULL);

      vars[v] = SCIPvarGetProbvar(vars[v]);
      assert(vars[v] != NULL);
   }
}

/** gets corresponding active, fixed, or multi-aggregated problem variable of a variable */
SCIP_VAR* SCIPvarGetProbvar(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_VAR* retvar;

   assert(var != NULL);

   retvar = var;

   SCIPdebugMessage("get problem variable of <%s>\n", var->name);

   while( TRUE ) /*lint !e716 */
   {
      assert(retvar != NULL);

      switch( SCIPvarGetStatus(retvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
	 if( retvar->data.original.transvar == NULL )
	 {
	    SCIPerrorMessage("original variable has no transformed variable attached\n");
	    SCIPABORT();
	    return NULL; /*lint !e527 */
	 }
	 retvar = retvar->data.original.transvar;
	 break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_FIXED:
	 return retvar;

      case SCIP_VARSTATUS_MULTAGGR:
	 /* handle multi-aggregated variables depending on one variable only (possibly caused by SCIPvarFlattenAggregationGraph()) */
	 if ( retvar->data.multaggr.nvars == 1 )
	    retvar = retvar->data.multaggr.vars[0];
	 else
	    return retvar;
	 break;

      case SCIP_VARSTATUS_AGGREGATED:
	 retvar = retvar->data.aggregate.var;
	 break;

      case SCIP_VARSTATUS_NEGATED:
	 retvar = retvar->negatedvar;
	 break;

      default:
	 SCIPerrorMessage("unknown variable status\n");
	 SCIPABORT();
	 return NULL; /*lint !e527*/
      }
   }
}

/** gets corresponding active, fixed, or multi-aggregated problem variables of binary variables and updates the given
 *  negation status of each variable
 */
SCIP_RETCODE SCIPvarsGetProbvarBinary(
   SCIP_VAR***           vars,               /**< pointer to binary problem variables */
   SCIP_Bool**           negatedarr,         /**< pointer to corresponding array to update the negation status */
   int                   nvars               /**< number of variables and values in vars and negated array */
   )
{
   SCIP_VAR** var;
   SCIP_Bool* negated;
   int v;

   assert(vars != NULL);
   assert(*vars != NULL || nvars == 0);
   assert(negatedarr != NULL);
   assert(*negatedarr != NULL || nvars == 0);

   for( v = nvars - 1; v >= 0; --v )
   {
      var = &((*vars)[v]);
      negated = &((*negatedarr)[v]);

      /* get problem variable */
      SCIP_CALL( SCIPvarGetProbvarBinary(var, negated) );
   }

   return SCIP_OKAY;
}


/** gets corresponding active, fixed, or multi-aggregated problem variable of a binary variable and updates the given
 *  negation status (this means you have to assign a value to SCIP_Bool negated before calling this method, usually
 *  FALSE is used)
 */
SCIP_RETCODE SCIPvarGetProbvarBinary(
   SCIP_VAR**            var,                /**< pointer to binary problem variable */
   SCIP_Bool*            negated             /**< pointer to update the negation status */
   )
{
   SCIP_Bool active = FALSE;
#ifndef NDEBUG
   SCIP_Real constant = 0.0;
   SCIP_Bool orignegated;
#endif

   assert(var != NULL);
   assert(*var != NULL);
   assert(negated != NULL);
   assert(SCIPvarIsBinary(*var));

#ifndef NDEBUG
   orignegated = *negated;
#endif

   while( !active && *var != NULL )
   {
      switch( SCIPvarGetStatus(*var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( (*var)->data.original.transvar == NULL )
            return SCIP_OKAY;
         *var = (*var)->data.original.transvar;
         break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_FIXED:
         active = TRUE;
         break;

      case SCIP_VARSTATUS_MULTAGGR:
         /* handle multi-aggregated variables depending on one variable only (possibly caused by SCIPvarFlattenAggregationGraph()) */
         if ( (*var)->data.multaggr.nvars == 1 )
         {
            assert( (*var)->data.multaggr.vars != NULL );
            assert( (*var)->data.multaggr.scalars != NULL );
            assert( SCIPvarIsBinary((*var)->data.multaggr.vars[0]) );
            assert(!EPSZ((*var)->data.multaggr.scalars[0], 1e-06));

            /* if not all variables were fully propagated, it might happen that a variable is multi-aggregated to
             * another variable which needs to be fixed
             *
             * e.g. x = y - 1 => (x = 0 && y = 1)
             * e.g. x = y + 1 => (x = 1 && y = 0)
             *
             * is this special case we need to return the muti-aggregation
             */
            if( EPSEQ((*var)->data.multaggr.constant, -1.0, 1e-06) || (EPSEQ((*var)->data.multaggr.constant, 1.0, 1e-06) && EPSEQ((*var)->data.multaggr.scalars[0], 1.0, 1e-06)) )
            {
               assert(EPSEQ((*var)->data.multaggr.scalars[0], 1.0, 1e-06));
            }
            else
            {
               /* @note due to fixations, a multi-aggregation can have a constant of zero and a negative scalar or even
                *       a scalar in absolute value unequal to one, in this case this aggregation variable needs to be
                *       fixed to zero, but this should be done by another enforcement; so not depending on the scalar,
                *       we will return the aggregated variable;
                */
               if( !EPSEQ(REALABS((*var)->data.multaggr.scalars[0]), 1.0, 1e-06) )
               {
                  active = TRUE;
                  break;
               }

               /* @note it may also happen that the constant is larger than 1 or smaller than 0, in that case the
                *       aggregation variable needs to be fixed to one, but this should be done by another enforcement;
                *       so if this is the case, we will return the aggregated variable
                */
               assert(EPSZ((*var)->data.multaggr.constant, 1e-06) || EPSEQ((*var)->data.multaggr.constant, 1.0, 1e-06)
                  || EPSZ((*var)->data.multaggr.constant + (*var)->data.multaggr.scalars[0], 1e-06)
                  || EPSEQ((*var)->data.multaggr.constant + (*var)->data.multaggr.scalars[0], 1.0, 1e-06));

               if( !EPSZ((*var)->data.multaggr.constant, 1e-06) && !EPSEQ((*var)->data.multaggr.constant, 1.0, 1e-06) )
               {
                  active = TRUE;
                  break;
               }

               assert(EPSEQ((*var)->data.multaggr.scalars[0], 1.0, 1e-06) || EPSEQ((*var)->data.multaggr.scalars[0], -1.0, 1e-06));

               if( EPSZ((*var)->data.multaggr.constant, 1e-06) )
               {
                  /* if the scalar is negative, either the aggregation variable is already fixed to zero or has at
                   * least one uplock (that hopefully will enforce this fixation to zero); can it happen that this
                   * variable itself is multi-aggregated again?
                   */
                  assert(EPSEQ((*var)->data.multaggr.scalars[0], -1.0, 1e-06) ?
                     ((SCIPvarGetUbGlobal((*var)->data.multaggr.vars[0]) < 0.5) ||
                        SCIPvarGetNLocksUp((*var)->data.multaggr.vars[0]) > 0) : TRUE);
               }
               else
               {
                  assert(EPSEQ((*var)->data.multaggr.scalars[0], -1.0, 1e-06));
#ifndef NDEBUG
                  constant += (*negated) != orignegated ? -1.0 : 1.0;
#endif

                  *negated = !(*negated);
               }
               *var = (*var)->data.multaggr.vars[0];
               break;
            }
         }
         active = TRUE;
         break;

      case SCIP_VARSTATUS_AGGREGATED:  /* x = a'*x' + c'  =>  a*x + c == (a*a')*x' + (a*c' + c) */
         assert((*var)->data.aggregate.var != NULL);
         assert(EPSEQ((*var)->data.aggregate.scalar, 1.0, 1e-06) || EPSEQ((*var)->data.aggregate.scalar, -1.0, 1e-06));
         assert(EPSLE((*var)->data.aggregate.var->glbdom.ub - (*var)->data.aggregate.var->glbdom.lb, 1.0, 1e-06));
#ifndef NDEBUG
         constant += (*negated) != orignegated ? -(*var)->data.aggregate.constant : (*var)->data.aggregate.constant;
#endif

         *negated = ((*var)->data.aggregate.scalar > 0.0) ? *negated : !(*negated);
         *var = (*var)->data.aggregate.var;
         break;

      case SCIP_VARSTATUS_NEGATED:     /* x =  - x' + c'  =>  a*x + c ==   (-a)*x' + (a*c' + c) */
         assert((*var)->negatedvar != NULL);
#ifndef NDEBUG
         constant += (*negated) != orignegated ? -1.0 : 1.0;
#endif

         *negated = !(*negated);
         *var = (*var)->negatedvar;
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }
   }
   assert(active == (*var != NULL));

   if( active )
   {
      assert(SCIPvarIsBinary(*var));
      assert(EPSZ(constant, 1e-06) || EPSEQ(constant, 1.0, 1e-06));
      assert(EPSZ(constant, 1e-06) == ((*negated) == orignegated));

      return SCIP_OKAY;
   }
   else
   {
      SCIPerrorMessage("active variable path leads to NULL pointer\n");
      return SCIP_INVALIDDATA;
   }
}

/** transforms given variable, boundtype and bound to the corresponding active, fixed, or multi-aggregated variable
 *  values
 */
SCIP_RETCODE SCIPvarGetProbvarBound(
   SCIP_VAR**            var,                /**< pointer to problem variable */
   SCIP_Real*            bound,              /**< pointer to bound value to transform */
   SCIP_BOUNDTYPE*       boundtype           /**< pointer to type of bound: lower or upper bound */
   )
{
   assert(var != NULL);
   assert(*var != NULL);
   assert(bound != NULL);
   assert(boundtype != NULL);

   SCIPdebugMessage("get probvar bound %g of type %d of variable <%s>\n", *bound, *boundtype, (*var)->name);

   switch( SCIPvarGetStatus(*var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( (*var)->data.original.transvar == NULL )
      {
         SCIPerrorMessage("original variable has no transformed variable attached\n");
         return SCIP_INVALIDDATA;
      }
      *var = (*var)->data.original.transvar;
      SCIP_CALL( SCIPvarGetProbvarBound(var, bound, boundtype) );
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      /* handle multi-aggregated variables depending on one variable only (possibly caused by SCIPvarFlattenAggregationGraph()) */
      if ( (*var)->data.multaggr.nvars == 1 )
      {
         assert( (*var)->data.multaggr.vars != NULL );
         assert( (*var)->data.multaggr.scalars != NULL );
         assert( (*var)->data.multaggr.scalars[0] != 0.0 );

         (*bound) /= (*var)->data.multaggr.scalars[0];
         (*bound) -= (*var)->data.multaggr.constant/(*var)->data.multaggr.scalars[0];
         if ( (*var)->data.multaggr.scalars[0] < 0.0 )
         {
            if ( *boundtype == SCIP_BOUNDTYPE_LOWER )
               *boundtype = SCIP_BOUNDTYPE_UPPER;
            else
               *boundtype = SCIP_BOUNDTYPE_LOWER;
         }
         *var = (*var)->data.multaggr.vars[0];
         SCIP_CALL( SCIPvarGetProbvarBound(var, bound, boundtype) );
      }
      break;

   case SCIP_VARSTATUS_AGGREGATED:  /* x = a*y + c  ->  y = x/a - c/a */
      assert((*var)->data.aggregate.var != NULL);
      assert((*var)->data.aggregate.scalar != 0.0);

      (*bound) /= (*var)->data.aggregate.scalar;
      (*bound) -= (*var)->data.aggregate.constant/(*var)->data.aggregate.scalar;
      if( (*var)->data.aggregate.scalar < 0.0 )
      {
         if( *boundtype == SCIP_BOUNDTYPE_LOWER )
            *boundtype = SCIP_BOUNDTYPE_UPPER;
         else
            *boundtype = SCIP_BOUNDTYPE_LOWER;
      }
      *var = (*var)->data.aggregate.var;
      SCIP_CALL( SCIPvarGetProbvarBound(var, bound, boundtype) );
      break;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert((*var)->negatedvar != NULL);
      assert(SCIPvarGetStatus((*var)->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert((*var)->negatedvar->negatedvar == *var);
      (*bound) = (*var)->data.negate.constant - *bound;
      if( *boundtype == SCIP_BOUNDTYPE_LOWER )
         *boundtype = SCIP_BOUNDTYPE_UPPER;
      else
         *boundtype = SCIP_BOUNDTYPE_LOWER;
      *var = (*var)->negatedvar;
      SCIP_CALL( SCIPvarGetProbvarBound(var, bound, boundtype) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** transforms given variable and domain hole to the corresponding active, fixed, or multi-aggregated variable
 *  values
 */
SCIP_RETCODE SCIPvarGetProbvarHole(
   SCIP_VAR**            var,                /**< pointer to problem variable */
   SCIP_Real*            left,               /**< pointer to left bound of open interval in hole to transform */
   SCIP_Real*            right               /**< pointer to right bound of open interval in hole to transform */
   )
{
   assert(var != NULL);
   assert(*var != NULL);
   assert(left != NULL);
   assert(right != NULL);

   SCIPdebugMessage("get probvar hole (%g,%g) of variable <%s>\n", *left, *right, (*var)->name);

   switch( SCIPvarGetStatus(*var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( (*var)->data.original.transvar == NULL )
      {
         SCIPerrorMessage("original variable has no transformed variable attached\n");
         return SCIP_INVALIDDATA;
      }
      *var = (*var)->data.original.transvar;
      SCIP_CALL( SCIPvarGetProbvarHole(var, left, right) );
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_MULTAGGR:
      break;

   case SCIP_VARSTATUS_AGGREGATED:  /* x = a*y + c  ->  y = x/a - c/a */
      assert((*var)->data.aggregate.var != NULL);
      assert((*var)->data.aggregate.scalar != 0.0);

      /* scale back */
      (*left) /= (*var)->data.aggregate.scalar;
      (*right) /= (*var)->data.aggregate.scalar;

      /* shift back */
      (*left) -= (*var)->data.aggregate.constant/(*var)->data.aggregate.scalar;
      (*right) -= (*var)->data.aggregate.constant/(*var)->data.aggregate.scalar;

      *var = (*var)->data.aggregate.var;

      /* check if the  interval bounds have to swapped */
      if( (*var)->data.aggregate.scalar < 0.0 )
      {
         SCIP_CALL( SCIPvarGetProbvarHole(var, right, left) );
      }
      else
      {
         SCIP_CALL( SCIPvarGetProbvarHole(var, left, right) );
      }
      break;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert((*var)->negatedvar != NULL);
      assert(SCIPvarGetStatus((*var)->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert((*var)->negatedvar->negatedvar == *var);

      /* shift and scale back */
      (*left) = (*var)->data.negate.constant - (*left);
      (*right) = (*var)->data.negate.constant - (*right);

      *var = (*var)->negatedvar;

      /* through the negated variable the left and right interval bound have to swapped */
      SCIP_CALL( SCIPvarGetProbvarHole(var, right, left) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** transforms given variable, scalar and constant to the corresponding active, fixed, or
 *  multi-aggregated variable, scalar and constant; if the variable resolves to a fixed variable,
 *  "scalar" will be 0.0 and the value of the sum will be stored in "constant"; a multi-aggregation
 *  with only one active variable (this can happen due to fixings after the multi-aggregation),
 *  is treated like an aggregation; if the multi-aggregation constant is infinite, "scalar" will be 0.0
 */
SCIP_RETCODE SCIPvarGetProbvarSum(
   SCIP_VAR**            var,                /**< pointer to problem variable x in sum a*x + c */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            scalar,             /**< pointer to scalar a in sum a*x + c */
   SCIP_Real*            constant            /**< pointer to constant c in sum a*x + c */
   )
{
   assert(var != NULL);
   assert(scalar != NULL);
   assert(constant != NULL);

   while( *var != NULL )
   {
      switch( SCIPvarGetStatus(*var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( (*var)->data.original.transvar == NULL )
         {
            SCIPerrorMessage("original variable has no transformed variable attached\n");
            return SCIP_INVALIDDATA;
         }
         *var = (*var)->data.original.transvar;
         break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
         return SCIP_OKAY;

      case SCIP_VARSTATUS_FIXED:       /* x = c'          =>  a*x + c ==             (a*c' + c) */
         if( !SCIPsetIsInfinity(set, (*constant)) && !SCIPsetIsInfinity(set, -(*constant)) )
         {
            if( SCIPsetIsInfinity(set, (*var)->glbdom.lb) || SCIPsetIsInfinity(set, -((*var)->glbdom.lb)) )
            {
               assert(*scalar != 0.0);
               if( (*scalar) * (*var)->glbdom.lb > 0.0 )
                  (*constant) = SCIPsetInfinity(set);
               else
                  (*constant) = -SCIPsetInfinity(set);
            }
            else
               (*constant) += *scalar * (*var)->glbdom.lb;
         }
#ifndef NDEBUG
         else
         {
            assert(!SCIPsetIsInfinity(set, (*constant)) || !((*scalar) * (*var)->glbdom.lb < 0.0 &&
                  (SCIPsetIsInfinity(set, (*var)->glbdom.lb) || SCIPsetIsInfinity(set, -((*var)->glbdom.lb)))));
            assert(!SCIPsetIsInfinity(set, -(*constant)) || !((*scalar) * (*var)->glbdom.lb > 0.0 &&
                  (SCIPsetIsInfinity(set, (*var)->glbdom.lb) || SCIPsetIsInfinity(set, -((*var)->glbdom.lb)))));
         }
#endif
         *scalar = 0.0;
         return SCIP_OKAY;

      case SCIP_VARSTATUS_MULTAGGR:
         /* handle multi-aggregated variables depending on one variable only (possibly caused by SCIPvarFlattenAggregationGraph()) */
         if ( (*var)->data.multaggr.nvars == 1 )
         {
            assert((*var)->data.multaggr.vars != NULL);
            assert((*var)->data.multaggr.scalars != NULL);
            assert((*var)->data.multaggr.vars[0] != NULL);
            if( !SCIPsetIsInfinity(set, (*constant)) && !SCIPsetIsInfinity(set, -(*constant)) )
            {
               /* the multi-aggregation constant can be infinite, if one of the multi-aggregation variables
                * was fixed to +/-infinity; ensure that the constant is set to +/-infinity, too, and the scalar
                * is set to 0.0, because the multi-aggregated variable can be seen as fixed, too
                */
               if( SCIPsetIsInfinity(set, (*var)->data.multaggr.constant)
                  || SCIPsetIsInfinity(set, -((*var)->data.multaggr.constant)) )
               {
                  if( (*scalar) * (*var)->data.multaggr.constant > 0 )
                  {
                     assert(!SCIPsetIsInfinity(set, -(*constant)));
                     (*constant) = SCIPsetInfinity(set);
                  }
                  else
                  {
                     assert(!SCIPsetIsInfinity(set, *constant));
                     (*constant) = -SCIPsetInfinity(set);
                  }
                  (*scalar) = 0.0;
               }
               else
                  (*constant) += *scalar * (*var)->data.multaggr.constant;
            }
            (*scalar) *= (*var)->data.multaggr.scalars[0];
            *var = (*var)->data.multaggr.vars[0];
            break;
         }
         return SCIP_OKAY;

      case SCIP_VARSTATUS_AGGREGATED:  /* x = a'*x' + c'  =>  a*x + c == (a*a')*x' + (a*c' + c) */
         assert((*var)->data.aggregate.var != NULL);
         assert(!SCIPsetIsInfinity(set, (*var)->data.aggregate.constant)
            && !SCIPsetIsInfinity(set, (*var)->data.aggregate.constant));
         if( !SCIPsetIsInfinity(set, (*constant)) && !SCIPsetIsInfinity(set, -(*constant)) )
            (*constant) += *scalar * (*var)->data.aggregate.constant;
         (*scalar) *= (*var)->data.aggregate.scalar;
         *var = (*var)->data.aggregate.var;
         break;

      case SCIP_VARSTATUS_NEGATED:     /* x =  - x' + c'  =>  a*x + c ==   (-a)*x' + (a*c' + c) */
         assert((*var)->negatedvar != NULL);
         assert(SCIPvarGetStatus((*var)->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert((*var)->negatedvar->negatedvar == *var);
         assert(!SCIPsetIsInfinity(set, (*var)->data.negate.constant)
            && !SCIPsetIsInfinity(set, (*var)->data.negate.constant));
         if( !SCIPsetIsInfinity(set, (*constant)) && !SCIPsetIsInfinity(set, -(*constant)) )
            (*constant) += *scalar * (*var)->data.negate.constant;
         (*scalar) *= -1.0;
         *var = (*var)->negatedvar;
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
	 SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }
   *scalar = 0.0;

   return SCIP_OKAY;
}

/** retransforms given variable, scalar and constant to the corresponding original variable, scalar
 *  and constant, if possible; if the retransformation is impossible, NULL is returned as variable
 */
SCIP_RETCODE SCIPvarGetOrigvarSum(
   SCIP_VAR**            var,                /**< pointer to problem variable x in sum a*x + c */
   SCIP_Real*            scalar,             /**< pointer to scalar a in sum a*x + c */
   SCIP_Real*            constant            /**< pointer to constant c in sum a*x + c */
   )
{
   SCIP_VAR* parentvar;

   assert(var != NULL);
   assert(*var != NULL);
   assert(scalar != NULL);
   assert(constant != NULL);

   while( !SCIPvarIsOriginal(*var) )
   {
      /* if the variable has no parent variables, it was generated during solving and has no corresponding original
       * var
       */
      if( (*var)->nparentvars == 0 )
      {
         /* negated variables do not need to have a parent variables, and negated variables can exist in original
          * space
          */
         if( SCIPvarGetStatus(*var) == SCIP_VARSTATUS_NEGATED &&
            ((*var)->negatedvar->nparentvars == 0 || (*var)->negatedvar->parentvars[0] != *var) )
         {
            *scalar *= -1.0;
            *constant -= (*var)->data.negate.constant * (*scalar);
            *var = (*var)->negatedvar;

            continue;
         }
         /* if the variables does not have any parent the variables was created during solving and has no original
          * counterpart
          */
         else
         {
            *var = NULL;

            return SCIP_OKAY;
         }
      }

      /* follow the link to the first parent variable */
      parentvar = (*var)->parentvars[0];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         break;

      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         SCIPerrorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         return SCIP_INVALIDDATA;

      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + b  ->  y = (x-b)/a,  s*y + c = (s/a)*x + c-b*s/a */
         assert(parentvar->data.aggregate.var == *var);
         assert(parentvar->data.aggregate.scalar != 0.0);
         *scalar /= parentvar->data.aggregate.scalar;
         *constant -= parentvar->data.aggregate.constant * (*scalar);
         break;

      case SCIP_VARSTATUS_NEGATED: /* x = b - y  ->  y = b - x,  s*y + c = -s*x + c+b*s */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         *scalar *= -1.0;
         *constant -= parentvar->data.negate.constant * (*scalar);
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }

      assert( parentvar != NULL );
      *var = parentvar;
   }

   return SCIP_OKAY;
}

/** returns whether the given variable is the direct counterpart of an original problem variable */
SCIP_Bool SCIPvarIsTransformedOrigvar(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_VAR* parentvar;
   assert(var != NULL);

   if( !SCIPvarIsTransformed(var) || var->nparentvars < 1 )
      return FALSE;

   assert(var->parentvars != NULL);
   parentvar = var->parentvars[0];
   assert(parentvar != NULL);

   /* we follow the aggregation tree to the root unless an original variable has been found - the first entries in the parentlist are candidates */
   while( parentvar->nparentvars >= 1 && SCIPvarGetStatus(parentvar) != SCIP_VARSTATUS_ORIGINAL )
      parentvar = parentvar->parentvars[0];
   assert( parentvar != NULL );

   return ( SCIPvarGetStatus(parentvar) == SCIP_VARSTATUS_ORIGINAL );
}

/** gets objective value of variable in current SCIP_LP; the value can be different from the objective value stored in
 *  the variable's own data due to diving, that operate only on the LP without updating the variables
 */
SCIP_Real SCIPvarGetObjLP(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   /* get bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      return SCIPvarGetObjLP(var->data.original.transvar);

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetObj(var->data.col);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      return var->obj;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      return var->data.aggregate.scalar * SCIPvarGetObjLP(var->data.aggregate.var);

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot get the objective value of a multiple aggregated variable\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return -SCIPvarGetObjLP(var->negatedvar);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** gets lower bound of variable in current SCIP_LP; the bound can be different from the bound stored in the variable's own
 *  data due to diving or conflict analysis, that operate only on the LP without updating the variables
 */
SCIP_Real SCIPvarGetLbLP(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* get bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      return SCIPvarGetLbLP(var->data.original.transvar, set);

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetLb(var->data.col);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      return var->locdom.lb;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( (var->data.aggregate.scalar > 0.0 && SCIPsetIsInfinity(set, -SCIPvarGetLbLP(var->data.aggregate.var, set)))
         || (var->data.aggregate.scalar < 0.0 && SCIPsetIsInfinity(set, SCIPvarGetUbLP(var->data.aggregate.var, set))) )
      {
         return -SCIPsetInfinity(set);
      }
      else if( var->data.aggregate.scalar > 0.0 )
      {
         /* a > 0 -> get lower bound of y */
         return var->data.aggregate.scalar * SCIPvarGetLbLP(var->data.aggregate.var, set) + var->data.aggregate.constant;
      }
      else if( var->data.aggregate.scalar < 0.0 )
      {
         /* a < 0 -> get upper bound of y */
         return var->data.aggregate.scalar * SCIPvarGetUbLP(var->data.aggregate.var, set) + var->data.aggregate.constant;
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         SCIPABORT();
         return SCIP_INVALID; /*lint !e527*/
      }

   case SCIP_VARSTATUS_MULTAGGR:
      /**@todo get the sides of the corresponding linear constraint */
      SCIPerrorMessage("getting the bounds of a multiple aggregated variable is not implemented yet\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetUbLP(var->negatedvar, set);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** gets upper bound of variable in current SCIP_LP; the bound can be different from the bound stored in the variable's own
 *  data due to diving or conflict analysis, that operate only on the LP without updating the variables
 */
SCIP_Real SCIPvarGetUbLP(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* get bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      return SCIPvarGetUbLP(var->data.original.transvar, set);

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetUb(var->data.col);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      return var->locdom.ub;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( (var->data.aggregate.scalar > 0.0 && SCIPsetIsInfinity(set, SCIPvarGetUbLP(var->data.aggregate.var, set)))
         || (var->data.aggregate.scalar < 0.0 && SCIPsetIsInfinity(set, -SCIPvarGetLbLP(var->data.aggregate.var, set))) )
      {
         return SCIPsetInfinity(set);
      }
      if( var->data.aggregate.scalar > 0.0 )
      {
         /* a > 0 -> get upper bound of y */
         return var->data.aggregate.scalar * SCIPvarGetUbLP(var->data.aggregate.var, set) + var->data.aggregate.constant;
      }
      else if( var->data.aggregate.scalar < 0.0 )
      {
         /* a < 0 -> get lower bound of y */
         return var->data.aggregate.scalar * SCIPvarGetLbLP(var->data.aggregate.var, set) + var->data.aggregate.constant;
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         SCIPABORT();
         return SCIP_INVALID; /*lint !e527*/
      }

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot get the bounds of a multi-aggregated variable.\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetLbLP(var->negatedvar, set);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** gets primal LP solution value of variable */
SCIP_Real SCIPvarGetLPSol_rec(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIP_INVALID;
      return SCIPvarGetLPSol(var->data.original.transvar);

   case SCIP_VARSTATUS_LOOSE:
      return SCIPvarGetBestBoundLocal(var);

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetPrimsol(var->data.col);

   case SCIP_VARSTATUS_FIXED:
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      return var->locdom.lb;

   case SCIP_VARSTATUS_AGGREGATED:
   {
      SCIP_Real lpsolval;

      assert(var->data.aggregate.var != NULL);
      lpsolval = SCIPvarGetLPSol(var->data.aggregate.var);

      /* a correct implementation would need to check the value of var->data.aggregate.var for infinity and return the
       * corresponding infinity value instead of performing an arithmetical transformation (compare method
       * SCIPvarGetLbLP()); however, we do not want to introduce a SCIP or SCIP_SET pointer to this method, since it is
       * (or is called by) a public interface method; instead, we only assert that values are finite
       * w.r.t. SCIP_DEFAULT_INFINITY, which seems to be true in our regression tests; note that this may yield false
       * positives and negatives if the parameter <numerics/infinity> is modified by the user
       */
      assert(lpsolval > -SCIP_DEFAULT_INFINITY);
      assert(lpsolval < +SCIP_DEFAULT_INFINITY);
      return var->data.aggregate.scalar * lpsolval + var->data.aggregate.constant;
   }
   case SCIP_VARSTATUS_MULTAGGR:
   {
      SCIP_Real primsol;
      int i;

      assert(!var->donotmultaggr);
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      /* Due to method SCIPvarFlattenAggregationGraph(), this assert is no longer correct
       * assert(var->data.multaggr.nvars >= 2); 
       */
      primsol = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         primsol += var->data.multaggr.scalars[i] * SCIPvarGetLPSol(var->data.multaggr.vars[i]);
      return primsol;
   }
   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetLPSol(var->negatedvar);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** gets primal NLP solution value of variable */
SCIP_Real SCIPvarGetNLPSol_rec(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real solval;
   int i;

   assert(var != NULL);

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      return SCIPvarGetNLPSol(var->data.original.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
         return var->nlpsol;

   case SCIP_VARSTATUS_FIXED:
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetUbGlobal(var));  /*lint !e777*/
      assert(SCIPvarGetLbLocal(var) == SCIPvarGetUbLocal(var));    /*lint !e777*/
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetLbLocal(var));   /*lint !e777*/
      return SCIPvarGetLbGlobal(var);

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      solval = SCIPvarGetNLPSol(var->data.aggregate.var);
      return var->data.aggregate.scalar * solval + var->data.aggregate.constant;

   case SCIP_VARSTATUS_MULTAGGR:
      solval = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         solval += var->data.multaggr.scalars[i] * SCIPvarGetNLPSol(var->data.multaggr.vars[i]);
      return solval;

   case SCIP_VARSTATUS_NEGATED:
      solval = SCIPvarGetNLPSol(var->negatedvar);
      return var->data.negate.constant - solval;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** gets pseudo solution value of variable at current node */
static
SCIP_Real SCIPvarGetPseudoSol_rec(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real pseudosol;
   int i;

   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIP_INVALID;
      return SCIPvarGetPseudoSol(var->data.original.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPvarGetBestBoundLocal(var);

   case SCIP_VARSTATUS_FIXED:
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      return var->locdom.lb;

   case SCIP_VARSTATUS_AGGREGATED:
   {
      SCIP_Real pseudosolval;
      assert(var->data.aggregate.var != NULL);
      /* a correct implementation would need to check the value of var->data.aggregate.var for infinity and return the
       * corresponding infinity value instead of performing an arithmetical transformation (compare method
       * SCIPvarGetLbLP()); however, we do not want to introduce a SCIP or SCIP_SET pointer to this method, since it is
       * (or is called by) a public interface method; instead, we only assert that values are finite
       * w.r.t. SCIP_DEFAULT_INFINITY, which seems to be true in our regression tests; note that this may yield false
       * positives and negatives if the parameter <numerics/infinity> is modified by the user
       */
      pseudosolval = SCIPvarGetPseudoSol(var->data.aggregate.var);
      assert(pseudosolval > -SCIP_DEFAULT_INFINITY);
      assert(pseudosolval < +SCIP_DEFAULT_INFINITY);
      return var->data.aggregate.scalar * pseudosolval + var->data.aggregate.constant;
   }
   case SCIP_VARSTATUS_MULTAGGR:
      assert(!var->donotmultaggr);
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      /* Due to method SCIPvarFlattenAggregationGraph(), this assert is no longer correct
       * assert(var->data.multaggr.nvars >= 2); 
       */
      pseudosol = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         pseudosol += var->data.multaggr.scalars[i] * SCIPvarGetPseudoSol(var->data.multaggr.vars[i]);
      return pseudosol;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetPseudoSol(var->negatedvar);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** gets current LP or pseudo solution value of variable */
SCIP_Real SCIPvarGetSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             getlpval            /**< should the LP solution value be returned? */
   )
{
   if( getlpval )
      return SCIPvarGetLPSol(var);
   else
      return SCIPvarGetPseudoSol(var);
}

/** remembers the current solution as root solution in the problem variables */
void SCIPvarStoreRootSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             roothaslp           /**< is the root solution from LP? */
   )
{
   assert(var != NULL);

   var->rootsol = SCIPvarGetSol(var, roothaslp);
}

/** updates the current solution as best root solution of the given variable if it is better */
void SCIPvarUpdateBestRootSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             rootsol,            /**< root solution value */
   SCIP_Real             rootredcost,        /**< root reduced cost */
   SCIP_Real             rootlpobjval        /**< objective value of the root LP */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* if reduced cost are zero nothing to update */
   if( SCIPsetIsDualfeasZero(set, rootredcost) )
      return;

   /* check if we have already a best combination stored */
   if( !SCIPsetIsDualfeasZero(set, var->bestrootredcost) )
   {
      SCIP_Real currcutoffbound;
      SCIP_Real cutoffbound;
      SCIP_Real bound;

      /* compute the cutoff bound which would improve the corresponding bound with the current stored root solution,
       * root reduced cost, and root LP objective value combination
       */
      if( var->bestrootredcost > 0.0 )
         bound = SCIPvarGetUbGlobal(var);
      else
         bound = SCIPvarGetLbGlobal(var);

      currcutoffbound = (bound - var->bestrootsol) * var->bestrootredcost + var->bestrootlpobjval;

      /* compute the cutoff bound which would improve the corresponding bound with new root solution, root reduced
       * cost, and root LP objective value combination
       */
      if( rootredcost > 0.0 )
         bound = SCIPvarGetUbGlobal(var);
      else
         bound = SCIPvarGetLbGlobal(var);

      cutoffbound = (bound - rootsol) * rootredcost + rootlpobjval;

      /* check if an improving root solution, root reduced cost, and root LP objective value is at hand */
      if( cutoffbound > currcutoffbound )
      {
         SCIPsetDebugMsg(set, "-> <%s> update potential cutoff bound <%g> -> <%g>\n",
            SCIPvarGetName(var), currcutoffbound, cutoffbound);

         var->bestrootsol = rootsol;
         var->bestrootredcost = rootredcost;
         var->bestrootlpobjval = rootlpobjval;
      }
   }
   else
   {
      SCIPsetDebugMsg(set, "-> <%s> initialize best root reduced cost information\n", SCIPvarGetName(var));
      SCIPsetDebugMsg(set, "   -> rootsol <%g>\n", rootsol);
      SCIPsetDebugMsg(set, "   -> rootredcost <%g>\n", rootredcost);
      SCIPsetDebugMsg(set, "   -> rootlpobjval <%g>\n", rootlpobjval);

      var->bestrootsol = rootsol;
      var->bestrootredcost = rootredcost;
      var->bestrootlpobjval = rootlpobjval;
   }
}

/** returns the solution of the variable in the last root node's relaxation, if the root relaxation is not yet
 *  completely solved, zero is returned
 */
SCIP_Real SCIPvarGetRootSol(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real rootsol;
   int i;

   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      return SCIPvarGetRootSol(var->data.original.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return var->rootsol;

   case SCIP_VARSTATUS_FIXED:
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      return var->locdom.lb;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      /* a correct implementation would need to check the value of var->data.aggregate.var for infinity and return the
       * corresponding infinity value instead of performing an arithmetical transformation (compare method
       * SCIPvarGetLbLP()); however, we do not want to introduce a SCIP or SCIP_SET pointer to this method, since it is
       * (or is called by) a public interface method; instead, we only assert that values are finite
       * w.r.t. SCIP_DEFAULT_INFINITY, which seems to be true in our regression tests; note that this may yield false
       * positives and negatives if the parameter <numerics/infinity> is modified by the user
       */
      assert(SCIPvarGetRootSol(var->data.aggregate.var) > -SCIP_DEFAULT_INFINITY);
      assert(SCIPvarGetRootSol(var->data.aggregate.var) < +SCIP_DEFAULT_INFINITY);
      return var->data.aggregate.scalar * SCIPvarGetRootSol(var->data.aggregate.var) + var->data.aggregate.constant;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(!var->donotmultaggr);
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      /* Due to method SCIPvarFlattenAggregationGraph(), this assert is no longer correct
       * assert(var->data.multaggr.nvars >= 2); 
       */
      rootsol = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         rootsol += var->data.multaggr.scalars[i] * SCIPvarGetRootSol(var->data.multaggr.vars[i]);
      return rootsol;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetRootSol(var->negatedvar);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** returns for given variable the reduced cost */
static
SCIP_Real getImplVarRedcost(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             varfixing,          /**< FALSE if for x == 0, TRUE for x == 1 */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
   {
      SCIP_COL* col;
      SCIP_Real primsol;
      SCIP_BASESTAT basestat;
      SCIP_Bool lpissolbasic;

      col = SCIPvarGetCol(var);
      assert(col != NULL);

      basestat = SCIPcolGetBasisStatus(col);
      lpissolbasic = SCIPlpIsSolBasic(lp);
      primsol = SCIPcolGetPrimsol(col);

      if( (lpissolbasic && (basestat == SCIP_BASESTAT_LOWER || basestat == SCIP_BASESTAT_UPPER)) ||
         (!lpissolbasic && (SCIPsetIsFeasEQ(set, SCIPvarGetLbLocal(var), primsol) || SCIPsetIsFeasEQ(set, SCIPvarGetUbLocal(var), primsol))) )
      {
         SCIP_Real redcost = SCIPcolGetRedcost(col, stat, lp);

         assert(((!lpissolbasic && SCIPsetIsFeasEQ(set, SCIPvarGetLbLocal(var), primsol)) ||
               (lpissolbasic && basestat == SCIP_BASESTAT_LOWER)) ? (!SCIPsetIsDualfeasNegative(set, redcost) ||
                  SCIPsetIsFeasEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))) : TRUE);
         assert(((!lpissolbasic && SCIPsetIsFeasEQ(set, SCIPvarGetUbLocal(var), primsol)) ||
               (lpissolbasic && basestat == SCIP_BASESTAT_UPPER)) ? (!SCIPsetIsDualfeasPositive(set, redcost) ||
                  SCIPsetIsFeasEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))) : TRUE);

         if( (varfixing && ((lpissolbasic && basestat == SCIP_BASESTAT_LOWER) ||
                  (!lpissolbasic && SCIPsetIsFeasEQ(set, SCIPvarGetLbLocal(var), primsol)))) ||
            (!varfixing && ((lpissolbasic && basestat == SCIP_BASESTAT_UPPER) ||
                  (!lpissolbasic && SCIPsetIsFeasEQ(set, SCIPvarGetUbLocal(var), primsol)))) )
            return redcost;
         else
            return 0.0;
      }

      return 0.0;
   }

   return 0.0;
}

#define MAX_CLIQUELENGTH 50
/** returns for the given binary variable the reduced cost which are given by the variable itself and its implication if
 *  the binary variable is fixed to the given value
 */
SCIP_Real SCIPvarGetImplRedcost(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             varfixing,          /**< FALSE if for x == 0, TRUE for x == 1 */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem, or NULL */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Real implredcost;
   int ncliques;
   int nvars;

   assert(SCIPvarIsBinary(var));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   /* get reduced cost of given variable */
   implredcost = getImplVarRedcost(var, set, varfixing, stat, lp);

#ifdef SCIP_MORE_DEBUG
   SCIPsetDebugMsg(set, "variable <%s> itself has reduced cost of %g\n", SCIPvarGetName(var), implredcost);
#endif

   /* the following algorithm is expensive */
   ncliques = SCIPvarGetNCliques(var, varfixing);

   if( ncliques > 0 )
   {
      SCIP_CLIQUE** cliques;
      SCIP_CLIQUE* clique;
      SCIP_VAR** clqvars;
      SCIP_VAR** probvars;
      SCIP_VAR* clqvar;
      SCIP_Bool* clqvalues;
      int* entries;
      int* ids;
      SCIP_Real redcost;
      SCIP_Bool cleanedup;
      int nclqvars;
      int nentries;
      int nids;
      int id;
      int c;
      int v;

      assert(prob != NULL);
      assert(SCIPprobIsTransformed(prob));

      nentries = SCIPprobGetNVars(prob) - SCIPprobGetNContVars(prob) + 1;

      SCIP_CALL_ABORT( SCIPsetAllocBufferArray(set, &ids, nentries) );
      nids = 0;
      SCIP_CALL_ABORT( SCIPsetAllocCleanBufferArray(set, &entries, nentries) );

      cliques = SCIPvarGetCliques(var, varfixing);
      assert(cliques != NULL);

      for( c = ncliques - 1; c >= 0; --c )
      {
         clique = cliques[c];
         assert(clique != NULL);
         nclqvars = SCIPcliqueGetNVars(clique);
         assert(nclqvars > 0);

         if( nclqvars > MAX_CLIQUELENGTH )
            continue;

         clqvars = SCIPcliqueGetVars(clique);
         clqvalues = SCIPcliqueGetValues(clique);
         assert(clqvars != NULL);
         assert(clqvalues != NULL);

         cleanedup = SCIPcliqueIsCleanedUp(clique);

         for( v = nclqvars - 1; v >= 0; --v )
         {
            clqvar = clqvars[v];
            assert(clqvar != NULL);

            /* ignore binary variable which are fixed */
            if( clqvar != var && (cleanedup || SCIPvarIsActive(clqvar)) &&
               (SCIPvarGetLbLocal(clqvar) < 0.5 && SCIPvarGetUbLocal(clqvar) > 0.5) )
            {
               int probindex = SCIPvarGetProbindex(clqvar) + 1;
               assert(0 < probindex && probindex < nentries);

#if 0
               /* check that the variable was not yet visited or does not appear with two contradicting implications, ->
                * can appear since there is no guarantee that all these infeasible bounds were found
                */
               assert(!entries[probindex] || entries[probindex] == (clqvalues[v] ? probindex : -probindex));
#endif
               if( entries[probindex] == 0 )
               {
                  ids[nids] = probindex;
                  ++nids;

                  /* mark variable as visited */
                  entries[probindex] = (clqvalues[v] ? probindex : -probindex);
               }
            }
         }
      }

      probvars = SCIPprobGetVars(prob);
      assert(probvars != NULL);

      /* add all implied reduced cost */
      for( v = nids - 1; v >= 0; --v )
      {
         id = ids[v];
         assert(0 < id && id < nentries);
         assert(entries[id] != 0);
         assert(probvars[id - 1] != NULL);
         assert(SCIPvarIsActive(probvars[id - 1]));
         assert(SCIPvarIsBinary(probvars[id - 1]));
         assert(SCIPvarGetLbLocal(probvars[id - 1]) < 0.5 && SCIPvarGetUbLocal(probvars[id - 1]) > 0.5);

         if( (entries[id] > 0) != varfixing )
            redcost = getImplVarRedcost(probvars[id - 1], set, (entries[id] < 0), stat, lp);
         else
            redcost = -getImplVarRedcost(probvars[id - 1], set, (entries[id] < 0), stat, lp);

         if( (varfixing && SCIPsetIsDualfeasPositive(set, redcost)) || (!varfixing && SCIPsetIsDualfeasNegative(set, redcost)) )
            implredcost += redcost;

         /* reset entries clear buffer array */
         entries[id] = 0;
      }

      SCIPsetFreeCleanBufferArray(set, &entries);
      SCIPsetFreeBufferArray(set, &ids);
   }

#ifdef SCIP_MORE_DEBUG
   SCIPsetDebugMsg(set, "variable <%s> incl. cliques (%d) has implied reduced cost of %g\n", SCIPvarGetName(var), ncliques,
      implredcost);
#endif

   /* collect non-binary implication information */
   nvars = SCIPimplicsGetNImpls(var->implics, varfixing);

   if( nvars > 0 )
   {
      SCIP_VAR** vars;
      SCIP_VAR* implvar;
      SCIP_COL* col;
      SCIP_Real* bounds;
      SCIP_BOUNDTYPE* boundtypes;
      SCIP_Real redcost;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool lpissolbasic;
      int v;

      vars =  SCIPimplicsGetVars(var->implics, varfixing);
      boundtypes = SCIPimplicsGetTypes(var->implics, varfixing);
      bounds = SCIPimplicsGetBounds(var->implics, varfixing);
      lpissolbasic = SCIPlpIsSolBasic(lp);

      for( v = nvars - 1; v >= 0; --v )
      {
         implvar = vars[v];
         assert(implvar != NULL);

         lb = SCIPvarGetLbLocal(implvar);
         ub = SCIPvarGetUbLocal(implvar);

         /* ignore binary variable which are fixed or not of column status */
         if( SCIPvarGetStatus(implvar) != SCIP_VARSTATUS_COLUMN || SCIPsetIsFeasEQ(set, lb, ub) )
            continue;

         col = SCIPvarGetCol(implvar);
         assert(col != NULL);
         redcost = 0.0;

         /* solved lp with basis information or not? */
         if( lpissolbasic )
         {
            SCIP_BASESTAT basestat = SCIPcolGetBasisStatus(col);

            /* check if the implication is not not yet applied */
            if( basestat == SCIP_BASESTAT_LOWER && boundtypes[v] == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasGT(set, bounds[v], lb) )
            {
               redcost = SCIPcolGetRedcost(col, stat, lp);
               assert(!SCIPsetIsDualfeasNegative(set, redcost));

               if( !varfixing )
                  redcost *= (lb - bounds[v]);
               else
                  redcost *= (bounds[v] - lb);
            }
            else if( basestat == SCIP_BASESTAT_UPPER && boundtypes[v] == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasLT(set, bounds[v], ub) )
            {
               redcost = SCIPcolGetRedcost(col, stat, lp);
               assert(!SCIPsetIsDualfeasPositive(set, redcost));

               if( varfixing )
                  redcost *= (bounds[v] - ub);
               else
                  redcost *= (ub - bounds[v]);
            }
         }
         else
         {
            SCIP_Real primsol = SCIPcolGetPrimsol(col);

            /* check if the implication is not not yet applied */
            if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasEQ(set, lb, primsol) && SCIPsetIsFeasGT(set, bounds[v], lb) )
            {
               redcost = SCIPcolGetRedcost(col, stat, lp);
               assert(!SCIPsetIsDualfeasNegative(set, redcost));

               if( varfixing )
                  redcost *= (lb - bounds[v]);
               else
                  redcost *= (bounds[v] - lb);
            }
            else if( boundtypes[v] == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasEQ(set, ub, primsol) && SCIPsetIsFeasLT(set, bounds[v], ub) )
            {
               redcost = SCIPcolGetRedcost(col, stat, lp);
               assert(!SCIPsetIsDualfeasPositive(set, redcost));

               if( varfixing )
                  redcost *= (bounds[v] - ub);
               else
                  redcost *= (ub - bounds[v]);
            }
         }

         /* improve implied reduced cost */
         if( (varfixing && SCIPsetIsDualfeasPositive(set, redcost)) || (!varfixing && SCIPsetIsDualfeasNegative(set, redcost)) )
            implredcost += redcost;
      }
   }

#ifdef SCIP_MORE_DEBUG
   SCIPsetDebugMsg(set, "variable <%s> incl. cliques (%d) and implications (%d) has implied reduced cost of %g\n",
      SCIPvarGetName(var), ncliques, nvars, implredcost);
#endif

   return implredcost;
}

/** returns the best solution (w.r.t. root reduced cost propagation) of the variable in the root node's relaxation, if
 *  the root relaxation is not yet completely solved, zero is returned
 */
SCIP_Real SCIPvarGetBestRootSol(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real rootsol;
   int i;

   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      return SCIPvarGetBestRootSol(var->data.original.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return var->bestrootsol;

   case SCIP_VARSTATUS_FIXED:
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      return var->locdom.lb;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      /* a correct implementation would need to check the value of var->data.aggregate.var for infinity and return the
       * corresponding infinity value instead of performing an arithmetical transformation (compare method
       * SCIPvarGetLbLP()); however, we do not want to introduce a SCIP or SCIP_SET pointer to this method, since it is
       * (or is called by) a public interface method; instead, we only assert that values are finite
       * w.r.t. SCIP_DEFAULT_INFINITY, which seems to be true in our regression tests; note that this may yield false
       * positives and negatives if the parameter <numerics/infinity> is modified by the user
       */
      assert(SCIPvarGetBestRootSol(var->data.aggregate.var) > -SCIP_DEFAULT_INFINITY);
      assert(SCIPvarGetBestRootSol(var->data.aggregate.var) < +SCIP_DEFAULT_INFINITY);
      return var->data.aggregate.scalar * SCIPvarGetBestRootSol(var->data.aggregate.var) + var->data.aggregate.constant;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(!var->donotmultaggr);
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      /* Due to method SCIPvarFlattenAggregationGraph(), this assert is no longer correct
       * assert(var->data.multaggr.nvars >= 2);
       */
      rootsol = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         rootsol += var->data.multaggr.scalars[i] * SCIPvarGetBestRootSol(var->data.multaggr.vars[i]);
      return rootsol;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetBestRootSol(var->negatedvar);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the best reduced costs (w.r.t. root reduced cost propagation) of the variable in the root node's relaxation,
 *  if the root relaxation is not yet completely solved, or the variable was no column of the root LP, SCIP_INVALID is
 *  returned
 */
SCIP_Real SCIPvarGetBestRootRedcost(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIP_INVALID;
      return SCIPvarGetBestRootRedcost(var->data.original.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return var->bestrootredcost;

   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_MULTAGGR:
   case SCIP_VARSTATUS_NEGATED:
      return 0.0;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the best objective value (w.r.t. root reduced cost propagation) of the root LP which belongs the root
 *  reduced cost which is accessible via SCIPvarGetRootRedcost() or the variable was no column of the root LP,
 *  SCIP_INVALID is returned
 */
SCIP_Real SCIPvarGetBestRootLPObjval(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIP_INVALID;
      return SCIPvarGetBestRootLPObjval(var->data.original.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return var->bestrootlpobjval;

   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_MULTAGGR:
   case SCIP_VARSTATUS_NEGATED:
      return SCIP_INVALID;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** set the given solution as the best root solution w.r.t. root reduced cost propagation in the variables */
void SCIPvarSetBestRootSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             rootsol,            /**< root solution value */
   SCIP_Real             rootredcost,        /**< root reduced cost */
   SCIP_Real             rootlpobjval        /**< objective value of the root LP */
   )
{
   assert(var != NULL);

   var->bestrootsol = rootsol;
   var->bestrootredcost = rootredcost;
   var->bestrootlpobjval = rootlpobjval;
}

/** stores the solution value as relaxation solution in the problem variable */
SCIP_RETCODE SCIPvarSetRelaxSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_Real             solval,             /**< solution value in the current relaxation solution */
   SCIP_Bool             updateobj           /**< should the objective value be updated? */
   )
{
   assert(var != NULL);
   assert(relaxation != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      SCIP_CALL( SCIPvarSetRelaxSol(var->data.original.transvar, set, relaxation, solval, updateobj) );
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      if( updateobj )
         SCIPrelaxationSolObjAdd(relaxation, var->obj * (solval - var->relaxsol));
      var->relaxsol = solval;
      break;

   case SCIP_VARSTATUS_FIXED:
      if( !SCIPsetIsEQ(set, solval, var->glbdom.lb) )
      {
         SCIPerrorMessage("cannot set relaxation solution value for variable <%s> fixed to %.15g to different value %.15g\n",
            SCIPvarGetName(var), var->glbdom.lb, solval);
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      SCIP_CALL( SCIPvarSetRelaxSol(var->data.aggregate.var, set, relaxation, 
            (solval - var->data.aggregate.constant)/var->data.aggregate.scalar, updateobj) );
      break;
   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot set solution value for multiple aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      SCIP_CALL( SCIPvarSetRelaxSol(var->negatedvar, set, relaxation, var->data.negate.constant - solval, updateobj) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** returns the solution value of the problem variable in the relaxation solution
 *
 *  @todo Inline this function - similar to SCIPvarGetLPSol_rec.
 */
SCIP_Real SCIPvarGetRelaxSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_Real solvalsum;
   SCIP_Real solval;
   int i;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      return SCIPvarGetRelaxSol(var->data.original.transvar, set);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return var->relaxsol;

   case SCIP_VARSTATUS_FIXED:
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetUbGlobal(var));  /*lint !e777*/
      assert(SCIPvarGetLbLocal(var) == SCIPvarGetUbLocal(var));    /*lint !e777*/
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetLbLocal(var));   /*lint !e777*/
      return SCIPvarGetLbGlobal(var);

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      solval = SCIPvarGetRelaxSol(var->data.aggregate.var, set);
      if( SCIPsetIsInfinity(set, solval) || SCIPsetIsInfinity(set, -solval) )
      {
         if( var->data.aggregate.scalar * solval > 0.0 )
            return SCIPsetInfinity(set);
         if( var->data.aggregate.scalar * solval < 0.0 )
            return -SCIPsetInfinity(set);
      }
      return var->data.aggregate.scalar * solval + var->data.aggregate.constant;

   case SCIP_VARSTATUS_MULTAGGR:
      solvalsum = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         solval = SCIPvarGetRelaxSol(var->data.multaggr.vars[i], set);
         if( SCIPsetIsInfinity(set, solval) || SCIPsetIsInfinity(set, -solval) )
         {
            if( var->data.multaggr.scalars[i] * solval > 0.0 )
               return SCIPsetInfinity(set);
            if( var->data.multaggr.scalars[i] * solval < 0.0 )
               return -SCIPsetInfinity(set);
         }
         solvalsum += var->data.multaggr.scalars[i] * solval;
      }
      return solvalsum;

   case SCIP_VARSTATUS_NEGATED:
      solval = SCIPvarGetRelaxSol(var->negatedvar, set);
      if( SCIPsetIsInfinity(set, solval) )
         return -SCIPsetInfinity(set);
      if( SCIPsetIsInfinity(set, -solval) )
         return SCIPsetInfinity(set);
      return var->data.negate.constant - solval;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** returns the solution value of the transformed problem variable in the relaxation solution */
SCIP_Real SCIPvarGetRelaxSolTransVar(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

   return var->relaxsol;
}

/** stores the solution value as NLP solution in the problem variable */
SCIP_RETCODE SCIPvarSetNLPSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             solval              /**< solution value in the current NLP solution */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      SCIP_CALL( SCIPvarSetNLPSol(var->data.original.transvar, set, solval) );
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      var->nlpsol = solval;
      break;

   case SCIP_VARSTATUS_FIXED:
      if( !SCIPsetIsEQ(set, solval, var->glbdom.lb) )
      {
         SCIPerrorMessage("cannot set NLP solution value for variable <%s> fixed to %.15g to different value %.15g\n",
            SCIPvarGetName(var), var->glbdom.lb, solval);
         SCIPABORT();
         return SCIP_INVALIDCALL; /*lint !e527*/
      }
      break;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      SCIP_CALL( SCIPvarSetNLPSol(var->data.aggregate.var, set, (solval - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot set solution value for multiple aggregated variable\n");
      SCIPABORT();
      return SCIP_INVALIDCALL; /*lint !e527*/

   case SCIP_VARSTATUS_NEGATED:
      SCIP_CALL( SCIPvarSetNLPSol(var->negatedvar, set, var->data.negate.constant - solval) );
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_ERROR; /*lint !e527*/
   }

   return SCIP_OKAY;
}

/** returns a weighted average solution value of the variable in all feasible primal solutions found so far */
SCIP_Real SCIPvarGetAvgSol(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real avgsol;
   int i;

   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      return SCIPvarGetAvgSol(var->data.original.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      avgsol = var->primsolavg;
      avgsol = MAX(avgsol, var->glbdom.lb);
      avgsol = MIN(avgsol, var->glbdom.ub);
      return avgsol;

   case SCIP_VARSTATUS_FIXED:
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      return var->locdom.lb;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      return var->data.aggregate.scalar * SCIPvarGetAvgSol(var->data.aggregate.var)
         + var->data.aggregate.constant;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(!var->donotmultaggr);
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      /* Due to method SCIPvarFlattenAggregationGraph(), this assert is no longer correct
       * assert(var->data.multaggr.nvars >= 2); 
       */
      avgsol = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         avgsol += var->data.multaggr.scalars[i] * SCIPvarGetAvgSol(var->data.multaggr.vars[i]);
      return avgsol;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetAvgSol(var->negatedvar);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns solution value and index of variable lower bound that is closest to the variable's value in the given primal solution
 *  or current LP solution if no primal solution is given; returns an index of -1 if no variable lower bound is available
 */
void SCIPvarGetClosestVlb(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            closestvlb,         /**< pointer to store the value of the closest variable lower bound */
   int*                  closestvlbidx       /**< pointer to store the index of the closest variable lower bound */
   )
{
   int nvlbs;

   assert(var != NULL);
   assert(stat != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(closestvlb != NULL);
   assert(closestvlbidx != NULL);

   *closestvlbidx = -1;
   *closestvlb = SCIP_REAL_MIN;

   nvlbs = SCIPvarGetNVlbs(var);
   if( nvlbs > 0 )
   {
      SCIP_VAR** vlbvars;
      SCIP_Real* vlbcoefs;
      SCIP_Real* vlbconsts;
      int i;

      vlbvars = SCIPvarGetVlbVars(var);
      vlbcoefs = SCIPvarGetVlbCoefs(var);
      vlbconsts = SCIPvarGetVlbConstants(var);

      /* check for cached values */
      if( var->closestvblpcount == stat->lpcount && var->closestvlbidx != -1 && sol == NULL)
      {
         i = var->closestvlbidx;
         assert(0 <= i && i < nvlbs);
         assert(SCIPvarIsActive(vlbvars[i]));
         *closestvlbidx = i;
         *closestvlb = vlbcoefs[i] * SCIPvarGetLPSol(vlbvars[i]) + vlbconsts[i];
      }
      else
      {
         /* search best VUB */
         for( i = 0; i < nvlbs; i++ )
         {
            if( SCIPvarIsActive(vlbvars[i]) )
            {
               SCIP_Real vlbsol;

               vlbsol = vlbcoefs[i] * (sol == NULL ? SCIPvarGetLPSol(vlbvars[i]) : SCIPsolGetVal(sol, set, stat, vlbvars[i])) + vlbconsts[i];
               if( vlbsol > *closestvlb )
               {
                  *closestvlb = vlbsol;
                  *closestvlbidx = i;
               }
            }
         }

         if( sol == NULL )
         {
            /* update cached value */
            if( var->closestvblpcount != stat->lpcount )
               var->closestvubidx = -1;
            var->closestvlbidx = *closestvlbidx;
            var->closestvblpcount = stat->lpcount;
         }
      }
   }
}

/** returns solution value and index of variable upper bound that is closest to the variable's value in the given primal solution;
 *  or current LP solution if no primal solution is given; returns an index of -1 if no variable upper bound is available
 */
void SCIPvarGetClosestVub(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            closestvub,         /**< pointer to store the value of the closest variable upper bound */
   int*                  closestvubidx       /**< pointer to store the index of the closest variable upper bound */
   )
{
   int nvubs;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(closestvub != NULL);
   assert(closestvubidx != NULL);

   *closestvubidx = -1;
   *closestvub = SCIP_REAL_MAX;

   nvubs = SCIPvarGetNVubs(var);
   if( nvubs > 0 )
   {
      SCIP_VAR** vubvars;
      SCIP_Real* vubcoefs;
      SCIP_Real* vubconsts;
      int i;

      vubvars = SCIPvarGetVubVars(var);
      vubcoefs = SCIPvarGetVubCoefs(var);
      vubconsts = SCIPvarGetVubConstants(var);

      /* check for cached values */
      if( var->closestvblpcount == stat->lpcount && var->closestvubidx != -1 && sol == NULL)
      {
         i = var->closestvubidx;
         assert(0 <= i && i < nvubs);
         assert(SCIPvarIsActive(vubvars[i]));
         *closestvubidx = i;
         *closestvub = vubcoefs[i] * SCIPvarGetLPSol(vubvars[i]) + vubconsts[i];
      }
      else
      {
         /* search best VUB */
         for( i = 0; i < nvubs; i++ )
         {
            if( SCIPvarIsActive(vubvars[i]) )
            {
               SCIP_Real vubsol;

               vubsol = vubcoefs[i] * (sol == NULL ? SCIPvarGetLPSol(vubvars[i]) : SCIPsolGetVal(sol, set, stat, vubvars[i])) + vubconsts[i];
               if( vubsol < *closestvub )
               {
                  *closestvub = vubsol;
                  *closestvubidx = i;
               }
            }
         }

         if( sol == NULL )
         {
            /* update cached value */
            if( var->closestvblpcount != stat->lpcount )
               var->closestvlbidx = -1;
            var->closestvubidx = *closestvubidx;
            var->closestvblpcount = stat->lpcount;
         }
      }
   }
}

/** resolves variable to columns and adds them with the coefficient to the row */
SCIP_RETCODE SCIPvarAddToRow(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   int i;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(row != NULL);
   assert(!SCIPsetIsInfinity(set, REALABS(val)));

   SCIPsetDebugMsg(set, "adding coefficient %g<%s> to row <%s>\n", val, var->name, row->name);

   if ( SCIPsetIsZero(set, val) )
      return SCIP_OKAY;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
      {
         SCIPerrorMessage("cannot add untransformed original variable <%s> to LP row <%s>\n", var->name, row->name);
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarAddToRow(var->data.original.transvar, blkmem, set, stat, eventqueue, prob, lp, row, val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
      /* add globally fixed variables as constant */
      if( SCIPsetIsEQ(set, var->glbdom.lb, var->glbdom.ub) )
      {
         SCIP_CALL( SCIProwAddConstant(row, blkmem, set, stat, eventqueue, lp, val * var->glbdom.lb) );
         return SCIP_OKAY;
      }
      /* convert loose variable into column */
      SCIP_CALL( SCIPvarColumn(var, blkmem, set, stat, prob, lp) );
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      /*lint -fallthrough*/

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      assert(var->data.col->var == var);
      SCIP_CALL( SCIProwIncCoef(row, blkmem, set, eventqueue, lp, var->data.col, val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      assert(var->glbdom.lb == var->glbdom.ub); /*lint !e777*/
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      assert(var->locdom.lb == var->glbdom.lb); /*lint !e777*/
      assert(!SCIPsetIsInfinity(set, REALABS(var->locdom.lb)));
      SCIP_CALL( SCIProwAddConstant(row, blkmem, set, stat, eventqueue, lp, val * var->locdom.lb) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      SCIP_CALL( SCIPvarAddToRow(var->data.aggregate.var, blkmem, set, stat, eventqueue, prob, lp,
            row, var->data.aggregate.scalar * val) );
      SCIP_CALL( SCIProwAddConstant(row, blkmem, set, stat, eventqueue, lp, var->data.aggregate.constant * val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(!var->donotmultaggr);
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      /* Due to method SCIPvarFlattenAggregationGraph(), this assert is no longer correct
       * assert(var->data.multaggr.nvars >= 2); 
       */
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         SCIP_CALL( SCIPvarAddToRow(var->data.multaggr.vars[i], blkmem, set, stat, eventqueue, prob, lp,
               row, var->data.multaggr.scalars[i] * val) );
      }
      SCIP_CALL( SCIProwAddConstant(row, blkmem, set, stat, eventqueue, lp, var->data.multaggr.constant * val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIP_CALL( SCIPvarAddToRow(var->negatedvar, blkmem, set, stat, eventqueue, prob, lp, row, -val) );
      SCIP_CALL( SCIProwAddConstant(row, blkmem, set, stat, eventqueue, lp, var->data.negate.constant * val) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/* optionally, define this compiler flag to write complete variable histories to a file */
#ifdef SCIP_HISTORYTOFILE
SCIP_Longint counter = 0l;
const char* historypath="."; /* allows for user-defined path; use '.' for calling directory of SCIP */
#include "scip/scip.h"
#endif

/** updates the pseudo costs of the given variable and the global pseudo costs after a change of
 *  "solvaldelta" in the variable's solution value and resulting change of "objdelta" in the in the LP's objective value
 */
SCIP_RETCODE SCIPvarUpdatePseudocost(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   SCIP_Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   SCIP_Real             weight              /**< weight in (0,1] of this update in pseudo cost sum */
   )
{
   SCIP_Real oldrootpseudocosts;
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(stat != NULL);

   /* check if history statistics should be collected for a variable */
   if( !stat->collectvarhistory )
      return SCIP_OKAY;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
      {
         SCIPerrorMessage("cannot update pseudo costs of original untransformed variable\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarUpdatePseudocost(var->data.original.transvar, set, stat, solvaldelta, objdelta, weight) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      /* store old pseudo-costs for root LP best-estimate update */
      oldrootpseudocosts = SCIPvarGetMinPseudocostScore(var, stat, set, SCIPvarGetRootSol(var));

      /* update history */
      SCIPhistoryUpdatePseudocost(var->history, set, solvaldelta, objdelta, weight);
      SCIPhistoryUpdatePseudocost(var->historycrun, set, solvaldelta, objdelta, weight);
      SCIPhistoryUpdatePseudocost(stat->glbhistory, set, solvaldelta, objdelta, weight);
      SCIPhistoryUpdatePseudocost(stat->glbhistorycrun, set, solvaldelta, objdelta, weight);

      /* update root LP best-estimate */
      SCIP_CALL( SCIPstatUpdateVarRootLPBestEstimate(stat, set, var, oldrootpseudocosts) );

      /* append history to file */
#ifdef SCIP_HISTORYTOFILE
   {
      FILE* f;
      char filename[256];
      SCIP_NODE* currentnode;
      SCIP_NODE* parentnode;
      currentnode = SCIPgetFocusNode(set->scip);
      parentnode = SCIPnodeGetParent(currentnode);

      sprintf(filename, "%s/%s.pse", historypath, SCIPgetProbName(set->scip));
      f = fopen(filename, "a");
      if( NULL != f )
      {
         fprintf(f, "%lld %s \t %lld \t %lld \t %lld \t %d \t %15.9f \t %.3f\n",
            ++counter,
            SCIPvarGetName(var),
            SCIPnodeGetNumber(currentnode),
            parentnode != NULL ? SCIPnodeGetNumber(parentnode) : -1,
            SCIPgetNLPIterations(set->scip),
            SCIPgetDepth(set->scip),
            objdelta,
            solvaldelta);
         fclose(f);
      }
   }
#endif
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot update pseudo cost values of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      SCIP_CALL( SCIPvarUpdatePseudocost(var->data.aggregate.var, set, stat,
            solvaldelta/var->data.aggregate.scalar, objdelta, weight) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot update pseudo cost values of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      SCIP_CALL( SCIPvarUpdatePseudocost(var->negatedvar, set, stat, -solvaldelta, objdelta, weight) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** gets the variable's pseudo cost value for the given step size "solvaldelta" in the variable's LP solution value */
SCIP_Real SCIPvarGetPseudocost(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP_BRANCHDIR dir;

   assert(var != NULL);
   assert(stat != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIPhistoryGetPseudocost(stat->glbhistory, solvaldelta);
      else
         return SCIPvarGetPseudocost(var->data.original.transvar, stat, solvaldelta);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      dir = (solvaldelta >= 0.0 ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS);

      return SCIPhistoryGetPseudocostCount(var->history, dir) > 0.0
         ? SCIPhistoryGetPseudocost(var->history, solvaldelta)
         : SCIPhistoryGetPseudocost(stat->glbhistory, solvaldelta);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      return SCIPvarGetPseudocost(var->data.aggregate.var, stat, var->data.aggregate.scalar * solvaldelta);

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetPseudocost(var->negatedvar, stat, -solvaldelta);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** gets the variable's pseudo cost value for the given step size "solvaldelta" in the variable's LP solution value,
 *  only using the pseudo cost information of the current run
 */
SCIP_Real SCIPvarGetPseudocostCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP_BRANCHDIR dir;

   assert(var != NULL);
   assert(stat != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIPhistoryGetPseudocost(stat->glbhistorycrun, solvaldelta);
      else
         return SCIPvarGetPseudocostCurrentRun(var->data.original.transvar, stat, solvaldelta);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      dir = (solvaldelta >= 0.0 ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS);

      return SCIPhistoryGetPseudocostCount(var->historycrun, dir) > 0.0
         ? SCIPhistoryGetPseudocost(var->historycrun, solvaldelta)
         : SCIPhistoryGetPseudocost(stat->glbhistorycrun, solvaldelta);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      return SCIPvarGetPseudocostCurrentRun(var->data.aggregate.var, stat, var->data.aggregate.scalar * solvaldelta);

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetPseudocostCurrentRun(var->negatedvar, stat, -solvaldelta);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction */
SCIP_Real SCIPvarGetPseudocostCount(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetPseudocostCount(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetPseudocostCount(var->history, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetPseudocostCount(var->data.aggregate.var, dir);
      else
         return SCIPvarGetPseudocostCount(var->data.aggregate.var, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetPseudocostCount(var->negatedvar, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 */
SCIP_Real SCIPvarGetPseudocostCountCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetPseudocostCountCurrentRun(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetPseudocostCount(var->historycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetPseudocostCountCurrentRun(var->data.aggregate.var, dir);
      else
         return SCIPvarGetPseudocostCountCurrentRun(var->data.aggregate.var, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetPseudocostCountCurrentRun(var->negatedvar, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** compares both possible directions for rounding the given solution value and returns the minimum pseudo-costs of the variable */
SCIP_Real SCIPvarGetMinPseudocostScore(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             solval              /**< solution value, e.g., LP solution value */
   )
{
   SCIP_Real upscore;
   SCIP_Real downscore;
   SCIP_Real solvaldeltaup;
   SCIP_Real solvaldeltadown;

   /* LP root estimate only works for variables with fractional LP root solution */
   if( SCIPsetIsFeasIntegral(set, solval) )
      return 0.0;

   /* no min pseudo-cost score is calculated as long as the variable was not initialized in a direction */
   if( SCIPvarGetPseudocostCount(var, SCIP_BRANCHDIR_DOWNWARDS) < 1.0 || SCIPvarGetPseudocostCount(var, SCIP_BRANCHDIR_UPWARDS) < 1.0 )
      return 0.0;

   /* compute delta's to ceil and floor of root LP solution value */
   solvaldeltaup = SCIPsetCeil(set, solval) - solval;
   solvaldeltadown = SCIPsetFloor(set, solval) - solval;


   upscore = SCIPvarGetPseudocost(var, stat, solvaldeltaup);
   downscore = SCIPvarGetPseudocost(var, stat, solvaldeltadown);

   return MIN(upscore, downscore);
}

/** gets the an estimate of the variable's pseudo cost variance in direction \p dir */
SCIP_Real SCIPvarGetPseudocostVariance(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Bool             onlycurrentrun      /**< return pseudo cost variance only for current branch and bound run */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetPseudocostVariance(var->data.original.transvar, dir, onlycurrentrun);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      if( onlycurrentrun )
         return SCIPhistoryGetPseudocostVariance(var->historycrun, dir);
      else
         return SCIPhistoryGetPseudocostVariance(var->history, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetPseudocostVariance(var->data.aggregate.var, dir, onlycurrentrun);
      else
         return SCIPvarGetPseudocostVariance(var->data.aggregate.var, SCIPbranchdirOpposite(dir), onlycurrentrun);

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetPseudocostVariance(var->negatedvar, SCIPbranchdirOpposite(dir), onlycurrentrun);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** calculates a confidence bound for this variable under the assumption of normally distributed pseudo costs
 *
 *  The confidence bound \f$ \theta \geq 0\f$ denotes the interval borders \f$ [X - \theta, \ X + \theta]\f$, which contains
 *  the true pseudo costs of the variable, i.e., the expected value of the normal distribution, with a probability
 *  of 2 * clevel - 1.
 *
 *  @return value of confidence bound for this variable
 */
SCIP_Real SCIPvarCalcPscostConfidenceBound(
   SCIP_VAR*             var,                /**< variable in question */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BRANCHDIR        dir,                /**< the branching direction for the confidence bound */
   SCIP_Bool             onlycurrentrun,     /**< should only the current run be taken into account */
   SCIP_CONFIDENCELEVEL  clevel              /**< confidence level for the interval */
   )
{
   SCIP_Real confidencebound;

   confidencebound = SCIPvarGetPseudocostVariance(var, dir, onlycurrentrun);
   if( SCIPsetIsFeasPositive(set, confidencebound) )
   {
      SCIP_Real count;

      if( onlycurrentrun )
         count = SCIPvarGetPseudocostCountCurrentRun(var, dir);
      else
         count = SCIPvarGetPseudocostCount(var, dir);
      /* assertion is valid because variance is positive */
      assert(count >= 1.9);

      confidencebound /= count; /*lint !e414 division by zero can obviously not occur */
      confidencebound = sqrt(confidencebound);

      /* the actual, underlying distribution of the mean is a student-t-distribution with degrees of freedom equal to
       * the number of pseudo cost evaluations of this variable in the respective direction. */
      confidencebound *= SCIPstudentTGetCriticalValue(clevel, (int)SCIPsetFloor(set, count) - 1);
   }
   else
      confidencebound = 0.0;

   return confidencebound;
}

/** check if the current pseudo cost relative error in a direction violates the given threshold. The Relative
 *  Error is calculated at a specific confidence level
 */
SCIP_Bool SCIPvarIsPscostRelerrorReliable(
   SCIP_VAR*             var,                /**< variable in question */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             threshold,          /**< threshold for relative errors to be considered reliable (enough) */
   SCIP_CONFIDENCELEVEL  clevel              /**< a given confidence level */
   )
{
   SCIP_Real downsize;
   SCIP_Real upsize;
   SCIP_Real size;
   SCIP_Real relerrorup;
   SCIP_Real relerrordown;
   SCIP_Real relerror;

   /* check, if the pseudo cost score of the variable is reliable */
   downsize = SCIPvarGetPseudocostCountCurrentRun(var, SCIP_BRANCHDIR_DOWNWARDS);
   upsize = SCIPvarGetPseudocostCountCurrentRun(var, SCIP_BRANCHDIR_UPWARDS);
   size = MIN(downsize, upsize);

   /* Pseudo costs relative error can only be reliable if both directions have been tried at least twice */
   if( size <= 1.9 )
      return FALSE;

   /* use the relative error between the current mean pseudo cost value of the candidate and its upper
    * confidence interval bound at confidence level of 95% for individual variable reliability.
    * this is only possible if we have at least 2 measurements and therefore a valid variance estimate.
    */
   if( downsize >= 1.9 )
   {
      SCIP_Real normval;

      relerrordown = SCIPvarCalcPscostConfidenceBound(var, set, SCIP_BRANCHDIR_DOWNWARDS, TRUE, clevel);
      normval = SCIPvarGetPseudocostCurrentRun(var, stat, -1.0);
      normval = MAX(1.0, normval);

      relerrordown /= normval;
   }
   else
      relerrordown = 0.0;

   if( upsize >= 1.9 )
   {
      SCIP_Real normval;

      relerrorup = SCIPvarCalcPscostConfidenceBound(var, set, SCIP_BRANCHDIR_UPWARDS, TRUE, clevel);
      normval = SCIPvarGetPseudocostCurrentRun(var, stat, +1.0);
      normval = MAX(1.0, normval);
      relerrorup /= normval;
   }
   else
      relerrorup = 0.0;

   /* consider the relative error threshold violated, if it is violated in at least one branching direction */
   relerror = MAX(relerrorup, relerrordown);

   return (relerror <= threshold);
}

/** check if variable pseudo-costs have a significant difference in location. The significance depends on
 *  the choice of \p clevel and on the kind of tested hypothesis. The one-sided hypothesis, which
 *  should be rejected, is that fracy * mu_y >= fracx * mu_x, where mu_y and mu_x denote the
 *  unknown location means of the underlying pseudo-cost distributions of x and y.
 *
 *  This method is applied best if variable x has a better pseudo-cost score than y. The method hypothesizes that y were actually
 *  better than x (despite the current information), meaning that y can be expected to yield branching
 *  decisions as least as good as x in the long run. If the method returns TRUE, the current history information is
 *  sufficient to safely rely on the alternative hypothesis that x yields indeed a better branching score (on average)
 *  than y.
 *
 *  @note The order of x and y matters for the one-sided hypothesis
 *
 *  @note set \p onesided to FALSE if you are not sure which variable is better. The hypothesis tested then reads
 *        fracy * mu_y == fracx * mu_x vs the alternative hypothesis fracy * mu_y != fracx * mu_x.
 *
 *  @return TRUE if the hypothesis can be safely rejected at the given confidence level
 */
SCIP_Bool SCIPvarSignificantPscostDifference(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             varx,               /**< variable x */
   SCIP_Real             fracx,              /**< the fractionality of variable x */
   SCIP_VAR*             vary,               /**< variable y */
   SCIP_Real             fracy,              /**< the fractionality of variable y */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_CONFIDENCELEVEL  clevel,             /**< confidence level for rejecting hypothesis */
   SCIP_Bool             onesided            /**< should a one-sided hypothesis y >= x be tested? */
   )
{
   SCIP_Real meanx;
   SCIP_Real meany;
   SCIP_Real variancex;
   SCIP_Real variancey;
   SCIP_Real countx;
   SCIP_Real county;
   SCIP_Real tresult;
   SCIP_Real realdirection;


   if( varx == vary )
      return FALSE;

   countx = SCIPvarGetPseudocostCount(varx, dir);
   county = SCIPvarGetPseudocostCount(vary, dir);

   /* if not at least 2 measurements were taken, return FALSE */
   if( countx <= 1.9 || county <= 1.9 )
      return FALSE;

   realdirection = (dir == SCIP_BRANCHDIR_DOWNWARDS ? -1.0 : 1.0);

   meanx = fracx * SCIPvarGetPseudocost(varx, stat, realdirection);
   meany = fracy * SCIPvarGetPseudocost(vary, stat, realdirection);

   variancex = SQR(fracx) * SCIPvarGetPseudocostVariance(varx, dir, FALSE);
   variancey = SQR(fracy) * SCIPvarGetPseudocostVariance(vary, dir, FALSE);

   /* if there is no variance, the means are taken from a constant distribution */
   if( SCIPsetIsFeasEQ(set, variancex + variancey, 0.0) )
      return (onesided ? SCIPsetIsFeasGT(set, meanx, meany) : !SCIPsetIsFeasEQ(set, meanx, meany));

   tresult = SCIPcomputeTwoSampleTTestValue(meanx, meany, variancex, variancey, countx, county);

   /* for the two-sided hypothesis, just take the absolute of t */
   if( !onesided )
      tresult = REALABS(tresult);

   return (tresult >= SCIPstudentTGetCriticalValue(clevel, (int)(countx + county - 2)));
}

/** tests at a given confidence level whether the variable pseudo-costs only have a small probability to
 *  exceed a \p threshold. This is useful to determine if past observations provide enough evidence
 *  to skip an expensive strong-branching step if there is already a candidate that has been proven to yield an improvement
 *  of at least \p threshold.
 *
 *  @note use \p clevel to adjust the level of confidence. For SCIP_CONFIDENCELEVEL_MIN, the method returns TRUE if
 *        the estimated probability to exceed \p threshold is less than 25 %.
 *
 *  @see  SCIP_Confidencelevel for a list of available levels. The used probability limits refer to the one-sided levels
 *        of confidence.
 *
 *  @return TRUE if the variable pseudo-cost probabilistic model is likely to be smaller than \p threshold
 *          at the given confidence level \p clevel.
 */
SCIP_Bool SCIPvarPscostThresholdProbabilityTest(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< variable x */
   SCIP_Real             frac,               /**< the fractionality of variable x */
   SCIP_Real             threshold,          /**< the threshold to test against */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_CONFIDENCELEVEL  clevel              /**< confidence level for rejecting hypothesis */
   )
{
   SCIP_Real mean;
   SCIP_Real variance;
   SCIP_Real count;
   SCIP_Real realdirection;
   SCIP_Real probability;
   SCIP_Real problimit;

   count = SCIPvarGetPseudocostCount(var, dir);

   /* if not at least 2 measurements were taken, return FALSE */
   if( count <= 1.9 )
      return FALSE;

   realdirection = (dir == SCIP_BRANCHDIR_DOWNWARDS ? -1.0 : 1.0);

   mean = frac * SCIPvarGetPseudocost(var, stat, realdirection);
   variance = SQR(frac) * SCIPvarGetPseudocostVariance(var, dir, FALSE);

   /* if mean is at least threshold, it has at least a 50% probability to exceed threshold, we therefore return FALSE */
   if( SCIPsetIsFeasGE(set, mean, threshold) )
      return FALSE;

   /* if there is no variance, the means are taken from a constant distribution */
   if( SCIPsetIsFeasEQ(set, variance, 0.0) )
      return SCIPsetIsFeasLT(set, mean, threshold);

   /* obtain probability of a normally distributed random variable at given mean and variance to yield at most threshold */
   probability = SCIPnormalCDF(mean, variance, threshold);

   /* determine a probability limit corresponding to the given confidence level */
   switch( clevel )
   {
      case SCIP_CONFIDENCELEVEL_MIN:
         problimit = 0.75;
         break;
      case SCIP_CONFIDENCELEVEL_LOW:
         problimit = 0.875;
         break;
      case SCIP_CONFIDENCELEVEL_MEDIUM:
         problimit = 0.9;
         break;
      case SCIP_CONFIDENCELEVEL_HIGH:
         problimit = 0.95;
         break;
      case SCIP_CONFIDENCELEVEL_MAX:
         problimit = 0.975;
         break;
      default:
         problimit = -1;
         SCIPerrorMessage("Confidence level set to unknown value <%d>", (int)clevel);
         SCIPABORT();
         break;
   }

   return (probability >= problimit);
}

/** find the corresponding history entry if already existing, otherwise create new entry */
static
SCIP_RETCODE findValuehistoryEntry(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_SET*             set,                /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_HISTORY**        history             /**< pointer to store the value based history, or NULL */
   )
{
   assert(var != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(history != NULL);

   (*history) = NULL;

   if( var->valuehistory == NULL )
   {
      SCIP_CALL( SCIPvaluehistoryCreate(&var->valuehistory, blkmem) );
   }

   SCIP_CALL( SCIPvaluehistoryFind(var->valuehistory, blkmem, set, value, history) );

   return SCIP_OKAY;
}

/** check if value based history should be used */
static
SCIP_Bool useValuehistory(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   SCIP_SET*             set                 /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   )
{
   /* check if the domain value is unknown (not specific) */
   if( value == SCIP_UNKNOWN ) /*lint !e777*/
      return FALSE;

   assert(set != NULL);

   /* check if value based history should be collected */
   if( !set->history_valuebased )
      return FALSE;

   /* value based history is not collected for binary variable since the standard history already contains all information */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      return FALSE;

   /* value based history is not collected for continuous variables */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      return FALSE;

   return TRUE;
}

/** increases VSIDS of the variable by the given weight */
SCIP_RETCODE SCIPvarIncVSIDS(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_SET*             set,                /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   SCIP_Real             weight              /**< weight of this update in VSIDS */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   /* check if history statistics should be collected for a variable */
   if( !stat->collectvarhistory )
      return SCIP_OKAY;

   if( SCIPsetIsZero(set, weight) )
      return SCIP_OKAY;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
      {
         SCIPerrorMessage("cannot update VSIDS of original untransformed variable\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarIncVSIDS(var->data.original.transvar, blkmem, set, stat, dir, value, weight) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   {
      SCIPhistoryIncVSIDS(var->history, dir, weight);
      SCIPhistoryIncVSIDS(var->historycrun, dir, weight);

      if( useValuehistory(var, value, set) )
      {
         SCIP_HISTORY* history;

         SCIP_CALL( findValuehistoryEntry(var, value, blkmem, set, &history) );
         assert(history != NULL);

         SCIPhistoryIncVSIDS(history, dir, weight);
         SCIPsetDebugMsg(set, "variable (<%s> %s %g) + <%g> = <%g>\n", SCIPvarGetName(var), dir == SCIP_BRANCHDIR_UPWARDS ? ">=" : "<=",
            value, weight, SCIPhistoryGetVSIDS(history, dir));
      }

      return SCIP_OKAY;
   }
   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot update VSIDS of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      value = (value - var->data.aggregate.constant)/var->data.aggregate.scalar;

      if( var->data.aggregate.scalar > 0.0 )
      {
         SCIP_CALL( SCIPvarIncVSIDS(var->data.aggregate.var, blkmem, set, stat, dir, value, weight) );
      }
      else
      {
         assert(var->data.aggregate.scalar < 0.0);
         SCIP_CALL( SCIPvarIncVSIDS(var->data.aggregate.var, blkmem, set, stat, SCIPbranchdirOpposite(dir), value, weight) );
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot update VSIDS of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      value = 1.0 - value;

      SCIP_CALL( SCIPvarIncVSIDS(var->negatedvar, blkmem, set, stat, SCIPbranchdirOpposite(dir), value, weight) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** scales the VSIDS of the variable by the given scalar */
SCIP_RETCODE SCIPvarScaleVSIDS(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             scalar              /**< scalar to multiply the VSIDSs with */
   )
{
   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
      {
         SCIPerrorMessage("cannot update VSIDS of original untransformed variable\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarScaleVSIDS(var->data.original.transvar, scalar) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   {
      SCIPhistoryScaleVSIDS(var->history, scalar);
      SCIPhistoryScaleVSIDS(var->historycrun, scalar);
      SCIPvaluehistoryScaleVSIDS(var->valuehistory, scalar);

      return SCIP_OKAY;
   }
   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot update VSIDS of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      SCIP_CALL( SCIPvarScaleVSIDS(var->data.aggregate.var, scalar) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot update VSIDS of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      SCIP_CALL( SCIPvarScaleVSIDS(var->negatedvar, scalar) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** increases the number of active conflicts by one and the overall length of the variable by the given length */
SCIP_RETCODE SCIPvarIncNActiveConflicts(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_SET*             set,                /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   SCIP_Real             length              /**< length of the conflict */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   /* check if history statistics should be collected for a variable */
   if( !stat->collectvarhistory )
      return SCIP_OKAY;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
      {
         SCIPerrorMessage("cannot update conflict score of original untransformed variable\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarIncNActiveConflicts(var->data.original.transvar, blkmem, set, stat, dir, value, length) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   {
      SCIPhistoryIncNActiveConflicts(var->history, dir, length);
      SCIPhistoryIncNActiveConflicts(var->historycrun, dir, length);

      if( useValuehistory(var, value, set) )
      {
         SCIP_HISTORY* history;

         SCIP_CALL( findValuehistoryEntry(var, value, blkmem, set, &history) );
         assert(history != NULL);

         SCIPhistoryIncNActiveConflicts(history, dir, length);
      }

      return SCIP_OKAY;
   }
   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot update conflict score of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      value = (value - var->data.aggregate.constant)/var->data.aggregate.scalar;

      if( var->data.aggregate.scalar > 0.0 )
      {
         SCIP_CALL( SCIPvarIncNActiveConflicts(var->data.aggregate.var, blkmem, set, stat, dir, value, length) );
      }
      else
      {
         assert(var->data.aggregate.scalar < 0.0);
         SCIP_CALL( SCIPvarIncNActiveConflicts(var->data.aggregate.var, blkmem, set, stat, SCIPbranchdirOpposite(dir), value, length) );
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot update conflict score of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      value = 1.0 - value;

      SCIP_CALL( SCIPvarIncNActiveConflicts(var->negatedvar, blkmem, set, stat, SCIPbranchdirOpposite(dir), value, length) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** gets the number of active conflicts containing this variable in given direction */
SCIP_Longint SCIPvarGetNActiveConflicts(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0;
      else
         return SCIPvarGetNActiveConflicts(var->data.original.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNActiveConflicts(var->history, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNActiveConflicts(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetNActiveConflicts(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetNActiveConflicts(var->negatedvar, stat, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }
}

/** gets the number of active conflicts containing this variable in given direction
 *  in the current run
 */
SCIP_Longint SCIPvarGetNActiveConflictsCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0;
      else
         return SCIPvarGetNActiveConflictsCurrentRun(var->data.original.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNActiveConflicts(var->historycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNActiveConflictsCurrentRun(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetNActiveConflictsCurrentRun(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetNActiveConflictsCurrentRun(var->negatedvar, stat, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }
}

/** gets the average conflict length in given direction due to branching on the variable */
SCIP_Real SCIPvarGetAvgConflictlength(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetAvgConflictlength(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
         return SCIPhistoryGetAvgConflictlength(var->history, dir);
   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgConflictlength(var->data.aggregate.var, dir);
      else
         return SCIPvarGetAvgConflictlength(var->data.aggregate.var, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgConflictlength(var->negatedvar, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/**  gets the average conflict length in given direction due to branching on the variable
 *   in the current run
 */
SCIP_Real SCIPvarGetAvgConflictlengthCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetAvgConflictlengthCurrentRun(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return  SCIPhistoryGetAvgConflictlength(var->historycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgConflictlengthCurrentRun(var->data.aggregate.var, dir);
      else
         return SCIPvarGetAvgConflictlengthCurrentRun(var->data.aggregate.var, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgConflictlengthCurrentRun(var->negatedvar, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** increases the number of branchings counter of the variable */
SCIP_RETCODE SCIPvarIncNBranchings(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_SET*             set,                /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   int                   depth               /**< depth at which the bound change took place */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   /* check if history statistics should be collected for a variable */
   if( !stat->collectvarhistory )
      return SCIP_OKAY;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
      {
         SCIPerrorMessage("cannot update branching counter of original untransformed variable\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarIncNBranchings(var->data.original.transvar, blkmem, set, stat, dir, value, depth) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   {
      SCIPhistoryIncNBranchings(var->history, dir, depth);
      SCIPhistoryIncNBranchings(var->historycrun, dir, depth);
      SCIPhistoryIncNBranchings(stat->glbhistory, dir, depth);
      SCIPhistoryIncNBranchings(stat->glbhistorycrun, dir, depth);

      if( useValuehistory(var, value, set) )
      {
         SCIP_HISTORY* history;

         SCIP_CALL( findValuehistoryEntry(var, value, blkmem, set, &history) );
         assert(history != NULL);

         SCIPhistoryIncNBranchings(history, dir, depth);
      }

      return SCIP_OKAY;
   }
   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot update branching counter of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      value = (value - var->data.aggregate.constant)/var->data.aggregate.scalar;

      if( var->data.aggregate.scalar > 0.0 )
      {
         SCIP_CALL( SCIPvarIncNBranchings(var->data.aggregate.var, blkmem, set, stat, dir, value, depth) );
      }
      else
      {
         assert(var->data.aggregate.scalar < 0.0);
         SCIP_CALL( SCIPvarIncNBranchings(var->data.aggregate.var, blkmem, set, stat, SCIPbranchdirOpposite(dir), value, depth) );
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot update branching counter of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      value = 1.0 - value;

      SCIP_CALL( SCIPvarIncNBranchings(var->negatedvar, blkmem, set, stat, SCIPbranchdirOpposite(dir), value, depth) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** increases the inference sum of the variable by the given weight */
SCIP_RETCODE SCIPvarIncInferenceSum(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_SET*             set,                /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   SCIP_Real             weight              /**< weight of this update in inference score */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   /* check if history statistics should be collected for a variable */
   if( !stat->collectvarhistory )
      return SCIP_OKAY;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
      {
         SCIPerrorMessage("cannot update inference counter of original untransformed variable\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarIncInferenceSum(var->data.original.transvar, blkmem, set, stat, dir, value, weight) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   {
      SCIPhistoryIncInferenceSum(var->history, dir, weight);
      SCIPhistoryIncInferenceSum(var->historycrun, dir, weight);
      SCIPhistoryIncInferenceSum(stat->glbhistory, dir, weight);
      SCIPhistoryIncInferenceSum(stat->glbhistorycrun, dir, weight);

      if( useValuehistory(var, value, set) )
      {
         SCIP_HISTORY* history;

         SCIP_CALL( findValuehistoryEntry(var, value, blkmem, set, &history) );
         assert(history != NULL);

         SCIPhistoryIncInferenceSum(history, dir, weight);
      }

      return SCIP_OKAY;
   }
   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot update inference counter of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      value = (value - var->data.aggregate.constant)/var->data.aggregate.scalar;

      if( var->data.aggregate.scalar > 0.0 )
      {
         SCIP_CALL( SCIPvarIncInferenceSum(var->data.aggregate.var, blkmem, set, stat, dir, value, weight) );
      }
      else
      {
         assert(var->data.aggregate.scalar < 0.0);
         SCIP_CALL( SCIPvarIncInferenceSum(var->data.aggregate.var, blkmem, set, stat, SCIPbranchdirOpposite(dir), value, weight) );
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot update inference counter of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      value = 1.0 - value;

      SCIP_CALL( SCIPvarIncInferenceSum(var->negatedvar, blkmem, set, stat, SCIPbranchdirOpposite(dir), value, weight) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** increases the cutoff sum of the variable by the given weight */
SCIP_RETCODE SCIPvarIncCutoffSum(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_SET*             set,                /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   SCIP_Real             weight              /**< weight of this update in cutoff score */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   /* check if history statistics should be collected for a variable */
   if( !stat->collectvarhistory )
      return SCIP_OKAY;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
      {
         SCIPerrorMessage("cannot update cutoff sum of original untransformed variable\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarIncCutoffSum(var->data.original.transvar, blkmem, set, stat, dir, value, weight) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   {
      SCIPhistoryIncCutoffSum(var->history, dir, weight);
      SCIPhistoryIncCutoffSum(var->historycrun, dir, weight);
      SCIPhistoryIncCutoffSum(stat->glbhistory, dir, weight);
      SCIPhistoryIncCutoffSum(stat->glbhistorycrun, dir, weight);

      if( useValuehistory(var, value, set) )
      {
         SCIP_HISTORY* history;

         SCIP_CALL( findValuehistoryEntry(var, value, blkmem, set, &history) );
         assert(history != NULL);

         SCIPhistoryIncCutoffSum(history, dir, weight);
      }

      return SCIP_OKAY;
   }
   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot update cutoff sum of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      value = (value - var->data.aggregate.constant)/var->data.aggregate.scalar;

      if( var->data.aggregate.scalar > 0.0 )
      {
         SCIP_CALL( SCIPvarIncCutoffSum(var->data.aggregate.var, blkmem, set, stat, dir, value, weight) );
      }
      else
      {
         assert(var->data.aggregate.scalar < 0.0);
         SCIP_CALL( SCIPvarIncCutoffSum(var->data.aggregate.var, blkmem, set, stat, SCIPbranchdirOpposite(dir), value, weight) );
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot update cutoff sum of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      value = 1.0 - value;

      SCIP_CALL( SCIPvarIncCutoffSum(var->negatedvar, blkmem, set, stat, SCIPbranchdirOpposite(dir), value, weight) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** returns the number of times, a bound of the variable was changed in given direction due to branching */
SCIP_Longint SCIPvarGetNBranchings(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0;
      else
         return SCIPvarGetNBranchings(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNBranchings(var->history, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNBranchings(var->data.aggregate.var, dir);
      else
         return SCIPvarGetNBranchings(var->data.aggregate.var, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetNBranchings(var->negatedvar, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }
}

/** returns the number of times, a bound of the variable was changed in given direction due to branching 
 *  in the current run
 */
SCIP_Longint SCIPvarGetNBranchingsCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0;
      else
         return SCIPvarGetNBranchingsCurrentRun(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNBranchings(var->historycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNBranchingsCurrentRun(var->data.aggregate.var, dir);
      else
         return SCIPvarGetNBranchingsCurrentRun(var->data.aggregate.var, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetNBranchingsCurrentRun(var->negatedvar, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }
}

/** returns the average depth of bound changes in given direction due to branching on the variable */
SCIP_Real SCIPvarGetAvgBranchdepth(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetAvgBranchdepth(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return  SCIPhistoryGetAvgBranchdepth(var->history, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgBranchdepth(var->data.aggregate.var, dir);
      else
         return SCIPvarGetAvgBranchdepth(var->data.aggregate.var, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgBranchdepth(var->negatedvar, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the average depth of bound changes in given direction due to branching on the variable
 *  in the current run
 */
SCIP_Real SCIPvarGetAvgBranchdepthCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetAvgBranchdepthCurrentRun(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return  SCIPhistoryGetAvgBranchdepth(var->historycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgBranchdepthCurrentRun(var->data.aggregate.var, dir);
      else
         return SCIPvarGetAvgBranchdepthCurrentRun(var->data.aggregate.var,
            dir == SCIP_BRANCHDIR_DOWNWARDS ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS);

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgBranchdepthCurrentRun(var->negatedvar,
         dir == SCIP_BRANCHDIR_DOWNWARDS ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the variable's VSIDS score */
SCIP_Real SCIPvarGetVSIDS_rec(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      return SCIPvarGetVSIDS(var->data.original.transvar, stat, dir);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetVSIDS(var->data.original.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE); /* column case already handled in if condition above */
      return SCIPhistoryGetVSIDS(var->history, dir)/stat->vsidsweight;

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetVSIDS(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetVSIDS(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetVSIDS(var->negatedvar, stat, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the variable's VSIDS score only using conflicts of the current run */
SCIP_Real SCIPvarGetVSIDSCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   if( dir != SCIP_BRANCHDIR_DOWNWARDS && dir != SCIP_BRANCHDIR_UPWARDS )
   {
      SCIPerrorMessage("invalid branching direction %d when asking for VSIDS value\n", dir);
      return SCIP_INVALID;
   }

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetVSIDSCurrentRun(var->data.original.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetVSIDS(var->historycrun, dir)/stat->vsidsweight;

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetVSIDSCurrentRun(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetVSIDSCurrentRun(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetVSIDSCurrentRun(var->negatedvar, stat, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the number of inferences branching on this variable in given direction triggered */
SCIP_Real SCIPvarGetInferenceSum(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetInferenceSum(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetInferenceSum(var->history, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetInferenceSum(var->data.aggregate.var, dir);
      else
         return SCIPvarGetInferenceSum(var->data.aggregate.var, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetInferenceSum(var->negatedvar, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the number of inferences branching on this variable in given direction triggered
 *  in the current run
 */
SCIP_Real SCIPvarGetInferenceSumCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetInferenceSumCurrentRun(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetInferenceSum(var->historycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetInferenceSumCurrentRun(var->data.aggregate.var, dir);
      else
         return SCIPvarGetInferenceSumCurrentRun(var->data.aggregate.var, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetInferenceSumCurrentRun(var->negatedvar, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the average number of inferences found after branching on the variable in given direction */
SCIP_Real SCIPvarGetAvgInferences(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIPhistoryGetAvgInferences(stat->glbhistory, dir);
      else
         return SCIPvarGetAvgInferences(var->data.original.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      if( SCIPhistoryGetNBranchings(var->history, dir) > 0 )
         return SCIPhistoryGetAvgInferences(var->history, dir);
      else
      {
         int nimpls;
         int ncliques;

         nimpls = SCIPvarGetNImpls(var, dir == SCIP_BRANCHDIR_UPWARDS);
         ncliques = SCIPvarGetNCliques(var, dir == SCIP_BRANCHDIR_UPWARDS);
         return nimpls + ncliques > 0 ? (SCIP_Real)(nimpls + 2*ncliques) : SCIPhistoryGetAvgInferences(stat->glbhistory, dir);  /*lint !e790*/
      }

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgInferences(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetAvgInferences(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgInferences(var->negatedvar, stat, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the average number of inferences found after branching on the variable in given direction
 *  in the current run
 */
SCIP_Real SCIPvarGetAvgInferencesCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIPhistoryGetAvgInferences(stat->glbhistorycrun, dir);
      else
         return SCIPvarGetAvgInferencesCurrentRun(var->data.original.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      if( SCIPhistoryGetNBranchings(var->historycrun, dir) > 0 )
         return SCIPhistoryGetAvgInferences(var->historycrun, dir);
      else
      {
         int nimpls;
         int ncliques;

         nimpls = SCIPvarGetNImpls(var, dir == SCIP_BRANCHDIR_UPWARDS);
         ncliques = SCIPvarGetNCliques(var, dir == SCIP_BRANCHDIR_UPWARDS);
         return nimpls + ncliques > 0 ? (SCIP_Real)(nimpls + 2*ncliques) : SCIPhistoryGetAvgInferences(stat->glbhistorycrun, dir);  /*lint !e790*/
      }

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgInferencesCurrentRun(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetAvgInferencesCurrentRun(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgInferencesCurrentRun(var->negatedvar, stat, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the number of cutoffs branching on this variable in given direction produced */
SCIP_Real SCIPvarGetCutoffSum(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0;
      else
         return SCIPvarGetCutoffSum(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetCutoffSum(var->history, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetCutoffSum(var->data.aggregate.var, dir);
      else
         return SCIPvarGetCutoffSum(var->data.aggregate.var, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetCutoffSum(var->negatedvar, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }
}

/** returns the number of cutoffs branching on this variable in given direction produced in the current run */
SCIP_Real SCIPvarGetCutoffSumCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return 0;
      else
         return SCIPvarGetCutoffSumCurrentRun(var->data.original.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetCutoffSum(var->historycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetCutoffSumCurrentRun(var->data.aggregate.var, dir);
      else
         return SCIPvarGetCutoffSumCurrentRun(var->data.aggregate.var, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetCutoffSumCurrentRun(var->negatedvar, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }
}

/** returns the average number of cutoffs found after branching on the variable in given direction */
SCIP_Real SCIPvarGetAvgCutoffs(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIPhistoryGetAvgCutoffs(stat->glbhistory, dir);
      else
         return SCIPvarGetAvgCutoffs(var->data.original.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNBranchings(var->history, dir) > 0
         ? SCIPhistoryGetAvgCutoffs(var->history, dir)
         : SCIPhistoryGetAvgCutoffs(stat->glbhistory, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgCutoffs(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetAvgCutoffs(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgCutoffs(var->negatedvar, stat, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the average number of cutoffs found after branching on the variable in given direction in the current run */
SCIP_Real SCIPvarGetAvgCutoffsCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIPhistoryGetAvgCutoffs(stat->glbhistorycrun, dir);
      else
         return SCIPvarGetAvgCutoffsCurrentRun(var->data.original.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNBranchings(var->historycrun, dir) > 0
         ? SCIPhistoryGetAvgCutoffs(var->historycrun, dir)
         : SCIPhistoryGetAvgCutoffs(stat->glbhistorycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgCutoffsCurrentRun(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetAvgCutoffsCurrentRun(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));

   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgCutoffsCurrentRun(var->negatedvar, stat, SCIPbranchdirOpposite(dir));

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}




/*
 * information methods for bound changes
 */

/** creates an artificial bound change information object with depth = INT_MAX and pos = -1 */
SCIP_RETCODE SCIPbdchginfoCreate(
   SCIP_BDCHGINFO**      bdchginfo,          /**< pointer to store bound change information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< active variable that changed the bounds */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound for var: lower or upper bound */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   assert(bdchginfo != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, bdchginfo) );
   (*bdchginfo)->oldbound = oldbound;
   (*bdchginfo)->newbound = newbound;
   (*bdchginfo)->var = var;
   (*bdchginfo)->inferencedata.var = var;
   (*bdchginfo)->inferencedata.reason.prop = NULL;
   (*bdchginfo)->inferencedata.info = 0;
   (*bdchginfo)->bdchgidx.depth = INT_MAX;
   (*bdchginfo)->bdchgidx.pos = -1;
   (*bdchginfo)->pos = 0;
   (*bdchginfo)->boundchgtype = SCIP_BOUNDCHGTYPE_BRANCHING; /*lint !e641*/
   (*bdchginfo)->boundtype = boundtype; /*lint !e641*/
   (*bdchginfo)->inferboundtype = boundtype; /*lint !e641*/
   (*bdchginfo)->redundant = FALSE;

   return SCIP_OKAY;
}

/** frees a bound change information object */
void SCIPbdchginfoFree(
   SCIP_BDCHGINFO**      bdchginfo,          /**< pointer to store bound change information */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(bdchginfo != NULL);

   BMSfreeBlockMemory(blkmem, bdchginfo);
}

/** returns the bound change information for the last lower bound change on given active problem variable before or
 *  after the bound change with the given index was applied;
 *  returns NULL, if no change to the lower bound was applied up to this point of time
 */
SCIP_BDCHGINFO* SCIPvarGetLbchgInfo(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   )
{
   int i;

   assert(var != NULL);
   assert(SCIPvarIsActive(var));

   /* search the correct bound change information for the given bound change index */
   if( after )
   {
      for( i = var->nlbchginfos-1; i >= 0; --i )
      {
         assert(var->lbchginfos[i].var == var);
         assert((SCIP_BOUNDTYPE)var->lbchginfos[i].boundtype == SCIP_BOUNDTYPE_LOWER);
         assert(var->lbchginfos[i].pos == i);

         /* if we reached the (due to global bounds) redundant bound changes, return NULL */
         if( var->lbchginfos[i].redundant )
            return NULL;
         assert(var->lbchginfos[i].oldbound < var->lbchginfos[i].newbound);

         /* if we reached the bound change index, return the current bound change info */
         if( !SCIPbdchgidxIsEarlier(bdchgidx, &var->lbchginfos[i].bdchgidx) )
            return &var->lbchginfos[i];
      }
   }
   else
   {
      for( i = var->nlbchginfos-1; i >= 0; --i )
      {
         assert(var->lbchginfos[i].var == var);
         assert((SCIP_BOUNDTYPE)var->lbchginfos[i].boundtype == SCIP_BOUNDTYPE_LOWER);
         assert(var->lbchginfos[i].pos == i);

         /* if we reached the (due to global bounds) redundant bound changes, return NULL */
         if( var->lbchginfos[i].redundant )
            return NULL;
         assert(var->lbchginfos[i].oldbound < var->lbchginfos[i].newbound);

         /* if we reached the bound change index, return the current bound change info */
         if( SCIPbdchgidxIsEarlier(&var->lbchginfos[i].bdchgidx, bdchgidx) )
            return &var->lbchginfos[i];
      }
   }

   return NULL;
}

/** returns the bound change information for the last upper bound change on given active problem variable before or
 *  after the bound change with the given index was applied;
 *  returns NULL, if no change to the upper bound was applied up to this point of time
 */
SCIP_BDCHGINFO* SCIPvarGetUbchgInfo(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   )
{
   int i;

   assert(var != NULL);
   assert(SCIPvarIsActive(var));

   /* search the correct bound change information for the given bound change index */
   if( after )
   {
      for( i = var->nubchginfos-1; i >= 0; --i )
      {
         assert(var->ubchginfos[i].var == var);
         assert((SCIP_BOUNDTYPE)var->ubchginfos[i].boundtype == SCIP_BOUNDTYPE_UPPER);
         assert(var->ubchginfos[i].pos == i);

         /* if we reached the (due to global bounds) redundant bound changes, return NULL */
         if( var->ubchginfos[i].redundant )
            return NULL;
         assert(var->ubchginfos[i].oldbound > var->ubchginfos[i].newbound);

         /* if we reached the bound change index, return the current bound change info */
         if( !SCIPbdchgidxIsEarlier(bdchgidx, &var->ubchginfos[i].bdchgidx) )
            return &var->ubchginfos[i];
      }
   }
   else
   {
      for( i = var->nubchginfos-1; i >= 0; --i )
      {
         assert(var->ubchginfos[i].var == var);
         assert((SCIP_BOUNDTYPE)var->ubchginfos[i].boundtype == SCIP_BOUNDTYPE_UPPER);
         assert(var->ubchginfos[i].pos == i);

         /* if we reached the (due to global bounds) redundant bound changes, return NULL */
         if( var->ubchginfos[i].redundant )
            return NULL;
         assert(var->ubchginfos[i].oldbound > var->ubchginfos[i].newbound);

         /* if we reached the bound change index, return the current bound change info */
         if( SCIPbdchgidxIsEarlier(&var->ubchginfos[i].bdchgidx, bdchgidx) )
            return &var->ubchginfos[i];
      }
   }

   return NULL;
}

/** returns the bound change information for the last lower or upper bound change on given active problem variable
 *  before or after the bound change with the given index was applied;
 *  returns NULL, if no change to the lower/upper bound was applied up to this point of time
 */
SCIP_BDCHGINFO* SCIPvarGetBdchgInfo(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   )
{
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
      return SCIPvarGetLbchgInfo(var, bdchgidx, after);
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      return SCIPvarGetUbchgInfo(var, bdchgidx, after);
   }
}

/** returns lower bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 *
 *  @deprecated Please use SCIPgetVarLbAtIndex()
 */
SCIP_Real SCIPvarGetLbAtIndex(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   )
{
   SCIP_VARSTATUS varstatus;
   assert(var != NULL);

   varstatus = SCIPvarGetStatus(var);

   if( varstatus == SCIP_VARSTATUS_COLUMN || varstatus == SCIP_VARSTATUS_LOOSE )
   {
      if( bdchgidx == NULL )
         return SCIPvarGetLbLocal(var);
      else
      {
         SCIP_BDCHGINFO* bdchginfo;

         bdchginfo = SCIPvarGetLbchgInfo(var, bdchgidx, after);
         if( bdchginfo != NULL )
            return SCIPbdchginfoGetNewbound(bdchginfo);
         else
            return var->glbdom.lb;
      }
   }

   /* get bounds of attached variables */
   switch( varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      return SCIPvarGetLbAtIndex(var->data.original.transvar, bdchgidx, after);

   case SCIP_VARSTATUS_FIXED:
      return var->glbdom.lb;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      /* a correct implementation would need to check the value of var->data.aggregate.var for infinity and return the
       * corresponding infinity value instead of performing an arithmetical transformation (compare method
       * SCIPvarGetLbLP()); however, we do not want to introduce a SCIP or SCIP_SET pointer to this method, since it is
       * (or is called by) a public interface method; instead, we only assert that values are finite
       * w.r.t. SCIP_DEFAULT_INFINITY, which seems to be true in our regression tests; note that this may yield false
       * positives and negatives if the parameter <numerics/infinity> is modified by the user
       */
      if( var->data.aggregate.scalar > 0.0 )
      {
         /* a > 0 -> get lower bound of y */
         assert(SCIPvarGetLbAtIndex(var->data.aggregate.var, bdchgidx, after) > -SCIP_DEFAULT_INFINITY);
         assert(SCIPvarGetLbAtIndex(var->data.aggregate.var, bdchgidx, after) < +SCIP_DEFAULT_INFINITY);
         return var->data.aggregate.scalar * SCIPvarGetLbAtIndex(var->data.aggregate.var, bdchgidx, after)
            + var->data.aggregate.constant;
      }
      else if( var->data.aggregate.scalar < 0.0 )
      {
         /* a < 0 -> get upper bound of y */
         assert(SCIPvarGetUbAtIndex(var->data.aggregate.var, bdchgidx, after) > -SCIP_DEFAULT_INFINITY);
         assert(SCIPvarGetUbAtIndex(var->data.aggregate.var, bdchgidx, after) < +SCIP_DEFAULT_INFINITY);
         return var->data.aggregate.scalar * SCIPvarGetUbAtIndex(var->data.aggregate.var, bdchgidx, after)
            + var->data.aggregate.constant;
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         SCIPABORT();
         return SCIP_INVALID; /*lint !e527*/
      }

   case SCIP_VARSTATUS_MULTAGGR:
      /* handle multi-aggregated variables depending on one variable only (possibly caused by SCIPvarFlattenAggregationGraph()) */
      if ( var->data.multaggr.nvars == 1 )
      {
         assert(var->data.multaggr.vars != NULL);
         assert(var->data.multaggr.scalars != NULL);
         assert(var->data.multaggr.vars[0] != NULL);

         if( var->data.multaggr.scalars[0] > 0.0 )
         {
            /* a > 0 -> get lower bound of y */
            assert(SCIPvarGetLbAtIndex(var->data.multaggr.vars[0], bdchgidx, after) > -SCIP_DEFAULT_INFINITY);
            assert(SCIPvarGetLbAtIndex(var->data.multaggr.vars[0], bdchgidx, after) < +SCIP_DEFAULT_INFINITY);
            return var->data.multaggr.scalars[0] * SCIPvarGetLbAtIndex(var->data.multaggr.vars[0], bdchgidx, after)
               + var->data.multaggr.constant;
         }
         else if( var->data.multaggr.scalars[0] < 0.0 )
         {
            /* a < 0 -> get upper bound of y */
            assert(SCIPvarGetUbAtIndex(var->data.multaggr.vars[0], bdchgidx, after) > -SCIP_DEFAULT_INFINITY);
            assert(SCIPvarGetUbAtIndex(var->data.multaggr.vars[0], bdchgidx, after) < +SCIP_DEFAULT_INFINITY);
            return var->data.multaggr.scalars[0] * SCIPvarGetUbAtIndex(var->data.multaggr.vars[0], bdchgidx, after)
               + var->data.multaggr.constant;
         }
         else
         {
            SCIPerrorMessage("scalar is zero in multi-aggregation\n");
            SCIPABORT();
            return SCIP_INVALID; /*lint !e527*/
         }
      }
      SCIPerrorMessage("cannot get the bounds of a multi-aggregated variable.\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetUbAtIndex(var->negatedvar, bdchgidx, after);

   case SCIP_VARSTATUS_COLUMN: /* for lint */
   case SCIP_VARSTATUS_LOOSE: /* for lint */
   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** returns upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 *
 *  @deprecated Please use SCIPgetVarUbAtIndex()
 */
SCIP_Real SCIPvarGetUbAtIndex(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   )
{
   SCIP_VARSTATUS varstatus;
   assert(var != NULL);

   varstatus = SCIPvarGetStatus(var);

   if( varstatus == SCIP_VARSTATUS_COLUMN || varstatus == SCIP_VARSTATUS_LOOSE )
   {
      if( bdchgidx == NULL )
         return SCIPvarGetUbLocal(var);
      else
      {
         SCIP_BDCHGINFO* bdchginfo;

         bdchginfo = SCIPvarGetUbchgInfo(var, bdchgidx, after);
         if( bdchginfo != NULL )
            return SCIPbdchginfoGetNewbound(bdchginfo);
         else
            return var->glbdom.ub;
      }
   }

   /* get bounds of attached variables */
   switch( varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      return SCIPvarGetUbAtIndex(var->data.original.transvar, bdchgidx, after);

   case SCIP_VARSTATUS_FIXED:
      return var->glbdom.ub;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      /* a correct implementation would need to check the value of var->data.aggregate.var for infinity and return the
       * corresponding infinity value instead of performing an arithmetical transformation (compare method
       * SCIPvarGetLbLP()); however, we do not want to introduce a SCIP or SCIP_SET pointer to this method, since it is
       * (or is called by) a public interface method; instead, we only assert that values are finite
       * w.r.t. SCIP_DEFAULT_INFINITY, which seems to be true in our regression tests; note that this may yield false
       * positives and negatives if the parameter <numerics/infinity> is modified by the user
       */
      if( var->data.aggregate.scalar > 0.0 )
      {
         /* a > 0 -> get lower bound of y */
         assert(SCIPvarGetUbAtIndex(var->data.aggregate.var, bdchgidx, after) > -SCIP_DEFAULT_INFINITY);
         assert(SCIPvarGetUbAtIndex(var->data.aggregate.var, bdchgidx, after) < +SCIP_DEFAULT_INFINITY);
         return var->data.aggregate.scalar * SCIPvarGetUbAtIndex(var->data.aggregate.var, bdchgidx, after)
            + var->data.aggregate.constant;
      }
      else if( var->data.aggregate.scalar < 0.0 )
      {
         /* a < 0 -> get upper bound of y */
         assert(SCIPvarGetLbAtIndex(var->data.aggregate.var, bdchgidx, after) > -SCIP_DEFAULT_INFINITY);
         assert(SCIPvarGetLbAtIndex(var->data.aggregate.var, bdchgidx, after) < +SCIP_DEFAULT_INFINITY);
         return var->data.aggregate.scalar * SCIPvarGetLbAtIndex(var->data.aggregate.var, bdchgidx, after)
            + var->data.aggregate.constant;
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         SCIPABORT();
         return SCIP_INVALID; /*lint !e527*/
      }

   case SCIP_VARSTATUS_MULTAGGR:
      /* handle multi-aggregated variables depending on one variable only (possibly caused by SCIPvarFlattenAggregationGraph()) */
      if ( var->data.multaggr.nvars == 1 )
      {
         assert(var->data.multaggr.vars != NULL);
         assert(var->data.multaggr.scalars != NULL);
         assert(var->data.multaggr.vars[0] != NULL);

         if( var->data.multaggr.scalars[0] > 0.0 )
         {
            /* a > 0 -> get lower bound of y */
            assert(SCIPvarGetUbAtIndex(var->data.multaggr.vars[0], bdchgidx, after) > -SCIP_DEFAULT_INFINITY);
            assert(SCIPvarGetUbAtIndex(var->data.multaggr.vars[0], bdchgidx, after) < +SCIP_DEFAULT_INFINITY);
            return var->data.multaggr.scalars[0] * SCIPvarGetUbAtIndex(var->data.multaggr.vars[0], bdchgidx, after)
               + var->data.multaggr.constant;
         }
         else if( var->data.multaggr.scalars[0] < 0.0 )
         {
            /* a < 0 -> get upper bound of y */
            assert(SCIPvarGetLbAtIndex(var->data.multaggr.vars[0], bdchgidx, after) > -SCIP_DEFAULT_INFINITY);
            assert(SCIPvarGetLbAtIndex(var->data.multaggr.vars[0], bdchgidx, after) < +SCIP_DEFAULT_INFINITY);
            return var->data.multaggr.scalars[0] * SCIPvarGetLbAtIndex(var->data.multaggr.vars[0], bdchgidx, after)
               + var->data.multaggr.constant;
         }
         else
         {
            SCIPerrorMessage("scalar is zero in multi-aggregation\n");
            SCIPABORT();
            return SCIP_INVALID; /*lint !e527*/
         }
      }
      SCIPerrorMessage("cannot get the bounds of a multiple aggregated variable.\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetLbAtIndex(var->negatedvar, bdchgidx, after);

   case SCIP_VARSTATUS_COLUMN: /* for lint */
   case SCIP_VARSTATUS_LOOSE: /* for lint */
   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** returns lower or upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 *
 *  @deprecated Please use SCIPgetVarBdAtIndex()
 */
SCIP_Real SCIPvarGetBdAtIndex(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   )
{
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
      return SCIPvarGetLbAtIndex(var, bdchgidx, after);
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      return SCIPvarGetUbAtIndex(var, bdchgidx, after);
   }
}

/** returns whether the binary variable was fixed at the time given by the bound change index
 *
 *  @deprecated Please use SCIPgetVarWasFixedAtIndex()
 */
SCIP_Bool SCIPvarWasFixedAtIndex(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsBinary(var));

   /* check the current bounds first in order to decide at which bound change information we have to look
    * (which is expensive because we have to follow the aggregation tree to the active variable)
    */
   return ((SCIPvarGetLbLocal(var) > 0.5 && SCIPvarGetLbAtIndex(var, bdchgidx, after) > 0.5)
      || (SCIPvarGetUbLocal(var) < 0.5 && SCIPvarGetUbAtIndex(var, bdchgidx, after) < 0.5));
}

/** bound change index representing the initial time before any bound changes took place */
static SCIP_BDCHGIDX initbdchgidx = {-2, 0};

/** bound change index representing the presolving stage */
static SCIP_BDCHGIDX presolvebdchgidx = {-1, 0};

/** returns the last bound change index, at which the bounds of the given variable were tightened */
SCIP_BDCHGIDX* SCIPvarGetLastBdchgIndex(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_BDCHGIDX* lbchgidx;
   SCIP_BDCHGIDX* ubchgidx;

   assert(var != NULL);

   var = SCIPvarGetProbvar(var);

   /* check, if variable is original without transformed variable */
   if( var == NULL )
      return &initbdchgidx;

   /* check, if variable was fixed in presolving */
   if( !SCIPvarIsActive(var) )
      return &presolvebdchgidx;

   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   /* get depths of last bound change information for the lower and upper bound */
   lbchgidx = (var->nlbchginfos > 0 && !var->lbchginfos[var->nlbchginfos-1].redundant
      ? &var->lbchginfos[var->nlbchginfos-1].bdchgidx : &initbdchgidx);
   ubchgidx = (var->nubchginfos > 0 && !var->ubchginfos[var->nubchginfos-1].redundant
      ? &var->ubchginfos[var->nubchginfos-1].bdchgidx : &initbdchgidx);

   if( SCIPbdchgidxIsEarlierNonNull(lbchgidx, ubchgidx) )
      return ubchgidx;
   else
      return lbchgidx;
}

/** returns the last depth level, at which the bounds of the given variable were tightened;
 *  returns -2, if the variable's bounds are still the global bounds
 *  returns -1, if the variable was fixed in presolving
 */
int SCIPvarGetLastBdchgDepth(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_BDCHGIDX* bdchgidx;

   bdchgidx = SCIPvarGetLastBdchgIndex(var);
   assert(bdchgidx != NULL);

   return bdchgidx->depth;
}

/** returns at which depth in the tree a bound change was applied to the variable that conflicts with the
 *  given bound; returns -1 if the bound does not conflict with the current local bounds of the variable
 */
int SCIPvarGetConflictingBdchgDepth(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BOUNDTYPE        boundtype,          /**< bound type of the conflicting bound */
   SCIP_Real             bound               /**< conflicting bound */
   )
{
   int i;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);

   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      /* check if the bound is in conflict with the current local bounds */
      if( SCIPsetIsLE(set, bound, var->locdom.ub) )
         return -1;

      /* check if the bound is in conflict with the global bound */
      if( SCIPsetIsGT(set, bound, var->glbdom.ub) )
         return 0;

      /* local bounds are in conflict with the given bound -> there must be at least one conflicting change! */
      assert(var->nubchginfos > 0);
      assert(SCIPsetIsGT(set, bound, var->ubchginfos[var->nubchginfos-1].newbound));

      /* search for the first conflicting bound change */
      for( i = var->nubchginfos-1; i > 0 && SCIPsetIsGT(set, bound, var->ubchginfos[i-1].newbound); --i )
      {
         assert(var->ubchginfos[i].var == var); /* perform sanity check on the search for the first conflicting bound */
         assert((SCIP_BOUNDTYPE)var->ubchginfos[i].boundtype == SCIP_BOUNDTYPE_UPPER);
      }
      assert(SCIPsetIsGT(set, bound, var->ubchginfos[i].newbound));             /* bound change i is conflicting */
      assert(i == 0 || SCIPsetIsLE(set, bound, var->ubchginfos[i-1].newbound)); /* bound change i-1 is not conflicting */

      /* return the depth at which the first conflicting bound change took place */
      return var->ubchginfos[i].bdchgidx.depth;
   }
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);

      /* check if the bound is in conflict with the current local bounds */
      if( SCIPsetIsGE(set, bound, var->locdom.lb) )
         return -1;

      /* check if the bound is in conflict with the global bound */
      if( SCIPsetIsLT(set, bound, var->glbdom.lb) )
         return 0;

      /* local bounds are in conflict with the given bound -> there must be at least one conflicting change! */
      assert(var->nlbchginfos > 0);
      assert(SCIPsetIsLT(set, bound, var->lbchginfos[var->nlbchginfos-1].newbound));

      /* search for the first conflicting bound change */
      for( i = var->nlbchginfos-1; i > 0 && SCIPsetIsLT(set, bound, var->lbchginfos[i-1].newbound); --i )
      {
         assert(var->lbchginfos[i].var == var); /* perform sanity check on the search for the first conflicting bound */
         assert((SCIP_BOUNDTYPE)var->lbchginfos[i].boundtype == SCIP_BOUNDTYPE_LOWER);
      }
      assert(SCIPsetIsLT(set, bound, var->lbchginfos[i].newbound));             /* bound change i is conflicting */
      assert(i == 0 || SCIPsetIsGE(set, bound, var->lbchginfos[i-1].newbound)); /* bound change i-1 is not conflicting */

      /* return the depth at which the first conflicting bound change took place */
      return var->lbchginfos[i].bdchgidx.depth;
   }
}

/** returns whether the first binary variable was fixed earlier than the second one;
 *  returns FALSE, if the first variable is not fixed, and returns TRUE, if the first variable is fixed, but the
 *  second one is not fixed
 */
SCIP_Bool SCIPvarWasFixedEarlier(
   SCIP_VAR*             var1,               /**< first binary variable */
   SCIP_VAR*             var2                /**< second binary variable */
   )
{
   SCIP_BDCHGIDX* bdchgidx1;
   SCIP_BDCHGIDX* bdchgidx2;

   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(SCIPvarIsBinary(var1));
   assert(SCIPvarIsBinary(var2));

   var1 = SCIPvarGetProbvar(var1);
   var2 = SCIPvarGetProbvar(var2);
   assert(var1 != NULL);
   assert(var2 != NULL);

   /* check, if variables are globally fixed */
   if( !SCIPvarIsActive(var2) || var2->glbdom.lb > 0.5 || var2->glbdom.ub < 0.5 )
      return FALSE;
   if( !SCIPvarIsActive(var1) || var1->glbdom.lb > 0.5 || var1->glbdom.ub < 0.5 )
      return TRUE;

   assert(SCIPvarGetStatus(var1) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var1) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetStatus(var2) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var2) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarIsBinary(var1));
   assert(SCIPvarIsBinary(var2));
   assert(var1->nlbchginfos + var1->nubchginfos <= 1);
   assert(var2->nlbchginfos + var2->nubchginfos <= 1);
   assert(var1->nlbchginfos == 0 || !var1->lbchginfos[0].redundant); /* otherwise, var would be globally fixed */
   assert(var1->nubchginfos == 0 || !var1->ubchginfos[0].redundant); /* otherwise, var would be globally fixed */
   assert(var2->nlbchginfos == 0 || !var2->lbchginfos[0].redundant); /* otherwise, var would be globally fixed */
   assert(var2->nubchginfos == 0 || !var2->ubchginfos[0].redundant); /* otherwise, var would be globally fixed */

   if( var1->nlbchginfos == 1 )
      bdchgidx1 = &var1->lbchginfos[0].bdchgidx;
   else if( var1->nubchginfos == 1 )
      bdchgidx1 = &var1->ubchginfos[0].bdchgidx;
   else
      bdchgidx1 = NULL;

   if( var2->nlbchginfos == 1 )
      bdchgidx2 = &var2->lbchginfos[0].bdchgidx;
   else if( var2->nubchginfos == 1 )
      bdchgidx2 = &var2->ubchginfos[0].bdchgidx;
   else
      bdchgidx2 = NULL;

   return SCIPbdchgidxIsEarlier(bdchgidx1, bdchgidx2);
}



/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given variable */
SCIP_DECL_HASHGETKEY(SCIPhashGetKeyVar)
{  /*lint --e{715}*/
   SCIP_VAR* var = (SCIP_VAR*)elem;

   assert(var != NULL);
   return var->name;
}




/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPboundchgGetNewbound
#undef SCIPboundchgGetVar
#undef SCIPboundchgGetBoundchgtype
#undef SCIPboundchgGetBoundtype
#undef SCIPboundchgIsRedundant
#undef SCIPdomchgGetNBoundchgs
#undef SCIPdomchgGetBoundchg
#undef SCIPholelistGetLeft
#undef SCIPholelistGetRight
#undef SCIPholelistGetNext
#undef SCIPvarGetName
#undef SCIPvarGetNUses
#undef SCIPvarGetData
#undef SCIPvarSetData
#undef SCIPvarSetDelorigData
#undef SCIPvarSetTransData
#undef SCIPvarSetDeltransData
#undef SCIPvarGetStatus
#undef SCIPvarIsOriginal
#undef SCIPvarIsTransformed
#undef SCIPvarIsNegated
#undef SCIPvarGetType
#undef SCIPvarIsBinary
#undef SCIPvarIsIntegral
#undef SCIPvarIsInitial
#undef SCIPvarIsRemovable
#undef SCIPvarIsDeleted
#undef SCIPvarIsDeletable
#undef SCIPvarMarkDeletable
#undef SCIPvarMarkNotDeletable
#undef SCIPvarIsActive
#undef SCIPvarGetIndex
#undef SCIPvarGetProbindex
#undef SCIPvarGetTransVar
#undef SCIPvarGetCol
#undef SCIPvarIsInLP
#undef SCIPvarGetAggrVar
#undef SCIPvarGetAggrScalar
#undef SCIPvarGetAggrConstant
#undef SCIPvarGetMultaggrNVars
#undef SCIPvarGetMultaggrVars
#undef SCIPvarGetMultaggrScalars
#undef SCIPvarGetMultaggrConstant
#undef SCIPvarGetNegatedVar
#undef SCIPvarGetNegationVar
#undef SCIPvarGetNegationConstant
#undef SCIPvarGetObj
#undef SCIPvarGetLbOriginal
#undef SCIPvarGetUbOriginal
#undef SCIPvarGetHolelistOriginal
#undef SCIPvarGetLbGlobal
#undef SCIPvarGetUbGlobal
#undef SCIPvarGetHolelistGlobal
#undef SCIPvarGetBestBoundGlobal
#undef SCIPvarGetWorstBoundGlobal
#undef SCIPvarGetLbLocal
#undef SCIPvarGetUbLocal
#undef SCIPvarGetHolelistLocal
#undef SCIPvarGetBestBoundLocal
#undef SCIPvarGetWorstBoundLocal
#undef SCIPvarGetBestBoundType
#undef SCIPvarGetWorstBoundType
#undef SCIPvarGetLbLazy
#undef SCIPvarGetUbLazy
#undef SCIPvarGetBranchFactor
#undef SCIPvarGetBranchPriority
#undef SCIPvarGetBranchDirection
#undef SCIPvarGetNVlbs
#undef SCIPvarGetVlbVars
#undef SCIPvarGetVlbCoefs
#undef SCIPvarGetVlbConstants
#undef SCIPvarGetNVubs
#undef SCIPvarGetVubVars
#undef SCIPvarGetVubCoefs
#undef SCIPvarGetVubConstants
#undef SCIPvarGetNImpls
#undef SCIPvarGetImplVars
#undef SCIPvarGetImplTypes
#undef SCIPvarGetImplBounds
#undef SCIPvarGetImplIds
#undef SCIPvarGetNCliques
#undef SCIPvarGetCliques
#undef SCIPvarGetLPSol
#undef SCIPvarGetNLPSol
#undef SCIPvarGetBdchgInfoLb
#undef SCIPvarGetNBdchgInfosLb
#undef SCIPvarGetBdchgInfoUb
#undef SCIPvarGetNBdchgInfosUb
#undef SCIPvarGetValuehistory
#undef SCIPvarGetPseudoSol
#undef SCIPvarCatchEvent
#undef SCIPvarDropEvent
#undef SCIPvarGetVSIDS
#undef SCIPvarGetCliqueComponentIdx
#undef SCIPbdchgidxGetPos
#undef SCIPbdchgidxIsEarlierNonNull
#undef SCIPbdchgidxIsEarlier
#undef SCIPbdchginfoGetOldbound
#undef SCIPbdchginfoGetNewbound
#undef SCIPbdchginfoGetVar
#undef SCIPbdchginfoGetChgtype
#undef SCIPbdchginfoGetBoundtype
#undef SCIPbdchginfoGetDepth
#undef SCIPbdchginfoGetPos
#undef SCIPbdchginfoGetIdx
#undef SCIPbdchginfoGetInferVar
#undef SCIPbdchginfoGetInferCons
#undef SCIPbdchginfoGetInferProp
#undef SCIPbdchginfoGetInferInfo
#undef SCIPbdchginfoGetInferBoundtype
#undef SCIPbdchginfoIsRedundant
#undef SCIPbdchginfoHasInferenceReason
#undef SCIPbdchginfoIsTighter


/** returns the new value of the bound in the bound change data */
SCIP_Real SCIPboundchgGetNewbound(
   SCIP_BOUNDCHG*        boundchg            /**< bound change data */
   )
{
   assert(boundchg != NULL);

   return boundchg->newbound;
}

/** returns the variable of the bound change in the bound change data */
SCIP_VAR* SCIPboundchgGetVar(
   SCIP_BOUNDCHG*        boundchg            /**< bound change data */
   )
{
   assert(boundchg != NULL);

   return boundchg->var;
}

/** returns the bound change type of the bound change in the bound change data */
SCIP_BOUNDCHGTYPE SCIPboundchgGetBoundchgtype(
   SCIP_BOUNDCHG*        boundchg            /**< bound change data */
   )
{
   assert(boundchg != NULL);

   return (SCIP_BOUNDCHGTYPE)(boundchg->boundchgtype);
}

/** returns the bound type of the bound change in the bound change data */
SCIP_BOUNDTYPE SCIPboundchgGetBoundtype(
   SCIP_BOUNDCHG*        boundchg            /**< bound change data */
   )
{
   assert(boundchg != NULL);

   return (SCIP_BOUNDTYPE)(boundchg->boundtype);
}

/** returns whether the bound change is redundant due to a more global bound that is at least as strong */
SCIP_Bool SCIPboundchgIsRedundant(
   SCIP_BOUNDCHG*        boundchg            /**< bound change data */
   )
{
   assert(boundchg != NULL);

   return boundchg->redundant;
}

/** returns the number of bound changes in the domain change data */
int SCIPdomchgGetNBoundchgs(
   SCIP_DOMCHG*          domchg              /**< domain change data */
   )
{
   return domchg != NULL ? domchg->domchgbound.nboundchgs : 0;
}

/** returns a particular bound change in the domain change data */
SCIP_BOUNDCHG* SCIPdomchgGetBoundchg(
   SCIP_DOMCHG*          domchg,             /**< domain change data */
   int                   pos                 /**< position of the bound change in the domain change data */
   )
{
   assert(domchg != NULL);
   assert(0 <= pos && pos < (int)domchg->domchgbound.nboundchgs);

   return &domchg->domchgbound.boundchgs[pos];
}

/** returns left bound of open interval in hole */
SCIP_Real SCIPholelistGetLeft(
   SCIP_HOLELIST*        holelist            /**< hole list pointer to hole of interest */
   )
{
   assert(holelist != NULL);

   return holelist->hole.left;
}

/** returns right bound of open interval in hole */
SCIP_Real SCIPholelistGetRight(
   SCIP_HOLELIST*        holelist            /**< hole list pointer to hole of interest */
   )
{
   assert(holelist != NULL);

   return holelist->hole.right;
}

/** returns next hole in list */
SCIP_HOLELIST* SCIPholelistGetNext(
   SCIP_HOLELIST*        holelist            /**< hole list pointer to hole of interest */
   )
{
   assert(holelist != NULL);

   return holelist->next;
}

/** returns the name of the variable
 *
 *  @note to change the name of a variable, use SCIPchgVarName() from scip.h
 */
const char* SCIPvarGetName(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->name;
}

/** gets number of times, the variable is currently captured */
int SCIPvarGetNUses(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->nuses;
}

/** returns the user data of the variable */
SCIP_VARDATA* SCIPvarGetData(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->vardata;
}

/** sets the user data for the variable */
void SCIPvarSetData(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_VARDATA*         vardata             /**< user variable data */
   )
{
   assert(var != NULL);

   var->vardata = vardata;
}

/** sets method to free user data for the original variable */
void SCIPvarSetDelorigData(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_DECL_VARDELORIG  ((*vardelorig))     /**< frees user data of original variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL);

   var->vardelorig = vardelorig;
}

/** sets method to transform user data of the variable */
void SCIPvarSetTransData(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_DECL_VARTRANS    ((*vartrans))       /**< creates transformed user data by transforming original user data */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL);

   var->vartrans = vartrans;
}

/** sets method to free transformed user data for the variable */
void SCIPvarSetDeltransData(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_DECL_VARDELTRANS ((*vardeltrans))    /**< frees user data of transformed variable */
   )
{
   assert(var != NULL);

   var->vardeltrans = vardeltrans;
}

/** sets method to copy this variable into sub-SCIPs */
void SCIPvarSetCopyData(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_DECL_VARCOPY     ((*varcopy))        /**< copy method of the variable */
   )
{
   assert(var != NULL);

   var->varcopy = varcopy;
}

/** sets the initial flag of a variable; only possible for original or loose variables */
SCIP_RETCODE SCIPvarSetInitial(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             initial             /**< initial flag */
   )
{
   assert(var != NULL);

   if( (SCIP_VARSTATUS)var->varstatus != SCIP_VARSTATUS_ORIGINAL && (SCIP_VARSTATUS)var->varstatus != SCIP_VARSTATUS_LOOSE )
      return SCIP_INVALIDCALL;

   var->initial = initial;

   return SCIP_OKAY;
}

/** sets the removable flag of a variable; only possible for original or loose variables */
SCIP_RETCODE SCIPvarSetRemovable(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             removable           /**< removable flag */
   )
{
   assert(var != NULL);

   if( (SCIP_VARSTATUS)var->varstatus != SCIP_VARSTATUS_ORIGINAL && (SCIP_VARSTATUS)var->varstatus != SCIP_VARSTATUS_LOOSE )
      return SCIP_INVALIDCALL;

   var->removable = removable;

   return SCIP_OKAY;
}

/** gets status of variable */
SCIP_VARSTATUS SCIPvarGetStatus(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (SCIP_VARSTATUS)(var->varstatus);
}

/** returns whether the variable belongs to the original problem */
SCIP_Bool SCIPvarIsOriginal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_NEGATED || var->negatedvar != NULL);

   return (SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
      || (SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED
         && SCIPvarGetStatus(var->negatedvar) == SCIP_VARSTATUS_ORIGINAL));
}

/** returns whether the variable belongs to the transformed problem */
SCIP_Bool SCIPvarIsTransformed(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_NEGATED || var->negatedvar != NULL);

   return (SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL
      && (SCIPvarGetStatus(var) != SCIP_VARSTATUS_NEGATED
         || SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_ORIGINAL));
}

/** returns whether the variable was created by negation of a different variable */
SCIP_Bool SCIPvarIsNegated(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
}

/** gets type of variable */
SCIP_VARTYPE SCIPvarGetType(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (SCIP_VARTYPE)(var->vartype);
}

/** returns TRUE if the variable is of binary type; this is the case if:
 *  (1) variable type is binary
 *  (2) variable type is integer or implicit integer and 
 *      (i)  the lazy lower bound or the global lower bound is greater than or equal to zero
 *      (ii) the lazy upper bound or the global upper bound is less than or equal to one
 */
SCIP_Bool SCIPvarIsBinary(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || 
      (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && MAX(var->glbdom.lb, var->lazylb) >= 0.0 && MIN(var->glbdom.ub, var->lazyub) <= 1.0));
}

/** returns whether variable is of integral type (binary, integer, or implicit integer) */
SCIP_Bool SCIPvarIsIntegral(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
}

/** returns whether variable's column should be present in the initial root LP */
SCIP_Bool SCIPvarIsInitial(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->initial;
}

/** returns whether variable's column is removable from the LP (due to aging or cleanup) */
SCIP_Bool SCIPvarIsRemovable(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->removable;
}

/** returns whether the variable was deleted from the problem */
SCIP_Bool SCIPvarIsDeleted(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->deleted;
}

/** marks the variable to be deletable, i.e., it may be deleted completely from the problem;
 *  method can only be called before the variable is added to the problem by SCIPaddVar() or SCIPaddPricedVar()
 */
void SCIPvarMarkDeletable(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->probindex == -1);

   var->deletable = TRUE;
}

/** marks the variable to be not deletable from the problem */
void SCIPvarMarkNotDeletable(
   SCIP_VAR*             var
   )
{
   assert(var != NULL);

   var->deletable = FALSE;
}

/** marks variable to be deleted from global structures (cliques etc.) when cleaning up
 *
 *  @note: this is not equivalent to marking the variable itself for deletion, this is done by using SCIPvarMarkDeletable()
 */
void SCIPvarMarkDeleteGlobalStructures(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   var->delglobalstructs = TRUE;
}

/** returns whether the variable was flagged for deletion from global structures (cliques etc.) */
SCIP_Bool SCIPvarIsMarkedDeleteGlobalStructures(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->delglobalstructs;
}

/** returns whether variable is allowed to be deleted completely from the problem */
SCIP_Bool SCIPvarIsDeletable(
   SCIP_VAR*             var
   )
{
   assert(var != NULL);

   return var->deletable;
}

/** returns whether variable is an active (neither fixed nor aggregated) variable */
SCIP_Bool SCIPvarIsActive(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (var->probindex >= 0);
}

/** gets unique index of variable */
int SCIPvarGetIndex(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->index;
}

/** gets position of variable in problem, or -1 if variable is not active */
int SCIPvarGetProbindex(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->probindex;
}

/** gets transformed variable of ORIGINAL variable */
SCIP_VAR* SCIPvarGetTransVar(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL);

   return var->data.original.transvar;
}

/** gets column of COLUMN variable */
SCIP_COL* SCIPvarGetCol(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   return var->data.col;
}

/** returns whether the variable is a COLUMN variable that is member of the current LP */
SCIP_Bool SCIPvarIsInLP(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN && SCIPcolIsInLP(var->data.col));
}

/** gets aggregation variable y of an aggregated variable x = a*y + c */
SCIP_VAR* SCIPvarGetAggrVar(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.var;
}

/** gets aggregation scalar a of an aggregated variable x = a*y + c */
SCIP_Real SCIPvarGetAggrScalar(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.scalar;
}

/** gets aggregation constant c of an aggregated variable x = a*y + c */
SCIP_Real SCIPvarGetAggrConstant(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.constant;
}

/** gets number n of aggregation variables of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
int SCIPvarGetMultaggrNVars(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
   assert(!var->donotmultaggr);

   return var->data.multaggr.nvars;
}

/** gets vector of aggregation variables y of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
SCIP_VAR** SCIPvarGetMultaggrVars(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
   assert(!var->donotmultaggr);

   return var->data.multaggr.vars;
}

/** gets vector of aggregation scalars a of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
SCIP_Real* SCIPvarGetMultaggrScalars(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
   assert(!var->donotmultaggr);

   return var->data.multaggr.scalars;
}

/** gets aggregation constant c of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
SCIP_Real SCIPvarGetMultaggrConstant(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
   assert(!var->donotmultaggr);

   return var->data.multaggr.constant;
}

/** gets the negation of the given variable; may return NULL, if no negation is existing yet */
SCIP_VAR* SCIPvarGetNegatedVar(
   SCIP_VAR*             var                 /**< negated problem variable */
   )
{
   assert(var != NULL);

   return var->negatedvar;
}

/** gets the negation variable x of a negated variable x' = offset - x */
SCIP_VAR* SCIPvarGetNegationVar(
   SCIP_VAR*             var                 /**< negated problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);

   return var->negatedvar;
}

/** gets the negation offset of a negated variable x' = offset - x */
SCIP_Real SCIPvarGetNegationConstant(
   SCIP_VAR*             var                 /**< negated problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);

   return var->data.negate.constant;
}

/** gets objective function value of variable */
SCIP_Real SCIPvarGetObj(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->obj;
}

/** gets the unchanged objective function value of a variable (ignoring temproray changes performed in probing mode) */
SCIP_Real SCIPvarGetUnchangedObj(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->unchangedobj;
}

/** gets corresponding objective value of active, fixed, or multi-aggregated problem variable of given variable
 *  e.g. obj(x) = 1 this method returns for ~x the value -1
 */
SCIP_RETCODE SCIPvarGetAggregatedObj(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real*            aggrobj             /**< pointer to store the aggregated objective value */
   )
{
   SCIP_VAR* probvar = var;
   SCIP_Real mult = 1.0;

   assert(probvar != NULL);
   assert(aggrobj != NULL);

   while( probvar != NULL )
   {
      switch( SCIPvarGetStatus(probvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
	 (*aggrobj) = mult * SCIPvarGetObj(probvar);
	 return SCIP_OKAY;

      case SCIP_VARSTATUS_FIXED:
	 assert(SCIPvarGetObj(probvar) == 0.0);
	 (*aggrobj) = 0.0;
         return SCIP_OKAY;

      case SCIP_VARSTATUS_MULTAGGR:
         /* handle multi-aggregated variables depending on one variable only (possibly caused by SCIPvarFlattenAggregationGraph()) */
         if ( probvar->data.multaggr.nvars == 1 )
         {
            assert( probvar->data.multaggr.vars != NULL );
            assert( probvar->data.multaggr.scalars != NULL );
            assert( probvar->data.multaggr.vars[0] != NULL );
            mult *= probvar->data.multaggr.scalars[0];
            probvar = probvar->data.multaggr.vars[0];
            break;
         }
	 else
	 {
	    SCIP_Real tmpobj;
	    int v;

	    (*aggrobj) = 0.0;

	    for( v = probvar->data.multaggr.nvars - 1; v >= 0; --v )
	    {
	       SCIP_CALL( SCIPvarGetAggregatedObj(probvar->data.multaggr.vars[v], &tmpobj) );
	       (*aggrobj) += probvar->data.multaggr.scalars[v] * tmpobj;
	    }
	    return SCIP_OKAY;
	 }

      case SCIP_VARSTATUS_AGGREGATED:  /* x = a'*x' + c'  =>  a*x + c == (a*a')*x' + (a*c' + c) */
         assert(probvar->data.aggregate.var != NULL);
         mult *= probvar->data.aggregate.scalar;
         probvar = probvar->data.aggregate.var;
         break;

      case SCIP_VARSTATUS_NEGATED:     /* x =  - x' + c'  =>  a*x + c ==   (-a)*x' + (a*c' + c) */
         assert(probvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(probvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(probvar->negatedvar->negatedvar == probvar);
         mult *= -1.0;
         probvar = probvar->negatedvar;
         break;

      default:
	 SCIPABORT();
	 return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   return SCIP_INVALIDDATA;
}

/** gets original lower bound of original problem variable (i.e. the bound set in problem creation) */
SCIP_Real SCIPvarGetLbOriginal(
   SCIP_VAR*             var                 /**< original problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsOriginal(var));

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      return var->data.original.origdom.lb;
   else
   {
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) == SCIP_VARSTATUS_ORIGINAL);

      return var->data.negate.constant - var->negatedvar->data.original.origdom.ub;
   }
}

/** gets original upper bound of original problem variable (i.e. the bound set in problem creation) */
SCIP_Real SCIPvarGetUbOriginal(
   SCIP_VAR*             var                 /**< original problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsOriginal(var));

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      return var->data.original.origdom.ub;
   else
   {
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) == SCIP_VARSTATUS_ORIGINAL);

      return var->data.negate.constant - var->negatedvar->data.original.origdom.lb;
   }
}

/** gets the original hole list of an original variable */
SCIP_HOLELIST* SCIPvarGetHolelistOriginal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsOriginal(var));

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      return var->data.original.origdom.holelist;

   return NULL;
}

/** gets global lower bound of variable */
SCIP_Real SCIPvarGetLbGlobal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->glbdom.lb;
}

/** gets global upper bound of variable */
SCIP_Real SCIPvarGetUbGlobal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->glbdom.ub;
}

/** gets the global hole list of an active variable */
SCIP_HOLELIST* SCIPvarGetHolelistGlobal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->glbdom.holelist;
}

/** gets best global bound of variable with respect to the objective function */
SCIP_Real SCIPvarGetBestBoundGlobal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return var->glbdom.lb;
   else
      return var->glbdom.ub;
}

/** gets worst global bound of variable with respect to the objective function */
SCIP_Real SCIPvarGetWorstBoundGlobal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return var->glbdom.ub;
   else
      return var->glbdom.lb;
}

/** gets current lower bound of variable */
SCIP_Real SCIPvarGetLbLocal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->locdom.lb;
}

/** gets current upper bound of variable */
SCIP_Real SCIPvarGetUbLocal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->locdom.ub;
}

/** gets the current hole list of an active variable */
SCIP_HOLELIST* SCIPvarGetHolelistLocal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->locdom.holelist;
}

/** gets best local bound of variable with respect to the objective function */
SCIP_Real SCIPvarGetBestBoundLocal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return var->locdom.lb;
   else
      return var->locdom.ub;
}

/** gets worst local bound of variable with respect to the objective function */
SCIP_Real SCIPvarGetWorstBoundLocal(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return var->locdom.ub;
   else
      return var->locdom.lb;
}

/** gets type (lower or upper) of best bound of variable with respect to the objective function */
SCIP_BOUNDTYPE SCIPvarGetBestBoundType(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return SCIP_BOUNDTYPE_LOWER;
   else
      return SCIP_BOUNDTYPE_UPPER;
}

/** gets type (lower or upper) of worst bound of variable with respect to the objective function */
SCIP_BOUNDTYPE SCIPvarGetWorstBoundType(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return SCIP_BOUNDTYPE_UPPER;
   else
      return SCIP_BOUNDTYPE_LOWER;
}

/** gets lazy lower bound of variable, returns -infinity if the variable has no lazy lower bound */
SCIP_Real SCIPvarGetLbLazy(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->lazylb;
}

/** gets lazy upper bound of variable, returns infinity if the variable has no lazy upper bound */
SCIP_Real SCIPvarGetUbLazy(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->lazyub;
}

/** gets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
SCIP_Real SCIPvarGetBranchFactor(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->branchfactor;
}

/** gets the branch priority of the variable; variables with higher priority should always be preferred to variables
 *  with lower priority
 */
int SCIPvarGetBranchPriority(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->branchpriority;
}

/** gets the preferred branch direction of the variable (downwards, upwards, or auto) */
SCIP_BRANCHDIR SCIPvarGetBranchDirection(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (SCIP_BRANCHDIR)var->branchdirection;
}

/** gets number of variable lower bounds x >= b_i*z_i + d_i of given variable x */
int SCIPvarGetNVlbs(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return SCIPvboundsGetNVbds(var->vlbs);
}

/** gets array with bounding variables z_i in variable lower bounds x >= b_i*z_i + d_i of given variable x;
 *  the variable bounds are sorted by increasing variable index of the bounding variable z_i (see SCIPvarGetIndex())
 */
SCIP_VAR** SCIPvarGetVlbVars(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return SCIPvboundsGetVars(var->vlbs);
}

/** gets array with bounding coefficients b_i in variable lower bounds x >= b_i*z_i + d_i of given variable x */
SCIP_Real* SCIPvarGetVlbCoefs(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return SCIPvboundsGetCoefs(var->vlbs);
}

/** gets array with bounding constants d_i in variable lower bounds x >= b_i*z_i + d_i of given variable x */
SCIP_Real* SCIPvarGetVlbConstants(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return SCIPvboundsGetConstants(var->vlbs);
}

/** gets number of variable upper bounds x <= b_i*z_i + d_i of given variable x */
int SCIPvarGetNVubs(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return SCIPvboundsGetNVbds(var->vubs);
}

/** gets array with bounding variables z_i in variable upper bounds x <= b_i*z_i + d_i of given variable x;
 *  the variable bounds are sorted by increasing variable index of the bounding variable z_i (see SCIPvarGetIndex())
 */
SCIP_VAR** SCIPvarGetVubVars(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return SCIPvboundsGetVars(var->vubs);
}

/** gets array with bounding coefficients b_i in variable upper bounds x <= b_i*z_i + d_i of given variable x */
SCIP_Real* SCIPvarGetVubCoefs(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return SCIPvboundsGetCoefs(var->vubs);
}

/** gets array with bounding constants d_i in variable upper bounds x <= b_i*z_i + d_i of given variable x */
SCIP_Real* SCIPvarGetVubConstants(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return SCIPvboundsGetConstants(var->vubs);
}

/** gets number of implications  y <= b or y >= b for x == 0 or x == 1 of given active problem variable x, 
 *  there are no implications for nonbinary variable x
 */
int SCIPvarGetNImpls(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_Bool             varfixing           /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsActive(var));

   return SCIPimplicsGetNImpls(var->implics, varfixing);
}

/** gets array with implication variables y of implications  y <= b or y >= b for x == 0 or x == 1 of given active
 *  problem variable x, there are no implications for nonbinary variable x;
 *  the implications are sorted such that implications with binary implied variables precede the ones with non-binary
 *  implied variables, and as a second criteria, the implied variables are sorted by increasing variable index
 *  (see SCIPvarGetIndex())
 */
SCIP_VAR** SCIPvarGetImplVars(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_Bool             varfixing           /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsActive(var));

   return SCIPimplicsGetVars(var->implics, varfixing);
}

/** gets array with implication types of implications  y <= b or y >= b for x == 0 or x == 1 of given active problem
 *  variable x (SCIP_BOUNDTYPE_UPPER if y <= b, SCIP_BOUNDTYPE_LOWER if y >= b), 
 *  there are no implications for nonbinary variable x
 */
SCIP_BOUNDTYPE* SCIPvarGetImplTypes(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_Bool             varfixing           /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsActive(var));

   return SCIPimplicsGetTypes(var->implics, varfixing);
}

/** gets array with implication bounds b of implications  y <= b or y >= b for x == 0 or x == 1 of given active problem
 *  variable x, there are no implications for nonbinary variable x
 */
SCIP_Real* SCIPvarGetImplBounds(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_Bool             varfixing           /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsActive(var));

   return SCIPimplicsGetBounds(var->implics, varfixing);
}

/** Gets array with unique ids of implications  y <= b or y >= b for x == 0 or x == 1 of given active problem variable x,
 *  there are no implications for nonbinary variable x.
 *  If an implication is a shortcut, i.e., it was added as part of the transitive closure of another implication,
 *  its id is negative, otherwise it is nonnegative.
 */
int* SCIPvarGetImplIds(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_Bool             varfixing           /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsActive(var));

   return SCIPimplicsGetIds(var->implics, varfixing);
}

/** gets number of cliques, the active variable is contained in */
int SCIPvarGetNCliques(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_Bool             varfixing           /**< FALSE for cliques containing x == 0, TRUE for x == 1 */
   )
{
   assert(var != NULL);

   return SCIPcliquelistGetNCliques(var->cliquelist, varfixing);
}

/** gets array of cliques, the active variable is contained in */
SCIP_CLIQUE** SCIPvarGetCliques(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_Bool             varfixing           /**< FALSE for cliques containing x == 0, TRUE for x == 1 */
   )
{
   assert(var != NULL);

   return SCIPcliquelistGetCliques(var->cliquelist, varfixing);
}

/** gets primal LP solution value of variable */
SCIP_Real SCIPvarGetLPSol(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      return SCIPcolGetPrimsol(var->data.col);
   else
      return SCIPvarGetLPSol_rec(var);
}

/** gets primal NLP solution value of variable */
SCIP_Real SCIPvarGetNLPSol(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( (SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE) )
      return var->nlpsol;
   else
      return SCIPvarGetNLPSol_rec(var);
}

/** return lower bound change info at requested position */
SCIP_BDCHGINFO* SCIPvarGetBdchgInfoLb(
   SCIP_VAR*             var,                /**< problem variable */
   int                   pos                 /**< requested position */
   )
{
   assert(pos >= 0);
   assert(pos < var->nlbchginfos);

   return &var->lbchginfos[pos];
} 

/** gets the number of lower bound change info array */
int SCIPvarGetNBdchgInfosLb(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   return var->nlbchginfos;
} 

/** return upper bound change info at requested position */
SCIP_BDCHGINFO* SCIPvarGetBdchgInfoUb(
   SCIP_VAR*             var,                /**< problem variable */
   int                   pos                 /**< requested position */
   )
{
   assert(pos >= 0);
   assert(pos < var->nubchginfos);

   return &var->ubchginfos[pos];
} 

/** gets the number upper bound change info array */
int SCIPvarGetNBdchgInfosUb(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->nubchginfos;
} 

/** returns the value based history for the variable */
SCIP_VALUEHISTORY* SCIPvarGetValuehistory(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->valuehistory;
}

/** gets pseudo solution value of variable */
SCIP_Real SCIPvarGetPseudoSol(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      return SCIPvarGetBestBoundLocal(var);
   else
      return SCIPvarGetPseudoSol_rec(var);
}

/** returns the variable's VSIDS score */
SCIP_Real SCIPvarGetVSIDS(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      return SCIPhistoryGetVSIDS(var->history, dir)/stat->vsidsweight;
   else
      return SCIPvarGetVSIDS_rec(var, stat, dir);
}

/** includes event handler with given data in variable's event filter */
SCIP_RETCODE SCIPvarCatchEvent(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(var->eventfilter != NULL);
   assert((eventtype & ~SCIP_EVENTTYPE_VARCHANGED) == 0);
   assert((eventtype & SCIP_EVENTTYPE_VARCHANGED) != 0);
   assert(SCIPvarIsTransformed(var));

   SCIPsetDebugMsg(set, "catch event of type 0x%" SCIP_EVENTTYPE_FORMAT " of variable <%s> with handler %p and data %p\n",
      eventtype, var->name, (void*)eventhdlr, (void*)eventdata);

   SCIP_CALL( SCIPeventfilterAdd(var->eventfilter, blkmem, set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** deletes event handler with given data from variable's event filter */
SCIP_RETCODE SCIPvarDropEvent(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int                   filterpos           /**< position of event filter entry returned by SCIPvarCatchEvent(), or -1 */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(var->eventfilter != NULL);
   assert(SCIPvarIsTransformed(var));

   SCIPsetDebugMsg(set, "drop event of variable <%s> with handler %p and data %p\n", var->name, (void*)eventhdlr,
         (void*)eventdata);

   SCIP_CALL( SCIPeventfilterDel(var->eventfilter, blkmem, set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** returns the position of the bound change index */
int SCIPbdchgidxGetPos(
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index */
   )
{
   assert(bdchgidx != NULL);

   return bdchgidx->pos;
}

/** returns whether first bound change index belongs to an earlier applied bound change than second one */
SCIP_Bool SCIPbdchgidxIsEarlierNonNull(
   SCIP_BDCHGIDX*        bdchgidx1,          /**< first bound change index */
   SCIP_BDCHGIDX*        bdchgidx2           /**< second bound change index */
   )
{
   assert(bdchgidx1 != NULL);
   assert(bdchgidx1->depth >= -2);
   assert(bdchgidx1->pos >= 0);
   assert(bdchgidx2 != NULL);
   assert(bdchgidx2->depth >= -2);
   assert(bdchgidx2->pos >= 0);

   return (bdchgidx1->depth < bdchgidx2->depth)
      || (bdchgidx1->depth == bdchgidx2->depth && (bdchgidx1->pos < bdchgidx2->pos));
}

/** returns whether first bound change index belongs to an earlier applied bound change than second one;
 *  if a bound change index is NULL, the bound change index represents the current time, i.e. the time after the
 *  last bound change was applied to the current node
 */
SCIP_Bool SCIPbdchgidxIsEarlier(
   SCIP_BDCHGIDX*        bdchgidx1,          /**< first bound change index, or NULL */
   SCIP_BDCHGIDX*        bdchgidx2           /**< second bound change index, or NULL */
   )
{
   assert(bdchgidx1 == NULL || bdchgidx1->depth >= -2);
   assert(bdchgidx1 == NULL || bdchgidx1->pos >= 0);
   assert(bdchgidx2 == NULL || bdchgidx2->depth >= -2);
   assert(bdchgidx2 == NULL || bdchgidx2->pos >= 0);

   if( bdchgidx1 == NULL )
      return FALSE;
   else if( bdchgidx2 == NULL )
      return TRUE;
   else
      return (bdchgidx1->depth < bdchgidx2->depth)
         || (bdchgidx1->depth == bdchgidx2->depth && (bdchgidx1->pos < bdchgidx2->pos));
}

/** returns old bound that was overwritten for given bound change information */
SCIP_Real SCIPbdchginfoGetOldbound(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return bdchginfo->oldbound;
}

/** returns new bound installed for given bound change information */
SCIP_Real SCIPbdchginfoGetNewbound(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return bdchginfo->newbound;
}

/** returns variable that belongs to the given bound change information */
SCIP_VAR* SCIPbdchginfoGetVar(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return bdchginfo->var;
}

/** returns whether the bound change information belongs to a branching decision or a deduction */
SCIP_BOUNDCHGTYPE SCIPbdchginfoGetChgtype(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return (SCIP_BOUNDCHGTYPE)(bdchginfo->boundchgtype);
}

/** returns whether the bound change information belongs to a lower or upper bound change */
SCIP_BOUNDTYPE SCIPbdchginfoGetBoundtype(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return (SCIP_BOUNDTYPE)(bdchginfo->boundtype);
}

/** returns depth level of given bound change information */
int SCIPbdchginfoGetDepth(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return bdchginfo->bdchgidx.depth;
}

/** returns bound change position in its depth level of given bound change information */
int SCIPbdchginfoGetPos(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return bdchginfo->bdchgidx.pos;
}

/** returns bound change index of given bound change information */
SCIP_BDCHGIDX* SCIPbdchginfoGetIdx(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return &bdchginfo->bdchgidx;
}

/** returns inference variable of given bound change information */
SCIP_VAR* SCIPbdchginfoGetInferVar(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert((SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER
      || (SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER);

   return bdchginfo->inferencedata.var;
}

/** returns inference constraint of given bound change information */
SCIP_CONS* SCIPbdchginfoGetInferCons(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert((SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER);
   assert(bdchginfo->inferencedata.reason.cons != NULL);

   return bdchginfo->inferencedata.reason.cons;
}

/** returns inference propagator of given bound change information, or NULL if no propagator was responsible */
SCIP_PROP* SCIPbdchginfoGetInferProp(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert((SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER);

   return bdchginfo->inferencedata.reason.prop;
}

/** returns inference user information of given bound change information */
int SCIPbdchginfoGetInferInfo(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert((SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER
      || (SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER);

   return bdchginfo->inferencedata.info;
}

/** returns inference bound of inference variable of given bound change information */
SCIP_BOUNDTYPE SCIPbdchginfoGetInferBoundtype(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert((SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER
      || (SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER);

   return (SCIP_BOUNDTYPE)(bdchginfo->inferboundtype);
}

/** returns the relaxed bound change type */
SCIP_Real SCIPbdchginfoGetRelaxedBound(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to add to the conflict set */
   )
{
   return ((SCIP_BOUNDTYPE)(bdchginfo->boundtype) == SCIP_BOUNDTYPE_LOWER ? bdchginfo->var->conflictrelaxedlb : bdchginfo->var->conflictrelaxedub);
}


/** returns whether the bound change information belongs to a redundant bound change */
SCIP_Bool SCIPbdchginfoIsRedundant(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert(bdchginfo->redundant == (bdchginfo->oldbound == bdchginfo->newbound)); /*lint !e777*/

   return bdchginfo->redundant;
}

/** returns whether the bound change has an inference reason (constraint or propagator), that can be resolved */
SCIP_Bool SCIPbdchginfoHasInferenceReason(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return ((SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER)
      || ((SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER
         && bdchginfo->inferencedata.reason.prop != NULL);
}

/** for two bound change informations belonging to the same variable and bound, returns whether the first bound change
 *  has a tighter new bound as the second bound change
 */
SCIP_Bool SCIPbdchginfoIsTighter(
   SCIP_BDCHGINFO*       bdchginfo1,         /**< first bound change information */
   SCIP_BDCHGINFO*       bdchginfo2          /**< second bound change information */
   )
{
   assert(bdchginfo1 != NULL);
   assert(bdchginfo2 != NULL);
   assert(bdchginfo1->var == bdchginfo2->var);
   assert(bdchginfo1->boundtype == bdchginfo2->boundtype);

   return (SCIPbdchginfoGetBoundtype(bdchginfo1) == SCIP_BOUNDTYPE_LOWER
      ? bdchginfo1->newbound > bdchginfo2->newbound
      : bdchginfo1->newbound < bdchginfo2->newbound);
}
