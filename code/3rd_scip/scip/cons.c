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

/**@file   cons.c
 * @brief  methods for constraints and constraint handlers
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/tree.h"
#include "scip/scip.h"
#include "scip/sepastore.h"
#include "scip/cons.h"
#include "scip/branch.h"
#include "scip/reopt.h"
#include "scip/pub_misc.h"

#ifndef NDEBUG
#include "scip/struct_cons.h"
#endif


#define AGERESETAVG_INIT         1000.0 /**< initial value of the exponentially decaying weighted sum for ages */
#define AGERESETAVG_MIN          100.0  /**< minimal value to use for weighted sum of ages */
#define AGERESETAVG_DECAY        0.0005 /**< weight of a new addend in the exponentially decyaing sum */
#define AGERESETAVG_AGELIMIT     2.0    /**< in dynamic setting, a constraint is deleted if its age exceeds the
                                         *   average reset age by this factor */
#define AGERESETAVG_OBSOLETEAGE  1.8    /**< in dynamic setting, a constraint is marked obsolete if its age exceeds the
                                         *   average reset age by this factor */


/* #define CHECKCONSARRAYS */


/*
 * dynamic memory arrays
 */


/** resizes conss array to be able to store at least num constraints */
static
SCIP_RETCODE conshdlrEnsureConssMem(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->consssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conshdlr->conss, newsize) );
      conshdlr->consssize = newsize;
   }
   assert(num <= conshdlr->consssize);

   return SCIP_OKAY;
}

/** resizes initconss array to be able to store at least num constraints */
static
SCIP_RETCODE conshdlrEnsureInitconssMem(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->initconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conshdlr->initconss, newsize) );
      conshdlr->initconsssize = newsize;
   }
   assert(num <= conshdlr->initconsssize);

   return SCIP_OKAY;
}

/** resizes sepaconss array to be able to store at least num constraints */
static
SCIP_RETCODE conshdlrEnsureSepaconssMem(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->sepaconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conshdlr->sepaconss, newsize) );
      conshdlr->sepaconsssize = newsize;
   }
   assert(num <= conshdlr->sepaconsssize);

   return SCIP_OKAY;
}

/** resizes enfoconss array to be able to store at least num constraints */
static
SCIP_RETCODE conshdlrEnsureEnfoconssMem(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->enfoconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conshdlr->enfoconss, newsize) );
      conshdlr->enfoconsssize = newsize;
   }
   assert(num <= conshdlr->enfoconsssize);

   return SCIP_OKAY;
}

/** resizes checkconss array to be able to store at least num constraints */
static
SCIP_RETCODE conshdlrEnsureCheckconssMem(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->checkconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conshdlr->checkconss, newsize) );
      conshdlr->checkconsssize = newsize;
   }
   assert(num <= conshdlr->checkconsssize);

   return SCIP_OKAY;
}

/** resizes propconss array to be able to store at least num constraints */
static
SCIP_RETCODE conshdlrEnsurePropconssMem(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->propconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conshdlr->propconss, newsize) );
      conshdlr->propconsssize = newsize;
   }
   assert(num <= conshdlr->propconsssize);

   return SCIP_OKAY;
}

/** resizes updateconss array to be able to store at least num constraints */
static
SCIP_RETCODE conshdlrEnsureUpdateconssMem(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->updateconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conshdlr->updateconss, newsize) );
      conshdlr->updateconsssize = newsize;
   }
   assert(num <= conshdlr->updateconsssize);

   return SCIP_OKAY;
}




/*
 * Constraint handler methods
 */

#define checkConssArrays(conshdlr) /**/
#ifndef NDEBUG
#ifdef CHECKCONSARRAYS
#undef checkConssArrays
/** sanity check for the constraint arrays of the constraint handler (only in debug mode) */
static
void checkConssArrays(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   int c;

   assert(conshdlr != NULL);
   assert(0 <= conshdlr->nactiveconss && conshdlr->nactiveconss <= conshdlr->nconss);

   for( c = 0; c < conshdlr->nconss; ++c )
   {
      assert(conshdlr->conss[c] != NULL);
      assert(!conshdlr->conss[c]->original);
      assert(conshdlr->conss[c]->active == (c < conshdlr->nactiveconss));
      assert(conshdlr->conss[c]->consspos == c);
   }

   for( c = 0; c < conshdlr->ninitconss; ++c )
   {
      assert(conshdlr->initconss[c] != NULL);
      assert(!conshdlr->initconss[c]->original);
      assert(c < conshdlr->ninitconsskept || conshdlr->initconss[c]->active);
      assert(conshdlr->initconss[c]->initial);
   }

   for( c = 0; c < conshdlr->nsepaconss; ++c )
   {
      assert(conshdlr->sepaconss[c] != NULL);
      assert(!conshdlr->sepaconss[c]->original);
      assert(conshdlr->sepaconss[c]->active);
      assert(conshdlr->sepaconss[c]->separate);
      assert(conshdlr->sepaconss[c]->sepaenabled);
      assert(conshdlr->sepaconss[c]->obsolete == (c >= conshdlr->nusefulsepaconss));
   }

   for( c = 0; c < conshdlr->nenfoconss; ++c )
   {
      assert(conshdlr->enfoconss[c] != NULL);
      assert(!conshdlr->enfoconss[c]->original);
      assert(conshdlr->enfoconss[c]->active);
      assert(conshdlr->enfoconss[c]->enforce);
      assert(conshdlr->enfoconss[c]->obsolete == (c >= conshdlr->nusefulenfoconss));
   }

   for( c = 0; c < conshdlr->ncheckconss; ++c )
   {
      assert(conshdlr->checkconss[c] != NULL);
      assert(!conshdlr->checkconss[c]->original);
      assert(conshdlr->checkconss[c]->active);
      assert(conshdlr->checkconss[c]->check);
      assert(conshdlr->checkconss[c]->obsolete == (c >= conshdlr->nusefulcheckconss));
   }

   for( c = 0; c < conshdlr->npropconss; ++c )
   {
      assert(conshdlr->propconss[c] != NULL);
      assert(!conshdlr->propconss[c]->original);
      assert(conshdlr->propconss[c]->active);
      assert(conshdlr->propconss[c]->propagate);
      assert(conshdlr->propconss[c]->propenabled);
      assert(conshdlr->propconss[c]->markpropagate == (c < conshdlr->nmarkedpropconss));
      assert(conshdlr->propconss[c]->markpropagate || (conshdlr->propconss[c]->obsolete == (c >= conshdlr->nusefulpropconss)));
   }
   assert(conshdlr->nmarkedpropconss <= conshdlr->npropconss);
}
#endif
#endif

/** returns whether the constraint updates of the constraint handler are currently delayed */
static
SCIP_Bool conshdlrAreUpdatesDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   return (conshdlr->delayupdatecount > 0);
}

/** returns the exponentially decaying weighted age average for age resets */
static
SCIP_Real conshdlrGetAgeresetavg(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return MAX(conshdlr->ageresetavg, AGERESETAVG_MIN);
}

/** updates the exponentially decaying weighted age average for age resets after a constraint age was reset */
static
void conshdlrUpdateAgeresetavg(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_Real             age                 /**< age of the constraint that is reset to zero */
   )
{
   assert(conshdlr != NULL);

   conshdlr->ageresetavg *= (1.0-AGERESETAVG_DECAY);
   conshdlr->ageresetavg += AGERESETAVG_DECAY * age;
}

/** returns whether the constraint's age exceeds the age limit */
static
SCIP_Bool consExceedsAgelimit(
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(set != NULL);

   return (cons->dynamic
      && ((set->cons_agelimit > 0 && cons->age > set->cons_agelimit)
         || (set->cons_agelimit == 0 && cons->age > AGERESETAVG_AGELIMIT * conshdlrGetAgeresetavg(cons->conshdlr))));
}

/** returns whether the constraint's age exceeds the obsolete age limit */
static
SCIP_Bool consExceedsObsoleteage(
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(set != NULL);

   return (cons->dynamic
      && ((set->cons_obsoleteage > 0 && cons->age > set->cons_obsoleteage)
         || (set->cons_obsoleteage == 0 && cons->age > AGERESETAVG_OBSOLETEAGE * conshdlrGetAgeresetavg(cons->conshdlr))));
}

/** marks constraint to be obsolete; it will be moved to the last part of the constraint arrays, such that
 *  it is checked, enforced, separated, and propagated after the useful constraints
 */
static
SCIP_RETCODE conshdlrMarkConsObsolete(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to be marked obsolete */
   )
{
   SCIP_CONS* tmpcons;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(!cons->original);
   assert(!cons->obsolete);
   assert(!conshdlrAreUpdatesDelayed(conshdlr));

   cons->obsolete = TRUE;

   if( cons->active )
   {
      if( cons->check )
      {
         assert(0 <= cons->checkconsspos && cons->checkconsspos < conshdlr->nusefulcheckconss);

         /* switch the last useful (non-obsolete) check constraint with this constraint */
         tmpcons = conshdlr->checkconss[conshdlr->nusefulcheckconss-1];
         assert(tmpcons->checkconsspos == conshdlr->nusefulcheckconss-1);

         conshdlr->checkconss[conshdlr->nusefulcheckconss-1] = cons;
         conshdlr->checkconss[cons->checkconsspos] = tmpcons;
         tmpcons->checkconsspos = cons->checkconsspos;
         cons->checkconsspos = conshdlr->nusefulcheckconss-1;

         conshdlr->nusefulcheckconss--;
      }
   }
   if( cons->enabled )
   {
      if( cons->separate && cons->sepaenabled )
      {
         assert(0 <= cons->sepaconsspos && cons->sepaconsspos < conshdlr->nusefulsepaconss);

         if( cons->sepaconsspos < conshdlr->lastnusefulsepaconss )
            conshdlr->lastnusefulsepaconss--;

         /* switch the last useful (non-obsolete) sepa constraint with this constraint */
         tmpcons = conshdlr->sepaconss[conshdlr->nusefulsepaconss-1];
         assert(tmpcons->sepaconsspos == conshdlr->nusefulsepaconss-1);

         conshdlr->sepaconss[conshdlr->nusefulsepaconss-1] = cons;
         conshdlr->sepaconss[cons->sepaconsspos] = tmpcons;
         tmpcons->sepaconsspos = cons->sepaconsspos;
         cons->sepaconsspos = conshdlr->nusefulsepaconss-1;

         conshdlr->nusefulsepaconss--;
      }
      if( cons->enforce )
      {
         assert(0 <= cons->enfoconsspos && cons->enfoconsspos < conshdlr->nusefulenfoconss);

         if( cons->enfoconsspos < conshdlr->lastnusefulenfoconss )
            conshdlr->lastnusefulenfoconss--;
         else
         {
            /* the constraint that becomes obsolete is not yet enforced on the current solution:
             * we have to make sure that it will be enforced the next time; this is not done, if the current
             * solution was already enforced and only enforcement on the additional constraints is performed
             * (because in this case, only the new useful constraints are enforced);
             * thus, we have to reset the enforcement counters in order to enforce all constraints again, especially
             * the now obsolete one;
             * this case should occur almost never, because a constraint that was not enforced in the last enforcement
             * is a newly added one, and it is very unlikely that this constraint will become obsolete before the next
             * enforcement call;
             * this reset is not performed for separation and propagation, because they are not vital for correctness
             */
            conshdlr->lastenfolplpcount = -1;
            conshdlr->lastenfolpdomchgcount = -1;
            conshdlr->lastenfopsdomchgcount = -1;
            conshdlr->lastenforelaxdomchgcount = -1;
            conshdlr->lastenforelaxrelaxcount = -1;
            conshdlr->lastenfolpnode = -1;
            conshdlr->lastenfopsnode = -1;
         }

         /* switch the last useful (non-obsolete) enfo constraint with this constraint */
         tmpcons = conshdlr->enfoconss[conshdlr->nusefulenfoconss-1];
         assert(tmpcons->enfoconsspos == conshdlr->nusefulenfoconss-1);

         conshdlr->enfoconss[conshdlr->nusefulenfoconss-1] = cons;
         conshdlr->enfoconss[cons->enfoconsspos] = tmpcons;
         tmpcons->enfoconsspos = cons->enfoconsspos;
         cons->enfoconsspos = conshdlr->nusefulenfoconss-1;

         conshdlr->nusefulenfoconss--;
      }
      /* in case the constraint is marked to be propagated, we do not move it in the propconss array since the first
       * part of the array contains all marked constraints independently of their age
       */
      assert((!cons->markpropagate) == (cons->propconsspos < conshdlr->nmarkedpropconss));
      if( cons->propagate && cons->propenabled && !cons->markpropagate )
      {
         assert(0 <= cons->propconsspos && cons->propconsspos < conshdlr->nusefulpropconss);

         if( cons->propconsspos < conshdlr->lastnusefulpropconss )
            conshdlr->lastnusefulpropconss--;

         /* switch the last useful (non-obsolete) prop constraint with this constraint */
         tmpcons = conshdlr->propconss[conshdlr->nusefulpropconss-1];
         assert(tmpcons->propconsspos == conshdlr->nusefulpropconss-1);

         conshdlr->propconss[conshdlr->nusefulpropconss-1] = cons;
         conshdlr->propconss[cons->propconsspos] = tmpcons;
         tmpcons->propconsspos = cons->propconsspos;
         cons->propconsspos = conshdlr->nusefulpropconss-1;

         conshdlr->nusefulpropconss--;
      }
   }

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** marks obsolete constraint to be not obsolete anymore;
 *  it will be moved to the first part of the constraint arrays, such that it is checked, enforced, separated,
 *  and propagated before the obsolete constraints
 */
static
SCIP_RETCODE conshdlrMarkConsUseful(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to be marked obsolete */
   )
{
   SCIP_CONS* tmpcons;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(!cons->original);
   assert(cons->obsolete);
   assert(!conshdlrAreUpdatesDelayed(conshdlr));

   cons->obsolete = FALSE;

   if( cons->active )
   {
      if( cons->check )
      {
         assert(conshdlr->nusefulcheckconss <= cons->checkconsspos && cons->checkconsspos < conshdlr->ncheckconss);

         /* switch the first obsolete check constraint with this constraint */
         tmpcons = conshdlr->checkconss[conshdlr->nusefulcheckconss];
         assert(tmpcons->checkconsspos == conshdlr->nusefulcheckconss);

         conshdlr->checkconss[conshdlr->nusefulcheckconss] = cons;
         conshdlr->checkconss[cons->checkconsspos] = tmpcons;
         tmpcons->checkconsspos = cons->checkconsspos;
         cons->checkconsspos = conshdlr->nusefulcheckconss;

         conshdlr->nusefulcheckconss++;
      }
   }
   if( cons->enabled )
   {
      if( cons->separate && cons->sepaenabled )
      {
         assert(conshdlr->nusefulsepaconss <= cons->sepaconsspos && cons->sepaconsspos < conshdlr->nsepaconss);

         /* switch the first obsolete sepa constraint with this constraint */
         tmpcons = conshdlr->sepaconss[conshdlr->nusefulsepaconss];
         assert(tmpcons->sepaconsspos == conshdlr->nusefulsepaconss);

         conshdlr->sepaconss[conshdlr->nusefulsepaconss] = cons;
         conshdlr->sepaconss[cons->sepaconsspos] = tmpcons;
         tmpcons->sepaconsspos = cons->sepaconsspos;
         cons->sepaconsspos = conshdlr->nusefulsepaconss;

         conshdlr->nusefulsepaconss++;
      }
      if( cons->enforce )
      {
         assert(conshdlr->nusefulenfoconss <= cons->enfoconsspos && cons->enfoconsspos < conshdlr->nenfoconss);

         /* switch the first obsolete enfo constraint with this constraint */
         tmpcons = conshdlr->enfoconss[conshdlr->nusefulenfoconss];
         assert(tmpcons->enfoconsspos == conshdlr->nusefulenfoconss);

         conshdlr->enfoconss[conshdlr->nusefulenfoconss] = cons;
         conshdlr->enfoconss[cons->enfoconsspos] = tmpcons;
         tmpcons->enfoconsspos = cons->enfoconsspos;
         cons->enfoconsspos = conshdlr->nusefulenfoconss;

         conshdlr->nusefulenfoconss++;
      }
      /* in case the constraint is marked to be propagated, we do not move it in the propconss array since the first
       * part of the array contains all marked constraints independently of their age
       */
      assert((!cons->markpropagate) == (cons->propconsspos < conshdlr->nmarkedpropconss));
      if( cons->propagate && cons->propenabled && !cons->markpropagate)
      {
         assert(conshdlr->nusefulpropconss <= cons->propconsspos && cons->propconsspos < conshdlr->npropconss);

         /* switch the first obsolete prop constraint with this constraint */
         tmpcons = conshdlr->propconss[conshdlr->nusefulpropconss];
         assert(tmpcons->propconsspos == conshdlr->nusefulpropconss);

         conshdlr->propconss[conshdlr->nusefulpropconss] = cons;
         conshdlr->propconss[cons->propconsspos] = tmpcons;
         tmpcons->propconsspos = cons->propconsspos;
         cons->propconsspos = conshdlr->nusefulpropconss;

         conshdlr->nusefulpropconss++;
      }
   }

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** marks constraint to be propagated in the next propagation round;
 *
 *  @note the propagation array is divided into three parts in contrast to the other constraint arrays;
 *        the first part contains constraints which were marked to be propagated (independently of its age)
 *        the second part contains the useful (non-obsolete) constraints which are not marked to be propagated
 *        finally, the third part contains obsolete constraints which are not marked to be propagated
 *
 *  @note if a constraint gets marked for propagation we put it into the first part regardless of its age
 */
static
void conshdlrMarkConsPropagate(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to be marked obsolete */
   )
{
   SCIP_CONS* tmpcons;

   assert(conshdlr != NULL);
   assert(conshdlr->nmarkedpropconss <= conshdlr->nusefulpropconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(!cons->original);

   /* it may happen that the constraint is deleted while updates are delayed: in this case we just return */
   if( !cons->enabled )
      return;

   if( cons->markpropagate )
      return;

   cons->markpropagate = TRUE;

   /* propagation of the constraint is globally or locally disabled, so we do not have to move the constraint in the
    * propconss array
    */
   if( !cons->propagate || !cons->propenabled )
   {
      assert(cons->propconsspos == -1);
      return;
   }
   assert(cons->propconsspos >= conshdlr->nmarkedpropconss);

   /* if the constraint is obsolete, we need to move it first to the non-obsolete part of the array */
   if( cons->obsolete )
   {
      assert(conshdlr->nusefulpropconss <= cons->propconsspos && cons->propconsspos < conshdlr->npropconss);

      /* switch the first obsolete prop constraint with this constraint */
      tmpcons = conshdlr->propconss[conshdlr->nusefulpropconss];
      assert(tmpcons->propconsspos == conshdlr->nusefulpropconss);

      conshdlr->propconss[conshdlr->nusefulpropconss] = cons;
      conshdlr->propconss[cons->propconsspos] = tmpcons;
      tmpcons->propconsspos = cons->propconsspos;
      cons->propconsspos = conshdlr->nusefulpropconss;

      conshdlr->nusefulpropconss++;
   }
   assert(conshdlr->nmarkedpropconss <= cons->propconsspos && cons->propconsspos < conshdlr->nusefulpropconss);

   /* switch the first useful prop constraint with this constraint */
   tmpcons = conshdlr->propconss[conshdlr->nmarkedpropconss];
   assert(tmpcons->propconsspos == conshdlr->nmarkedpropconss);

   conshdlr->propconss[conshdlr->nmarkedpropconss] = cons;
   conshdlr->propconss[cons->propconsspos] = tmpcons;
   tmpcons->propconsspos = cons->propconsspos;
   cons->propconsspos = conshdlr->nmarkedpropconss;

   conshdlr->nmarkedpropconss++;
   assert(conshdlr->nmarkedpropconss <= conshdlr->npropconss);

   checkConssArrays(conshdlr);
}

/** unmarks constraint to be propagated in the next propagation round;
 *
 *  @note the propagation array is divided into three parts in contrast to the other constraint arrays;
 *        the first part contains constraints which were marked to be propagated (independently of its age)
 *        the second part contains the useful (non-obsolete) constraints which are not marked to be propagated
 *        finally, the third part contains obsolete constraints which are not marked to be propagated
 *
 *  @note if a constraint gets unmarked for propagation, it is put into the right part depending on its age
 */
static
void conshdlrUnmarkConsPropagate(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to be marked obsolete */
   )
{
   SCIP_CONS* tmpcons;

   assert(conshdlr != NULL);
   assert(conshdlr->nmarkedpropconss <= conshdlr->nusefulpropconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(!cons->original);

   /* it may happen that the constraint is deleted while updates are delayed: in this case we just return */
   if( !cons->enabled )
      return;

   if( !cons->markpropagate )
      return;

   cons->markpropagate = FALSE;

   /* propagation of the constraint is globally or locally disabled, so we do not have to move the constraint in the
    * propconss array
    */
   if( !cons->propagate || !cons->propenabled )
   {
      assert(cons->propconsspos == -1);
      return;
   }
   assert(cons->propconsspos >= 0);
   assert(cons->propconsspos < conshdlr->nmarkedpropconss);

   /* first, move the constraint out of the first part to the second part of the constraint array */
   if( cons->propconsspos < conshdlr->nmarkedpropconss - 1 )
   {
      conshdlr->nmarkedpropconss--;

      /* switch the last marked prop constraint with this constraint */
      tmpcons = conshdlr->propconss[conshdlr->nmarkedpropconss];
      assert(tmpcons->propconsspos == conshdlr->nmarkedpropconss);

      conshdlr->propconss[conshdlr->nmarkedpropconss] = cons;
      conshdlr->propconss[cons->propconsspos] = tmpcons;
      tmpcons->propconsspos = cons->propconsspos;
      cons->propconsspos = conshdlr->nmarkedpropconss;
   }
   else if( cons->propconsspos == conshdlr->nmarkedpropconss - 1 )
      conshdlr->nmarkedpropconss--;
   assert(cons->propconsspos == conshdlr->nmarkedpropconss);

   /* if the constraint is obsolete, move it to the last part of the constraint array */
   if( cons->obsolete )
   {
      conshdlr->nusefulpropconss--;

      /* switch the last useful prop constraint with this constraint */
      tmpcons = conshdlr->propconss[conshdlr->nusefulpropconss];
      assert(tmpcons->propconsspos == conshdlr->nusefulpropconss);

      conshdlr->propconss[conshdlr->nusefulpropconss] = cons;
      conshdlr->propconss[cons->propconsspos] = tmpcons;
      tmpcons->propconsspos = cons->propconsspos;
      cons->propconsspos = conshdlr->nusefulpropconss;
   }

   checkConssArrays(conshdlr);
}


/** adds constraint to the conss array of constraint handler */
static
SCIP_RETCODE conshdlrAddCons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(!cons->active);
   assert(cons->consspos == -1);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   /* insert the constraint as inactive constraint into the transformed constraints array */
   SCIP_CALL( conshdlrEnsureConssMem(conshdlr, set, conshdlr->nconss+1) );
   conshdlr->conss[conshdlr->nconss] = cons;
   cons->consspos = conshdlr->nconss;
   conshdlr->nconss++;

   return SCIP_OKAY;
}

/** deletes constraint from the conss array of constraint handler */
static
void conshdlrDelCons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to remove */
   )
{
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(!cons->active);
   assert(conshdlr->nactiveconss <= cons->consspos && cons->consspos < conshdlr->nconss);

   conshdlr->conss[cons->consspos] = conshdlr->conss[conshdlr->nconss-1];
   conshdlr->conss[cons->consspos]->consspos = cons->consspos;
   conshdlr->nconss--;
   cons->consspos = -1;
}

/** adds constraint to the initconss array of constraint handler */
static
SCIP_RETCODE conshdlrAddInitcons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->active);
   assert(cons->initial);
   assert(cons->initconsspos == -1 || cons->initconsspos < conshdlr->ninitconsskept);

   SCIP_CALL( conshdlrEnsureInitconssMem(conshdlr, set, conshdlr->ninitconss+1) );

   insertpos = conshdlr->ninitconss;

   conshdlr->initconss[insertpos] = cons;
   conshdlr->ninitconss++;
   stat->ninitconssadded++;

   /* if the constraint is kept, we keep the stored position at the beginning of the array */
   if( cons->initconsspos == -1 )
      cons->initconsspos = insertpos;

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** deletes constraint from the initconss array of constraint handler */
static
void conshdlrDelInitcons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to remove */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(0 <= cons->initconsspos && cons->initconsspos < conshdlr->ninitconss);

   delpos = cons->initconsspos;
   if( delpos < conshdlr->ninitconsskept )
   {
      conshdlr->ninitconsskept--;
      conshdlr->initconss[delpos] = conshdlr->initconss[conshdlr->ninitconsskept];
      conshdlr->initconss[delpos]->initconsspos = delpos;
      delpos = conshdlr->ninitconsskept;
   }

   if( delpos < conshdlr->ninitconss-1 )
   {
      conshdlr->initconss[delpos] = conshdlr->initconss[conshdlr->ninitconss-1];
      conshdlr->initconss[delpos]->initconsspos = delpos;
   }
   conshdlr->ninitconss--;
   cons->initconsspos = -1;

   checkConssArrays(conshdlr);
}

/** adds constraint to the sepaconss array of constraint handler */
static
SCIP_RETCODE conshdlrAddSepacons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->active);
   assert(cons->separate);
   assert(cons->sepaenabled);
   assert(cons->sepaconsspos == -1);
   assert(set != NULL);
   assert(cons->scip == set->scip);
   assert(!conshdlrAreUpdatesDelayed(conshdlr) || !conshdlr->duringsepa);

   SCIP_CALL( conshdlrEnsureSepaconssMem(conshdlr, set, conshdlr->nsepaconss+1) );
   insertpos = conshdlr->nsepaconss;
   if( !cons->obsolete )
   {
      if( conshdlr->nusefulsepaconss < conshdlr->nsepaconss )
      {
         conshdlr->sepaconss[conshdlr->nsepaconss] = conshdlr->sepaconss[conshdlr->nusefulsepaconss];
         conshdlr->sepaconss[conshdlr->nsepaconss]->sepaconsspos = conshdlr->nsepaconss;
         insertpos = conshdlr->nusefulsepaconss;
      }
      conshdlr->nusefulsepaconss++;
   }
   conshdlr->sepaconss[insertpos] = cons;
   cons->sepaconsspos = insertpos;
   conshdlr->nsepaconss++;

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** deletes constraint from the sepaconss array of constraint handler */
static
void conshdlrDelSepacons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to remove */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->separate);
   assert(cons->sepaenabled);
   assert(cons->sepaconsspos != -1);
   assert(!conshdlrAreUpdatesDelayed(conshdlr) || !conshdlr->duringsepa);

   delpos = cons->sepaconsspos;
   if( !cons->obsolete )
   {
      assert(0 <= delpos && delpos < conshdlr->nusefulsepaconss);

      if( delpos < conshdlr->lastnusefulsepaconss )
         conshdlr->lastnusefulsepaconss--;

      conshdlr->sepaconss[delpos] = conshdlr->sepaconss[conshdlr->nusefulsepaconss-1];
      conshdlr->sepaconss[delpos]->sepaconsspos = delpos;
      delpos = conshdlr->nusefulsepaconss-1;
      conshdlr->nusefulsepaconss--;
      assert(conshdlr->nusefulsepaconss >= 0);
      assert(conshdlr->lastnusefulsepaconss >= 0);
   }
   assert(conshdlr->nusefulsepaconss <= delpos && delpos < conshdlr->nsepaconss);
   if( delpos < conshdlr->nsepaconss-1 )
   {
      conshdlr->sepaconss[delpos] = conshdlr->sepaconss[conshdlr->nsepaconss-1];
      conshdlr->sepaconss[delpos]->sepaconsspos = delpos;
   }
   conshdlr->nsepaconss--;
   cons->sepaconsspos = -1;

   checkConssArrays(conshdlr);
}

/** adds constraint to the enfoconss array of constraint handler */
static
SCIP_RETCODE conshdlrAddEnfocons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->active);
   assert(cons->enforce);
   assert(cons->enfoconsspos == -1);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   SCIP_CALL( conshdlrEnsureEnfoconssMem(conshdlr, set, conshdlr->nenfoconss+1) );
   insertpos = conshdlr->nenfoconss;
   if( !cons->obsolete )
   {
      if( conshdlr->nusefulenfoconss < conshdlr->nenfoconss )
      {
         conshdlr->enfoconss[conshdlr->nenfoconss] = conshdlr->enfoconss[conshdlr->nusefulenfoconss];
         conshdlr->enfoconss[conshdlr->nenfoconss]->enfoconsspos = conshdlr->nenfoconss;
         insertpos = conshdlr->nusefulenfoconss;
      }
      conshdlr->nusefulenfoconss++;
   }
   else
   {
      /* we have to make sure that even this obsolete constraint is enforced in the next enforcement call;
       * if the same LP or pseudo solution is enforced again, only the newly added useful constraints are
       * enforced; thus, we have to reset the enforcement counters and force all constraints to be 
       * enforced again; this is not needed for separation and propagation, because they are not vital for correctness
       */
      conshdlr->lastenfolplpcount = -1;
      conshdlr->lastenfolpdomchgcount = -1;
      conshdlr->lastenfopsdomchgcount = -1;
      conshdlr->lastenforelaxdomchgcount = -1;
      conshdlr->lastenforelaxrelaxcount = -1;
      conshdlr->lastenfolpnode = -1;
      conshdlr->lastenfopsnode = -1;
   }
   conshdlr->enfoconss[insertpos] = cons;
   cons->enfoconsspos = insertpos;
   conshdlr->nenfoconss++;

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** deletes constraint from the enfoconss array of constraint handler */
static
void conshdlrDelEnfocons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to remove */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->enforce);
   assert(cons->enfoconsspos != -1);

   delpos = cons->enfoconsspos;
   if( !cons->obsolete )
   {
      assert(0 <= delpos && delpos < conshdlr->nusefulenfoconss);

      if( delpos < conshdlr->lastnusefulenfoconss )
         conshdlr->lastnusefulenfoconss--;

      conshdlr->enfoconss[delpos] = conshdlr->enfoconss[conshdlr->nusefulenfoconss-1];
      conshdlr->enfoconss[delpos]->enfoconsspos = delpos;
      delpos = conshdlr->nusefulenfoconss-1;
      conshdlr->nusefulenfoconss--;

      /* if the constraint that moved to the free position was a newly added constraint and not enforced in the last
       * enforcement, we have to make sure it will be enforced in the next run;
       * this check is not performed for separation and propagation, because they are not vital for correctness
       */
      if( delpos >= conshdlr->lastnusefulenfoconss )
         conshdlr->lastnusefulenfoconss = cons->enfoconsspos;

      assert(conshdlr->nusefulenfoconss >= 0);
      assert(conshdlr->lastnusefulenfoconss >= 0);
   }
   assert(conshdlr->nusefulenfoconss <= delpos && delpos < conshdlr->nenfoconss);
   if( delpos < conshdlr->nenfoconss-1 )
   {
      conshdlr->enfoconss[delpos] = conshdlr->enfoconss[conshdlr->nenfoconss-1];
      conshdlr->enfoconss[delpos]->enfoconsspos = delpos;
   }
   conshdlr->nenfoconss--;
   cons->enfoconsspos = -1;

   checkConssArrays(conshdlr);
}

/** adds constraint to the checkconss array of constraint handler */
static
SCIP_RETCODE conshdlrAddCheckcons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->active);
   assert(cons->check);
   assert(cons->checkconsspos == -1);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   SCIP_CALL( conshdlrEnsureCheckconssMem(conshdlr, set, conshdlr->ncheckconss+1) );
   insertpos = conshdlr->ncheckconss;
   if( !cons->obsolete )
   {
      if( conshdlr->nusefulcheckconss < conshdlr->ncheckconss )
      {
         assert(conshdlr->checkconss[conshdlr->nusefulcheckconss] != NULL);
         conshdlr->checkconss[conshdlr->ncheckconss] = conshdlr->checkconss[conshdlr->nusefulcheckconss];
         conshdlr->checkconss[conshdlr->ncheckconss]->checkconsspos = conshdlr->ncheckconss;
         insertpos = conshdlr->nusefulcheckconss;
      }
      conshdlr->nusefulcheckconss++;
   }
   assert(0 <= insertpos && insertpos <= conshdlr->ncheckconss);
   conshdlr->checkconss[insertpos] = cons;
   cons->checkconsspos = insertpos;
   conshdlr->ncheckconss++;

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** deletes constraint from the checkconss array of constraint handler */
static
void conshdlrDelCheckcons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->active);
   assert(cons->check);
   assert(cons->checkconsspos != -1);

   delpos = cons->checkconsspos;
   if( !cons->obsolete )
   {
      assert(0 <= delpos && delpos < conshdlr->nusefulcheckconss);
      conshdlr->checkconss[delpos] = conshdlr->checkconss[conshdlr->nusefulcheckconss-1];
      conshdlr->checkconss[delpos]->checkconsspos = delpos;
      delpos = conshdlr->nusefulcheckconss-1;
      conshdlr->nusefulcheckconss--;
   }
   assert(conshdlr->nusefulcheckconss <= delpos && delpos < conshdlr->ncheckconss);
   if( delpos < conshdlr->ncheckconss-1 )
   {
      conshdlr->checkconss[delpos] = conshdlr->checkconss[conshdlr->ncheckconss-1];
      conshdlr->checkconss[delpos]->checkconsspos = delpos;
   }
   conshdlr->ncheckconss--;
   cons->checkconsspos = -1;

   checkConssArrays(conshdlr);
}

/** adds constraint to the propconss array of constraint handler */
static
SCIP_RETCODE conshdlrAddPropcons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->active);
   assert(cons->enabled);
   assert(cons->propagate);
   assert(cons->propenabled);
   assert(cons->propconsspos == -1);
   assert(set != NULL);
   assert(cons->scip == set->scip);
   assert(!conshdlrAreUpdatesDelayed(conshdlr) || !conshdlr->duringprop);

   /* add constraint to the propagation array */
   SCIP_CALL( conshdlrEnsurePropconssMem(conshdlr, set, conshdlr->npropconss+1) );
   insertpos = conshdlr->npropconss;
   if( !cons->obsolete )
   {
      if( conshdlr->nusefulpropconss < conshdlr->npropconss )
      {
         conshdlr->propconss[conshdlr->npropconss] = conshdlr->propconss[conshdlr->nusefulpropconss];
         conshdlr->propconss[conshdlr->npropconss]->propconsspos = conshdlr->npropconss;
         insertpos = conshdlr->nusefulpropconss;
      }
      conshdlr->nusefulpropconss++;
   }
   conshdlr->propconss[insertpos] = cons;
   cons->propconsspos = insertpos;
   conshdlr->npropconss++;

   /* if the constraint is marked to be propagated, we have to move it to the first part of the array */
   if( cons->markpropagate )
   {
      /* temporarily unmark the constraint to be propagated, such that we can use the method below */
      cons->markpropagate = FALSE;

      conshdlrMarkConsPropagate(cons->conshdlr, cons);
      assert(cons->markpropagate);
   }

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** deletes constraint from the propconss array of constraint handler */
static
void conshdlrDelPropcons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to remove */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->propagate);
   assert(cons->propenabled);
   assert(cons->propconsspos != -1);
   assert(!conshdlrAreUpdatesDelayed(conshdlr) || !conshdlr->duringprop);

   /* unmark constraint to be propagated; this will move the constraint to the obsolete or non-obsolete part of the
    * array, depending on its age
    */
   if( cons->markpropagate )
   {
      conshdlrUnmarkConsPropagate(cons->conshdlr, cons);
      assert(!cons->markpropagate);
   }

   /* delete constraint from the propagation array */
   delpos = cons->propconsspos;
   assert(delpos >= conshdlr->nmarkedpropconss);
   if( !cons->obsolete )
   {
      assert(0 <= delpos && delpos < conshdlr->nusefulpropconss);

      if( delpos < conshdlr->lastnusefulpropconss )
         conshdlr->lastnusefulpropconss--;

      conshdlr->propconss[delpos] = conshdlr->propconss[conshdlr->nusefulpropconss-1];
      conshdlr->propconss[delpos]->propconsspos = delpos;
      delpos = conshdlr->nusefulpropconss-1;
      conshdlr->nusefulpropconss--;
      assert(conshdlr->nusefulpropconss >= 0);
      assert(conshdlr->lastnusefulpropconss >= 0);
   }
   assert(conshdlr->nusefulpropconss <= delpos && delpos < conshdlr->npropconss);

   if( delpos < conshdlr->npropconss-1 )
   {
      conshdlr->propconss[delpos] = conshdlr->propconss[conshdlr->npropconss-1];
      conshdlr->propconss[delpos]->propconsspos = delpos;
   }
   conshdlr->npropconss--;
   cons->propconsspos = -1;
   assert(conshdlr->nmarkedpropconss <= conshdlr->npropconss);

   checkConssArrays(conshdlr);
}

/** enables separation of constraint */
static
SCIP_RETCODE conshdlrEnableConsSeparation(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->sepaenabled);
   assert(cons->sepaconsspos == -1);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   SCIPsetDebugMsg(set, "enable separation of constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* enable separation of constraint */
   cons->sepaenabled = TRUE;

   /* add constraint to the separation array */
   if( cons->enabled && cons->separate )
   {
      SCIP_CALL( conshdlrAddSepacons(conshdlr, set, cons) );
   }

   return SCIP_OKAY;
}

/** disables separation of constraint */
static
SCIP_RETCODE conshdlrDisableConsSeparation(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to remove */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->sepaenabled);
   assert((cons->separate && cons->enabled) == (cons->sepaconsspos != -1));

   SCIPdebugMessage("disable separation of constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* delete constraint from the separation array */
   if( cons->separate && cons->enabled )
   {
      conshdlrDelSepacons(conshdlr, cons);
   }
   assert(cons->sepaconsspos == -1);

   /* disable separation of constraint */
   cons->sepaenabled = FALSE;

   return SCIP_OKAY;
}

/** enables propagation of constraint */
static
SCIP_RETCODE conshdlrEnableConsPropagation(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->propenabled);
   assert(cons->propconsspos == -1);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   SCIPsetDebugMsg(set, "enable propagation of constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* enable propagation of constraint */
   cons->propenabled = TRUE;

   /* add constraint to the propagation array */
   if( cons->enabled && cons->propagate )
   {
      SCIP_CALL( conshdlrAddPropcons(conshdlr, set, cons) );
   }

   return SCIP_OKAY;
}

/** disables propagation of constraint */
static
SCIP_RETCODE conshdlrDisableConsPropagation(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint to remove */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->propenabled);
   assert((cons->propagate && cons->enabled) == (cons->propconsspos != -1));

   SCIPdebugMessage("disable propagation of constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* delete constraint from the propagation array */
   if( cons->propagate && cons->enabled )
   {
      conshdlrDelPropcons(conshdlr, cons);
   }
   assert(cons->propconsspos == -1);

   /* disable propagation of constraint */
   cons->propenabled = FALSE;

   return SCIP_OKAY;
}

/** enables separation, enforcement, and propagation of constraint */
static
SCIP_RETCODE conshdlrEnableCons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(stat != NULL);
   assert(cons != NULL);
   assert(cons->scip == set->scip);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->active);
   assert(!cons->enabled);
   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->propconsspos == -1);

   SCIPsetDebugMsg(set, "enable constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* enable constraint */
   cons->enabled = TRUE;
   conshdlr->nenabledconss++;
   stat->nenabledconss++;

   /* add constraint to the separation array */
   if( cons->separate && cons->sepaenabled )
   {
      SCIP_CALL( conshdlrAddSepacons(conshdlr, set, cons) );
   }

   /* add constraint to the enforcement array */
   if( cons->enforce )
   {
      SCIP_CALL( conshdlrAddEnfocons(conshdlr, set, cons) );
   }

   /* add constraint to the propagation array */
   if( cons->propagate && cons->propenabled )
   {
      SCIP_CALL( conshdlrAddPropcons(conshdlr, set, cons) );
   }

   /* call constraint handler's enabling notification method */
   if( conshdlr->consenable != NULL )
   {
      SCIP_CALL( conshdlr->consenable(set->scip, conshdlr, cons) );
   }

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** disables separation, enforcement, and propagation of constraint */
static
SCIP_RETCODE conshdlrDisableCons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_CONS*            cons                /**< constraint to remove */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(stat != NULL);
   assert(cons != NULL);
   assert(cons->scip == set->scip);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->active);
   assert(cons->enabled);
   assert((cons->separate && cons->sepaenabled) == (cons->sepaconsspos != -1));
   assert(cons->enforce == (cons->enfoconsspos != -1));
   assert((cons->propagate && cons->propenabled) == (cons->propconsspos != -1));

   SCIPsetDebugMsg(set, "disable constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* call constraint handler's disabling notification method */
   if( conshdlr->consdisable != NULL )
   {
      SCIP_CALL( conshdlr->consdisable(set->scip, conshdlr, cons) );
   }

   /* delete constraint from the separation array */
   if( cons->separate && cons->sepaenabled )
   {
      conshdlrDelSepacons(conshdlr, cons);
   }

   /* delete constraint from the enforcement array */
   if( cons->enforce )
   {
      conshdlrDelEnfocons(conshdlr, cons);
   }

   /* delete constraint from the propagation array */
   if( cons->propagate && cons->propenabled )
   {
      conshdlrDelPropcons(conshdlr, cons);
   }

   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->propconsspos == -1);

   /* disable constraint */
   cons->enabled = FALSE;
   conshdlr->nenabledconss--;
   stat->nenabledconss--;

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** activates and adds constraint to constraint handler's constraint arrays */
static
SCIP_RETCODE conshdlrActivateCons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_CONS*            cons,               /**< constraint to add */
   int                   depth,              /**< depth in the tree where the activation takes place, or -1 for global problem */
   SCIP_Bool             focusnode           /**< does the constraint activation take place at the focus node? */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(stat != NULL);
   assert(cons != NULL);
   assert(cons->scip == set->scip);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(!cons->active);
   assert(!cons->enabled);
   assert(conshdlr->nactiveconss <= cons->consspos && cons->consspos < conshdlr->nconss);
   assert(conshdlr->conss[cons->consspos] == cons);
   assert(cons->initconsspos < conshdlr->ninitconsskept);
   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->checkconsspos == -1);
   assert(cons->propconsspos == -1);
   assert(depth >= -1);

   SCIPsetDebugMsg(set, "activate constraint <%s> in constraint handler <%s> (depth %d, focus=%u)\n",
      cons->name, conshdlr->name, depth, focusnode);

   /* activate constraint, switch positions with first inactive constraint */
   cons->active = TRUE;
   cons->activedepth = depth;
   conshdlr->conss[cons->consspos] = conshdlr->conss[conshdlr->nactiveconss];
   conshdlr->conss[cons->consspos]->consspos = cons->consspos;
   conshdlr->conss[conshdlr->nactiveconss] = cons;
   cons->consspos = conshdlr->nactiveconss;
   conshdlr->nactiveconss++;
   conshdlr->maxnactiveconss = MAX(conshdlr->maxnactiveconss, conshdlr->nactiveconss);
   stat->nactiveconss++;

   /* add constraint to the check array */
   if( cons->check )
   {
      SCIP_CALL( conshdlrAddCheckcons(conshdlr, set, cons) );
   }

   /* add constraint to the initconss array if the constraint is initial and added to the focus node */
   if( cons->initial )
   {
      SCIP_CALL( conshdlrAddInitcons(conshdlr, set, stat, cons) );
   }

   /* call constraint handler's activation notification method */
   if( conshdlr->consactive != NULL )
   {
      SCIP_CALL( conshdlr->consactive(set->scip, conshdlr, cons) );
   }

   /* enable separation, enforcement, and propagation of constraint */
   SCIP_CALL( conshdlrEnableCons(conshdlr, set, stat, cons) );

   assert(0 <= cons->consspos && cons->consspos < conshdlr->nactiveconss);

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** deactivates and removes constraint from constraint handler's conss array */
static
SCIP_RETCODE conshdlrDeactivateCons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_CONS*            cons                /**< constraint to remove */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(stat != NULL);
   assert(cons != NULL);
   assert(cons->scip == set->scip);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->original);
   assert(cons->active);
   assert(0 <= cons->consspos && cons->consspos < conshdlr->nactiveconss);
   assert(conshdlr->conss[cons->consspos] == cons);
   assert(cons->check == (cons->checkconsspos != -1));

   SCIPsetDebugMsg(set, "deactivate constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* disable constraint */
   if( cons->enabled )
   {
      SCIP_CALL( conshdlrDisableCons(conshdlr, set, stat, cons) );
   }
   assert(!cons->enabled);

   /* call constraint handler's deactivation notification method */
   if( conshdlr->consdeactive != NULL )
   {
      SCIP_CALL( conshdlr->consdeactive(set->scip, conshdlr, cons) );
   }

   /* delete constraint from the initconss array */
   if( cons->initconsspos >= 0 )
   {
      conshdlrDelInitcons(conshdlr, cons);
   }

   /* delete constraint from the check array */
   if( cons->check )
   {
      conshdlrDelCheckcons(conshdlr, cons);
   }

   /* switch constraint with the last active constraint in the conss array */
   conshdlr->conss[cons->consspos] = conshdlr->conss[conshdlr->nactiveconss-1];
   conshdlr->conss[cons->consspos]->consspos = cons->consspos;
   conshdlr->conss[conshdlr->nactiveconss-1] = cons;
   cons->consspos = conshdlr->nactiveconss-1;
   conshdlr->nactiveconss--;
   cons->active = FALSE;
   cons->activedepth = -2;
   stat->nactiveconss--;

   assert(conshdlr->nactiveconss <= cons->consspos && cons->consspos < conshdlr->nconss);
   assert(cons->initconsspos == -1);
   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->checkconsspos == -1);
   assert(cons->propconsspos == -1);

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** processes all delayed updates of constraints:
 *  recently (de)activated constraints will be (de)activated;
 *  recently en/disabled constraints will be en/disabled;
 *  recent obsolete non-check constraints will be globally deleted;
 *  recent obsolete check constraints will be moved to the last positions in the sepa-, enfo-, check-, and prop-arrays;
 *  recent useful constraints will be moved to the first positions in the sepa-, enfo-, check-, and prop-arrays;
 *  constraints which were recently marked to be propagated are moved to the first positions in the prop-array;
 *  no longer used constraints will be freed and removed from the conss array
 */
static
SCIP_RETCODE conshdlrProcessUpdates(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   SCIP_CONS* cons;
   int i;

   assert(conshdlr != NULL);
   assert(!conshdlrAreUpdatesDelayed(conshdlr));
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);

   SCIPsetDebugMsg(set, "processing %d constraints that have to be updated in constraint handler <%s>\n",
      conshdlr->nupdateconss, conshdlr->name);

   for( i = conshdlr->nupdateconss - 1; i >= 0; --i )
   {
      cons = conshdlr->updateconss[i];
      assert(cons != NULL);
      assert(cons->conshdlr == conshdlr);
      assert(cons->update);
      assert(cons->updateinsert || cons->updateactivate || cons->updatedeactivate
         || cons->updateenable || cons->updatedisable
         || cons->updatesepaenable || cons->updatesepadisable
         || cons->updatepropenable || cons->updatepropdisable
         || cons->updateobsolete || cons->updatefree
         || cons->updatemarkpropagate || cons->updateunmarkpropagate);

      SCIPsetDebugMsg(set, " -> constraint <%s>: insert=%u, activate=%u, deactivate=%u, enable=%u, disable=%u, sepaenable=%u, sepadisable=%u, propenable=%u, propdisable=%u, obsolete=%u, free=%u (consdata=%p)\n",
         cons->name, cons->updateinsert, cons->updateactivate, cons->updatedeactivate, 
         cons->updateenable, cons->updatedisable,
         cons->updatesepaenable, cons->updatesepadisable, 
         cons->updatepropenable, cons->updatepropdisable, 
         cons->updateobsolete, cons->updatefree, (void*)cons->consdata);

      if( cons->updateinsert )
      {
         SCIP_CALL( conshdlrAddCons(conshdlr, set, cons) );
         cons->updateinsert = FALSE;
      }

      if( cons->updateactivate )
      {
         assert(!cons->active);
         assert(!cons->updatedeactivate);
         assert(!cons->updateenable);
         assert(!cons->updatedisable);
         assert(!cons->updateobsolete);
         assert(!cons->updatefree);

         /* the activation depth was already stored in SCIPconsActivate() */
         SCIP_CALL( conshdlrActivateCons(conshdlr, set, stat, cons, cons->activedepth, cons->updateactfocus) );
         assert(cons->active);
         cons->updateactivate = FALSE;
      }
      else if( cons->updatedeactivate )
      {
         assert(cons->active);

         SCIP_CALL( conshdlrDeactivateCons(conshdlr, set, stat, cons) );
         assert(!cons->active);
         cons->updatedeactivate = FALSE;
         cons->updateenable = FALSE;
         cons->updatedisable = FALSE;
         cons->obsolete = consExceedsObsoleteage(cons, set);
         cons->updateobsolete = FALSE;
      }
      else if( cons->updateenable )
      {
         assert(!cons->enabled);
         assert(!cons->updatedisable);

         SCIP_CALL( conshdlrEnableCons(conshdlr, set, stat, cons) );
         assert(cons->enabled);
         cons->updateenable = FALSE;
      }
      else if( cons->updatedisable )
      {
         assert(cons->enabled);

         SCIP_CALL( conshdlrDisableCons(conshdlr, set, stat, cons) );
         assert(!cons->enabled);
         cons->updatedisable = FALSE;
      }

      if( cons->updatesepaenable )
      {
         assert(!cons->updatesepadisable);
         if( !cons->sepaenabled )
         {
            SCIP_CALL( conshdlrEnableConsSeparation(conshdlr, set, cons) );
            assert(cons->sepaenabled);
         }
         cons->updatesepaenable = FALSE;
      }
      else if( cons->updatesepadisable )
      {
         if( cons->sepaenabled )
         {         
            SCIP_CALL( conshdlrDisableConsSeparation(conshdlr, cons) );
            assert(!cons->sepaenabled);
         }
         cons->updatesepadisable = FALSE;
      }

      if( cons->updatepropenable )
      {
         assert(!cons->updatepropdisable);
         if( !cons->propenabled )
         {
            SCIP_CALL( conshdlrEnableConsPropagation(conshdlr, set, cons) );
            assert(cons->propenabled);
         }
         cons->updatepropenable = FALSE;
      }
      else if( cons->updatepropdisable )
      {
         if( cons->propenabled )
         {         
            SCIP_CALL( conshdlrDisableConsPropagation(conshdlr, cons) );
            assert(!cons->propenabled);
         }
         cons->updatepropdisable = FALSE;
      }

      if( cons->updatefree )
      {
         /* nothing to do here: the constraint is freed, when it is released from the updateconss array */
         assert(cons->nuses == 1); /* it only exists in the updateconss array */
         cons->updatefree = FALSE;
         cons->updateobsolete = FALSE;
      }
      else
      {
         if( cons->updateobsolete )
         {
            if( !cons->obsolete && consExceedsObsoleteage(cons, set) )
            {
               /* the constraint's status must be switched to obsolete */
               SCIP_CALL( conshdlrMarkConsObsolete(conshdlr, cons) );
            }
            else if( cons->obsolete && !consExceedsObsoleteage(cons, set) )
            {
               /* the constraint's status must be switched to useful */
               SCIP_CALL( conshdlrMarkConsUseful(conshdlr, cons) );
            }
            cons->updateobsolete = FALSE;
         }

         if( cons->updatemarkpropagate )
         {
            /* the constraint must be marked to be propagated */
            conshdlrMarkConsPropagate(conshdlr, cons);
            cons->updatemarkpropagate = FALSE;
         }
         else if( cons->updateunmarkpropagate )
         {
            /* the constraint must be unmarked to be propagated */
            conshdlrUnmarkConsPropagate(conshdlr, cons);
            cons->updateunmarkpropagate = FALSE;
         }
      }

      assert(!cons->updateinsert);
      assert(!cons->updateactivate);
      assert(!cons->updatedeactivate);
      assert(!cons->updateenable);
      assert(!cons->updatedisable);
      assert(!cons->updatesepaenable);
      assert(!cons->updatesepadisable);
      assert(!cons->updatepropenable);
      assert(!cons->updatepropdisable);
      assert(!cons->updateobsolete);
      assert(!cons->updatemarkpropagate);
      assert(!cons->updateunmarkpropagate);
      assert(!cons->updatefree);
      cons->update = FALSE;

      /* release the constraint */
      SCIP_CALL( SCIPconsRelease(&conshdlr->updateconss[i], blkmem, set) );
   }

   conshdlr->nupdateconss = 0;

   return SCIP_OKAY;
}

/** marks constraint handler to delay all constraint updates until the next conshdlrProcessUpdates() call */
static
void conshdlrDelayUpdates(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   SCIPdebugMessage("constraint updates of constraint handler <%s> will be delayed (count:%d)\n",
      conshdlr->name, conshdlr->delayupdatecount+1);

   conshdlr->delayupdatecount++;
}

/** marks constraint handler to perform all constraint updates immediately;
 *  all delayed constraint updates will be processed
 */
static
SCIP_RETCODE conshdlrForceUpdates(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlrAreUpdatesDelayed(conshdlr));

   SCIPsetDebugMsg(set, "constraint updates of constraint handler <%s> will be processed immediately (count:%d)\n",
      conshdlr->name, conshdlr->delayupdatecount);
   conshdlr->delayupdatecount--;

   /* only run the update if all delays are taken away (reference counting) */
   if( !conshdlrAreUpdatesDelayed(conshdlr) )
   {
      SCIP_CALL( conshdlrProcessUpdates(conshdlr, blkmem, set, stat) );
      assert(conshdlr->nupdateconss == 0);
   }

   return SCIP_OKAY;
}

/** adds constraint to constraint handler's update constraint array and captures it */
static
SCIP_RETCODE conshdlrAddUpdateCons(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);

   if( !cons->update )
   {
      SCIPsetDebugMsg(set, "constraint <%s> of age %g has to be updated in constraint handler <%s> (consdata=%p)\n",
         cons->name, cons->age, conshdlr->name, (void*)cons->consdata);

      /* add constraint to the updateconss array */
      SCIP_CALL( conshdlrEnsureUpdateconssMem(conshdlr, set, conshdlr->nupdateconss+1) );
      conshdlr->updateconss[conshdlr->nupdateconss] = cons;
      conshdlr->nupdateconss++;

      /* capture constraint */
      SCIPconsCapture(cons);

      cons->update = TRUE;
   }

   return SCIP_OKAY;
}

/** compares two constraint handlers w. r. to their separation priority */
SCIP_DECL_SORTPTRCOMP(SCIPconshdlrCompSepa)
{  /*lint --e{715}*/
   return ((SCIP_CONSHDLR*)elem2)->sepapriority - ((SCIP_CONSHDLR*)elem1)->sepapriority;
}

/** compares two constraint handlers w. r. to their enforcing priority */
SCIP_DECL_SORTPTRCOMP(SCIPconshdlrCompEnfo)
{  /*lint --e{715}*/
   return ((SCIP_CONSHDLR*)elem2)->enfopriority - ((SCIP_CONSHDLR*)elem1)->enfopriority;
}

/** compares two constraint handlers w. r. to their feasibility check priority */
SCIP_DECL_SORTPTRCOMP(SCIPconshdlrCompCheck)
{  /*lint --e{715}*/
   return ((SCIP_CONSHDLR*)elem2)->checkpriority - ((SCIP_CONSHDLR*)elem1)->checkpriority;
}

/** copies the given constraint handler to a new scip */
SCIP_RETCODE SCIPconshdlrCopyInclude(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< SCIP_SET of SCIP to copy to */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(valid != NULL);
   assert(set->scip != NULL);

   if( conshdlr->conshdlrcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including constraint handler %s in subscip %p\n", SCIPconshdlrGetName(conshdlr), (void*)set->scip);
      SCIP_CALL( conshdlr->conshdlrcopy(set->scip, conshdlr, valid) );
   }

   return SCIP_OKAY;
}

/** creates a constraint handler */
SCIP_RETCODE SCIPconshdlrCreate(
   SCIP_CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of constraint handler */
   const char*           desc,               /**< description of constraint handler */
   int                   sepapriority,       /**< priority of the constraint handler for separation */
   int                   enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int                   checkpriority,      /**< priority of the constraint handler for checking feasibility (and propagation) */
   int                   sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   int                   eagerfreq,          /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
   int                   maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   SCIP_Bool             delaysepa,          /**< should separation method be delayed, if other separators found cuts? */
   SCIP_Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   SCIP_Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   SCIP_PROPTIMING       proptiming,         /**< positions in the node solving loop where propagation method of constraint handlers should be executed */
   SCIP_PRESOLTIMING     presoltiming,       /**< timing mask of the constraint handler's presolving method */
   SCIP_DECL_CONSHDLRCOPY((*conshdlrcopy)),  /**< copy method of constraint handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   SCIP_DECL_CONSINIT    ((*consinit)),      /**< initialize constraint handler */
   SCIP_DECL_CONSEXIT    ((*consexit)),      /**< deinitialize constraint handler */
   SCIP_DECL_CONSINITPRE ((*consinitpre)),   /**< presolving initialization method of constraint handler */
   SCIP_DECL_CONSEXITPRE ((*consexitpre)),   /**< presolving deinitialization method of constraint handler */
   SCIP_DECL_CONSINITSOL ((*consinitsol)),   /**< solving process initialization method of constraint handler */
   SCIP_DECL_CONSEXITSOL ((*consexitsol)),   /**< solving process deinitialization method of constraint handler */
   SCIP_DECL_CONSDELETE  ((*consdelete)),    /**< free specific constraint data */
   SCIP_DECL_CONSTRANS   ((*constrans)),     /**< transform constraint data into data belonging to the transformed problem */
   SCIP_DECL_CONSINITLP  ((*consinitlp)),    /**< initialize LP with relaxations of "initial" constraints */
   SCIP_DECL_CONSSEPALP  ((*conssepalp)),    /**< separate cutting planes for LP solution */
   SCIP_DECL_CONSSEPASOL ((*conssepasol)),   /**< separate cutting planes for arbitrary primal solution */
   SCIP_DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   SCIP_DECL_CONSENFORELAX ((*consenforelax)), /**< enforcing constraints for relaxation solutions */
   SCIP_DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   SCIP_DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   SCIP_DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   SCIP_DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method */
   SCIP_DECL_CONSRESPROP ((*consresprop)),   /**< propagation conflict resolving method */
   SCIP_DECL_CONSLOCK    ((*conslock)),      /**< variable rounding lock method */
   SCIP_DECL_CONSACTIVE  ((*consactive)),    /**< activation notification method */
   SCIP_DECL_CONSDEACTIVE((*consdeactive)),  /**< deactivation notification method */
   SCIP_DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   SCIP_DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
   SCIP_DECL_CONSDELVARS ((*consdelvars)),   /**< variable deletion method */
   SCIP_DECL_CONSPRINT   ((*consprint)),     /**< constraint display method */
   SCIP_DECL_CONSCOPY    ((*conscopy)),      /**< constraint copying method */
   SCIP_DECL_CONSPARSE   ((*consparse)),     /**< constraint parsing method */
   SCIP_DECL_CONSGETVARS ((*consgetvars)),   /**< constraint get variables method */
   SCIP_DECL_CONSGETNVARS((*consgetnvars)),  /**< constraint get number of variable method */
   SCIP_DECL_CONSGETDIVEBDCHGS((*consgetdivebdchgs)), /**< constraint handler diving solution enforcement method */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(conshdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(conssepalp != NULL || conssepasol != NULL || sepafreq == -1);
   assert(consprop != NULL || propfreq == -1);
   assert(eagerfreq >= -1);
   assert(!needscons || ((conshdlrcopy == NULL) == (conscopy == NULL)));

   /* the interface change from delay flags to timings cannot be recognized at compile time: Exit with an appropriate
    * error message
    */
   if( presoltiming < SCIP_PRESOLTIMING_NONE || presoltiming > SCIP_PRESOLTIMING_MAX )
   {
      SCIPmessagePrintError("ERROR: 'PRESOLDELAY'-flag no longer available since SCIP 3.2, use an appropriate "
         "'SCIP_PRESOLTIMING' for <%s> constraint handler instead.\n", name);

      return SCIP_PARAMETERWRONGVAL;
   }

   /* both callbacks have to exist or not exist */
   assert((consgetvars != NULL) == (consgetnvars != NULL));

   SCIP_ALLOC( BMSallocMemory(conshdlr) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*conshdlr)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*conshdlr)->desc, desc, strlen(desc)+1) );
   (*conshdlr)->sepapriority = sepapriority;
   (*conshdlr)->enfopriority = enfopriority;
   (*conshdlr)->checkpriority = checkpriority;
   (*conshdlr)->sepafreq = sepafreq;
   (*conshdlr)->propfreq = propfreq;
   (*conshdlr)->eagerfreq = eagerfreq;
   (*conshdlr)->maxprerounds = maxprerounds;
   (*conshdlr)->conshdlrcopy = conshdlrcopy;
   (*conshdlr)->consfree = consfree;
   (*conshdlr)->consinit = consinit;
   (*conshdlr)->consexit = consexit;
   (*conshdlr)->consinitpre = consinitpre;
   (*conshdlr)->consexitpre = consexitpre;
   (*conshdlr)->consinitsol = consinitsol;
   (*conshdlr)->consexitsol = consexitsol;
   (*conshdlr)->consdelete = consdelete;
   (*conshdlr)->constrans = constrans;
   (*conshdlr)->consinitlp = consinitlp;
   (*conshdlr)->conssepalp = conssepalp;
   (*conshdlr)->conssepasol = conssepasol;
   (*conshdlr)->consenfolp = consenfolp;
   (*conshdlr)->consenforelax = consenforelax;
   (*conshdlr)->consenfops = consenfops;
   (*conshdlr)->conscheck = conscheck;
   (*conshdlr)->consprop = consprop;
   (*conshdlr)->conspresol = conspresol;
   (*conshdlr)->consresprop = consresprop;
   (*conshdlr)->conslock = conslock;
   (*conshdlr)->consactive = consactive;
   (*conshdlr)->consdeactive = consdeactive;
   (*conshdlr)->consenable = consenable;
   (*conshdlr)->consdisable = consdisable;
   (*conshdlr)->consprint = consprint;
   (*conshdlr)->consdelvars = consdelvars;
   (*conshdlr)->conscopy = conscopy;
   (*conshdlr)->consparse = consparse;
   (*conshdlr)->consgetvars = consgetvars;
   (*conshdlr)->consgetnvars = consgetnvars;
   (*conshdlr)->conshdlrdata = conshdlrdata;
   (*conshdlr)->consgetdivebdchgs = NULL;
   (*conshdlr)->conss = NULL;
   (*conshdlr)->consssize = 0;
   (*conshdlr)->nconss = 0;
   (*conshdlr)->nactiveconss = 0;
   (*conshdlr)->maxnactiveconss = 0;
   (*conshdlr)->startnactiveconss = 0;
   (*conshdlr)->initconss = NULL;
   (*conshdlr)->initconsssize = 0;
   (*conshdlr)->ninitconss = 0;
   (*conshdlr)->ninitconsskept = 0;
   (*conshdlr)->sepaconss = NULL;
   (*conshdlr)->sepaconsssize = 0;
   (*conshdlr)->nsepaconss = 0;
   (*conshdlr)->nusefulsepaconss = 0;
   (*conshdlr)->enfoconss = NULL;
   (*conshdlr)->enfoconsssize = 0;
   (*conshdlr)->nenfoconss = 0;
   (*conshdlr)->nusefulenfoconss = 0;
   (*conshdlr)->checkconss = NULL;
   (*conshdlr)->checkconsssize = 0;
   (*conshdlr)->ncheckconss = 0;
   (*conshdlr)->nusefulcheckconss = 0;
   (*conshdlr)->propconss = NULL;
   (*conshdlr)->propconsssize = 0;
   (*conshdlr)->npropconss = 0;
   (*conshdlr)->nusefulpropconss = 0;
   (*conshdlr)->nmarkedpropconss = 0;
   (*conshdlr)->updateconss = NULL;
   (*conshdlr)->updateconsssize = 0;
   (*conshdlr)->nupdateconss = 0;
   (*conshdlr)->nenabledconss = 0;
   (*conshdlr)->lastnusefulpropconss = 0;
   (*conshdlr)->lastnusefulsepaconss = 0;
   (*conshdlr)->lastnusefulenfoconss = 0;

   (*conshdlr)->storedpropconss = NULL;
   (*conshdlr)->storedpropconsssize = 0;
   (*conshdlr)->storednmarkedpropconss = 0;
   (*conshdlr)->storedpropdomchgcount = 0;

   SCIP_CALL( SCIPclockCreate(&(*conshdlr)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conshdlr)->presoltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conshdlr)->sepatime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conshdlr)->enfolptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conshdlr)->enfopstime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conshdlr)->enforelaxtime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conshdlr)->proptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conshdlr)->sbproptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conshdlr)->checktime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conshdlr)->resproptime, SCIP_CLOCKTYPE_DEFAULT) );

   (*conshdlr)->nsepacalls = 0;
   (*conshdlr)->nenfolpcalls = 0;
   (*conshdlr)->nenfopscalls = 0;
   (*conshdlr)->nenforelaxcalls = 0;
   (*conshdlr)->npropcalls = 0;
   (*conshdlr)->ncheckcalls = 0;
   (*conshdlr)->nrespropcalls = 0;
   (*conshdlr)->ncutoffs = 0;
   (*conshdlr)->ncutsfound = 0;
   (*conshdlr)->ncutsapplied = 0;
   (*conshdlr)->nconssfound = 0;
   (*conshdlr)->ndomredsfound = 0;
   (*conshdlr)->nchildren = 0;
   (*conshdlr)->lastpropdomchgcount = -1;
   (*conshdlr)->lastsepalpcount = -1;
   (*conshdlr)->lastenfolplpcount = -1;
   (*conshdlr)->lastenfolpdomchgcount = -1;
   (*conshdlr)->lastenfopsdomchgcount = -1;
   (*conshdlr)->lastenforelaxdomchgcount = -1;
   (*conshdlr)->lastenforelaxrelaxcount = -1;
   (*conshdlr)->lastenfolpnode = -1;
   (*conshdlr)->lastenfopsnode = -1;
   (*conshdlr)->lastenfolpresult = SCIP_DIDNOTRUN;
   (*conshdlr)->lastenfopsresult = SCIP_DIDNOTRUN;
   (*conshdlr)->lastenforelaxresult = SCIP_DIDNOTRUN;
   (*conshdlr)->lastnfixedvars = 0;
   (*conshdlr)->lastnaggrvars = 0;
   (*conshdlr)->lastnchgvartypes = 0;
   (*conshdlr)->lastnchgbds = 0;
   (*conshdlr)->lastnaddholes = 0;
   (*conshdlr)->lastndelconss = 0;
   (*conshdlr)->lastnaddconss = 0;
   (*conshdlr)->lastnupgdconss = 0;
   (*conshdlr)->lastnchgcoefs = 0;
   (*conshdlr)->lastnchgsides = 0;
   (*conshdlr)->nfixedvars = 0;
   (*conshdlr)->naggrvars = 0;
   (*conshdlr)->nchgvartypes = 0;
   (*conshdlr)->nchgbds = 0;
   (*conshdlr)->naddholes = 0;
   (*conshdlr)->ndelconss = 0;
   (*conshdlr)->naddconss = 0;
   (*conshdlr)->nupgdconss = 0;
   (*conshdlr)->nchgcoefs = 0;
   (*conshdlr)->nchgsides = 0;
   (*conshdlr)->npresolcalls = 0;
   (*conshdlr)->delayupdatecount = 0;
   (*conshdlr)->ageresetavg = AGERESETAVG_INIT;
   (*conshdlr)->needscons = needscons;
   (*conshdlr)->sepalpwasdelayed = FALSE;
   (*conshdlr)->sepasolwasdelayed = FALSE;
   (*conshdlr)->propwasdelayed = FALSE;
   (*conshdlr)->duringsepa = FALSE;
   (*conshdlr)->duringprop = FALSE;
   (*conshdlr)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/sepafreq", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "frequency for separating cuts (-1: never, 0: only in root node)",
         &(*conshdlr)->sepafreq, FALSE, sepafreq, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/propfreq", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "frequency for propagating domains (-1: never, 0: only in root node)",
         &(*conshdlr)->propfreq, FALSE, propfreq, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/proptiming", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "timing when constraint propagation should be called (%u:BEFORELP, %u:DURINGLPLOOP, %u:AFTERLPLOOP, %u:ALWAYS)", SCIP_PROPTIMING_BEFORELP, SCIP_PROPTIMING_DURINGLPLOOP, SCIP_PROPTIMING_AFTERLPLOOP, SCIP_PROPTIMING_ALWAYS);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         (int*)(&(*conshdlr)->proptiming), TRUE, proptiming, (int) SCIP_PROPTIMING_BEFORELP, (int) SCIP_PROPTIMING_ALWAYS, NULL, NULL) ); /*lint !e713*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/eagerfreq", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "frequency for using all instead of only the useful constraints in separation, propagation and enforcement (-1: never, 0: only in first evaluation)",
         &(*conshdlr)->eagerfreq, TRUE, eagerfreq, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/maxprerounds", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "maximal number of presolving rounds the constraint handler participates in (-1: no limit)",
         &(*conshdlr)->maxprerounds, TRUE, maxprerounds, -1, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/delaysepa", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should separation method be delayed, if other separators found cuts?",
         &(*conshdlr)->delaysepa, TRUE, delaysepa, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/delayprop", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should propagation method be delayed, if other propagators found reductions?",
         &(*conshdlr)->delayprop, TRUE, delayprop, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/presoltiming", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "timing mask of the constraint handler's presolving method (%u:FAST, %u:MEDIUM, %u:EXHAUSTIVE, %u:FINAL)",
      SCIP_PRESOLTIMING_FAST, SCIP_PRESOLTIMING_MEDIUM, SCIP_PRESOLTIMING_EXHAUSTIVE, SCIP_PRESOLTIMING_FINAL);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         (int*)&(*conshdlr)->presoltiming, TRUE, presoltiming, (int) SCIP_PRESOLTIMING_FAST, (int) SCIP_PRESOLTIMING_MAX, NULL, NULL) ); /*lint !e740 !e713*/

   return SCIP_OKAY;
} /*lint !e715*/

/** calls destructor and frees memory of constraint handler */
SCIP_RETCODE SCIPconshdlrFree(
   SCIP_CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conshdlr != NULL);
   assert(*conshdlr != NULL);
   assert(!(*conshdlr)->initialized);
   assert((*conshdlr)->nconss == 0);
   assert(set != NULL);

   /* call destructor of constraint handler */
   if( (*conshdlr)->consfree != NULL )
   {
      SCIP_CALL( (*conshdlr)->consfree(set->scip, *conshdlr) );
   }

   SCIPclockFree(&(*conshdlr)->resproptime);
   SCIPclockFree(&(*conshdlr)->checktime);
   SCIPclockFree(&(*conshdlr)->sbproptime);
   SCIPclockFree(&(*conshdlr)->proptime);
   SCIPclockFree(&(*conshdlr)->enforelaxtime);
   SCIPclockFree(&(*conshdlr)->enfopstime);
   SCIPclockFree(&(*conshdlr)->enfolptime);
   SCIPclockFree(&(*conshdlr)->sepatime);
   SCIPclockFree(&(*conshdlr)->presoltime);
   SCIPclockFree(&(*conshdlr)->setuptime);

   BMSfreeMemoryArray(&(*conshdlr)->name);
   BMSfreeMemoryArray(&(*conshdlr)->desc);
   BMSfreeMemoryArrayNull(&(*conshdlr)->conss);
   BMSfreeMemoryArrayNull(&(*conshdlr)->initconss);
   BMSfreeMemoryArrayNull(&(*conshdlr)->sepaconss);
   BMSfreeMemoryArrayNull(&(*conshdlr)->enfoconss);
   BMSfreeMemoryArrayNull(&(*conshdlr)->checkconss);
   BMSfreeMemoryArrayNull(&(*conshdlr)->propconss);
   BMSfreeMemoryArrayNull(&(*conshdlr)->updateconss);
   BMSfreeMemoryArrayNull(&(*conshdlr)->storedpropconss);
   BMSfreeMemory(conshdlr);

   return SCIP_OKAY;
}

/** calls initialization method of constraint handler */
SCIP_RETCODE SCIPconshdlrInit(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( conshdlr->initialized )
   {
      SCIPerrorMessage("constraint handler <%s> already initialized\n", conshdlr->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(conshdlr->setuptime);
      SCIPclockReset(conshdlr->presoltime);
      SCIPclockReset(conshdlr->sepatime);
      SCIPclockReset(conshdlr->enfolptime);
      SCIPclockReset(conshdlr->enfopstime);
      SCIPclockReset(conshdlr->enforelaxtime);
      SCIPclockReset(conshdlr->proptime);
      SCIPclockReset(conshdlr->sbproptime);
      SCIPclockReset(conshdlr->checktime);
      SCIPclockReset(conshdlr->resproptime);

      conshdlr->nsepacalls = 0;
      conshdlr->nenfolpcalls = 0;
      conshdlr->nenfopscalls = 0;
      conshdlr->nenforelaxcalls = 0;
      conshdlr->npropcalls = 0;
      conshdlr->ncheckcalls = 0;
      conshdlr->nrespropcalls = 0;
      conshdlr->ncutoffs = 0;
      conshdlr->ncutsfound = 0;
      conshdlr->ncutsapplied = 0;
      conshdlr->nconssfound = 0;
      conshdlr->ndomredsfound = 0;
      conshdlr->nchildren = 0;
      conshdlr->lastpropdomchgcount = -1;
      conshdlr->lastenfolpdomchgcount = -1;
      conshdlr->lastenfopsdomchgcount = -1;
      conshdlr->lastenforelaxdomchgcount = -1;
      conshdlr->lastenforelaxrelaxcount = -1;
      conshdlr->lastenfolpnode = -1;
      conshdlr->lastenfopsnode = -1;
      conshdlr->lastenfolpresult = SCIP_DIDNOTRUN;
      conshdlr->lastenfopsresult = SCIP_DIDNOTRUN;
      conshdlr->maxnactiveconss = conshdlr->nactiveconss;
      conshdlr->startnactiveconss = 0;
      conshdlr->lastsepalpcount = -1;
      conshdlr->lastenfolplpcount = -1;
      conshdlr->lastnusefulpropconss = 0;
      conshdlr->lastnusefulsepaconss = 0;
      conshdlr->lastnusefulenfoconss = 0;
      conshdlr->lastnfixedvars = 0;
      conshdlr->lastnaggrvars = 0;
      conshdlr->lastnchgvartypes = 0;
      conshdlr->lastnchgbds = 0;
      conshdlr->lastnaddholes = 0;
      conshdlr->lastndelconss = 0;
      conshdlr->lastnaddconss = 0;
      conshdlr->lastnupgdconss = 0;
      conshdlr->lastnchgcoefs = 0;
      conshdlr->lastnchgsides = 0;
      conshdlr->nfixedvars = 0;
      conshdlr->naggrvars = 0;
      conshdlr->nchgvartypes = 0;
      conshdlr->nchgbds = 0;
      conshdlr->naddholes = 0;
      conshdlr->ndelconss = 0;
      conshdlr->naddconss = 0;
      conshdlr->nupgdconss = 0;
      conshdlr->nchgcoefs = 0;
      conshdlr->nchgsides = 0;
      conshdlr->npresolcalls = 0;
      conshdlr->ageresetavg = AGERESETAVG_INIT;
      conshdlr->sepalpwasdelayed = FALSE;
      conshdlr->sepasolwasdelayed = FALSE;
      conshdlr->propwasdelayed = FALSE;
   }

   /* call initialization method of constraint handler */
   if( conshdlr->consinit != NULL )
   {
      /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* start timing */
      SCIPclockStart(conshdlr->setuptime, set);

      /* call external method */
      SCIP_CALL( conshdlr->consinit(set->scip, conshdlr, conshdlr->conss, conshdlr->nconss) );

      /* stop timing */
      SCIPclockStop(conshdlr->setuptime, set);

      /* perform the cached constraint updates */
      SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );
   }
   conshdlr->initialized = TRUE;
   assert(!conshdlrAreUpdatesDelayed(conshdlr));

   return SCIP_OKAY;
}

/** calls exit method of constraint handler */
SCIP_RETCODE SCIPconshdlrExit(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( !conshdlr->initialized )
   {
      SCIPerrorMessage("constraint handler <%s> not initialized\n", conshdlr->name);
      return SCIP_INVALIDCALL;
   }

   /* call deinitialization method of constraint handler */
   if( conshdlr->consexit != NULL )
   {
      /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* start timing */
      SCIPclockStart(conshdlr->setuptime, set);

      /* call external method */
      SCIP_CALL( conshdlr->consexit(set->scip, conshdlr, conshdlr->conss, conshdlr->nconss) );

      /* stop timing */
      SCIPclockStop(conshdlr->setuptime, set);

      /* perform the cached constraint updates */
      SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );
   }
   conshdlr->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs constraint handler that the presolving process is being started */
SCIP_RETCODE SCIPconshdlrInitpre(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   /* reset conshdlr last presolved data in case of a restart */
   conshdlr->lastpropdomchgcount = -1;
   conshdlr->lastenfolpdomchgcount = -1;
   conshdlr->lastenfopsdomchgcount = -1;
   conshdlr->lastenforelaxdomchgcount = -1;
   conshdlr->lastenforelaxrelaxcount = -1;
   conshdlr->lastenfolpnode = -1;
   conshdlr->lastenfopsnode = -1;
   conshdlr->lastenfolpresult = SCIP_DIDNOTRUN;
   conshdlr->lastenfopsresult = SCIP_DIDNOTRUN;
   conshdlr->lastenforelaxresult = SCIP_DIDNOTRUN;
   conshdlr->maxnactiveconss = conshdlr->nactiveconss;
   conshdlr->startnactiveconss = 0;
   conshdlr->lastsepalpcount = -1;
   conshdlr->lastenfolplpcount = -1;
   conshdlr->lastnusefulpropconss = 0;
   conshdlr->lastnusefulsepaconss = 0;
   conshdlr->lastnusefulenfoconss = 0;
   conshdlr->lastnfixedvars = 0;
   conshdlr->lastnaggrvars = 0;
   conshdlr->lastnchgvartypes = 0;
   conshdlr->lastnchgbds = 0;
   conshdlr->lastnaddholes = 0;
   conshdlr->lastndelconss = 0;
   conshdlr->lastnaddconss = 0;
   conshdlr->lastnupgdconss = 0;
   conshdlr->lastnchgcoefs = 0;
   conshdlr->lastnchgsides = 0;
   conshdlr->propwasdelayed = FALSE;

   /* call presolving initialization method of constraint handler */
   if( conshdlr->consinitpre != NULL )
   {
      /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* start timing */
      SCIPclockStart(conshdlr->setuptime, set);

      /* call external method */
      SCIP_CALL( conshdlr->consinitpre(set->scip, conshdlr, conshdlr->conss, conshdlr->nconss) );

      /* stop timing */
      SCIPclockStop(conshdlr->setuptime, set);

      /* perform the cached constraint updates */
      SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );
   }

   /* after a restart the LP is empty but the initial constraints are not included in the initialconss array anymore;
    * we have to put them back into this array in order to obtain the correct initial root relaxation
    */
   if( stat->nruns >= 2 )
   {
      int c;

      for( c = 0; c < conshdlr->nconss; ++c )
      {

         /**@todo should only active constraints be added to the initconss array? at least cons->active is asserted in
          *       conshdlrAddInitcons(conshdlr, set, conshdlr->conss[c])
          */
         if( conshdlr->conss[c]->addarraypos >= 0 && !conshdlr->conss[c]->deleted &&
            conshdlr->conss[c]->initial && conshdlr->conss[c]->initconsspos == -1 )
         {
            SCIP_CALL( conshdlrAddInitcons(conshdlr, set, stat, conshdlr->conss[c]) );
         }
      }
   }

#if 0
   /* check if all initial constraints are included in the initconss array */
   {
      int c;

      for( c = 0; c < conshdlr->nconss; ++c )
      {
         assert(conshdlr->conss[c]->deleted || conshdlr->conss[c]->initial == (conshdlr->conss[c]->initconsspos >= 0));
         assert(conshdlr->conss[c]->deleted || !conshdlr->conss[c]->initial
            || conshdlr->initconss[conshdlr->conss[c]->initconsspos] == conshdlr->conss[c]);
      }
   }
#endif

   return SCIP_OKAY;
}

/** informs constraint handler that the presolving is finished */
SCIP_RETCODE SCIPconshdlrExitpre(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   /* call presolving deinitialization method of constraint handler */
   if( conshdlr->consexitpre != NULL )
   {
      /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* start timing */
      SCIPclockStart(conshdlr->setuptime, set);

      /* call external method */
      SCIP_CALL( conshdlr->consexitpre(set->scip, conshdlr, conshdlr->conss, conshdlr->nconss) );

      /* stop timing */
      SCIPclockStop(conshdlr->setuptime, set);

      /* perform the cached constraint updates */
      SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );
   }

   /* update statistics */
   conshdlr->maxnactiveconss = conshdlr->nactiveconss;
   conshdlr->startnactiveconss = conshdlr->nactiveconss;

   return SCIP_OKAY;
}

/** informs constraint handler that the branch and bound process is being started */
SCIP_RETCODE SCIPconshdlrInitsol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(stat != NULL);

   conshdlr->sepalpwasdelayed = FALSE;
   conshdlr->sepasolwasdelayed = FALSE;

   /* call solving process initialization method of constraint handler */
   if( conshdlr->consinitsol != NULL )
   {
      /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* start timing */
      SCIPclockStart(conshdlr->setuptime, set);

      /* call external method */
      SCIP_CALL( conshdlr->consinitsol(set->scip, conshdlr, conshdlr->conss, conshdlr->nconss) );

      /* stop timing */
      SCIPclockStop(conshdlr->setuptime, set);

      /* perform the cached constraint updates */
      SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );
   }

   return SCIP_OKAY;
}

/** informs constraint handler that the branch and bound process data is being freed */
SCIP_RETCODE SCIPconshdlrExitsol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_Bool             restart             /**< was this exit solve call triggered by a restart? */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of constraint handler */
   if( conshdlr->consexitsol != NULL )
   {
      /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* start timing */
      SCIPclockStart(conshdlr->setuptime, set);

      /* call external method */
      SCIP_CALL( conshdlr->consexitsol(set->scip, conshdlr, conshdlr->conss, conshdlr->nconss, restart) );

      /* stop timing */
      SCIPclockStop(conshdlr->setuptime, set);

      /* perform the cached constraint updates */
      SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );
   }

   return SCIP_OKAY;
}

/** calls LP initialization method of constraint handler to separate all initial active constraints */
SCIP_RETCODE SCIPconshdlrInitLP(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             initkeptconss,      /**< Also initialize constraints which are valid at a more global node,
                                              *   but were not activated there? Should be FALSE for repeated calls at
                                              *   one node or if the current focusnode is a child of the former one */
   SCIP_Bool*            cutoff              /**< pointer to store whether infeasibility was detected while building the LP */
   )
{
   assert(conshdlr != NULL);
   assert(cutoff != NULL);
#ifdef MORE_DEBUG
   assert(stat->nnodes > 1 || conshdlr->ninitconsskept == 0 || SCIPtreeProbing(tree));
#endif

   *cutoff = FALSE;

   if( conshdlr->consinitlp != NULL )
   {
      int currentdepth;
      int oldninitconss;
      int c;

      SCIPsetDebugMsg(set, "initializing LP with %d initial constraints of handler <%s> (ninitconss=%d, kept=%d, initkept=%u)\n",
         initkeptconss ? conshdlr->ninitconss : conshdlr->ninitconss - conshdlr->ninitconsskept, conshdlr->name,
         conshdlr->ninitconss, conshdlr->ninitconsskept, initkeptconss);

      /* no constraints to initialize (or only kept constraints which do not need to be initialized this time) -> return */
      if( conshdlr->needscons && (conshdlr->ninitconss == 0 || (!initkeptconss && conshdlr->ninitconss == conshdlr->ninitconsskept)) )
         return SCIP_OKAY;

      /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      oldninitconss = conshdlr->ninitconss;

      /* start timing */
      SCIPclockStart(conshdlr->sepatime, set);

      if( initkeptconss )
      {
         /* add all kept initial constraints which are currently active to the second part of the initconss array */
         /* @todo keep track of where a constraint was already initialized (e.g., in the conssetchg)? */
         for( c = 0; c < conshdlr->ninitconsskept; ++c )
         {
            assert(conshdlr->initconss[c]->initconsspos == c);

            if( SCIPconsIsActive(conshdlr->initconss[c]) )
            {
               SCIP_CALL( conshdlrAddInitcons(conshdlr, set, stat, conshdlr->initconss[c]) );
            }
         }
      }

      /* call external method */
      SCIP_CALL( conshdlr->consinitlp(set->scip, conshdlr, &conshdlr->initconss[conshdlr->ninitconsskept],
            conshdlr->ninitconss - conshdlr->ninitconsskept, cutoff) );

      /* stop timing */
      SCIPclockStop(conshdlr->sepatime, set);

      /* perform the cached constraint updates */
      SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );

      currentdepth = SCIPtreeGetCurrentDepth(tree);
      assert(currentdepth >= 0);

      /* clear the initconss array */
      for( c = conshdlr->ninitconsskept; c < oldninitconss; ++c )
      {
         assert(SCIPconsGetActiveDepth(conshdlr->initconss[c]) >= -1);
         assert(SCIPconsGetActiveDepth(conshdlr->initconss[c]) <= currentdepth);

         /* if the constraint was not initialized at its valid node, we keep it */
         if( currentdepth > 0 ? SCIPconsGetActiveDepth(conshdlr->initconss[c]) != currentdepth :
            SCIPconsGetActiveDepth(conshdlr->initconss[c]) > 0 )
         {
            conshdlr->initconss[conshdlr->ninitconsskept] = conshdlr->initconss[c];
            conshdlr->initconss[conshdlr->ninitconsskept]->initconsspos = conshdlr->ninitconsskept;
            ++(conshdlr->ninitconsskept);
         }
         else
            conshdlr->initconss[c]->initconsspos = -1;
      }
#ifndef NDEBUG
      for( ; c < conshdlr->ninitconss; ++c )
         assert(conshdlr->initconss[c]->initconsspos < conshdlr->ninitconsskept);
#endif
      conshdlr->ninitconss = conshdlr->ninitconsskept;

      if( conshdlr->ninitconss == 0 )
      {
         BMSfreeMemoryArrayNull(&conshdlr->initconss);
         conshdlr->initconsssize = 0;
      }
   }

   return SCIP_OKAY;
}

/** calls separator method of constraint handler to separate LP solution */
SCIP_RETCODE SCIPconshdlrSeparateLP(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             execdelayed,        /**< execute separation method even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(stat != NULL);
   assert(conshdlr->lastsepalpcount != stat->lpcount
      || (0 <= conshdlr->lastnusefulsepaconss && conshdlr->lastnusefulsepaconss <= conshdlr->nusefulsepaconss));
   assert(set != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( conshdlr->conssepalp != NULL
      && ((depth == 0 && conshdlr->sepafreq == 0)
         || (conshdlr->sepafreq > 0 && depth % conshdlr->sepafreq == 0)
         || conshdlr->sepalpwasdelayed) )
   {
      /* check, if separation method should be delayed */
      if( !conshdlr->delaysepa || execdelayed )
      {
         int nconss;
         int nusefulconss;
         int firstcons;

         /* check, if this LP solution was already separated */
         if( conshdlr->lastsepalpcount == stat->lpcount )
         {
            /* all constraints that were not yet separated on the new LP solution must be useful constraints, which means,
             * that the new constraints are the last constraints of the useful ones
             */
            nconss = conshdlr->nusefulsepaconss - conshdlr->lastnusefulsepaconss;
            nusefulconss = nconss;
            firstcons = conshdlr->lastnusefulsepaconss;
         }
         else
         {
            /* on a new LP solution, we want to separate all constraints */
            nconss = conshdlr->nsepaconss;
            nusefulconss = conshdlr->nusefulsepaconss;
            firstcons = 0;
         }
         assert(firstcons >= 0);
         assert(firstcons + nconss <= conshdlr->nsepaconss);
         assert(nusefulconss <= nconss);

         /* constraint handlers without constraints should only be called once */
         if( nconss > 0 || (!conshdlr->needscons && conshdlr->lastsepalpcount != stat->lpcount) )
         {
            SCIP_CONS** conss;
            SCIP_Longint oldndomchgs;
            SCIP_Longint oldnprobdomchgs;
            SCIP_Longint lastsepalpcount;
            int oldncuts;
            int oldnactiveconss;
            int lastnusefulsepaconss;

            SCIPsetDebugMsg(set, "separating constraints %d to %d of %d constraints of handler <%s> (%s LP solution)\n",
               firstcons, firstcons + nconss - 1, conshdlr->nsepaconss, conshdlr->name,
               conshdlr->lastsepalpcount == stat->lpcount ? "old" : "new");

            /* remember the number of processed constraints on the current LP solution */
            lastsepalpcount = stat->lpcount;
            lastnusefulsepaconss = conshdlr->nusefulsepaconss;

            /* get the array of the constraints to be processed */
            conss = &(conshdlr->sepaconss[firstcons]);

            oldndomchgs = stat->nboundchgs + stat->nholechgs;
            oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;
            oldncuts = SCIPsepastoreGetNCuts(sepastore);
            oldnactiveconss = stat->nactiveconss;

            /* check, if we want to use eager evaluation */
            if( (conshdlr->eagerfreq == 0 && conshdlr->nsepacalls == 0)
               || (conshdlr->eagerfreq > 0 && conshdlr->nsepacalls % conshdlr->eagerfreq == 0) )
               nusefulconss = nconss;

            /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
             * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
             * external method; to avoid this, these changes will be buffered and processed after the method call
             */
            conshdlrDelayUpdates(conshdlr);
            conshdlr->duringsepa = TRUE;

            /* start timing */
            SCIPclockStart(conshdlr->sepatime, set);

            /* call external method */
            SCIP_CALL( conshdlr->conssepalp(set->scip, conshdlr, conss, nconss, nusefulconss, result) );
            SCIPsetDebugMsg(set, " -> separating LP returned result <%d>\n", *result);

            /* stop timing */
            SCIPclockStop(conshdlr->sepatime, set);

            /* perform the cached constraint updates */
            conshdlr->duringsepa = FALSE;
            SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );

            /* update statistics */
            if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
            {
               conshdlr->lastsepalpcount = lastsepalpcount;
               conshdlr->lastnusefulsepaconss = MIN(lastnusefulsepaconss, conshdlr->nusefulsepaconss);
               conshdlr->nsepacalls++;
            }
            if( *result == SCIP_CUTOFF )
               conshdlr->ncutoffs++;
            conshdlr->ncutsfound += SCIPsepastoreGetNCuts(sepastore) - oldncuts; /*lint !e776*/
            conshdlr->nconssfound += MAX(stat->nactiveconss - oldnactiveconss, 0); /*lint !e776*/

            /* update domain reductions; therefore remove the domain
             * reduction counts which were generated in probing mode */
            conshdlr->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
            conshdlr->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

            /* evaluate result */
            if( *result != SCIP_CUTOFF
               && *result != SCIP_CONSADDED
               && *result != SCIP_REDUCEDDOM
               && *result != SCIP_SEPARATED
               && *result != SCIP_NEWROUND
               && *result != SCIP_DIDNOTFIND
               && *result != SCIP_DIDNOTRUN
               && *result != SCIP_DELAYED )
            {
               SCIPerrorMessage("LP separation method of constraint handler <%s> returned invalid result <%d>\n", 
                  conshdlr->name, *result);
               return SCIP_INVALIDRESULT;
            }
         }
      }
      else
      {
         SCIPsetDebugMsg(set, "LP separation method of constraint handler <%s> was delayed\n", conshdlr->name);
         *result = SCIP_DELAYED;
      }

      /* remember whether separation method was delayed */
      conshdlr->sepalpwasdelayed = (*result == SCIP_DELAYED);
   }

   return SCIP_OKAY;
}

/** calls separator method of constraint handler to separate given primal solution */
SCIP_RETCODE SCIPconshdlrSeparateSol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             execdelayed,        /**< execute separation method even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(stat != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( conshdlr->conssepasol != NULL
      && ((depth == 0 && conshdlr->sepafreq == 0)
         || (conshdlr->sepafreq > 0 && depth % conshdlr->sepafreq == 0)
         || conshdlr->sepasolwasdelayed) )
   {
      /* check, if separation method should be delayed */
      if( !conshdlr->delaysepa || execdelayed )
      {
         int nconss;
         int nusefulconss;

         /* always separate all constraints */
         nconss = conshdlr->nsepaconss;
         nusefulconss = conshdlr->nusefulsepaconss;
         assert(nusefulconss <= nconss);

         if( nconss > 0 || !conshdlr->needscons )
         {
            SCIP_CONS** conss;
            SCIP_Longint oldndomchgs;
            SCIP_Longint oldnprobdomchgs;
            int oldncuts;
            int oldnactiveconss;

            SCIPsetDebugMsg(set, "separating %d constraints of handler <%s> (primal solution %p)\n",
               nconss, conshdlr->name, (void*)sol);

            /* get the array of the constraints to be processed */
            conss = conshdlr->sepaconss;

            oldndomchgs = stat->nboundchgs + stat->nholechgs;
            oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;
            oldncuts = SCIPsepastoreGetNCuts(sepastore);
            oldnactiveconss = stat->nactiveconss;

            /* check, if we want to use eager evaluation */
            if( (conshdlr->eagerfreq == 0 && conshdlr->nsepacalls == 0)
               || (conshdlr->eagerfreq > 0 && conshdlr->nsepacalls % conshdlr->eagerfreq == 0) )
               nusefulconss = nconss;

            /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
             * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
             * external method; to avoid this, these changes will be buffered and processed after the method call
             */
            conshdlrDelayUpdates(conshdlr);
            conshdlr->duringsepa = TRUE;

            /* start timing */
            SCIPclockStart(conshdlr->sepatime, set);

            /* call external method */
            SCIP_CALL( conshdlr->conssepasol(set->scip, conshdlr, conss, nconss, nusefulconss, sol, result) );
            SCIPsetDebugMsg(set, " -> separating sol returned result <%d>\n", *result);

            /* stop timing */
            SCIPclockStop(conshdlr->sepatime, set);

            /* perform the cached constraint updates */
            conshdlr->duringsepa = FALSE;
            SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );

            /* update statistics */
            if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
               conshdlr->nsepacalls++;
            if( *result == SCIP_CUTOFF )
               conshdlr->ncutoffs++;
            conshdlr->ncutsfound += SCIPsepastoreGetNCuts(sepastore) - oldncuts; /*lint !e776*/
            conshdlr->nconssfound += MAX(stat->nactiveconss - oldnactiveconss, 0); /*lint !e776*/

            /* update domain reductions; therefore remove the domain
             * reduction counts which were generated in probing mode */
            conshdlr->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
            conshdlr->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

            /* evaluate result */
            if( *result != SCIP_CUTOFF
               && *result != SCIP_CONSADDED
               && *result != SCIP_REDUCEDDOM
               && *result != SCIP_SEPARATED
               && *result != SCIP_NEWROUND
               && *result != SCIP_DIDNOTFIND
               && *result != SCIP_DIDNOTRUN
               && *result != SCIP_DELAYED )
            {
               SCIPerrorMessage("SOL separation method of constraint handler <%s> returned invalid result <%d>\n", 
                  conshdlr->name, *result);
               return SCIP_INVALIDRESULT;
            }
         }
      }
      else
      {
         SCIPsetDebugMsg(set, "SOL separation method of constraint handler <%s> was delayed\n", conshdlr->name);
         *result = SCIP_DELAYED;
      }

      /* remember whether separation method was delayed */
      conshdlr->sepasolwasdelayed = (*result == SCIP_DELAYED);
   }

   return SCIP_OKAY;
}

/** calls enforcing method of constraint handler for a relaxation solution for all constraints added after last
 *  conshdlrResetEnfo() call
 */
SCIP_RETCODE SCIPconshdlrEnforceRelaxSol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             relaxsol,           /**< solution to be enforced */
   SCIP_Bool             solinfeasible,      /**< was the solution already found out to be infeasible? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   int nconss;
   int nusefulconss;
   int firstcons;
   SCIP_Bool relaxchanged;
   SCIP_Bool lastinfeasible;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(stat != NULL);
   assert(conshdlr->lastenfopsdomchgcount != stat->domchgcount
      || conshdlr->lastenfopsnode != stat->nnodes
      || (0 <= conshdlr->lastnusefulenfoconss && conshdlr->lastnusefulenfoconss <= conshdlr->nusefulenfoconss));
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(relaxsol != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* check, if this relaxation solution was already enforced at this node
    * the integrality constraint handler always needs to be enforced for all constraints since external branching
    * candidates are cleared before each resolve
    */
   if( conshdlr->lastenforelaxrelaxcount == stat->relaxcount
      && conshdlr->lastenforelaxdomchgcount == stat->domchgcount
      && conshdlr->lastenforelaxnode == stat->nnodes
      && conshdlr->lastenforelaxresult != SCIP_CONSADDED
      && conshdlr->lastenforelaxresult != SCIP_SOLVELP
      && ( strcmp(conshdlr->name, "integral") != 0 )
      )
   {
      assert(conshdlr->lastenforelaxresult != SCIP_CUTOFF);
      assert(conshdlr->lastenforelaxresult != SCIP_BRANCHED);
      assert(conshdlr->lastenforelaxresult != SCIP_REDUCEDDOM);
      assert(conshdlr->lastenforelaxresult != SCIP_DIDNOTRUN);

      /* if we already enforced the same relaxation solution at this node, we will only enforce new constraints in the
       * following; however, the result of the last call for the old constraint is still valid and we have to ensure
       * that an infeasibility in the last call is not lost because we only enforce new constraints
       */
      if( conshdlr->lastenforelaxresult == SCIP_INFEASIBLE )
      {
         *result = SCIP_INFEASIBLE;
         lastinfeasible = TRUE;
      }
      else
         lastinfeasible = FALSE;

      /* all constraints that were not yet enforced on the new relaxation solution must be useful constraints, which means,
       * that the new constraints are the last constraints of the useful ones
       */
      nconss = conshdlr->nusefulenfoconss - conshdlr->lastnusefulenfoconss;
      nusefulconss = nconss;
      firstcons = conshdlr->lastnusefulenfoconss;
      relaxchanged = FALSE;
   }
   else
   {
      /* on a new relaxation solution or a new node, we want to enforce all constraints */
      nconss = conshdlr->nenfoconss;
      nusefulconss = conshdlr->nusefulenfoconss;
      firstcons = 0;
      relaxchanged = TRUE;
      lastinfeasible = FALSE;
   }
   assert(firstcons >= 0);
   assert(firstcons + nconss <= conshdlr->nenfoconss);
   assert(nusefulconss <= nconss);

   /* constraint handlers without constraints should only be called once */
   if( nconss > 0 || (!conshdlr->needscons && relaxchanged) )
   {
      SCIP_CONS** conss;
      SCIP_Longint oldndomchgs;
      SCIP_Longint oldnprobdomchgs;
      int oldncuts;
      int oldnactiveconss;

      assert(conshdlr->consenforelax != NULL);

      SCIPdebugMessage("enforcing constraints %d to %d of %d constraints of handler <%s> (%s relaxation solution)\n",
         firstcons, firstcons + nconss - 1, conshdlr->nenfoconss, conshdlr->name, relaxchanged ? "new" : "old");

      /* remember the number of processed constraints on the current relaxation solution */
      conshdlr->lastenforelaxrelaxcount = stat->relaxcount;
      conshdlr->lastenforelaxdomchgcount = stat->domchgcount;
      conshdlr->lastenforelaxnode = stat->nnodes;
      conshdlr->lastnusefulenfoconss = conshdlr->nusefulenfoconss;

      /* get the array of the constraints to be processed */
      conss = &(conshdlr->enfoconss[firstcons]);

      oldncuts = SCIPsepastoreGetNCuts(sepastore);
      oldnactiveconss = stat->nactiveconss;
      oldndomchgs = stat->nboundchgs + stat->nholechgs;
      oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;

      /* check, if we want to use eager evaluation */
      if( (conshdlr->eagerfreq == 0 && conshdlr->nenforelaxcalls == 0)
         || (conshdlr->eagerfreq > 0 && conshdlr->nenforelaxcalls % conshdlr->eagerfreq == 0) )
         nusefulconss = nconss;

      /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* start timing */
      SCIPclockStart(conshdlr->enforelaxtime, set);

      /* call external method */
      SCIP_CALL( conshdlr->consenforelax(set->scip, relaxsol, conshdlr, conss, nconss, nusefulconss, solinfeasible, result) );
      SCIPdebugMessage(" -> enforcing returned result <%d>\n", *result);

      /* stop timing */
      SCIPclockStop(conshdlr->enforelaxtime, set);

      /* perform the cached constraint updates */
      SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );

      /* update statistics */
      if( *result != SCIP_DIDNOTRUN )
         conshdlr->nenforelaxcalls++;
      else
      {
         SCIPerrorMessage("enforcing method of constraint handler <%s> for relaxation returned an invalid result %d\n",
            conshdlr->name, *result);
         conshdlr->lastenforelaxresult = *result;

         return SCIP_INVALIDRESULT;
      }

      if( *result == SCIP_CUTOFF )
         conshdlr->ncutoffs++;
      conshdlr->ncutsfound += SCIPsepastoreGetNCuts(sepastore) - oldncuts; /*lint !e776*/
      conshdlr->nconssfound += MAX(stat->nactiveconss - oldnactiveconss, 0); /*lint !e776*/

      if( *result != SCIP_BRANCHED )
      {
         assert(tree->nchildren == 0);

         /* update domain reductions; therefore remove the domain
          * reduction counts which were generated in probing mode */
         conshdlr->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
         conshdlr->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);
      }
      else
         conshdlr->nchildren += tree->nchildren;

      /* remember the result of the enforcement call */
      conshdlr->lastenforelaxresult = *result;

      /* evaluate result */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_CONSADDED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_SEPARATED
         && *result != SCIP_BRANCHED
         && *result != SCIP_SOLVELP
         && *result != SCIP_INFEASIBLE
         && *result != SCIP_FEASIBLE
         && *result != SCIP_DIDNOTRUN )
      {
         SCIPerrorMessage("enforcing method of constraint handler <%s> for relaxation solutions returned invalid result <%d>\n",
            conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }

      /* if the same relaxation solution was already enforced at this node, we only enforced new constraints this time;
       * if the enforelax call returns feasible now, the solution is only feasible w.r.t. the new constraints, if the
       * last call detected infeasibility for the old constraints, we have to change the result to infeasible
       */
      if( lastinfeasible && *result == SCIP_FEASIBLE )
         *result = SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}

/** calls enforcing method of constraint handler for LP solution for all constraints added after last
 *  conshdlrResetEnfo() call
 */
SCIP_RETCODE SCIPconshdlrEnforceLPSol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_Bool             solinfeasible,      /**< was the solution already found out to be infeasible? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(stat != NULL);
   assert(conshdlr->lastenfolplpcount != stat->lpcount
      || conshdlr->lastenfolpdomchgcount != stat->domchgcount
      || conshdlr->lastenfolpnode != stat->nnodes
      || (0 <= conshdlr->lastnusefulenfoconss && conshdlr->lastnusefulenfoconss <= conshdlr->nusefulenfoconss));
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( conshdlr->consenfolp != NULL )
   {
      int nconss;
      int nusefulconss;
      int firstcons;
      SCIP_Bool lpchanged;
      SCIP_Bool lastinfeasible;

      /* check, if this LP solution was already enforced at this node */
      if( conshdlr->lastenfolplpcount == stat->lpcount
         && conshdlr->lastenfolpdomchgcount == stat->domchgcount
         && conshdlr->lastenfolpnode == stat->nnodes
         && conshdlr->lastenfolpresult != SCIP_CONSADDED )
      {
         assert(conshdlr->lastenfolpresult == SCIP_FEASIBLE || conshdlr->lastenfolpresult == SCIP_INFEASIBLE
            || conshdlr->lastenfolpresult == SCIP_SEPARATED );

         /* if we already enforced the same pseudo solution at this node, we will only enforce new constraints in the
          * following; however, the result of the last call for the old constraint is still valid and we have to ensure
          * that an infeasibility in the last call is not lost because we only enforce new constraints
          */
         if( conshdlr->lastenfolpresult == SCIP_FEASIBLE )
            lastinfeasible = FALSE;
         else
         {
            assert(conshdlr->lastenfolpresult == SCIP_INFEASIBLE || conshdlr->lastenfolpresult == SCIP_SEPARATED);
            *result = SCIP_INFEASIBLE;
            lastinfeasible = TRUE;
         }

         /* all constraints that were not yet enforced on the new LP solution must be useful constraints, which means,
          * that the new constraints are the last constraints of the useful ones
          */
         nconss = conshdlr->nusefulenfoconss - conshdlr->lastnusefulenfoconss;
         nusefulconss = nconss;
         firstcons = conshdlr->lastnusefulenfoconss;
         lpchanged = FALSE;
      }
      else
      {
         /* on a new LP solution or a new node, we want to enforce all constraints */
         nconss = conshdlr->nenfoconss;
         nusefulconss = conshdlr->nusefulenfoconss;
         firstcons = 0;
         lpchanged = TRUE;
         lastinfeasible = FALSE;
      }
      assert(firstcons >= 0);
      assert(firstcons + nconss <= conshdlr->nenfoconss);
      assert(nusefulconss <= nconss);

      /* constraint handlers without constraints should only be called once */
      if( nconss > 0 || (!conshdlr->needscons && lpchanged) )
      {
         SCIP_CONS** conss;
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;
         int oldncuts;
         int oldnactiveconss;

         SCIPsetDebugMsg(set, "enforcing constraints %d to %d of %d constraints of handler <%s> (%s LP solution)\n",
            firstcons, firstcons + nconss - 1, conshdlr->nenfoconss, conshdlr->name, lpchanged ? "new" : "old");

         /* remember the number of processed constraints on the current LP solution */
         conshdlr->lastenfolplpcount = stat->lpcount;
         conshdlr->lastenfolpdomchgcount = stat->domchgcount;
         conshdlr->lastenfolpnode = stat->nnodes;
         conshdlr->lastnusefulenfoconss = conshdlr->nusefulenfoconss;

         /* get the array of the constraints to be processed */
         conss = &(conshdlr->enfoconss[firstcons]);

         oldncuts = SCIPsepastoreGetNCuts(sepastore);
         oldnactiveconss = stat->nactiveconss;
         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;

         /* check, if we want to use eager evaluation */
         if( (conshdlr->eagerfreq == 0 && conshdlr->nenfolpcalls == 0)
            || (conshdlr->eagerfreq > 0 && conshdlr->nenfolpcalls % conshdlr->eagerfreq == 0) )
            nusefulconss = nconss;

         /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
          * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
          * external method; to avoid this, these changes will be buffered and processed after the method call
          */
         conshdlrDelayUpdates(conshdlr);

         /* start timing */
         SCIPclockStart(conshdlr->enfolptime, set);

         /* call external method */
         SCIP_CALL( conshdlr->consenfolp(set->scip, conshdlr, conss, nconss, nusefulconss, solinfeasible, result) );
         SCIPsetDebugMsg(set, " -> enforcing returned result <%d>\n", *result);

         /* stop timing */
         SCIPclockStop(conshdlr->enfolptime, set);

         /* perform the cached constraint updates */
         SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );

         /* remember the result of the enforcement call */
         conshdlr->lastenfolpresult = *result;

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN )
            conshdlr->nenfolpcalls++;
         if( *result == SCIP_CUTOFF )
            conshdlr->ncutoffs++;
         conshdlr->ncutsfound += SCIPsepastoreGetNCuts(sepastore) - oldncuts; /*lint !e776*/
         conshdlr->nconssfound += MAX(stat->nactiveconss - oldnactiveconss, 0); /*lint !e776*/
         if( *result != SCIP_BRANCHED )
         {
            assert(tree->nchildren == 0);

            /* update domain reductions; therefore remove the domain
             * reduction counts which were generated in probing mode */
            conshdlr->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
            conshdlr->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);
         }
         else
            conshdlr->nchildren += tree->nchildren;

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_BRANCHED
            && *result != SCIP_INFEASIBLE
            && *result != SCIP_FEASIBLE )
         {
            SCIPerrorMessage("enforcing method of constraint handler <%s> for LP solutions returned invalid result <%d>\n", 
               conshdlr->name, *result);
            return SCIP_INVALIDRESULT;
         }

         /* if the same LP solution was already enforced at this node, we only enforced new constraints this time;
          * if the enfolp call returns feasible now, the solution is only feasible w.r.t. the new constraints, if the
          * last call detected infeasibility for the old constraints, we have to change the result to infeasible
          */
         if( lastinfeasible && *result == SCIP_FEASIBLE )
            *result = SCIP_INFEASIBLE;
      }
   }

   return SCIP_OKAY;
}

/** calls diving solution enforcement callback of constraint handler, if it exists */
SCIP_RETCODE SCIPconshdlrGetDiveBoundChanges(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DIVESET*         diveset,            /**< diving settings to control scoring */
   SCIP_SOL*             sol,                /**< current solution of diving mode */
   SCIP_Bool*            success,            /**< pointer to store whether constraint handler successfully found a variable */
   SCIP_Bool*            infeasible          /**< pointer to store whether the current node was detected to be infeasible */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(diveset != NULL);
   assert(sol != NULL);
   assert(success != NULL);
   assert(infeasible != NULL);

   if( conshdlr->consgetdivebdchgs != NULL )
   {
      SCIP_CALL( conshdlr->consgetdivebdchgs(set->scip, conshdlr, diveset, sol, success, infeasible) );
   }

   return SCIP_OKAY;
}

/** calls enforcing method of constraint handler for pseudo solution for all constraints added after last
 *  conshdlrResetEnfo() call
 */
SCIP_RETCODE SCIPconshdlrEnforcePseudoSol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_Bool             solinfeasible,      /**< was the solution already found out to be infeasible? */
   SCIP_Bool             objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
   SCIP_Bool             forced,             /**< should enforcement of pseudo solution be forced? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(stat != NULL);
   assert(conshdlr->lastenfopsdomchgcount != stat->domchgcount
      || conshdlr->lastenfopsnode != stat->nnodes
      || (0 <= conshdlr->lastnusefulenfoconss && conshdlr->lastnusefulenfoconss <= conshdlr->nusefulenfoconss));
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(result != NULL);

   /* no enforcing of pseudo solution */ 
   if( set->cons_disableenfops && SCIPbranchcandGetNPseudoCands(branchcand) > 0 ) 
   {
      *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_FEASIBLE;
   if( conshdlr->consenfops != NULL )
   {
      int nconss;
      int nusefulconss;
      int firstcons;
      SCIP_Bool pschanged;
      SCIP_Bool lastinfeasible;

      /* check, if this pseudo solution was already enforced at this node */
      if( !forced && conshdlr->lastenfopsdomchgcount == stat->domchgcount
         && conshdlr->lastenfopsnode == stat->nnodes
         && conshdlr->lastenfopsresult != SCIP_CONSADDED
         && conshdlr->lastenfopsresult != SCIP_SOLVELP
         )
      {
         assert(conshdlr->lastenfopsresult != SCIP_CUTOFF);
         assert(conshdlr->lastenfopsresult != SCIP_BRANCHED);
         assert(conshdlr->lastenfopsresult != SCIP_REDUCEDDOM);
         assert(conshdlr->lastenfopsresult != SCIP_DIDNOTRUN || objinfeasible);

         /* if we already enforced the same pseudo solution at this node, we will only enforce new constraints in the
          * following; however, the result of the last call for the old constraint is still valid and we have to ensure
          * that an infeasibility in the last call is not lost because we only enforce new constraints
          */
         if( conshdlr->lastenfopsresult == SCIP_INFEASIBLE )
         {
            *result = SCIP_INFEASIBLE;
            lastinfeasible = TRUE;
         }
         else
            lastinfeasible = FALSE;

         /* all constraints that were not yet enforced on the new LP solution must be useful constraints, which means,
          * that the new constraints are the last constraints of the useful ones
          */
         nconss = conshdlr->nusefulenfoconss - conshdlr->lastnusefulenfoconss;
         nusefulconss = nconss;
         firstcons = conshdlr->lastnusefulenfoconss;
         pschanged = FALSE;
      }
      else
      {
         /* on a new pseudo solution or a new node, we want to enforce all constraints */
         nconss = conshdlr->nenfoconss;
         nusefulconss = conshdlr->nusefulenfoconss;
         firstcons = 0;
         pschanged = TRUE;
         lastinfeasible = FALSE;
      }
      assert(firstcons >= 0);
      assert(firstcons + nconss <= conshdlr->nenfoconss);
      assert(nusefulconss <= nconss);

      /* constraint handlers without constraints should only be called once */
      if( nconss > 0 || (!conshdlr->needscons && pschanged) ) 
      {
         SCIP_CONS** conss;
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;

         SCIPsetDebugMsg(set, "enforcing constraints %d to %d of %d constraints of handler <%s> (%s pseudo solution, objinfeasible=%u)\n",
            firstcons, firstcons + nconss - 1, conshdlr->nenfoconss, conshdlr->name, pschanged ? "new" : "old", objinfeasible);

         /* remember the number of processed constraints on the current pseudo solution */
         conshdlr->lastenfopsdomchgcount = stat->domchgcount;
         conshdlr->lastenfopsnode = stat->nnodes;
         conshdlr->lastnusefulenfoconss = conshdlr->nusefulenfoconss;

         /* get the array of the constraints to be processed */
         conss = &(conshdlr->enfoconss[firstcons]);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;

         /* check, if we want to use eager evaluation */
         if( (conshdlr->eagerfreq == 0 && conshdlr->nenfopscalls == 0)
            || (conshdlr->eagerfreq > 0 && conshdlr->nenfopscalls % conshdlr->eagerfreq == 0) )
            nusefulconss = nconss;

         /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
          * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
          * external method; to avoid this, these changes will be buffered and processed after the method call
          */
         conshdlrDelayUpdates(conshdlr);

         /* start timing */
         SCIPclockStart(conshdlr->enfopstime, set);

         /* call external method */
         SCIP_CALL( conshdlr->consenfops(set->scip, conshdlr, conss, nconss, nusefulconss, solinfeasible, objinfeasible, result) );
         SCIPsetDebugMsg(set, " -> enforcing returned result <%d>\n", *result);

         /* stop timing */
         SCIPclockStop(conshdlr->enfopstime, set);

         /* perform the cached constraint updates */
         SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN )
            conshdlr->nenfopscalls++;
         else if( !objinfeasible )
         {
            SCIPerrorMessage("enforcing method of constraint handler <%s> for pseudo solutions was skipped, even though the solution was not objective-infeasible\n", 
               conshdlr->name);
            conshdlr->lastenfopsresult = *result;

            return SCIP_INVALIDRESULT;
         }
         /* A constraint handler might return SCIP_DIDNOTRUN and not check any constraints in case objinfeasible was
          * TRUE; we change the result pointer to SCIP_INFEASIBLE in this case.
          */
         else
            *result = SCIP_INFEASIBLE;

         if( *result == SCIP_CUTOFF )
            conshdlr->ncutoffs++;

         if( *result != SCIP_BRANCHED )
         {
            assert(tree->nchildren == 0);

            /* update domain reductions; therefore remove the domain
             * reduction counts which were generated in probing mode */
            conshdlr->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
            conshdlr->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);
         }
         else
            conshdlr->nchildren += tree->nchildren;

         /* remember the result of the enforcement call */
         conshdlr->lastenfopsresult = *result;

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_BRANCHED
            && *result != SCIP_SOLVELP
            && *result != SCIP_INFEASIBLE
            && *result != SCIP_FEASIBLE
            && *result != SCIP_DIDNOTRUN )
         {
            SCIPerrorMessage("enforcing method of constraint handler <%s> for pseudo solutions returned invalid result <%d>\n", 
               conshdlr->name, *result);
            return SCIP_INVALIDRESULT;
         }

         /* if the same pseudo solution was already enforced at this node, we only enforced new constraints this time;
          * if the enfops call returns feasible now, the solution is only feasible w.r.t. the new constraints, if the
          * last call detected infeasibility for the old constraints, we have to change the result to infeasible
          */
         if( lastinfeasible && *result == SCIP_FEASIBLE )
            *result = SCIP_INFEASIBLE;
      }
   }

   return SCIP_OKAY;
}

/** calls feasibility check method of constraint handler */
SCIP_RETCODE SCIPconshdlrCheck(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             printreason,        /**< Should the reason for the violation be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( conshdlr->conscheck != NULL && (!conshdlr->needscons || conshdlr->ncheckconss > 0) )
   {
      SCIPsetDebugMsg(set, "checking %d constraints of handler <%s>\n", conshdlr->ncheckconss, conshdlr->name);

      /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* start timing */
      SCIPclockStart(conshdlr->checktime, set);

      /* call external method */
      SCIP_CALL( conshdlr->conscheck(set->scip, conshdlr, conshdlr->checkconss, conshdlr->ncheckconss,
            sol, checkintegrality, checklprows, printreason, completely, result) );
      SCIPsetDebugMsg(set, " -> checking returned result <%d>\n", *result);

      /* stop timing */
      SCIPclockStop(conshdlr->checktime, set);

      /* update statistics */
      conshdlr->ncheckcalls++;

      /* perform the cached constraint updates */
      SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );

      /* evaluate result */
      if( *result != SCIP_INFEASIBLE && *result != SCIP_FEASIBLE )
      {
         SCIPerrorMessage("feasibility check of constraint handler <%s> returned invalid result <%d>\n", conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** calls propagation method of constraint handler */
SCIP_RETCODE SCIPconshdlrPropagate(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             fullpropagation,    /**< should all constraints be propagated (or only new ones)? */
   SCIP_Bool             execdelayed,        /**< execute propagation method even if it is marked to be delayed */
   SCIP_Bool             instrongbranching,  /**< are we currently doing strong branching? */
   SCIP_PROPTIMING       proptiming,         /**< current point in the node solving process */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(stat != NULL);
   assert(conshdlr->lastpropdomchgcount != stat->domchgcount
      || (0 <= conshdlr->lastnusefulpropconss && conshdlr->lastnusefulpropconss <= conshdlr->nusefulpropconss));
   assert(set != NULL);
   assert(depth >= 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( conshdlr->consprop != NULL
      && (!conshdlr->needscons || conshdlr->npropconss > 0)
      && ((depth == 0 && conshdlr->propfreq == 0)
         || (conshdlr->propfreq > 0 && depth % conshdlr->propfreq == 0)
         || conshdlr->propwasdelayed) )
   {
      /* check, if propagation method should be delayed */
      if( !conshdlr->delayprop || execdelayed )
      {
         int nconss;
         int nusefulconss;
         int nmarkedpropconss;
         int firstcons;

         /* check, if the current domains were already propagated */
         if( !fullpropagation && conshdlr->lastpropdomchgcount == stat->domchgcount && conshdlr->nmarkedpropconss == 0 )
         {
            /* all constraints that were not yet propagated on the new domains must be useful constraints, which means,
             * that the new constraints are the last constraints of the useful ones
             */
            nconss = conshdlr->nusefulpropconss - conshdlr->lastnusefulpropconss;
            nusefulconss = nconss;
            firstcons = conshdlr->lastnusefulpropconss;
         }
         else
         {
            /* on new domains, we want to propagate all constraints */
            nconss = conshdlr->npropconss;
            nusefulconss = conshdlr->nusefulpropconss;
            firstcons = 0;
         }
         assert(firstcons >= 0);
         assert(firstcons + nconss <= conshdlr->npropconss);
         assert(nusefulconss <= nconss);

         nmarkedpropconss = conshdlr->nmarkedpropconss;

         /* constraint handlers without constraints should only be called once */
         if( nconss > 0 || fullpropagation
            || (!conshdlr->needscons && conshdlr->lastpropdomchgcount != stat->domchgcount) )
         {
            SCIP_CONS** conss;
            SCIP_Longint oldndomchgs;
            SCIP_Longint oldnprobdomchgs;
            SCIP_Longint lastpropdomchgcount;
            int lastnusefulpropconss;

            SCIPsetDebugMsg(set, "propagating constraints %d to %d of %d constraints of handler <%s> (%s pseudo solution, %d useful)\n",
               firstcons, firstcons + nconss - 1, conshdlr->npropconss, conshdlr->name,
               !fullpropagation && conshdlr->lastpropdomchgcount == stat->domchgcount ? "old" : "new", nusefulconss);

            /* remember the number of processed constraints on the current domains */
            lastpropdomchgcount = stat->domchgcount;
            lastnusefulpropconss = conshdlr->nusefulpropconss;

            /* get the array of the constraints to be processed */
            conss = &(conshdlr->propconss[firstcons]);

            oldndomchgs = stat->nboundchgs + stat->nholechgs;
            oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;

            /* check, if we want to use eager evaluation */
            if( (conshdlr->eagerfreq == 0 && conshdlr->npropcalls == 0)
               || (conshdlr->eagerfreq > 0 && conshdlr->npropcalls % conshdlr->eagerfreq == 0) )
               nusefulconss = nconss;

            /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
             * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
             * external method; to avoid this, these changes will be buffered and processed after the method call
             */
            conshdlrDelayUpdates(conshdlr);
            conshdlr->duringprop = TRUE;

            /* start timing */
            if( instrongbranching )
               SCIPclockStart(conshdlr->sbproptime, set);
            else
               SCIPclockStart(conshdlr->proptime, set);

            assert(nusefulconss <= nconss);
            assert(nmarkedpropconss <= nconss);

            /* call external method */
            SCIP_CALL( conshdlr->consprop(set->scip, conshdlr, conss, nconss, nusefulconss, nmarkedpropconss, proptiming, result) );
            SCIPsetDebugMsg(set, " -> propagation returned result <%d>\n", *result);

            /* stop timing */
            if( instrongbranching )
               SCIPclockStop(conshdlr->sbproptime, set);
            else
               SCIPclockStop(conshdlr->proptime, set);

            /* perform the cached constraint updates */
            conshdlr->duringprop = FALSE;
            SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );

            /* update statistics */
            if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
            {
               conshdlr->lastpropdomchgcount = lastpropdomchgcount;
               conshdlr->lastnusefulpropconss = MIN(conshdlr->nusefulpropconss, lastnusefulpropconss);
               conshdlr->npropcalls++;
            }
            else
            {
               assert(lastpropdomchgcount == stat->domchgcount);
               assert(lastnusefulpropconss == conshdlr->nusefulpropconss);
            }
            if( *result == SCIP_CUTOFF )
               conshdlr->ncutoffs++;

            /* update domain reductions; therefore remove the domain
             * reduction counts which were generated in probing mode */
            conshdlr->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
            conshdlr->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

            /* check result code of callback method */
            if( *result != SCIP_CUTOFF
               && *result != SCIP_REDUCEDDOM
               && *result != SCIP_DIDNOTFIND
               && *result != SCIP_DIDNOTRUN
               && *result != SCIP_DELAYED
               && *result != SCIP_DELAYNODE )
            {
               SCIPerrorMessage("propagation method of constraint handler <%s> returned invalid result <%d>\n",
                  conshdlr->name, *result);
               return SCIP_INVALIDRESULT;
            }
         }
      }
      else
      {
         SCIPsetDebugMsg(set, "propagation method of constraint handler <%s> was delayed\n", conshdlr->name);
         *result = SCIP_DELAYED;
      }

      /* remember whether propagation method was delayed */
      conshdlr->propwasdelayed = (*result == SCIP_DELAYED);
   }

   return SCIP_OKAY;
}

/** calls presolving method of constraint handler */
SCIP_RETCODE SCIPconshdlrPresolve(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PRESOLTIMING     timing,             /**< current presolving timing */
   int                   nrounds,            /**< number of presolving rounds already done */
   int*                  nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*                  naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*                  nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*                  nchgbds,            /**< pointer to total number of variable bounds tightened of all presolvers */
   int*                  naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*                  ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*                  naddconss,          /**< pointer to total number of added constraints of all presolvers */
   int*                  nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*                  nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*                  nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nchgvartypes != NULL);
   assert(nchgbds != NULL);
   assert(naddholes != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);
   assert(nupgdconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( conshdlr->conspresol != NULL
      && (!conshdlr->needscons || conshdlr->nactiveconss > 0)
      && (conshdlr->maxprerounds == -1 || nrounds < conshdlr->maxprerounds ) )
   {
      SCIPsetDebugMsg(set, "presolving %d constraints of handler <%s>\n", conshdlr->nactiveconss, conshdlr->name);

      /* check, if presolving method should be executed for the current timing */
      if( timing & conshdlr->presoltiming )
      {
         int nnewfixedvars;
         int nnewaggrvars;
         int nnewchgvartypes;
         int nnewchgbds;
         int nnewholes;
         int nnewdelconss;
         int nnewaddconss;
         int nnewupgdconss;
         int nnewchgcoefs;
         int nnewchgsides;

         /* calculate the number of changes since last call */
         nnewfixedvars = *nfixedvars - conshdlr->lastnfixedvars;
         nnewaggrvars = *naggrvars - conshdlr->lastnaggrvars;
         nnewchgvartypes = *nchgvartypes - conshdlr->lastnchgvartypes;
         nnewchgbds = *nchgbds - conshdlr->lastnchgbds;
         nnewholes = *naddholes - conshdlr->lastnaddholes;
         nnewdelconss = *ndelconss - conshdlr->lastndelconss;
         nnewaddconss = *naddconss - conshdlr->lastnaddconss;
         nnewupgdconss = *nupgdconss - conshdlr->lastnupgdconss;
         nnewchgcoefs = *nchgcoefs - conshdlr->lastnchgcoefs;
         nnewchgsides = *nchgsides - conshdlr->lastnchgsides;

         /* remember the old number of changes */
         conshdlr->lastnfixedvars = *nfixedvars;
         conshdlr->lastnaggrvars = *naggrvars;
         conshdlr->lastnchgvartypes = *nchgvartypes;
         conshdlr->lastnchgbds = *nchgbds;
         conshdlr->lastnaddholes = *naddholes;
         conshdlr->lastndelconss = *ndelconss;
         conshdlr->lastnaddconss = *naddconss;
         conshdlr->lastnupgdconss = *nupgdconss;
         conshdlr->lastnchgcoefs = *nchgcoefs;
         conshdlr->lastnchgsides = *nchgsides;

         /* because during constraint processing, constraints of this handler may be deleted, activated, deactivated,
          * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
          * external method; to avoid this, these changes will be buffered and processed after the method call
          */
         conshdlrDelayUpdates(conshdlr);

         /* start timing */
         SCIPclockStart(conshdlr->presoltime, set);

         /* call external method */
         SCIP_CALL( conshdlr->conspresol(set->scip, conshdlr, conshdlr->conss, conshdlr->nactiveconss, nrounds, timing,
               nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
               nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides,
               nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes,
               ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

         /* stop timing */
         SCIPclockStop(conshdlr->presoltime, set);

         /* perform the cached constraint updates */
         SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );

         /* count the new changes */
         conshdlr->nfixedvars += *nfixedvars - conshdlr->lastnfixedvars;
         conshdlr->naggrvars += *naggrvars - conshdlr->lastnaggrvars;
         conshdlr->nchgvartypes += *nchgvartypes - conshdlr->lastnchgvartypes;
         conshdlr->nchgbds += *nchgbds - conshdlr->lastnchgbds;
         conshdlr->naddholes += *naddholes - conshdlr->lastnaddholes;
         conshdlr->ndelconss += *ndelconss - conshdlr->lastndelconss;
         conshdlr->naddconss += *naddconss - conshdlr->lastnaddconss;
         conshdlr->nupgdconss += *nupgdconss - conshdlr->lastnupgdconss;
         conshdlr->nchgcoefs += *nchgcoefs - conshdlr->lastnchgcoefs;
         conshdlr->nchgsides += *nchgsides - conshdlr->lastnchgsides;

         /* check result code of callback method */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_UNBOUNDED
            && *result != SCIP_SUCCESS
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_DELAYED )
         {
            SCIPerrorMessage("presolving method of constraint handler <%s> returned invalid result <%d>\n", 
               conshdlr->name, *result);
            return SCIP_INVALIDRESULT;
         }

         /* increase the number of calls, if the presolving method tried to find reductions */
         if( *result != SCIP_DIDNOTRUN )
            ++(conshdlr->npresolcalls);
      }

      SCIPsetDebugMsg(set, "after presolving %d constraints left of handler <%s>\n", conshdlr->nactiveconss, conshdlr->name);
   }

   return SCIP_OKAY;
}

/** calls variable deletion method of constraint handler */
SCIP_RETCODE SCIPconshdlrDelVars(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( conshdlr->consdelvars != NULL )
   {
      SCIPsetDebugMsg(set, "deleting variables in constraints of handler <%s>\n", conshdlr->name);

      /* during constraint processing, constraints of this handler may be deleted, activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* call external method */
      SCIP_CALL( conshdlr->consdelvars(set->scip, conshdlr, conshdlr->conss, conshdlr->nconss) );

      /* perform the cached constraint updates */
      SCIP_CALL( conshdlrForceUpdates(conshdlr, blkmem, set, stat) );
   }

   return SCIP_OKAY;
}

/** locks rounding of variables involved in the given constraint constraint handler that doesn't need constraints */
SCIP_RETCODE SCIPconshdlrLockVars(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->conslock != NULL);
   assert(!conshdlr->needscons);

   SCIP_CALL( conshdlr->conslock(set->scip, conshdlr, NULL, +1, 0) );

   return SCIP_OKAY;
}

/** unlocks rounding of variables involved in the given constraint constraint handler that doesn't need constraints */
SCIP_RETCODE SCIPconshdlrUnlockVars(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->conslock != NULL);
   assert(!conshdlr->needscons);

   SCIP_CALL( conshdlr->conslock(set->scip, conshdlr, NULL, -1, 0) );

   return SCIP_OKAY;
}

/** gets name of constraint handler */
const char* SCIPconshdlrGetName(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->name;
}

/** gets description of constraint handler */
const char* SCIPconshdlrGetDesc(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->desc;
}

/** gets user data of constraint handler */
SCIP_CONSHDLRDATA* SCIPconshdlrGetData(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->conshdlrdata;
}

/** sets user data of constraint handler; user has to free old data in advance! */
void SCIPconshdlrSetData(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< new constraint handler user data */
   )
{
   assert(conshdlr != NULL);

   conshdlr->conshdlrdata = conshdlrdata;
}

/** sets all separation related callbacks of the constraint handler */
void SCIPconshdlrSetSepa(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSSEPALP  ((*conssepalp)),    /**< separate cutting planes for LP solution */
   SCIP_DECL_CONSSEPASOL ((*conssepasol)),   /**< separate cutting planes for arbitrary primal solution */
   int                   sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   sepapriority,       /**< priority of the constraint handler for separation */
   SCIP_Bool             delaysepa           /**< should separation method be delayed, if other separators found cuts? */
   )
{
   assert(conshdlr != NULL);

   assert(conssepalp != NULL || conssepasol != NULL || sepafreq == -1);

   conshdlr->conssepalp = conssepalp;
   conshdlr->conssepasol = conssepasol;
   conshdlr->sepafreq = sepafreq;
   conshdlr->sepapriority = sepapriority;
   conshdlr->delaysepa = delaysepa;
}

/** sets both the propagation callback and the propagation frequency of the constraint handler */
void SCIPconshdlrSetProp(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   int                   propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   SCIP_Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       timingmask          /**< positions in the node solving loop where propagators should be executed */
   )
{
   assert(conshdlr != NULL);

   assert(consprop != NULL || propfreq == -1);

   conshdlr->consprop = consprop;
   conshdlr->propfreq = propfreq;
   conshdlr->delayprop = delayprop;
   conshdlr->proptiming = timingmask;
}

/** sets copy method of both the constraint handler and each associated constraint */
void SCIPconshdlrSetEnforelax(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSENFORELAX ((*consenforelax)) /**< constraint copying method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consenforelax = consenforelax;
}

/** sets copy method of both the constraint handler and each associated constraint */
void SCIPconshdlrSetCopy(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSHDLRCOPY((*conshdlrcopy)),  /**< copy method of constraint handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONSCOPY    ((*conscopy))       /**< constraint copying method */
   )
{
   assert(conshdlr != NULL);

   assert(!conshdlr->needscons || (conshdlrcopy == NULL) == (conscopy == NULL));

   conshdlr->conshdlrcopy = conshdlrcopy;
   conshdlr->conscopy = conscopy;
}

/** sets destructor method of constraint handler */
void SCIPconshdlrSetFree(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSFREE    ((*consfree))       /**< destructor of constraint handler */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consfree = consfree;
}

/** sets initialization method of constraint handler */
void SCIPconshdlrSetInit(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSINIT    ((*consinit))       /**< initialize constraint handler */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consinit = consinit;
}

/** sets deinitialization method of constraint handler */
void SCIPconshdlrSetExit(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSEXIT    ((*consexit))       /**< deinitialize constraint handler */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consexit = consexit;
}

/** sets solving process initialization method of constraint handler */
void SCIPconshdlrSetInitsol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSINITSOL((*consinitsol))     /**< solving process initialization method of constraint handler */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consinitsol = consinitsol;
}

/** sets solving process deinitialization method of constraint handler */
void SCIPconshdlrSetExitsol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSEXITSOL ((*consexitsol))    /**< solving process deinitialization method of constraint handler */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consexitsol = consexitsol;
}

/** sets preprocessing initialization method of constraint handler */
void SCIPconshdlrSetInitpre(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSINITPRE((*consinitpre))     /**< preprocessing initialization method of constraint handler */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consinitpre = consinitpre;
}

/** sets preprocessing deinitialization method of constraint handler */
void SCIPconshdlrSetExitpre(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSEXITPRE((*consexitpre))     /**< preprocessing deinitialization method of constraint handler */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consexitpre = consexitpre;
}

/** sets presolving method of constraint handler */
SCIP_RETCODE SCIPconshdlrSetPresol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method of constraint handler */
   int                   maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   SCIP_PRESOLTIMING     presoltiming        /**< timing mask of the constraint handler's presolving method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->conspresol = conspresol;
   conshdlr->maxprerounds = maxprerounds;

   /* the interface change from delay flags to timings cannot be recognized at compile time: Exit with an appropriate
    * error message
    */
   if( presoltiming < SCIP_PRESOLTIMING_FAST || presoltiming > SCIP_PRESOLTIMING_MAX )
   {
      SCIPmessagePrintError("ERROR: 'PRESOLDELAY'-flag no longer available since SCIP 3.2, use an appropriate "
         "'SCIP_PRESOLTIMING' for <%s> constraint handler instead.\n", conshdlr->name);

      return SCIP_PARAMETERWRONGVAL;
   }

   conshdlr->presoltiming = presoltiming;

   return SCIP_OKAY;
}

/** sets method of constraint handler to free specific constraint data */
void SCIPconshdlrSetDelete(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSDELETE  ((*consdelete))     /**< free specific constraint data */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consdelete = consdelete;
}

/** sets method of constraint handler to transform constraint data into data belonging to the transformed problem */
void SCIPconshdlrSetTrans(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSTRANS   ((*constrans))      /**< transform constraint data into data belonging to the transformed problem */
   )
{
   assert(conshdlr != NULL);

   conshdlr->constrans = constrans;
}

/** sets method of constraint handler to initialize LP with relaxations of "initial" constraints */
void SCIPconshdlrSetInitlp(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSINITLP  ((*consinitlp))     /**< initialize LP with relaxations of "initial" constraints */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consinitlp = consinitlp;
}

/** sets propagation conflict resolving method of constraint handler */
void SCIPconshdlrSetResprop(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSRESPROP ((*consresprop))    /**< propagation conflict resolving method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consresprop = consresprop;
}

/** sets activation notification method of constraint handler */
void SCIPconshdlrSetActive(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSACTIVE  ((*consactive))     /**< activation notification method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consactive = consactive;
}

/** sets deactivation notification method of constraint handler */
void SCIPconshdlrSetDeactive(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSDEACTIVE((*consdeactive))   /**< deactivation notification method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consdeactive = consdeactive;
}

/** sets enabling notification method of constraint handler */
void SCIPconshdlrSetEnable(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSENABLE  ((*consenable))     /**< enabling notification method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consenable = consenable;
}

/** sets disabling notification method of constraint handler */
void SCIPconshdlrSetDisable(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSDISABLE ((*consdisable))    /**< disabling notification method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consdisable = consdisable;
}

/** sets variable deletion method of constraint handler */
void SCIPconshdlrSetDelvars(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSDELVARS ((*consdelvars))    /**< variable deletion method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consdelvars = consdelvars;
}

/** sets constraint display method of constraint handler */
void SCIPconshdlrSetPrint(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSPRINT   ((*consprint))      /**< constraint display method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consprint = consprint;
}

/** sets constraint parsing method of constraint handler */
void SCIPconshdlrSetParse(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSPARSE   ((*consparse))      /**< constraint parsing method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consparse = consparse;
}

/** sets constraint variable getter method of constraint handler */
void SCIPconshdlrSetGetVars(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSGETVARS ((*consgetvars))    /**< constraint variable getter method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consgetvars = consgetvars;
}

/** sets constraint variable number getter method of constraint handler */
void SCIPconshdlrSetGetNVars(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSGETNVARS((*consgetnvars))   /**< constraint variable number getter method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consgetnvars = consgetnvars;
}

/** sets diving enforcement method of constraint handler */
void SCIPconshdlrSetGetDiveBdChgs(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSGETDIVEBDCHGS((*consgetdivebdchgs)) /**< constraint handler diving solution enforcement method */
   )
{
   assert(conshdlr != NULL);

   conshdlr->consgetdivebdchgs = consgetdivebdchgs;
}

/** gets array with constraints of constraint handler; the first SCIPconshdlrGetNActiveConss() entries are the active
 *  constraints, the last SCIPconshdlrGetNConss() - SCIPconshdlrGetNActiveConss() constraints are deactivated
 *
 *  @note A constraint is active if it is global and was not removed or it was added locally (in that case the local
 *        flag is TRUE) and the current node belongs to the corresponding sub tree.
 */
SCIP_CONS** SCIPconshdlrGetConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->conss;
}

/** gets array with enforced constraints of constraint handler; this is local information */
SCIP_CONS** SCIPconshdlrGetEnfoConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->enfoconss;
}

/** gets array with checked constraints of constraint handler; this is local information */
SCIP_CONS** SCIPconshdlrGetCheckConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->checkconss;
}

/** gets total number of existing transformed constraints of constraint handler */
int SCIPconshdlrGetNConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nconss;
}

/** gets number of enforced constraints of constraint handler; this is local information */
int SCIPconshdlrGetNEnfoConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nenfoconss;
}

/** gets number of checked constraints of constraint handler; this is local information */
int SCIPconshdlrGetNCheckConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ncheckconss;
}

/** gets number of active constraints of constraint handler
 *
 *  @note A constraint is active if it is global and was not removed or it was added locally (in that case the local
 *        flag is TRUE) and the current node belongs to the corresponding sub tree.
 */
int SCIPconshdlrGetNActiveConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nactiveconss;
}

/** gets number of enabled constraints of constraint handler */
int SCIPconshdlrGetNEnabledConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nenabledconss;
}

/** enables or disables all clocks of \p conshdlr, depending on the value of the flag */
void SCIPconshdlrEnableOrDisableClocks(
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the constraint handler be enabled? */
   )
{
   assert(conshdlr != NULL);

   SCIPclockEnableOrDisable(conshdlr->setuptime, enable);
   SCIPclockEnableOrDisable(conshdlr->checktime, enable);
   SCIPclockEnableOrDisable(conshdlr->enfolptime, enable);
   SCIPclockEnableOrDisable(conshdlr->enfopstime, enable);
   SCIPclockEnableOrDisable(conshdlr->enforelaxtime, enable);
   SCIPclockEnableOrDisable(conshdlr->presoltime, enable);
   SCIPclockEnableOrDisable(conshdlr->proptime, enable);
   SCIPclockEnableOrDisable(conshdlr->resproptime, enable);
   SCIPclockEnableOrDisable(conshdlr->sbproptime, enable);
   SCIPclockEnableOrDisable(conshdlr->sepatime, enable);
}

/** gets time in seconds used for setting up this constraint handler for new stages */
SCIP_Real SCIPconshdlrGetSetupTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->setuptime);
}

/** gets time in seconds used for presolving in this constraint handler */
SCIP_Real SCIPconshdlrGetPresolTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->presoltime);
}

/** gets time in seconds used for separation in this constraint handler */
SCIP_Real SCIPconshdlrGetSepaTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->sepatime);
}

/** gets time in seconds used for LP enforcement in this constraint handler */
SCIP_Real SCIPconshdlrGetEnfoLPTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->enfolptime);
}

/** gets time in seconds used for pseudo enforcement in this constraint handler */
SCIP_Real SCIPconshdlrGetEnfoPSTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->enfopstime);
}

/** gets time in seconds used for relaxation enforcement in this constraint handler */
SCIP_Real SCIPconshdlrGetEnfoRelaxTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->enforelaxtime);
}

/** gets time in seconds used for propagation in this constraint handler */
SCIP_Real SCIPconshdlrGetPropTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->proptime);
}

/** gets time in seconds used for propagation in this constraint handler during strong branching */
SCIP_Real SCIPconshdlrGetStrongBranchPropTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->sbproptime);
}

/** gets time in seconds used for feasibility checking in this constraint handler */
SCIP_Real SCIPconshdlrGetCheckTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->checktime);
}

/** gets time in seconds used for resolving propagation in this constraint handler */
SCIP_Real SCIPconshdlrGetRespropTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->resproptime);
}

/** gets number of calls to the constraint handler's separation method */
SCIP_Longint SCIPconshdlrGetNSepaCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nsepacalls;
}

/** gets number of calls to the constraint handler's LP enforcing method */
SCIP_Longint SCIPconshdlrGetNEnfoLPCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nenfolpcalls;
}

/** gets number of calls to the constraint handler's pseudo enforcing method */
SCIP_Longint SCIPconshdlrGetNEnfoPSCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nenfopscalls;
}

/** gets number of calls to the constraint handler's relaxation enforcing method */
SCIP_Longint SCIPconshdlrGetNEnfoRelaxCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nenforelaxcalls;
}

/** gets number of calls to the constraint handler's propagation method */
SCIP_Longint SCIPconshdlrGetNPropCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->npropcalls;
}

/** gets number of calls to the constraint handler's checking method */
SCIP_Longint SCIPconshdlrGetNCheckCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ncheckcalls;
}

/** gets number of calls to the constraint handler's resolve propagation method */
SCIP_Longint SCIPconshdlrGetNRespropCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nrespropcalls;
}

/** gets total number of times, this constraint handler detected a cutoff */
SCIP_Longint SCIPconshdlrGetNCutoffs(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ncutoffs;
}

/** gets total number of cuts found by this constraint handler */
SCIP_Longint SCIPconshdlrGetNCutsFound(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ncutsfound;
}

/** gets total number of cuts found by this constraint handler applied to lp */
SCIP_Longint SCIPconshdlrGetNCutsApplied(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ncutsapplied;
}

/** increase count of applied cuts */
void SCIPconshdlrIncNAppliedCuts(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   ++conshdlr->ncutsapplied;
}

/** increase count of found cuts */
void SCIPconshdlrIncNCutsFound(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   ++conshdlr->ncutsfound;
}

/** gets total number of additional constraints added by this constraint handler */
SCIP_Longint SCIPconshdlrGetNConssFound(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nconssfound;
}

/** gets total number of domain reductions found by this constraint handler */
SCIP_Longint SCIPconshdlrGetNDomredsFound(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ndomredsfound;
}

/** gets number of children created by this constraint handler */
SCIP_Longint SCIPconshdlrGetNChildren(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nchildren;
}

/** gets maximum number of active constraints of constraint handler existing at the same time */
int SCIPconshdlrGetMaxNActiveConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->maxnactiveconss;
}

/** gets initial number of active constraints of constraint handler */
int SCIPconshdlrGetStartNActiveConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->startnactiveconss;
}

/** gets number of variables fixed in presolving method of constraint handler */
int SCIPconshdlrGetNFixedVars(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nfixedvars;
}

/** gets number of variables aggregated in presolving method of constraint handler */
int SCIPconshdlrGetNAggrVars(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->naggrvars;
}

/** gets number of variable types changed in presolving method of constraint handler */
int SCIPconshdlrGetNChgVarTypes(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nchgvartypes;
}

/** gets number of bounds changed in presolving method of constraint handler */
int SCIPconshdlrGetNChgBds(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nchgbds;
}

/** gets number of holes added to domains of variables in presolving method of constraint handler */
int SCIPconshdlrGetNAddHoles(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->naddholes;
}

/** gets number of constraints deleted in presolving method of constraint handler */
int SCIPconshdlrGetNDelConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ndelconss;
}

/** gets number of constraints added in presolving method of constraint handler */
int SCIPconshdlrGetNAddConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->naddconss;
}

/** gets number of constraints upgraded in presolving method of constraint handler */
int SCIPconshdlrGetNUpgdConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nupgdconss;
}

/** gets number of coefficients changed in presolving method of constraint handler */
int SCIPconshdlrGetNChgCoefs(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nchgcoefs;
}

/** gets number of constraint sides changed in presolving method of constraint handler */
int SCIPconshdlrGetNChgSides(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nchgsides;
}

/** gets number of times the presolving method of the constraint handler was called and tried to find reductions */
int SCIPconshdlrGetNPresolCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->npresolcalls;
}

/** gets separation priority of constraint handler */
int SCIPconshdlrGetSepaPriority(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->sepapriority;
}

/** gets enforcing priority of constraint handler */
int SCIPconshdlrGetEnfoPriority(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->enfopriority;
}

/** gets checking priority of constraint handler */
int SCIPconshdlrGetCheckPriority(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->checkpriority;
}

/** gets separation frequency of constraint handler */
int SCIPconshdlrGetSepaFreq(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->sepafreq;
}

/** gets propagation frequency of constraint handler */
int SCIPconshdlrGetPropFreq(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->propfreq;
}

/** gets frequency of constraint handler for eager evaluations in separation, propagation and enforcement */
int SCIPconshdlrGetEagerFreq(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->eagerfreq;
}

/** needs constraint handler a constraint to be called? */
SCIP_Bool SCIPconshdlrNeedsCons(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->needscons;
}

/** does the constraint handler perform presolving? */
SCIP_Bool SCIPconshdlrDoesPresolve(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return (conshdlr->conspresol != NULL);
}

/** should separation method be delayed, if other separators found cuts? */
SCIP_Bool SCIPconshdlrIsSeparationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->delaysepa;
}

/** should propagation method be delayed, if other propagators found reductions? */
SCIP_Bool SCIPconshdlrIsPropagationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->delayprop;
}

/** was LP separation method delayed at the last call? */
SCIP_Bool SCIPconshdlrWasLPSeparationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->sepalpwasdelayed;
}

/** was primal solution separation method delayed at the last call? */
SCIP_Bool SCIPconshdlrWasSolSeparationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->sepasolwasdelayed;
}

/** was propagation method delayed at the last call? */
SCIP_Bool SCIPconshdlrWasPropagationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->propwasdelayed;
}

/** is constraint handler initialized? */
SCIP_Bool SCIPconshdlrIsInitialized(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->initialized;
}

/** does the constraint handler have a copy function? */
SCIP_Bool SCIPconshdlrIsClonable(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return (conshdlr->conshdlrcopy != NULL);
}

/** returns the timing mask of the propagation method of the constraint handler */
SCIP_PROPTIMING SCIPconshdlrGetPropTiming(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->proptiming;
}

/** sets the timing mask of the propagation method of the constraint handler */
void SCIPconshdlrSetPropTiming(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_PROPTIMING       proptiming          /**< timing mask to be set */
   )
{
   assert(conshdlr != NULL);

   conshdlr->proptiming = proptiming;
}


/** returns the timing mask of the presolving method of the constraint handler */
SCIP_PRESOLTIMING SCIPconshdlrGetPresolTiming(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->presoltiming;
}

/** sets the timing mask of the presolving method of the constraint handler */
void SCIPconshdlrSetPresolTiming(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_PRESOLTIMING     presoltiming        /** timing mask to be set */
   )
{
   assert(conshdlr != NULL);

   conshdlr->presoltiming = presoltiming;
}


/*
 * Constraint set change methods
 */

/** creates empty constraint set change data */
static
SCIP_RETCODE conssetchgCreate(
   SCIP_CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(conssetchg != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, conssetchg) );
   (*conssetchg)->addedconss = NULL;
   (*conssetchg)->disabledconss = NULL;
   (*conssetchg)->addedconsssize = 0;
   (*conssetchg)->naddedconss = 0;
   (*conssetchg)->disabledconsssize = 0;
   (*conssetchg)->ndisabledconss = 0;

   return SCIP_OKAY;
}

/** releases all constraints of the constraint set change data */
static
SCIP_RETCODE conssetchgRelease(
   SCIP_CONSSETCHG*      conssetchg,         /**< constraint set change data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(conssetchg != NULL);

   /* release constraints */
   for( i = 0; i < conssetchg->naddedconss; ++i )
   {
      if( conssetchg->addedconss[i] != NULL )
      {
         SCIP_CALL( SCIPconsRelease(&conssetchg->addedconss[i], blkmem, set) );
      }
   }
   for( i = 0; i < conssetchg->ndisabledconss; ++i )
   {
      if( conssetchg->disabledconss[i] != NULL )
      {
         SCIP_CALL( SCIPconsRelease(&conssetchg->disabledconss[i], blkmem, set) );
      }
   }

   return SCIP_OKAY;
}

/** frees constraint set change data and releases all included constraints */
SCIP_RETCODE SCIPconssetchgFree(
   SCIP_CONSSETCHG**     conssetchg,         /**< pointer to constraint set change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conssetchg != NULL);
   assert(blkmem != NULL);

   if( *conssetchg != NULL )
   {
      /* release constraints */
      SCIP_CALL( conssetchgRelease(*conssetchg, blkmem, set) );

      /* free memory */
      BMSfreeBlockMemoryArrayNull(blkmem, &(*conssetchg)->addedconss, (*conssetchg)->addedconsssize);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*conssetchg)->disabledconss, (*conssetchg)->disabledconsssize);
      BMSfreeBlockMemory(blkmem, conssetchg);
   }

   return SCIP_OKAY;
}

/** ensures, that addedconss array can store at least num entries */
static
SCIP_RETCODE conssetchgEnsureAddedconssSize(
   SCIP_CONSSETCHG*      conssetchg,         /**< constraint set change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(conssetchg != NULL);

   if( num > conssetchg->addedconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conssetchg->addedconss, conssetchg->addedconsssize, newsize) );
      conssetchg->addedconsssize = newsize;
   }
   assert(num <= conssetchg->addedconsssize);

   return SCIP_OKAY;
}

/** ensures, that disabledconss array can store at least num entries */
static
SCIP_RETCODE conssetchgEnsureDisabledconssSize(
   SCIP_CONSSETCHG*      conssetchg,         /**< constraint set change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(conssetchg != NULL);

   if( num > conssetchg->disabledconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conssetchg->disabledconss, conssetchg->disabledconsssize, newsize) );
      conssetchg->disabledconsssize = newsize;
   }
   assert(num <= conssetchg->disabledconsssize);

   return SCIP_OKAY;
}

/** adds constraint addition to constraint set changes, and captures constraint; activates constraint if the
 *  constraint set change data is currently active
 */
SCIP_RETCODE SCIPconssetchgAddAddedCons(
   SCIP_CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_CONS*            cons,               /**< added constraint */
   int                   depth,              /**< depth of constraint set change's node */
   SCIP_Bool             focusnode,          /**< does the constraint set change belong to the focus node? */
   SCIP_Bool             active              /**< is the constraint set change currently active? */
   )
{
   assert(conssetchg != NULL);
   assert(cons != NULL);

   /* if constraint set change doesn't exist, create it */
   if( *conssetchg == NULL )
   {
      SCIP_CALL( conssetchgCreate(conssetchg, blkmem) );
   }

   /* add constraint to the addedconss array */
   SCIP_CALL( conssetchgEnsureAddedconssSize(*conssetchg, blkmem, set, (*conssetchg)->naddedconss+1) );
   (*conssetchg)->addedconss[(*conssetchg)->naddedconss] = cons;
   (*conssetchg)->naddedconss++;

   /* undelete constraint, if it was globally deleted in the past */
   cons->deleted = FALSE;

   /* capture constraint */
   SCIPconsCapture(cons);

   /* activate constraint, if node is active */
   if( active && !SCIPconsIsActive(cons) )
   {
      SCIP_CALL( SCIPconsActivate(cons, set, stat, depth, focusnode) );
      assert(SCIPconsIsActive(cons));

      /* remember, that this constraint set change data was responsible for the constraint's addition */
      cons->addconssetchg = *conssetchg;
      cons->addarraypos = (*conssetchg)->naddedconss-1;
   }

   return SCIP_OKAY;
}

/** adds constraint disabling to constraint set changes, and captures constraint */
SCIP_RETCODE SCIPconssetchgAddDisabledCons(
   SCIP_CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS*            cons                /**< disabled constraint */
   )
{
   assert(conssetchg != NULL);
   assert(cons != NULL);

   /* if constraint set change doesn't exist, create it */
   if( *conssetchg == NULL )
   {
      SCIP_CALL( conssetchgCreate(conssetchg, blkmem) );
   }

   /* add constraint to the disabledconss array */
   SCIP_CALL( conssetchgEnsureDisabledconssSize(*conssetchg, blkmem, set, (*conssetchg)->ndisabledconss+1) );
   (*conssetchg)->disabledconss[(*conssetchg)->ndisabledconss] = cons;
   (*conssetchg)->ndisabledconss++;

   /* capture constraint */
   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

/** deactivates, deletes, and releases constraint from the addedconss array of the constraint set change data */
static
SCIP_RETCODE conssetchgDelAddedCons(
   SCIP_CONSSETCHG*      conssetchg,         /**< constraint set change to delete constraint from */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   arraypos            /**< position of constraint in disabledconss array */
   )
{
   SCIP_CONS* cons;

   assert(conssetchg != NULL);
   assert(conssetchg->addedconss != NULL);
   assert(0 <= arraypos && arraypos < conssetchg->naddedconss);

   cons = conssetchg->addedconss[arraypos];
   assert(cons != NULL);

   SCIPsetDebugMsg(set, "delete added constraint <%s> at position %d from constraint set change data\n", cons->name, arraypos);

   /* remove the link to the constraint set change data */
   if( cons->addconssetchg == conssetchg )
   {
      cons->addconssetchg = NULL;
      cons->addarraypos = -1;
   }

   /* release constraint */
   SCIP_CALL( SCIPconsRelease(&conssetchg->addedconss[arraypos], blkmem, set) );

   /* we want to keep the order of the constraint additions: move all subsequent constraints one slot to the front */
   for( ; arraypos < conssetchg->naddedconss-1; ++arraypos )
   {
      conssetchg->addedconss[arraypos] = conssetchg->addedconss[arraypos+1];
      assert(conssetchg->addedconss[arraypos] != NULL);
      if( conssetchg->addedconss[arraypos]->addconssetchg == conssetchg )
      {
         assert(conssetchg->addedconss[arraypos]->addarraypos == arraypos+1);
         conssetchg->addedconss[arraypos]->addarraypos = arraypos;
      }
   }
   conssetchg->naddedconss--;

   return SCIP_OKAY;
}

/** deletes and releases deactivated constraint from the disabledconss array of the constraint set change data */
static
SCIP_RETCODE conssetchgDelDisabledCons(
   SCIP_CONSSETCHG*      conssetchg,         /**< constraint set change to apply */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   arraypos            /**< position of constraint in disabledconss array */
   )
{
   assert(conssetchg != NULL);
   assert(0 <= arraypos && arraypos < conssetchg->ndisabledconss);
   assert(conssetchg->disabledconss[arraypos] != NULL);

   SCIPsetDebugMsg(set, "delete disabled constraint <%s> at position %d from constraint set change data\n",
      conssetchg->disabledconss[arraypos]->name, arraypos);

   /* release constraint */
   SCIP_CALL( SCIPconsRelease(&conssetchg->disabledconss[arraypos], blkmem, set) );

   /* we want to keep the order of the constraint disablings: move all subsequent constraints one slot to the front */
   for( ; arraypos < conssetchg->ndisabledconss-1; ++arraypos )
   {
      conssetchg->disabledconss[arraypos] = conssetchg->disabledconss[arraypos+1];
      assert(conssetchg->disabledconss[arraypos] != NULL);
   }
   conssetchg->ndisabledconss--;

   return SCIP_OKAY;
}

/** gets added constraints data for a constraint set change */
void SCIPconssetchgGetAddedConsData(
   SCIP_CONSSETCHG*      conssetchg,         /**< constraint set change to get data from */
   SCIP_CONS***          conss,              /**< reference to constraints array added in the conssetchg, or NULL */
   int*                  nconss              /**< reference to store the size of the constraints array, or NULL */
   )
{
   assert(conssetchg != NULL);
   if( conss != NULL )
      *conss = conssetchg->addedconss;
   if( nconss != NULL )
      *nconss = conssetchg->naddedconss;
}

/** applies constraint set change */
SCIP_RETCODE SCIPconssetchgApply(
   SCIP_CONSSETCHG*      conssetchg,         /**< constraint set change to apply */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of constraint set change's node */
   SCIP_Bool             focusnode           /**< does the constraint set change belong to the focus node? */
   )
{
   SCIP_CONS* cons;
   int i;

   if( conssetchg == NULL )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "applying constraint set changes at %p: %d constraint additions, %d constraint disablings\n",
      (void*)conssetchg, conssetchg->naddedconss, conssetchg->ndisabledconss);

   /* apply constraint additions */
   i = 0;
   while( i < conssetchg->naddedconss )
   {
      cons = conssetchg->addedconss[i];
      assert(cons != NULL);
      assert(!cons->update);

      /* if constraint is already active, or if constraint is globally deleted, it can be removed from addedconss array */
      if( cons->active || cons->deleted )
      {
         /* delete constraint from addedcons array, the empty slot is now used by the next constraint,
          * and naddedconss was decreased, so do not increase i
          */
         SCIP_CALL( conssetchgDelAddedCons(conssetchg, blkmem, set, i) );
      }
      else
      {
         assert(cons->addconssetchg == NULL);
         assert(cons->addarraypos == -1);

         /* activate constraint */
         SCIP_CALL( SCIPconsActivate(cons, set, stat, depth, focusnode) );
         assert(cons->active);
         assert(!cons->update);

         /* remember, that this constraint set change data was responsible for the constraint's addition */
         cons->addconssetchg = conssetchg;
         cons->addarraypos = i;

         ++i; /* handle the next constraint */
      }
   }

   /* apply constraint disablings */
   i = 0;
   while( i < conssetchg->ndisabledconss )
   {
      cons = conssetchg->disabledconss[i];
      assert(cons != NULL);
      assert(!cons->update);

      /* if the constraint is disabled, we can permanently remove it from the disabledconss array */
      if( !cons->enabled )
      {
         SCIPsetDebugMsg(set, "constraint <%s> of handler <%s> was deactivated -> remove it from disabledconss array\n",
            cons->name, cons->conshdlr->name);

         /* release and remove constraint from the disabledconss array, the empty slot is now used by the next constraint
          * and ndisabledconss was decreased, so do not increase i
          */
         SCIP_CALL( conssetchgDelDisabledCons(conssetchg, blkmem, set, i) );
      }
      else
      {
         assert(cons->addarraypos >= 0);
         assert(!cons->deleted); /* deleted constraints must not be enabled! */
         SCIP_CALL( SCIPconsDisable(conssetchg->disabledconss[i], set, stat) );
         assert(!cons->update);
         assert(!cons->enabled);

         ++i; /* handle the next constraint */
      }
   }

   return SCIP_OKAY;
}

/** undoes constraint set change */
SCIP_RETCODE SCIPconssetchgUndo(
   SCIP_CONSSETCHG*      conssetchg,         /**< constraint set change to undo */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   SCIP_CONS* cons;
   int i;

   if( conssetchg == NULL )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "undoing constraint set changes at %p: %d constraint additions, %d constraint disablings\n",
      (void*)conssetchg, conssetchg->naddedconss, conssetchg->ndisabledconss);

   /* undo constraint disablings */
   for( i = conssetchg->ndisabledconss-1; i >= 0; --i )
   {
      cons = conssetchg->disabledconss[i];
      assert(cons != NULL);
      assert(!cons->update);

      /* If the constraint is inactive, we can permanently remove it from the disabledconss array. It was deactivated
       * in the subtree of the current node but not reactivated on the switching way back to the current node, which
       * means, the deactivation was more global (i.e. valid on a higher level node) than the current node and the
       * disabling at the current node doesn't have any effect anymore.
       * If the constraint is already enabled, we need not to do anything. This may happen on a path A -> B,
       * if the constraint is disabled at node B, and while processing the subtree of B, it is also disabled at
       * the more global node A. Then on the switching path back to A, the constraint is enabled at node B (which is
       * actually wrong, since it now should be disabled in the whole subtree of A, but we cannot know this), and
       * again enabled at node A (where enabling is ignored). If afterwards, a subnode of B is processed, the
       * switching disables the constraint in node A, and the disabling is then removed from node B.
       */
      if( !cons->active )
      {
         SCIPsetDebugMsg(set, "constraint <%s> of handler <%s> was deactivated -> remove it from disabledconss array\n",
            cons->name, cons->conshdlr->name);

         /* release and remove constraint from the disabledconss array */
         SCIP_CALL( conssetchgDelDisabledCons(conssetchg, blkmem, set, i) );
      }
      else if( !cons->enabled )
      {
         assert(cons->addarraypos >= 0);
         assert(!cons->deleted); /* deleted constraints must not be active! */
         SCIP_CALL( SCIPconsEnable(cons, set, stat) );
         assert(!cons->update);
         assert(!cons->active || cons->enabled);
      }
   }

   /* undo constraint additions */
   for( i = conssetchg->naddedconss-1; i >= 0; --i )
   {
      cons = conssetchg->addedconss[i];
      assert(cons != NULL);
      assert(!cons->update);

      /* If the constraint is already deactivated, we need not to do anything. This may happen on a path A -> B,
       * if the constraint is added at node B, and while processing the subtree of B, it is also added at
       * the more global node A. Then on the switching path back to A, the node is deactivated at node B (which is
       * actually wrong, since it now should be active in the whole subtree of A, but we cannot know this), and
       * again deactivated at node A (where deactivation is ignored). If afterwards, a subnode of B is processed, the
       * switching activates the constraint in node A, and the activation is then removed from node B.
       */
      if( cons->active )
      {
         assert(cons->addconssetchg == conssetchg);
         assert(cons->addarraypos == i);

         /* deactivate constraint */
         SCIP_CALL( SCIPconsDeactivate(cons, set, stat) );

         /* unlink the constraint and the constraint set change */
         cons->addconssetchg = NULL;
         cons->addarraypos = -1;
      }
      assert(!cons->active);
      assert(!cons->update);
   }

   return SCIP_OKAY;
}

/** applies constraint set change to the global problem and deletes the constraint set change data */
SCIP_RETCODE SCIPconssetchgMakeGlobal(
   SCIP_CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   )
{
   SCIP_CONS* cons;
   int i;

   assert(conssetchg != NULL);

   /* nothing to do on empty constraint set change data */
   if( *conssetchg == NULL )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "moving constraint set changes at %p to global problem: %d constraint additions, %d constraint disablings\n",
      (void*)*conssetchg, (*conssetchg)->naddedconss, (*conssetchg)->ndisabledconss);

   /* apply constraint additions to the global problem (loop backwards, because then conssetchgDelAddedCons() is
    * more efficient)
    */
   for( i = (*conssetchg)->naddedconss-1; i >= 0; --i )
   {
      cons = (*conssetchg)->addedconss[i];
      assert(cons != NULL);
      assert(!cons->update);

      /* only move constraints that are not sticking at the current node */
      if( !SCIPconsIsStickingAtNode(cons) )
      {
         /* because we first have to delete the constraint, we have to capture it in order to not loose it */
         SCIPconsCapture(cons);

         /* delete constraint addition from constraint set change data */
         SCIP_CALL( conssetchgDelAddedCons(*conssetchg, blkmem, set, i) );

         /* don't move deleted constraints to the global problem */
         if( !cons->deleted )
         {
            SCIP_CALL( SCIPprobAddCons(prob, set, stat, cons) );
         }

         /* release constraint */
         SCIP_CALL( SCIPconsRelease(&cons, blkmem, set) );
      }
   }

   /* apply constraint disablings to the global problem (loop backwards, because then conssetchgDelDisabledCons() is
    * more efficient)
    */
   for( i = (*conssetchg)->ndisabledconss-1; i >= 0; --i )
   {
      cons = (*conssetchg)->disabledconss[i];
      assert(cons != NULL);
      assert(!cons->update);

      /* only delete constraints that are not sticking at the current node */
      if( !SCIPconsIsStickingAtNode(cons) )
      {
         /* globally delete constraint */
         if( !cons->deleted )
         {
            SCIP_CALL( SCIPconsDelete(cons, blkmem, set, stat, prob, reopt) );
         }

         /* release and remove constraint from the disabledconss array */
         SCIP_CALL( conssetchgDelDisabledCons(*conssetchg, blkmem, set, i) );
      }
   }

   if( (*conssetchg)->naddedconss == 0 && (*conssetchg)->ndisabledconss == 0 )
   {
      /* free empty constraint set change data */
      SCIP_CALL( SCIPconssetchgFree(conssetchg, blkmem, set) );
   }

   return SCIP_OKAY;
}




/*
 * Constraint methods
 */

/** creates and captures a constraint, and inserts it into the conss array of its constraint handler
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution
 *  may be declared feasible even if it violates this particular constraint.
 *  This constellation should only be used, if no LP or pseudo solution can violate the constraint -- e.g. if a
 *  local constraint is redundant due to the variable's local bounds.
 */
SCIP_RETCODE SCIPconsCreate(
   SCIP_CONS**           cons,               /**< pointer to constraint */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of constraint */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   SCIP_CONSDATA*        consdata,           /**< data for this specific constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   SCIP_Bool             original,           /**< is constraint belonging to the original problem? */
   SCIP_Bool             deleteconsdata      /**< has the constraint data to be deleted if constraint is freed? */
   )
{
   assert(cons != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(name != NULL);
   assert(conshdlr != NULL);
   assert(!original || deleteconsdata);

   /* create constraint data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, cons) );
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*cons)->name, name, strlen(name)+1) );
#ifndef NDEBUG
   (*cons)->scip = set->scip;
#endif
   (*cons)->conshdlr = conshdlr;
   (*cons)->consdata = consdata;
   (*cons)->transorigcons = NULL;
   (*cons)->addconssetchg = NULL;
   (*cons)->addarraypos = -1;
   (*cons)->consspos = -1;
   (*cons)->initconsspos = -1;
   (*cons)->sepaconsspos = -1;
   (*cons)->enfoconsspos = -1;
   (*cons)->checkconsspos = -1;
   (*cons)->propconsspos = -1;
   (*cons)->activedepth = -2;
   (*cons)->validdepth = (local ? -1 : 0);
   (*cons)->age = 0.0;
   (*cons)->nlockspos = 0;
   (*cons)->nlocksneg = 0;
   (*cons)->nuses = 0;
   (*cons)->nupgradelocks = 0;
   (*cons)->initial = initial;
   (*cons)->separate = separate;
   (*cons)->enforce = enforce;
   (*cons)->check = check;
   (*cons)->propagate = propagate;
   (*cons)->sepaenabled = separate;
   (*cons)->propenabled = propagate;
   (*cons)->local = local;
   (*cons)->modifiable = modifiable;
   (*cons)->dynamic = dynamic;
   (*cons)->removable = removable;
   (*cons)->stickingatnode = stickingatnode;
   (*cons)->original = original;
   (*cons)->deleteconsdata = deleteconsdata;
   (*cons)->active = FALSE;
   (*cons)->conflict = FALSE;
   (*cons)->enabled = FALSE;
   (*cons)->obsolete = FALSE;
   (*cons)->markpropagate = TRUE;
   (*cons)->deleted = FALSE;
   (*cons)->update = FALSE;
   (*cons)->updateinsert = FALSE;
   (*cons)->updateactivate = FALSE;
   (*cons)->updatedeactivate = FALSE;
   (*cons)->updateenable = FALSE;
   (*cons)->updatedisable = FALSE;
   (*cons)->updatesepaenable = FALSE;
   (*cons)->updatesepadisable = FALSE;
   (*cons)->updatepropenable = FALSE;
   (*cons)->updatepropdisable = FALSE;
   (*cons)->updateobsolete = FALSE;
   (*cons)->updatemarkpropagate = FALSE;
   (*cons)->updateunmarkpropagate = FALSE;
   (*cons)->updatefree = FALSE;
   (*cons)->updateactfocus = FALSE;

   /* capture constraint */
   SCIPconsCapture(*cons);

   /* insert the constraint as inactive constraint into the transformed constraints array */
   if( !original )
   {
      /* check, if inserting constraint should be delayed */
      if( conshdlrAreUpdatesDelayed(conshdlr) )
      {
         SCIPsetDebugMsg(set, " -> delaying insertion of constraint <%s>\n", (*cons)->name);
         (*cons)->updateinsert = TRUE;
         SCIP_CALL( conshdlrAddUpdateCons((*cons)->conshdlr, set, *cons) );
         assert((*cons)->update);
         assert((*cons)->nuses == 2);
      }
      else
      {
         SCIP_CALL( conshdlrAddCons(conshdlr, set, *cons) );
      }
   }

   checkConssArrays(conshdlr);

   return SCIP_OKAY;
}

/** copies source constraint of source SCIP into the target constraint for the target SCIP, using the variable map for
 *  mapping the variables of the source SCIP to the variables of the target SCIP; if the copying process was successful
 *  a constraint is created and captured;
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, an LP or pseudo solution
 *  may be declared feasible even if it violates this particular constraint.
 *  This constellation should only be used, if no LP or pseudo solution can violate the constraint -- e.g. if a
 *  local constraint is redundant due to the variable's local bounds.
 */
SCIP_RETCODE SCIPconsCopy(
   SCIP_CONS**           cons,               /**< pointer to store the created target constraint */
   SCIP_SET*             set,                /**< global SCIP settings of the target SCIP */
   const char*           name,               /**< name of constraint, or NULL if the name of the source constraint should be used */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP_CONSHDLR*        sourceconshdlr,     /**< source constraint handler for this constraint */
   SCIP_CONS*            sourcecons,         /**< source constraint of the source SCIP */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to corresponding
                                              *   variables of the target SCIP */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, must not be NULL! */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid or not */
   )
{
   assert(cons != NULL);
   assert(set != NULL);
   assert(sourcescip != NULL);
   assert(sourceconshdlr != NULL);
   assert(sourcecons != NULL);
   assert(varmap != NULL);
   assert(consmap != NULL);
   assert(valid != NULL);

   /* if constraint handler does not support copying, success will return false. Constraints handlers have to actively set this to true. */
   (*valid) = FALSE;

   if( sourceconshdlr->conscopy != NULL )
   {
      SCIP_CALL( sourceconshdlr->conscopy(set->scip, cons, name, sourcescip, sourceconshdlr, sourcecons, varmap, consmap,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );
   }

   return SCIP_OKAY;
}


/** parses constraint information (in cip format) out of a string; if the parsing process was successful a constraint is
 *  created, captured, and inserted into the conss array of its constraint handler.
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, an LP or pseudo solution
 *  may be declared feasible even if it violates this particular constraint.
 *  This constellation should only be used, if no LP or pseudo solution can violate the constraint -- e.g. if a
 *  local constraint is redundant due to the variable's local bounds.
 */
SCIP_RETCODE SCIPconsParse(
   SCIP_CONS**           cons,               /**< pointer to constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler of target SCIP */
   const char*           str,                /**< string to parse for constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   )
{
   SCIP_CONSHDLR* conshdlr;
   char conshdlrname[SCIP_MAXSTRLEN];
   char consname[SCIP_MAXSTRLEN];
   char* endptr;

   assert(cons != NULL);
   assert(set != NULL);

   (*success) = FALSE;

   /* scan constraint handler name */
   assert(str != NULL);
   SCIPstrCopySection(str, '[', ']', conshdlrname, SCIP_MAXSTRLEN, &endptr);
   if ( endptr == NULL || endptr == str )
   {
      SCIPmessagePrintWarning(messagehdlr, "Syntax error: Could not find constraint handler name.\n");
      return SCIP_OKAY;
   }
   assert(endptr != NULL);
   SCIPsetDebugMsg(set, "constraint handler name <%s>\n", conshdlrname);

   /* scan constraint name */
   SCIPstrCopySection(endptr, '<', '>', consname, SCIP_MAXSTRLEN, &endptr);
   if ( endptr == NULL || endptr == str )
   {
      SCIPmessagePrintWarning(messagehdlr, "Syntax error: Could not find constraint name.\n");
      return SCIP_OKAY;
   }
   assert(endptr != NULL);
   SCIPsetDebugMsg(set, "constraint name <%s>\n", consname);

   str = endptr;

   /* skip white space */
   while ( isspace((unsigned char)* str) )
      ++str;

   /* check for colon */
   if( *str != ':' )
   {
      SCIPmessagePrintWarning(messagehdlr, "Syntax error: Could not find colon ':' after constraint name.\n");
      return SCIP_OKAY;
   }

   /* skip colon */
   ++str;

   /* skip white space */
   while ( isspace((unsigned char)* str) )
      ++str;

   /* check if a constraint handler with parsed name exists */
   conshdlr = SCIPsetFindConshdlr(set, conshdlrname);

   if( conshdlr == NULL )
   {
      SCIPmessagePrintWarning(messagehdlr, "constraint handler <%s> doesn't exist in SCIP data structure\n", conshdlrname);
   }
   else
   {
      assert( conshdlr != NULL );
      if ( conshdlr->consparse == NULL )
      {
         SCIPmessagePrintWarning(messagehdlr, "constraint handler <%s> does not support parsing constraints\n", conshdlrname);
      }
      else
      {
         SCIP_CALL( conshdlr->consparse(set->scip, conshdlr, cons, consname, str,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, success) );
      }
   }

   return SCIP_OKAY;
}

/** change name of given constraint */
SCIP_RETCODE SCIPconsChgName(
   SCIP_CONS*            cons,               /**< problem constraint */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   const char*           name                /**< new name of constraint */
   )
{
   assert(cons != NULL);
   assert(cons->name != NULL);

   /* free old constraint name */
   BMSfreeBlockMemoryArray(blkmem, &cons->name, strlen(cons->name)+1);

   /* copy new constraint name */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &cons->name, name, strlen(name)+1) );

   return SCIP_OKAY;
}


/** frees a constraint and removes it from the conss array of its constraint handler */
SCIP_RETCODE SCIPconsFree(
   SCIP_CONS**           cons,               /**< constraint to free */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(*cons != NULL);
   assert((*cons)->conshdlr != NULL);
   assert((*cons)->nuses == 0);
   assert(!(*cons)->active);
   assert(!(*cons)->update);
   assert(!(*cons)->original || (*cons)->transorigcons == NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert((*cons)->scip == set->scip);

   SCIPsetDebugMsg(set, "freeing constraint <%s> at conss pos %d of handler <%s>\n",
      (*cons)->name, (*cons)->consspos, (*cons)->conshdlr->name);

   /* free constraint data */
   if( (*cons)->conshdlr->consdelete != NULL && (*cons)->consdata != NULL && (*cons)->deleteconsdata )
   {
      SCIP_CALL( (*cons)->conshdlr->consdelete(set->scip, (*cons)->conshdlr, *cons, &(*cons)->consdata) );
   }
   else if( !(*cons)->deleteconsdata )
      (*cons)->consdata = NULL;
   assert((*cons)->consdata == NULL);

   /* unlink transformed and original constraint */
   if( (*cons)->transorigcons != NULL )
   {
      assert(!(*cons)->original);
      assert((*cons)->transorigcons->original);
      assert((*cons)->transorigcons->transorigcons == *cons);

      (*cons)->transorigcons->transorigcons = NULL;
   }

   /* remove constraint from the transformed constraints array */
   if( !(*cons)->original )
   {
      conshdlrDelCons((*cons)->conshdlr, *cons);
      checkConssArrays((*cons)->conshdlr);
   }
   assert((*cons)->consspos == -1);

   /* free constraint */
   BMSfreeBlockMemoryArray(blkmem, &(*cons)->name, strlen((*cons)->name)+1);
   BMSfreeBlockMemory(blkmem, cons);

   return SCIP_OKAY;
}

/** increases usage counter of constraint */
void SCIPconsCapture(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(cons->nuses >= 0);

   SCIPdebugMessage("capture constraint <%s> with nuses=%d, cons pointer %p\n", cons->name, cons->nuses, (void*)cons);
   cons->nuses++;
}

/** decreases usage counter of constraint, and frees memory if necessary */
SCIP_RETCODE SCIPconsRelease(
   SCIP_CONS**           cons,               /**< pointer to constraint */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(blkmem != NULL);
   assert(cons != NULL);
   assert(*cons != NULL);
   assert((*cons)->conshdlr != NULL);
   assert((*cons)->nuses >= 1);
   assert(set != NULL);
   assert((*cons)->scip == set->scip);

   SCIPsetDebugMsg(set, "release constraint <%s> with nuses=%d, cons pointer %p\n", (*cons)->name, (*cons)->nuses, (void*)(*cons));
   (*cons)->nuses--;
   if( (*cons)->nuses == 0 )
   {
      assert(!(*cons)->active || (*cons)->updatedeactivate);

      /* check, if freeing constraint should be delayed */
      if( conshdlrAreUpdatesDelayed((*cons)->conshdlr) )
      {
         SCIPsetDebugMsg(set, " -> delaying freeing constraint <%s>\n", (*cons)->name);
         (*cons)->updatefree = TRUE;
         SCIP_CALL( conshdlrAddUpdateCons((*cons)->conshdlr, set, *cons) );
         assert((*cons)->update);
         assert((*cons)->nuses == 1);
      }
      else
      {
         SCIP_CALL( SCIPconsFree(cons, blkmem, set) );
      }
   }
   *cons = NULL;

   return SCIP_OKAY;
}

/** outputs constraint information to file stream */
SCIP_RETCODE SCIPconsPrint(
   SCIP_CONS*            cons,               /**< constraint to print */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   SCIPmessageFPrintInfo(messagehdlr, file, "  [%s] <%s>: ", conshdlr->name, cons->name);

   if( conshdlr->consprint != NULL )
   {
      SCIP_CALL( conshdlr->consprint(set->scip, conshdlr, cons, file) );
   }
   else
      SCIPmessageFPrintInfo(messagehdlr, file, "constraint handler <%s> doesn't support printing constraint", conshdlr->name);

   return SCIP_OKAY;
}

/** method to collect the variables of a constraint
 *
 *  If the number of variables is greater than the available slots in the variable array, nothing happens except that
 *  the success point is set to FALSE. With the method SCIPconsGetNVars() it is possible to get the number of variables
 *  a constraint has in its scope.
 *
 *  @note The success pointer indicates if all variables were copied into the vars arrray.
 *
 *  @note It might be that a constraint handler does not support this functionality, in that case the success pointer is
 *        set to FALSE.
 */
SCIP_RETCODE SCIPconsGetVars(
   SCIP_CONS*            cons,               /**< constraint to print */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< array to store the involved variable of the constraint */
   int                   varssize,           /**< available slots in vars array which is needed to check if the array is large enough */
   SCIP_Bool*            success             /**< pointer to store whether the variables are successfully copied */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   if( conshdlr->consgetvars != NULL )
   {
      SCIP_CALL( conshdlr->consgetvars(set->scip, conshdlr, cons, vars, varssize, success) );
   }
   else
   {
      (*success) = FALSE;
   }

   return SCIP_OKAY;
}

/** methed to collect the number of variables of a constraint
 *
 *  @note The success pointer indicates if the contraint handler was able to return the number of variables
 *
 *  @note It might be that a constraint handler does not support this functionality, in that case the success pointer is
 *        set to FALSE
 */
SCIP_RETCODE SCIPconsGetNVars(
   SCIP_CONS*            cons,               /**< constraint to print */
   SCIP_SET*             set,                /**< global SCIP settings */
   int*                  nvars,              /**< pointer to store the number of variables */
   SCIP_Bool*            success             /**< pointer to store whether the constraint successfully returned the number of variables */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   if( conshdlr->consgetnvars != NULL )
   {
      SCIP_CALL( conshdlr->consgetnvars(set->scip, conshdlr, cons, nvars, success) );
   }
   else
   {
      (*nvars) = 0;
      (*success) = FALSE;
   }

   return SCIP_OKAY;
}

/** globally removes constraint from all subproblems; removes constraint from the constraint set change data of the
 *  node, where it was created, or from the problem, if it was a problem constraint
 */
SCIP_RETCODE SCIPconsDelete(
   SCIP_CONS*            cons,               /**< constraint to delete */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(!cons->active || cons->updatedeactivate || cons->addarraypos >= 0);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   SCIPsetDebugMsg(set, "globally deleting constraint <%s> (delay updates: %d)\n",
      cons->name, cons->conshdlr->delayupdatecount);

   /* mark constraint deleted */
   cons->deleted = TRUE;

   /* deactivate constraint, if it is currently active */
   if( cons->active && !cons->updatedeactivate )
   {
      SCIP_CALL( SCIPconsDeactivate(cons, set, stat) );
   }
   else
      cons->updateactivate = FALSE;

   if( set->reopt_enable && !SCIPreoptConsCanBeDeleted(reopt, cons) )
      return SCIP_OKAY;

   assert(!cons->active || cons->updatedeactivate);
   assert(!cons->enabled || cons->updatedeactivate);

   /* remove formerly active constraint from the conssetchg's addedconss / prob's conss array */
   if( cons->addarraypos >= 0 )
   {
      if( cons->addconssetchg == NULL )
      {
         /* remove problem constraint from the problem */
         SCIP_CALL( SCIPprobDelCons(prob, blkmem, set, stat, cons) );
      }
      else
      {
         assert(cons->addconssetchg->addedconss != NULL);
         assert(0 <= cons->addarraypos && cons->addarraypos < cons->addconssetchg->naddedconss);
         assert(cons->addconssetchg->addedconss[cons->addarraypos] == cons);

         /* remove constraint from the constraint set change addedconss array */
         SCIP_CALL( conssetchgDelAddedCons(cons->addconssetchg, blkmem, set, cons->addarraypos) );
      }
   }

   return SCIP_OKAY;
}

/** gets and captures transformed constraint of a given original constraint; if the constraint is not yet transformed,
 *  a new transformed constraint for this constraint is created
 */
SCIP_RETCODE SCIPconsTransform(
   SCIP_CONS*            origcons,           /**< original constraint */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS**           transcons           /**< pointer to store the transformed constraint */
   )
{
   assert(origcons != NULL);
   assert(set != NULL);
   assert(origcons->scip == set->scip);
   assert(origcons->conshdlr != NULL);
   assert(origcons->original);
   assert(transcons != NULL);

   /* check, if the constraint is already transformed */
   if( origcons->transorigcons != NULL )
   {
      *transcons = origcons->transorigcons;
      SCIPconsCapture(*transcons);
   }
   else
   {
      /* create transformed constraint */
      if( origcons->conshdlr->constrans != NULL )
      {
         /* use constraint handler's own method to transform constraint */
         SCIP_CALL( origcons->conshdlr->constrans(set->scip, origcons->conshdlr, origcons, transcons) );
      }
      else
      {
         /* create new constraint with a pointer copy of the constraint data */
         SCIP_CALL( SCIPconsCreate(transcons, blkmem, set, origcons->name, origcons->conshdlr, origcons->consdata, origcons->initial,
               origcons->separate, origcons->enforce, origcons->check, origcons->propagate, 
               origcons->local, origcons->modifiable, origcons->dynamic, origcons->removable, origcons->stickingatnode,
               FALSE, FALSE) );
      }

      /* link original and transformed constraint */
      origcons->transorigcons = *transcons;
      (*transcons)->transorigcons = origcons;

      /* copy the number of upgradelocks */
      (*transcons)->nupgradelocks = origcons->nupgradelocks; /*lint !e732*/
   }
   assert(*transcons != NULL);

   return SCIP_OKAY;
}

/** sets the initial flag of the given constraint */
SCIP_RETCODE SCIPconsSetInitial(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_Bool             initial             /**< new value */
   )
{
   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( cons->initial != initial )
   {
      cons->initial = initial;
      if( !cons->original )
      {
         if( cons->initial )
         {
            SCIP_CALL( conshdlrAddInitcons(SCIPconsGetHdlr(cons), set, stat, cons) );
         }
         else
         {
            if( cons->initconsspos >= 0 )
            {
               conshdlrDelInitcons(SCIPconsGetHdlr(cons), cons);
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** sets the separate flag of the given constraint */
SCIP_RETCODE SCIPconsSetSeparated(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             separate            /**< new value */
   )
{
   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( cons->separate != separate )
   {
      if( SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM )
      {
         cons->separate = separate;
      }
      else if( cons->enabled && cons->sepaenabled )
      {
         if( separate )
         {
            cons->separate = separate;
            SCIP_CALL( conshdlrAddSepacons(cons->conshdlr, set, cons) );
         }
         else
         {
            conshdlrDelSepacons(cons->conshdlr, cons);
            cons->separate = separate;
         }
      }
   }

   return SCIP_OKAY;
}

/** sets the enforce flag of the given constraint */
SCIP_RETCODE SCIPconsSetEnforced(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             enforce             /**< new value */
   )
{
   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( cons->enforce != enforce )
   {
      if( SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM )
      {
         cons->enforce = enforce;
      }
      else if( cons->enabled )
      {
         if( enforce )
         {
            cons->enforce = enforce;
            SCIP_CALL( conshdlrAddEnfocons(cons->conshdlr, set, cons) );
         }
         else
         {
            conshdlrDelEnfocons(cons->conshdlr, cons);
            cons->enforce = enforce;
         }
      }
   }

   return SCIP_OKAY;
}

/** sets the check flag of the given constraint */
SCIP_RETCODE SCIPconsSetChecked(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             check               /**< new value */
   )
{
   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( cons->check != check )
   {
      cons->check = check;

      if( !cons->original )
      {
         /* if constraint is a problem constraint, update variable roundings locks */
         if( cons->addconssetchg == NULL && cons->addarraypos >= 0 )
         {
            if( cons->check )
            {
               SCIP_CALL( SCIPconsAddLocks(cons, set, +1, 0) );
            }
            else
            {
               SCIP_CALL( SCIPconsAddLocks(cons, set, -1, 0) );
            }
         }

         /* if constraint is active, update the checkconss array of the constraint handler */
         if( cons->active )
         {
            if( cons->check )
            {
               SCIP_CALL( conshdlrAddCheckcons(cons->conshdlr, set, cons) );
            }
            else
            {
               conshdlrDelCheckcons(cons->conshdlr, cons);
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** sets the propagate flag of the given constraint */
SCIP_RETCODE SCIPconsSetPropagated(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             propagate           /**< new value */
   )
{
   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( cons->propagate != propagate )
   {
      if( SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM )
      {
         cons->propagate = propagate;
      }
      else if( cons->enabled && cons->propenabled )
      {
         if( propagate )
         {
            cons->propagate = propagate;
            SCIP_CALL( conshdlrAddPropcons(cons->conshdlr, set, cons) );
         }
         else
         {
            conshdlrDelPropcons(cons->conshdlr, cons);
            cons->propagate = propagate;
         }
      }
   }

   return SCIP_OKAY;
}

/** sets the local flag of the given constraint */
void SCIPconsSetLocal(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             local               /**< new value */
   )
{
   assert(cons != NULL);

   cons->local = local;
   if( !local )
      cons->validdepth = 0;
}

/** sets the modifiable flag of the given constraint */
void SCIPconsSetModifiable(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             modifiable          /**< new value */
   )
{
   assert(cons != NULL);

   cons->modifiable = modifiable;
}

/** sets the dynamic flag of the given constraint */
void SCIPconsSetDynamic(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             dynamic             /**< new value */
   )
{
   assert(cons != NULL);

   cons->dynamic = dynamic;
}

/** sets the removable flag of the given constraint */
void SCIPconsSetRemovable(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             removable           /**< new value */
   )
{
   assert(cons != NULL);

   cons->removable = removable;
}

/** sets the stickingatnode flag of the given constraint */
void SCIPconsSetStickingAtNode(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             stickingatnode      /**< new value */
   )
{
   assert(cons != NULL);

   cons->stickingatnode = stickingatnode;
}

/** gives the constraint a new name; ATTENTION: to old pointer is over written that might
 *  result in a memory leakage */
void SCIPconsSetNamePointer(
   SCIP_CONS*            cons,               /**< constraint */
   const char*           name                /**< new name of constraint */
   )
{
   assert( cons != NULL );
   assert( name != NULL );

   cons->name = (char*)name;
}

/** gets associated transformed constraint of an original constraint, or NULL if no associated transformed constraint
 *  exists
 */
SCIP_CONS* SCIPconsGetTransformed(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons->original);

   return cons->transorigcons;
}

/** activates constraint or marks constraint to be activated in next update */
SCIP_RETCODE SCIPconsActivate(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth in the tree where the constraint activation takes place, or -1 for global problem */
   SCIP_Bool             focusnode           /**< does the constraint activation take place at the focus node? */
   )
{
   assert(cons != NULL);
   assert(!cons->original);
   assert(!cons->active);
   assert(!cons->updateactivate);
   assert(!cons->updatedeactivate);
   assert(!cons->updateenable);
   assert(!cons->updatedisable);
   assert(!cons->updateobsolete);
   assert(!cons->updatefree);
   assert(cons->activedepth == -2);
   assert(cons->conshdlr != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
   {
      SCIPsetDebugMsg(set, "delayed activation of constraint <%s> in constraint handler <%s> (depth %d)\n",
         cons->name, cons->conshdlr->name, depth);
      cons->updateactivate = TRUE;
      cons->activedepth = depth;
      cons->updateactfocus = focusnode;
      SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      SCIP_CALL( conshdlrActivateCons(cons->conshdlr, set, stat, cons, depth, focusnode) );
      assert(cons->active);
   }

   return SCIP_OKAY;
}

/** deactivates constraint or marks constraint to be deactivated in next update */
SCIP_RETCODE SCIPconsDeactivate(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(cons != NULL);
   assert(!cons->original);
   assert(cons->active);
   assert(!cons->updateactivate);
   assert(!cons->updatedeactivate);
   assert(cons->activedepth >= -1);
   assert(cons->conshdlr != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
   {
      SCIPsetDebugMsg(set, "delayed deactivation of constraint <%s> in constraint handler <%s>\n",
         cons->name, cons->conshdlr->name);
      cons->updatedeactivate = TRUE;
      cons->activedepth = -2;
      SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      SCIP_CALL( conshdlrDeactivateCons(cons->conshdlr, set, stat, cons) );
      assert(!cons->active);
   }

   return SCIP_OKAY;
}

/** enables constraint's separation, enforcing, and propagation capabilities or marks them to be enabled in next update */
SCIP_RETCODE SCIPconsEnable(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(cons != NULL);
   assert(!cons->original);
   assert(cons->conshdlr != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( !cons->active || cons->updatedeactivate || cons->updateenable || (cons->enabled && !cons->updatedisable) )
      return SCIP_OKAY;

   assert(!cons->updateactivate);

   if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
   {
      cons->updateenable = TRUE;
      SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      SCIP_CALL( conshdlrEnableCons(cons->conshdlr, set, stat, cons) );
      assert(cons->enabled);
   }

   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities or marks them to be disabled in next update */
SCIP_RETCODE SCIPconsDisable(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(cons != NULL);
   assert(!cons->original);
   assert(cons->conshdlr != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( cons->updatedisable || (!cons->enabled && !cons->updateenable) )
      return SCIP_OKAY;

   assert(cons->active);
   assert(!cons->updateactivate);

   if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
   {
      cons->updatedisable = TRUE;
      SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      SCIP_CALL( conshdlrDisableCons(cons->conshdlr, set, stat, cons) );
      assert(!cons->enabled);
   }

   return SCIP_OKAY;
}

/** enables constraint's separation capabilities or marks them to be enabled in next update */
SCIP_RETCODE SCIPconsEnableSeparation(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( cons->updatesepaenable || (cons->sepaenabled && !cons->updatesepadisable) )
      return SCIP_OKAY;

   if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
   {
      cons->updatesepadisable = FALSE;
      cons->updatesepaenable = TRUE;
      SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      SCIP_CALL( conshdlrEnableConsSeparation(cons->conshdlr, set, cons) );
      assert(cons->sepaenabled);
   }

   return SCIP_OKAY;
}

/** disables constraint's separation capabilities or marks them to be disabled in next update */
SCIP_RETCODE SCIPconsDisableSeparation(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);

   if( cons->updatesepadisable || (!cons->sepaenabled && !cons->updatesepaenable) )
      return SCIP_OKAY;

   if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
   {
      cons->updatesepaenable = FALSE;
      cons->updatesepadisable = TRUE;
      SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      SCIP_CALL( conshdlrDisableConsSeparation(cons->conshdlr, cons) );
      assert(!cons->sepaenabled);
   }

   return SCIP_OKAY;
}

/** enables constraint's propagation capabilities or marks them to be enabled in next update */
SCIP_RETCODE SCIPconsEnablePropagation(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( cons->updatepropenable || (cons->propenabled && !cons->updatepropdisable) )
      return SCIP_OKAY;

   if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
   {
      cons->updatepropdisable = FALSE;
      cons->updatepropenable = TRUE;
      SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      SCIP_CALL( conshdlrEnableConsPropagation(cons->conshdlr, set, cons) );
      assert(cons->propenabled);
   }

   return SCIP_OKAY;
}

/** disables constraint's propagation capabilities or marks them to be disabled in next update */
SCIP_RETCODE SCIPconsDisablePropagation(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( cons->updatepropdisable || (!cons->propenabled && !cons->updatepropenable) )
      return SCIP_OKAY;

   if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
   {
      cons->updatepropenable = FALSE;
      cons->updatepropdisable = TRUE;
      SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      SCIP_CALL( conshdlrDisableConsPropagation(cons->conshdlr, cons) );
      assert(!cons->propenabled);
   }

   return SCIP_OKAY;
}

/** marks the constraint to be a conflict */
void SCIPconsMarkConflict(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   cons->conflict = TRUE;
}

/** marks the constraint to be propagated (update might be delayed) */
SCIP_RETCODE SCIPconsMarkPropagate(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( cons->updatemarkpropagate || (cons->markpropagate && !cons->updateunmarkpropagate) )
         return SCIP_OKAY;

   if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
   {
      cons->updateunmarkpropagate = FALSE;
      cons->updatemarkpropagate = TRUE;
      SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      conshdlrMarkConsPropagate(cons->conshdlr, cons);
      assert(cons->markpropagate || !cons->enabled);
   }

   return SCIP_OKAY;
}

/** unmarks the constraint to be propagated (update might be delayed) */
SCIP_RETCODE SCIPconsUnmarkPropagate(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   if( cons->updateunmarkpropagate || (!cons->markpropagate && !cons->updatemarkpropagate) )
      return SCIP_OKAY;

   if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
   {
      cons->updatemarkpropagate = FALSE;
      cons->updateunmarkpropagate = TRUE;
      SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      conshdlrUnmarkConsPropagate(cons->conshdlr, cons);
      assert(!cons->markpropagate || !cons->enabled);
   }

   return SCIP_OKAY;

}

/** adds given value to age of constraint, but age can never become negative;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *  if it's age exceeds the constraint age limit, makes constraint obsolete or marks constraint to be made obsolete
 *  in next update
 */
SCIP_RETCODE SCIPconsAddAge(
   SCIP_CONS*            cons,               /**< constraint */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             deltaage,           /**< value to add to the constraint's age */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(!cons->updateactivate);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   /* no aging in presolving */
   if( set->stage == SCIP_STAGE_PRESOLVING )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "adding %g to age (%g) of constraint <%s> of handler <%s>\n",
      deltaage, cons->age, cons->name, cons->conshdlr->name);

   cons->age += deltaage;
   cons->age = MAX(cons->age, 0.0);

   if( !cons->original )
   {
      if( !cons->check && consExceedsAgelimit(cons, set) )
      {
         SCIP_CALL( SCIPconsDelete(cons, blkmem, set, stat, prob, reopt) );
      }
      else if( !cons->obsolete && consExceedsObsoleteage(cons, set) )
      {
         if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
         {
            cons->updateobsolete = TRUE;
            SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
            assert(cons->update);
         }
         else
         {
            SCIP_CALL( conshdlrMarkConsObsolete(cons->conshdlr, cons) );
            assert(cons->obsolete);
         }
      }
   }

   return SCIP_OKAY;
}

/** increases age of constraint by 1.0;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *  if it's age exceeds the constraint age limit, makes constraint obsolete or marks constraint to be made obsolete
 *  in next update
 */
SCIP_RETCODE SCIPconsIncAge(
   SCIP_CONS*            cons,               /**< constraint */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   )
{
   SCIP_CALL( SCIPconsAddAge(cons, blkmem, set, stat, prob, 1.0, reopt) );

   return SCIP_OKAY;
}

/** resets age of constraint to zero;
 *  should be called
 *   - in constraint separation, if a cut was found for this constraint,
 *   - in constraint enforcing, if the constraint was violated, and
 *   - in constraint propagation, if a domain reduction was deduced;
 *  if it was obsolete, makes constraint useful again or marks constraint to be made useful again in next update
 */
SCIP_RETCODE SCIPconsResetAge(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(!cons->updateactivate);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   SCIPsetDebugMsg(set, "resetting age %g of constraint <%s> of handler <%s>\n", cons->age, cons->name, cons->conshdlr->name);

   conshdlrUpdateAgeresetavg(cons->conshdlr, cons->age);
   cons->age = 0.0;

   if( cons->obsolete )
   {
      assert(!cons->original);
      if( conshdlrAreUpdatesDelayed(cons->conshdlr) )
      {
         cons->updateobsolete = TRUE;
         SCIP_CALL( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
         assert(cons->update);
      }
      else
      {
         SCIP_CALL( conshdlrMarkConsUseful(cons->conshdlr, cons) );
         assert(!cons->obsolete);
      }
   }

   return SCIP_OKAY;
}

/** resolves the given conflicting bound, that was deduced by the given constraint, by putting all "reason" bounds
 *  leading to the deduction into the conflict queue with calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or SCIPaddConflictBinvar();
 *
 *  @note it is sufficient to explain the relaxed bound change
 */
SCIP_RETCODE SCIPconsResolvePropagation(
   SCIP_CONS*            cons,               /**< constraint that deduced the assignment */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             infervar,           /**< variable whose bound was deduced by the constraint */
   int                   inferinfo,          /**< user inference information attached to the bound change */
   SCIP_BOUNDTYPE        inferboundtype,     /**< bound that was deduced (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index, representing the point of time where change took place */
   SCIP_Real             relaxedbd,          /**< the relaxed bound */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(set != NULL);
   assert(cons != NULL);
   assert((inferboundtype == SCIP_BOUNDTYPE_LOWER
         && SCIPgetVarLbAtIndex(set->scip, infervar, bdchgidx, TRUE) > SCIPvarGetLbGlobal(infervar))
      || (inferboundtype == SCIP_BOUNDTYPE_UPPER
         && SCIPgetVarUbAtIndex(set->scip, infervar, bdchgidx, TRUE) < SCIPvarGetUbGlobal(infervar)));
   assert(result != NULL);
   assert(cons->scip == set->scip);

   *result = SCIP_DIDNOTRUN;

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   if( conshdlr->consresprop != NULL )
   {
      /* start timing */
      SCIPclockStart(conshdlr->resproptime, set);

      SCIP_CALL( conshdlr->consresprop(set->scip, conshdlr, cons, infervar, inferinfo, inferboundtype, bdchgidx,
            relaxedbd, result) );

      /* stop timing */
      SCIPclockStop(conshdlr->resproptime, set);

      /* update statistics */
      conshdlr->nrespropcalls++;

      /* check result code */
      if( *result != SCIP_SUCCESS && *result != SCIP_DIDNOTFIND )
      {
         SCIPerrorMessage("propagation conflict resolving method of constraint handler <%s> returned invalid result <%d>\n", 
            conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }
   else
   {
      SCIPerrorMessage("propagation conflict resolving method of constraint handler <%s> is not implemented\n", 
         conshdlr->name);
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** adds given values to lock status of the constraint and updates the rounding locks of the involved variables */
SCIP_RETCODE SCIPconsAddLocks(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nlockspos,          /**< increase in number of rounding locks for constraint */
   int                   nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   int oldnlockspos;
   int oldnlocksneg;
   int updlockpos;
   int updlockneg;

   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(cons->conshdlr->conslock != NULL);
   assert(cons->nlockspos >= 0);
   assert(cons->nlocksneg >= 0);
   assert(-2 <= nlockspos && nlockspos <= 2);
   assert(-2 <= nlocksneg && nlocksneg <= 2);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   /* update the rounding locks */
   oldnlockspos = cons->nlockspos;
   oldnlocksneg = cons->nlocksneg;
   cons->nlockspos += nlockspos;
   cons->nlocksneg += nlocksneg;
   assert(cons->nlockspos >= 0);
   assert(cons->nlocksneg >= 0);

   /* check, if the constraint switched from unlocked to locked, or from locked to unlocked */
   updlockpos = (int)(cons->nlockspos > 0) - (int)(oldnlockspos > 0);
   updlockneg = (int)(cons->nlocksneg > 0) - (int)(oldnlocksneg > 0);

   /* lock the variables, if the constraint switched from unlocked to locked or from locked to unlocked */
   if( updlockpos != 0 || updlockneg != 0 )
   {
      SCIP_CALL( cons->conshdlr->conslock(set->scip, cons->conshdlr, cons, updlockpos, updlockneg) );
   }

   return SCIP_OKAY;
}

/** checks single constraint for feasibility of the given solution */
SCIP_RETCODE SCIPconsCheck(
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             printreason,        /**< Should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);
   assert(result != NULL);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   assert(conshdlr->conscheck != NULL);

   SCIP_CALL( conshdlr->conscheck(set->scip, conshdlr, &cons, 1, sol, checkintegrality, checklprows, printreason,
         FALSE, result) );
   SCIPsetDebugMsg(set, " -> checking returned result <%d>\n", *result);

   if( *result != SCIP_INFEASIBLE && *result != SCIP_FEASIBLE )
   {
      SCIPerrorMessage("feasibility check of constraint handler <%s> on constraint <%s> returned invalid result <%d>\n",
         conshdlr->name, cons->name, *result);
      return SCIP_INVALIDRESULT;
   }

   return SCIP_OKAY;
}

/** enforces single constraint for a given pseudo solution */
SCIP_RETCODE SCIPconsEnfops(
   SCIP_CONS*            cons,               /**< constraint to enforce */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_Bool             objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);
   assert(result != NULL);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   assert(conshdlr->consenfops != NULL);

   SCIP_CALL( conshdlr->consenfops(set->scip, conshdlr, &cons, 1, 1, solinfeasible, objinfeasible, result) );
   SCIPsetDebugMsg(set, " -> enfops returned result <%d>\n", *result);

   if( *result != SCIP_CUTOFF
      && *result != SCIP_CONSADDED
      && *result != SCIP_REDUCEDDOM
      && *result != SCIP_BRANCHED
      && *result != SCIP_SOLVELP
      && *result != SCIP_INFEASIBLE
      && *result != SCIP_FEASIBLE
      && *result != SCIP_DIDNOTRUN )
   {
      SCIPerrorMessage("enforcing method of constraint handler <%s> for pseudo solutions returned invalid result <%d>\n",
         conshdlr->name, *result);
      return SCIP_INVALIDRESULT;
   }

   /* do not update statistics */

   return SCIP_OKAY;
}

/** enforces single constraint for a given LP solution */
SCIP_RETCODE SCIPconsEnfolp(
   SCIP_CONS*            cons,               /**< constraint to enforce */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);
   assert(result != NULL);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   assert(conshdlr->consenfolp != NULL);

   SCIP_CALL( conshdlr->consenfolp(set->scip, conshdlr, &cons, 1, 1, solinfeasible, result) );
   SCIPsetDebugMsg(set, " -> enfolp returned result <%d>\n", *result);

   if( *result != SCIP_CUTOFF
      && *result != SCIP_CONSADDED
      && *result != SCIP_REDUCEDDOM
      && *result != SCIP_BRANCHED
      && *result != SCIP_SEPARATED
      && *result != SCIP_INFEASIBLE
      && *result != SCIP_FEASIBLE)
   {
      SCIPerrorMessage("enforcing method of constraint handler <%s> for LP returned invalid result <%d>\n",
         conshdlr->name, *result);
      return SCIP_INVALIDRESULT;
   }

   /* do not update statistics */

   return SCIP_OKAY;
}

/** enforces single constraint for a given relaxation solution */
SCIP_RETCODE SCIPconsEnforelax(
   SCIP_CONS*            cons,               /**< constraint to enforce */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< solution to be enforced */
   SCIP_Bool             solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);
   assert(sol != NULL);
   assert(result != NULL);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   assert(conshdlr->consenfolp != NULL);

   SCIP_CALL( conshdlr->consenforelax(set->scip, sol, conshdlr, &cons, 1, 1, solinfeasible, result) );
   SCIPdebugMessage(" -> enforelax returned result <%d>\n", *result);

   if( *result != SCIP_CUTOFF
      && *result != SCIP_CONSADDED
      && *result != SCIP_REDUCEDDOM
      && *result != SCIP_BRANCHED
      && *result != SCIP_SEPARATED
      && *result != SCIP_INFEASIBLE
      && *result != SCIP_FEASIBLE)
   {
      SCIPerrorMessage("enforcing method of constraint handler <%s> for relaxation returned invalid result <%d>\n",
         conshdlr->name, *result);
      return SCIP_INVALIDRESULT;
   }

   /* do not update statistics */

   return SCIP_OKAY;
}

/** calls LP initialization method for single constraint */
SCIP_RETCODE SCIPconsInitlp(
   SCIP_CONS*            cons,               /**< constraint to initialize */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool*            infeasible          /**< pointer to store whether infeasibility was detected while building the LP */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(infeasible != NULL);
   assert(cons->scip == set->scip);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   if( conshdlr->consinitlp != NULL )
   {
      SCIP_CALL( conshdlr->consinitlp(set->scip, conshdlr, &cons, 1, infeasible) );
   }

   return SCIP_OKAY;
}

/** calls separation method of single constraint for LP solution */
SCIP_RETCODE SCIPconsSepalp(
   SCIP_CONS*            cons,               /**< constraint to separate */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);
   assert(result != NULL);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   if( conshdlr->conssepalp != NULL )
   {
      SCIP_CALL( conshdlr->conssepalp(set->scip, conshdlr, &cons, 1, 1, result) );
      SCIPsetDebugMsg(set, " -> sepalp returned result <%d>\n", *result);

      if( *result != SCIP_CUTOFF
         && *result != SCIP_CONSADDED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_SEPARATED
         && *result != SCIP_NEWROUND
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN
         && *result != SCIP_DELAYED )
      {
         SCIPerrorMessage("separation method of constraint handler <%s> returned invalid result <%d>\n", conshdlr->name,
            *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** calls separation method of single constraint for given primal solution */
SCIP_RETCODE SCIPconsSepasol(
   SCIP_CONS*            cons,               /**< constraint to separate */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);
   assert(sol != NULL);
   assert(result != NULL);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   if( conshdlr->conssepasol != NULL )
   {
      SCIP_CALL( conshdlr->conssepasol(set->scip, conshdlr, &cons, 1, 1, sol, result) );
      SCIPsetDebugMsg(set, " -> sepasol returned result <%d>\n", *result);

      if( *result != SCIP_CUTOFF
         && *result != SCIP_CONSADDED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_SEPARATED
         && *result != SCIP_NEWROUND
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN
         && *result != SCIP_DELAYED )
      {
         SCIPerrorMessage("separation method of constraint handler for arbitrary primal solution <%s> returned invalid result <%d>\n",
            conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** calls domain propagation method of single constraint */
SCIP_RETCODE SCIPconsProp(
   SCIP_CONS*            cons,               /**< constraint to propagate */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROPTIMING       proptiming,         /**< current point in the node solving loop */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);
   assert(result != NULL);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   if( conshdlr->consprop != NULL )
   {
      SCIP_CALL( conshdlr->consprop(set->scip, conshdlr, &cons, 1, 1, 1, proptiming, result) );
      SCIPsetDebugMsg(set, " -> prop returned result <%d>\n", *result);

      if( *result != SCIP_CUTOFF
         && *result != SCIP_CONSADDED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN
         && *result != SCIP_DELAYED )
      {
         SCIPerrorMessage("propagation method of constraint handler <%s> returned invalid result <%d>\n",
            conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** resolves propagation conflict of single constraint */
SCIP_RETCODE SCIPconsResprop(
   SCIP_CONS*            cons,               /**< constraint to resolve conflict for */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             infervar,           /**< the conflict variable whose bound change has to be resolved */
   int                   inferinfo,          /**< the user information passed to the corresponding SCIPinferVarLbCons() or SCIPinferVarUbCons() call */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Real             relaxedbd,          /**< the relaxed bound which is sufficient to be explained */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);
   assert(result != NULL);
   assert(infervar != NULL);
   assert(bdchgidx != NULL);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   if( conshdlr->consresprop != NULL )
   {
      SCIP_CALL( conshdlr->consresprop(set->scip, conshdlr, cons, infervar, inferinfo, boundtype, bdchgidx, relaxedbd, result) );
      SCIPsetDebugMsg(set, " -> resprop returned result <%d>\n", *result);

      if( *result != SCIP_SUCCESS
         && *result != SCIP_DIDNOTFIND )
      {
         SCIPerrorMessage("propagation conflict resolving method of constraint handler <%s> returned invalid result <%d>\n",
            conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** presolves single constraint */
SCIP_RETCODE SCIPconsPresol(
   SCIP_CONS*            cons,               /**< constraint to presolve */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nrounds,            /**< number of presolving rounds already done */
   SCIP_PRESOLTIMING     timing,             /**< current presolving timing */
   int                   nnewfixedvars,      /**< number of variables fixed since the last call to the presolving method */
   int                   nnewaggrvars,       /**< number of variables aggregated since the last call to the presolving method */
   int                   nnewchgvartypes,    /**< number of variable type changes since the last call to the presolving method */
   int                   nnewchgbds,         /**< number of variable bounds tightened since the last call to the presolving method */
   int                   nnewholes,          /**< number of domain holes added since the last call to the presolving method */
   int                   nnewdelconss,       /**< number of deleted constraints since the last call to the presolving method */
   int                   nnewaddconss,       /**< number of added constraints since the last call to the presolving method */
   int                   nnewupgdconss,      /**< number of upgraded constraints since the last call to the presolving method */
   int                   nnewchgcoefs,       /**< number of changed coefficients since the last call to the presolving method */
   int                   nnewchgsides,       /**< number of changed left or right hand sides since the last call to the presolving method */
   int*                  nfixedvars,         /**< pointer to count total number of variables fixed of all presolvers */
   int*                  naggrvars,          /**< pointer to count total number of variables aggregated of all presolvers */
   int*                  nchgvartypes,       /**< pointer to count total number of variable type changes of all presolvers */
   int*                  nchgbds,            /**< pointer to count total number of variable bounds tightened of all presolvers */
   int*                  naddholes,          /**< pointer to count total number of domain holes added of all presolvers */
   int*                  ndelconss,          /**< pointer to count total number of deleted constraints of all presolvers */
   int*                  naddconss,          /**< pointer to count total number of added constraints of all presolvers */
   int*                  nupgdconss,         /**< pointer to count total number of upgraded constraints of all presolvers */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients of all presolvers */
   int*                  nchgsides,          /**< pointer to count total number of changed left/right hand sides of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nchgvartypes != NULL);
   assert(nchgbds != NULL);
   assert(naddholes != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);
   assert(nupgdconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(result != NULL);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   if( conshdlr->conspresol != NULL )
   {
      SCIP_CALL( conshdlr->conspresol(set->scip, conshdlr, &cons, 1, nrounds,  timing,
            nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,  nnewdelconss, nnewaddconss,
            nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars, nchgvartypes,
            nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );
      SCIPsetDebugMsg(set, " -> presol returned result <%d>\n", *result);

      if( *result != SCIP_UNBOUNDED
         && *result != SCIP_CUTOFF
         && *result != SCIP_SUCCESS
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN
         && *result != SCIP_DELAYED )
      {
         SCIPerrorMessage("presolving method of constraint handler <%s> returned invalid result <%d>\n",
            conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** calls constraint activation notification method of single constraint */
SCIP_RETCODE SCIPconsActive(
   SCIP_CONS*            cons,               /**< constraint to notify */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   if( conshdlr->consactive != NULL )
   {
      SCIP_CALL( conshdlr->consactive(set->scip, conshdlr, cons) );
   }

   return SCIP_OKAY;
}

/** calls constraint deactivation notification method of single constraint */
SCIP_RETCODE SCIPconsDeactive(
   SCIP_CONS*            cons,               /**< constraint to notify */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);
   assert(cons->scip == set->scip);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   /* call external method */
   if( conshdlr->consdeactive != NULL )
   {
      SCIP_CALL( conshdlr->consdeactive(set->scip, conshdlr, cons) );
   }

   return SCIP_OKAY;
}



/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given constraint */
SCIP_DECL_HASHGETKEY(SCIPhashGetKeyCons)
{  /*lint --e{715}*/
   SCIP_CONS* cons = (SCIP_CONS*)elem;

   assert(cons != NULL);
   return cons->name;
}


/*
 * method for arrays of contraint handlers
 */

/** ensures size of storage for propagable constraints with a minimum size of num */
static
SCIP_RETCODE ensurePropagationStorage(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(set != NULL);
   assert(conshdlr != NULL);

   if( num > conshdlr->storedpropconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&(conshdlr->storedpropconss), newsize) );

      conshdlr->storedpropconsssize = newsize;
   }
   assert(num <= conshdlr->storedpropconsssize);

   return SCIP_OKAY;
}

/** stores all constraints marked for propagation away when probing is started */
SCIP_RETCODE SCIPconshdlrsStorePropagationStatus(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONSHDLR**       conshdlrs,          /**< all constraint handlers */
   int                   nconshdlrs          /**< number of contraint handlers */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int c;

   assert(set != NULL);
   assert(conshdlrs != NULL || nconshdlrs == 0);

   for( c = nconshdlrs - 1; c >= 0; --c )
   {
      conshdlr = conshdlrs[c]; /*lint !e613*/
      assert(conshdlr != NULL);
      assert(conshdlr->storednmarkedpropconss == 0);

      if( conshdlr->nmarkedpropconss > 0 )
      {
         int v;

         SCIP_CALL( ensurePropagationStorage(set, conshdlr, conshdlr->nmarkedpropconss) );
         BMScopyMemoryArray(conshdlr->storedpropconss, conshdlr->propconss, conshdlr->nmarkedpropconss);

         conshdlr->storednmarkedpropconss = conshdlr->nmarkedpropconss;
         conshdlr->storedpropdomchgcount = conshdlr->lastpropdomchgcount;

         for( v = conshdlr->storednmarkedpropconss - 1; v >= 0; --v )
         {
            SCIPconsCapture(conshdlr->storedpropconss[v]);
         }
         /* assert(conshdlr->nmarkedpropconss == 0); this assert does not hold if updates are delayed */
      }
   }

   return SCIP_OKAY;
}

/** reset all constraints marked for propagation when probing was finished */
SCIP_RETCODE SCIPconshdlrsResetPropagationStatus(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONSHDLR**       conshdlrs,          /**< all constraint handlers */
   int                   nconshdlrs          /**< number of contraint handlers */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int c;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(conshdlrs != NULL || nconshdlrs == 0);

   for( c = nconshdlrs - 1; c >= 0; --c )
   {
      conshdlr = conshdlrs[c]; /*lint !e613*/
      assert(conshdlr != NULL);

      if( conshdlr->storednmarkedpropconss > 0 )
      {
#ifndef NDEBUG
         int ndisabled = 0;
#endif
         int v;

         for( v = conshdlr->nmarkedpropconss - 1; v >= 0; --v )
         {
            SCIP_CALL( SCIPconsUnmarkPropagate(conshdlr->propconss[v], set) );
         }

         /* mark all previously marked constraint, which were marked before probing */
         for( v = 0; v < conshdlr->storednmarkedpropconss; ++v )
         {
            SCIP_CONS* cons = conshdlr->storedpropconss[v];
            assert(cons != NULL);

            if( cons->enabled && cons->propagate && cons->propenabled )
            {
               SCIP_CALL( SCIPconsMarkPropagate(cons, set) );
            }
#ifndef NDEBUG
            else
               ++ndisabled;
#endif
            SCIP_CALL( SCIPconsRelease(&cons, blkmem, set) );
         }

         assert(conshdlr->storednmarkedpropconss - ndisabled <= conshdlr->npropconss);
         assert(conshdlr->nmarkedpropconss + ndisabled >= conshdlr->storednmarkedpropconss || (conshdlrAreUpdatesDelayed(conshdlr) && conshdlr->nupdateconss + ndisabled >= conshdlr->storednmarkedpropconss));

         conshdlr->lastpropdomchgcount = conshdlr->storedpropdomchgcount;
         conshdlr->storednmarkedpropconss = 0;
      }
   }

   return SCIP_OKAY;
}

/** create linear constraint statistics */
SCIP_RETCODE SCIPlinConsStatsCreate(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_LINCONSSTATS**   linconsstats        /**< pointer to linear constraint classification statistics */
   )
{
   assert(linconsstats != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, linconsstats) );

   return SCIP_OKAY;
}

/** free linear constraint statistics */
void SCIPlinConsStatsFree(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_LINCONSSTATS**   linconsstats        /**< pointer to linear constraint classification statistics */
   )
{
   assert(linconsstats != NULL);
   assert(*linconsstats != NULL);

   SCIPfreeBlockMemory(scip, linconsstats);
}

/** resets linear constraint statistics */
void SCIPlinConsStatsReset(
   SCIP_LINCONSSTATS*    linconsstats        /**< linear constraint classification statistics */
   )
{
   BMSclearMemoryArray(linconsstats->counter, SCIP_NLINCONSTYPES);
   linconsstats->sum = 0;
}

/** returns the number of occurrences of a specific type of linear constraint */
int SCIPlinConsStatsGetTypeCount(
   SCIP_LINCONSSTATS*    linconsstats,       /**< linear constraint classification statistics */
   SCIP_LINCONSTYPE      linconstype         /**< linear constraint type */
   )
{
   assert(linconsstats != NULL);
   assert((int)linconstype < SCIP_NLINCONSTYPES);

   return linconsstats->counter[(int)linconstype];
}

/** returns the total number of classified constraints */
int SCIPlinConsStatsGetSum(
   SCIP_LINCONSSTATS*    linconsstats        /**< linear constraint classification statistics */
   )
{
   assert(linconsstats != NULL);

   return linconsstats->sum;
}

/** increases the number of occurrences of a specific type of linear constraint */
void SCIPlinConsStatsIncTypeCount(
   SCIP_LINCONSSTATS*    linconsstats,       /**< linear constraint classification statistics */
   SCIP_LINCONSTYPE      linconstype,        /**< linear constraint type */
   int                   increment           /**< positive increment */
   )
{
   assert(linconsstats != NULL);
   assert(increment >= 1);
   assert((int)linconstype < SCIP_NLINCONSTYPES);

   linconsstats->counter[(int)linconstype] += increment;
   linconsstats->sum += increment;
}

/** print linear constraint classification statistics */
void SCIPprintLinConsStats(
   SCIP*                 scip,               /**< scip data structure */
   FILE*                 file,               /**< file handle or NULL to print to standard out */
   SCIP_LINCONSSTATS*    linconsstats        /**< linear constraint classification statistics */
   )
{
   assert(scip != NULL);
   assert(linconsstats != NULL);

   /* print statistics */
   SCIPinfoMessage(scip, file, "\n");
   SCIPinfoMessage(scip, file, "%-19s : %10s\n", "Linear cons types", "count");
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "total", SCIPlinConsStatsGetSum(linconsstats));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "empty",        SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_EMPTY));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "free",         SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_FREE));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "singleton",    SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_SINGLETON));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "aggregation",  SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_AGGREGATION));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "precedence",   SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_PRECEDENCE));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "varbound",     SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_VARBOUND));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "setpartition", SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_SETPARTITION));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "setpacking",   SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_SETPACKING));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "setcovering",  SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_SETCOVERING));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "cardinality",  SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_CARDINALITY));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "invknapsack",  SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_INVKNAPSACK));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "eqknapsack",   SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_EQKNAPSACK));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "binpacking",   SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_BINPACKING));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "knapsack",     SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_KNAPSACK));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "intknapsack",  SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_INTKNAPSACK));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "mixedbinary",  SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_MIXEDBINARY));
   SCIPinfoMessage(scip, file, "  %-17s : %10d\n", "general",      SCIPlinConsStatsGetTypeCount(linconsstats, SCIP_LINCONSTYPE_GENERAL));
   SCIPinfoMessage(scip, file, "\n");
}

/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPconsGetName
#undef SCIPconsGetPos
#undef SCIPconsGetHdlr
#undef SCIPconsGetData
#undef SCIPconsGetNUses
#undef SCIPconsGetActiveDepth
#undef SCIPconsGetValidDepth
#undef SCIPconsIsActive
#undef SCIPconsIsEnabled
#undef SCIPconsIsSeparationEnabled
#undef SCIPconsIsPropagationEnabled
#undef SCIPconsIsDeleted
#undef SCIPconsIsObsolete
#undef SCIPconsIsConflict
#undef SCIPconsGetAge
#undef SCIPconsIsInitial
#undef SCIPconsIsSeparated
#undef SCIPconsIsEnforced
#undef SCIPconsIsChecked
#undef SCIPconsIsMarkedPropagate
#undef SCIPconsIsPropagated
#undef SCIPconsIsGlobal
#undef SCIPconsIsLocal
#undef SCIPconsIsModifiable
#undef SCIPconsIsDynamic
#undef SCIPconsIsRemovable
#undef SCIPconsIsStickingAtNode
#undef SCIPconsIsInProb
#undef SCIPconsIsOriginal
#undef SCIPconsIsTransformed
#undef SCIPconsIsLockedPos
#undef SCIPconsIsLockedNeg
#undef SCIPconsIsLocked
#undef SCIPconsGetNLocksPos
#undef SCIPconsGetNLocksNeg
#undef SCIPconsIsAdded
#undef SCIPconsGetNUpgradeLocks

/** returns the name of the constraint 
 *
 *  @note to change the name of a constraint, use SCIPchgConsName() from scip.h
 */
const char* SCIPconsGetName(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->name;
}

/** returns the position of constraint in the corresponding handler's conss array */
int SCIPconsGetPos(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->consspos;
}

/** returns the constraint handler of the constraint */
SCIP_CONSHDLR* SCIPconsGetHdlr(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->conshdlr;
}

/** returns the constraint data field of the constraint */
SCIP_CONSDATA* SCIPconsGetData(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->consdata;
}

/** gets number of times, the constraint is currently captured */
int SCIPconsGetNUses(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->nuses;
}

/** for an active constraint, returns the depth in the tree at which the constraint was activated */
int SCIPconsGetActiveDepth(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsIsActive(cons));

   return cons->activedepth;
}

/** returns TRUE iff constraint is active in the current node */
SCIP_Bool SCIPconsIsActive(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->updateactivate || (cons->active && !cons->updatedeactivate);
}

/** returns TRUE iff constraint is active in the current node */
SCIP_Bool SCIPconsIsUpdatedeactivate(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->updatedeactivate;
}

/** returns the depth in the tree at which the constraint is valid; returns INT_MAX, if the constraint is local
 *  and currently not active
 */
int SCIPconsGetValidDepth(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(cons->validdepth == 0 || cons->local);

   return (!cons->local ? 0
      : !SCIPconsIsActive(cons) ? INT_MAX
      : cons->validdepth == -1 ? SCIPconsGetActiveDepth(cons)
      : cons->validdepth);
}

/** returns TRUE iff constraint is enabled in the current node */
SCIP_Bool SCIPconsIsEnabled(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->updateenable || (cons->enabled && !cons->updatedisable);
}

/** returns TRUE iff constraint's separation is enabled in the current node */
SCIP_Bool SCIPconsIsSeparationEnabled(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return SCIPconsIsEnabled(cons)
      && (cons->updatesepaenable || (cons->sepaenabled && !cons->updatesepadisable));
}

/** returns TRUE iff constraint's propagation is enabled in the current node */
SCIP_Bool SCIPconsIsPropagationEnabled(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return SCIPconsIsEnabled(cons)
      && (cons->updatepropenable || (cons->propenabled && !cons->updatepropdisable));
}

/** returns TRUE iff constraint is deleted or marked to be deleted */
SCIP_Bool SCIPconsIsDeleted(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->deleted;
}

/** returns TRUE iff constraint is marked obsolete */
SCIP_Bool SCIPconsIsObsolete(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->updateobsolete || cons->obsolete;
}

/** returns TRUE iff constraint is marked as a conflict */
SCIP_Bool SCIPconsIsConflict(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->conflict;
}

/** gets age of constraint */
SCIP_Real SCIPconsGetAge(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->age;
}

/** returns TRUE iff the LP relaxation of constraint should be in the initial LP */
SCIP_Bool SCIPconsIsInitial(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->initial;
}

/** returns TRUE iff constraint should be separated during LP processing */
SCIP_Bool SCIPconsIsSeparated(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->separate;
}

/** returns TRUE iff constraint should be enforced during node processing */
SCIP_Bool SCIPconsIsEnforced(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->enforce;
}

/** returns TRUE iff constraint should be checked for feasibility */
SCIP_Bool SCIPconsIsChecked(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->check;
}

/** returns whether the constraint is marked for propagation */
SCIP_Bool SCIPconsIsMarkedPropagate(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return (cons->updatemarkpropagate || (cons->markpropagate && !cons->updateunmarkpropagate));
}

/** returns TRUE iff constraint should be propagated during node processing */
SCIP_Bool SCIPconsIsPropagated(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->propagate;
}

/** returns TRUE iff constraint is globally valid */
SCIP_Bool SCIPconsIsGlobal(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return !cons->local;
}

/** returns TRUE iff constraint is only locally valid or not added to any (sub)problem */
SCIP_Bool SCIPconsIsLocal(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->local;
}

/** returns TRUE iff constraint is modifiable (subject to column generation) */
SCIP_Bool SCIPconsIsModifiable(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->modifiable;
}

/** returns TRUE iff constraint is subject to aging */
SCIP_Bool SCIPconsIsDynamic(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->dynamic;
}

/** returns TRUE iff constraint's relaxation should be removed from the LP due to aging or cleanup */
SCIP_Bool SCIPconsIsRemovable(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->removable;
}

/** returns TRUE iff constraint's relaxation should be removed from the LP due to aging or cleanup */
SCIP_Bool SCIPconsIsStickingAtNode(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->stickingatnode;
}

/** returns TRUE iff constraint belongs to the global problem */
SCIP_Bool SCIPconsIsInProb(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return (cons->addconssetchg == NULL && cons->addarraypos >= 0);
}

/** returns TRUE iff constraint is belonging to original space */
SCIP_Bool SCIPconsIsOriginal(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->original;
}

/** returns TRUE iff constraint is belonging to transformed space */
SCIP_Bool SCIPconsIsTransformed(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return !cons->original;
}

/** returns TRUE iff roundings for variables in constraint are locked */
SCIP_Bool SCIPconsIsLockedPos(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return (cons->nlockspos > 0);
}

/** returns TRUE iff roundings for variables in constraint's negation are locked */
SCIP_Bool SCIPconsIsLockedNeg(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return (cons->nlocksneg > 0);
}

/** returns TRUE iff roundings for variables in constraint or in constraint's negation are locked */
SCIP_Bool SCIPconsIsLocked(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return (cons->nlockspos > 0 || cons->nlocksneg > 0);
}

/** get number of times the roundings for variables in constraint are locked */
int SCIPconsGetNLocksPos(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->nlockspos;
}

/** get number of times the roundings for variables in constraint's negation are locked */
int SCIPconsGetNLocksNeg(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->nlocksneg;
}

/** returns if the constraint was already added to a SCIP instance */
SCIP_Bool SCIPconsIsAdded(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return (cons->addarraypos >= 0);
}

/** adds locks to (dis-)allow upgrading of constraint */
void SCIPconsAddUpgradeLocks(
   SCIP_CONS*            cons,               /**< constraint to add locks */
   int                   nlocks              /**< number of locks to add */
   )
{
   assert(cons != NULL);

   assert(cons->nupgradelocks < (1 << 28) - nlocks); /*lint !e574*/
   cons->nupgradelocks += (unsigned int) nlocks;
}

/** gets number of locks against upgrading the constraint, 0 means this constraint can be upgraded */
int SCIPconsGetNUpgradeLocks(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->nupgradelocks;
}
