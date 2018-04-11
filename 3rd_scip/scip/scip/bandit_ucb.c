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

/**@file   bandit_ucb.c
 * @brief  methods for UCB bandit selection
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/bandit_ucb.h"
#include "blockmemshell/memory.h"

#define BANDIT_NAME "ucb"
#define NUMEPS 1e-6

/*
 * Data structures
 */

/** implementation specific data of UCB bandit algorithm */
struct SCIP_BanditData
{
   int                   nselections;        /**< counter for the number of selections */
   int*                  counter;            /**< array of counters how often every action has been chosen */
   int*                  startperm;          /**< indices for starting permutation */
   SCIP_Real*            meanscores;         /**< array of average scores for the actions */
   SCIP_Real             alpha;              /**< parameter to increase confidence width */
};


/*
 * Local methods
 */

/** data reset method */
static
SCIP_RETCODE dataReset(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_BANDIT*          ucb,                /**< ucb bandit algorithm */
   SCIP_BANDITDATA*      banditdata,         /**< UCB bandit data structure */
   SCIP_Real*            priorities,         /**< priorities for start permutation, or NULL */
   int                   nactions            /**< number of actions */
   )
{
   int i;
   SCIP_RANDNUMGEN* rng;

   assert(bufmem != NULL);
   assert(ucb != NULL);
   assert(nactions > 0);

   /* clear counters and scores */
   BMSclearMemoryArray(banditdata->counter, nactions);
   BMSclearMemoryArray(banditdata->meanscores, nactions);
   banditdata->nselections = 0;

   rng = SCIPbanditGetRandnumgen(ucb);
   assert(rng != NULL);

   /* initialize start permutation as identity */
   for( i = 0; i < nactions; ++i )
      banditdata->startperm[i] = i;

   /* prepare the start permutation in decreasing order of priority */
   if( priorities != NULL )
   {
      SCIP_Real* prioritycopy;

      SCIP_ALLOC( BMSduplicateBufferMemoryArray(bufmem, &prioritycopy, priorities, nactions) );

      /* randomly wiggle priorities a little bit to make them unique */
      for( i = 0; i < nactions; ++i )
         prioritycopy[i] += SCIPrandomGetReal(rng, -NUMEPS, NUMEPS);

      SCIPsortDownRealInt(prioritycopy, banditdata->startperm, nactions);

      BMSfreeBufferMemoryArray(bufmem, &prioritycopy);
   }
   else
   {
      /* use a random start permutation */
      SCIPrandomPermuteIntArray(rng, banditdata->startperm, 0, nactions);
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of bandit algorithm
 */

/** callback to free bandit specific data structures */
SCIP_DECL_BANDITFREE(SCIPbanditFreeUcb)
{  /*lint --e{715}*/

   SCIP_BANDITDATA* banditdata;
   int nactions;
   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   BMSfreeBlockMemoryArray(blkmem, &banditdata->counter, nactions);
   BMSfreeBlockMemoryArray(blkmem, &banditdata->startperm, nactions);
   BMSfreeBlockMemoryArray(blkmem, &banditdata->meanscores, nactions);
   BMSfreeBlockMemory(blkmem, &banditdata);

   SCIPbanditSetData(bandit, NULL);

   return SCIP_OKAY;
}

/** selection callback for bandit selector */
SCIP_DECL_BANDITSELECT(SCIPbanditSelectUcb)
{  /*lint --e{715}*/

   SCIP_BANDITDATA* banditdata;
   int nactions;
   int* counter;

   assert(bandit != NULL);
   assert(selection != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   counter = banditdata->counter;
   /* select the next uninitialized action from the start permutation */
   if( banditdata->nselections < nactions )
   {
      *selection = banditdata->startperm[banditdata->nselections];
      assert(counter[*selection] == 0);
   }
   else
   {
      /* select the action with the highest upper confidence bound */
      SCIP_Real* meanscores;
      SCIP_Real widthfactor;
      SCIP_Real maxucb;
      int i;
      SCIP_RANDNUMGEN* rng = SCIPbanditGetRandnumgen(bandit);
      meanscores = banditdata->meanscores;

      assert(rng != NULL);
      assert(meanscores != NULL);

      /* compute the confidence width factor that is common for all actions */
      /* cppcheck-suppress unpreciseMathCall */
      widthfactor = banditdata->alpha * LOG1P((SCIP_Real)banditdata->nselections);
      widthfactor = sqrt(widthfactor);
      maxucb = -1.0;

      /* loop over the actions and determine the maximum upper confidence bound.
       * The upper confidence bound of an action is the sum of its mean score
       * plus a confidence term that decreases with increasing number of observations of
       * this action.
       */
      for( i = 0; i < nactions; ++i )
      {
         SCIP_Real uppercb;
         SCIP_Real rootcount;
         assert(counter[i] > 0);

         /* compute the upper confidence bound for action i */
         uppercb = meanscores[i];
         rootcount = sqrt((SCIP_Real)counter[i]);
         uppercb += widthfactor / rootcount;
         assert(uppercb > 0);

         /* update maximum, breaking ties uniformly at random */
         if( EPSGT(uppercb, maxucb, NUMEPS) || (EPSEQ(uppercb, maxucb, NUMEPS) && SCIPrandomGetReal(rng, 0.0, 1.0) >= 0.5) )
         {
            maxucb = uppercb;
            *selection = i;
         }
      }
   }

   assert(*selection >= 0);
   assert(*selection < nactions);

   return SCIP_OKAY;
}

/** update callback for bandit algorithm */
SCIP_DECL_BANDITUPDATE(SCIPbanditUpdateUcb)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   SCIP_Real delta;

   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   assert(selection >= 0);
   assert(selection < SCIPbanditGetNActions(bandit));

   /* increase the mean by the incremental formula: A_n = A_n-1 + 1/n (a_n - A_n-1) */
   delta = score - banditdata->meanscores[selection];
   ++banditdata->counter[selection];
   banditdata->meanscores[selection] += delta / (SCIP_Real)banditdata->counter[selection];

   banditdata->nselections++;

   return SCIP_OKAY;
}

/** reset callback for bandit algorithm */
SCIP_DECL_BANDITRESET(SCIPbanditResetUcb)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   int nactions;

   assert(bufmem != NULL);
   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   /* call the data reset for the given priorities */
   SCIP_CALL( dataReset(bufmem, bandit, banditdata, priorities, nactions) );

   return SCIP_OKAY;
}

/*
 * bandit algorithm specific interface methods
 */

/** returns the upper confidence bound of a selected action */
SCIP_Real SCIPgetConfidenceBoundUcb(
   SCIP_BANDIT*          ucb,                /**< UCB bandit algorithm */
   int                   action              /**< index of the queried action */
   )
{
   SCIP_Real uppercb;
   SCIP_BANDITDATA* banditdata;
   int nactions;

   assert(ucb != NULL);
   banditdata = SCIPbanditGetData(ucb);
   nactions = SCIPbanditGetNActions(ucb);
   assert(action < nactions);

   /* since only scores between 0 and 1 are allowed, 1.0 is a sure upper confidence bound */
   if( banditdata->nselections < nactions )
      return 1.0;

   /* the bandit algorithm must have picked every action once */
   assert(banditdata->counter[action] > 0);
   uppercb = banditdata->meanscores[action];

   /* cppcheck-suppress unpreciseMathCall */
   uppercb += sqrt(banditdata->alpha * LOG1P((SCIP_Real)banditdata->nselections) / (SCIP_Real)banditdata->counter[action]);

   return uppercb;
}

/** return start permutation of the UCB bandit algorithm */
int* SCIPgetStartPermutationUcb(
   SCIP_BANDIT*          ucb                 /**< UCB bandit algorithm */
   )
{
   SCIP_BANDITDATA* banditdata = SCIPbanditGetData(ucb);

   assert(banditdata != NULL);

   return banditdata->startperm;
}

/** internal method to create and reset UCB bandit algorithm */
SCIP_RETCODE SCIPbanditCreateUcb(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_BANDITVTABLE*    vtable,             /**< virtual function table for UCB bandit algorithm */
   SCIP_BANDIT**         ucb,                /**< pointer to store bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   SCIP_Real             alpha,              /**< parameter to increase confidence width */
   int                   nactions,           /**< the positive number of actions for this bandit algorithm */
   unsigned int          initseed            /**< initial random seed */
   )
{
   SCIP_BANDITDATA* banditdata;

   if( alpha < 0.0 )
   {
      SCIPerrorMessage("UCB requires nonnegative alpha parameter, have %f\n", alpha);
      return SCIP_INVALIDDATA;
   }

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &banditdata) );
   assert(banditdata != NULL);

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &banditdata->counter, nactions) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &banditdata->startperm, nactions) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &banditdata->meanscores, nactions) );

   banditdata->alpha = alpha;

   SCIP_CALL( SCIPbanditCreate(ucb, vtable, blkmem, bufmem, priorities, nactions, initseed, banditdata) );

   return SCIP_OKAY;
}

/** create and reset UCB bandit algorithm */
SCIP_RETCODE SCIPcreateBanditUcb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         ucb,                /**< pointer to store bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   SCIP_Real             alpha,              /**< parameter to increase confidence width */
   int                   nactions,           /**< the positive number of actions for this bandit algorithm */
   unsigned int          initseed            /**< initial random number seed */
   )
{
   SCIP_BANDITVTABLE* vtable;

   vtable = SCIPfindBanditvtable(scip, BANDIT_NAME);
   if( vtable == NULL )
   {
      SCIPerrorMessage("Could not find virtual function table for %s bandit algorithm\n", BANDIT_NAME);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbanditCreateUcb(SCIPblkmem(scip), SCIPbuffer(scip), vtable, ucb,
         priorities, alpha, nactions, SCIPinitializeRandomSeed(scip, (int)(initseed % INT_MAX))) );

   return SCIP_OKAY;
}

/** include virtual function table for UCB bandit algorithms */
SCIP_RETCODE SCIPincludeBanditvtableUcb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BANDITVTABLE* vtable;

   SCIP_CALL( SCIPincludeBanditvtable(scip, &vtable, BANDIT_NAME,
         SCIPbanditFreeUcb, SCIPbanditSelectUcb, SCIPbanditUpdateUcb, SCIPbanditResetUcb) );
   assert(vtable != NULL);

   return SCIP_OKAY;
}
