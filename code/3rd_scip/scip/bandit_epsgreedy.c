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

/**@file   bandit_epsgreedy.c
 * @ingroup OTHER_CFILES
 * @brief  implementation of epsilon greedy bandit algorithm
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/bandit.h"
#include "scip/bandit_epsgreedy.h"
#include "scip/pub_bandit.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/scip_bandit.h"
#include "scip/scip_mem.h"
#include "scip/scip_randnumgen.h"

#define BANDIT_NAME           "eps-greedy"
#define EPSGREEDY_SMALL       1e-6

/*
 * Data structures
 */

/** private data structure of epsilon greedy bandit algorithm */
struct SCIP_BanditData
{
   SCIP_Real*            weights;            /**< weights for every action */
   SCIP_Real*            priorities;         /**< saved priorities for tie breaking */
   int*                  sels;               /**< individual number of selections per action */
   SCIP_Real             eps;                /**< epsilon parameter (between 0 and 1) to control epsilon greedy */
   SCIP_Real             decayfactor;        /**< the factor to reduce the weight of older observations if exponential decay is enabled */
   int                   avglim;             /**< nonnegative limit on observation number before the exponential decay starts,
                                               *  only relevant if exponential decay is enabled
                                               */
   int                   nselections;        /**< counter for the number of selection calls */
   SCIP_Bool             preferrecent;       /**< should the weights be updated in an exponentially decaying way? */
};

/*
 * Callback methods of bandit algorithm virtual function table
 */

/** callback to free bandit specific data structures */
SCIP_DECL_BANDITFREE(SCIPbanditFreeEpsgreedy)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   int nactions;

   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   assert(banditdata->weights != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   BMSfreeBlockMemoryArray(blkmem, &banditdata->weights, nactions);
   BMSfreeBlockMemoryArray(blkmem, &banditdata->priorities, nactions);
   BMSfreeBlockMemoryArray(blkmem, &banditdata->sels, nactions);
   BMSfreeBlockMemory(blkmem, &banditdata);

   SCIPbanditSetData(bandit, NULL);

   return SCIP_OKAY;
}

/** selection callback for bandit algorithm */
SCIP_DECL_BANDITSELECT(SCIPbanditSelectEpsgreedy)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   SCIP_Real randnr;
   SCIP_Real curreps;
   SCIP_RANDNUMGEN* rng;
   int nactions;
   assert(bandit != NULL);
   assert(selection != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   rng = SCIPbanditGetRandnumgen(bandit);
   assert(rng != NULL);

   nactions = SCIPbanditGetNActions(bandit);

   /* roll the dice to check if the best element should be picked, or an element at random */
   randnr = SCIPrandomGetReal(rng, 0.0, 1.0);

   /* make epsilon decrease with an increasing number of selections */
   banditdata->nselections++;
   curreps = banditdata->eps * sqrt((SCIP_Real)nactions/(SCIP_Real)banditdata->nselections);

   /* select the best action seen so far */
   if( randnr >= curreps )
   {
      SCIP_Real* weights = banditdata->weights;
      SCIP_Real* priorities = banditdata->priorities;
      int j;
      SCIP_Real maxweight;

      assert(weights != NULL);
      assert(priorities != NULL);

      /* pick the element with the largest reward */
      maxweight = weights[0];
      *selection = 0;

      /* determine reward for every element */
      for( j = 1; j < nactions; ++j )
      {
         SCIP_Real weight = weights[j];

         /* select the action that maximizes the reward, breaking ties by action priorities */
         if( maxweight < weight
               || (weight >= maxweight - EPSGREEDY_SMALL && priorities[j] > priorities[*selection] ) )
         {
            *selection = j;
            maxweight = weight;
         }
      }
   }
   else
   {
      /* play one of the actions at random */
      *selection = SCIPrandomGetInt(rng, 0, nactions - 1);
   }

   return SCIP_OKAY;
}

/** update callback for bandit algorithm */
SCIP_DECL_BANDITUPDATE(SCIPbanditUpdateEpsgreedy)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;

   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);

   /* increase the selection count */
   ++banditdata->sels[selection];

   /* the very first observation is directly stored as weight for both average or exponential decay */
   if( banditdata->sels[selection] == 1 )
      banditdata->weights[selection] = score;
   else
   {
      /* use exponentially decreasing weights for older observations */
      if( banditdata->preferrecent && banditdata->sels[selection] > banditdata->avglim )
      {
         /* decrease old weights by decay factor */
         banditdata->weights[selection] *= banditdata->decayfactor;
         banditdata->weights[selection] += (1.0 - banditdata->decayfactor) * score;
      }
      else
      {
         /* update average score */
         SCIP_Real diff = score - banditdata->weights[selection];
         banditdata->weights[selection] += diff / (SCIP_Real)(banditdata->sels[selection]);
      }
   }

   return SCIP_OKAY;
}

/** reset callback for bandit algorithm */
SCIP_DECL_BANDITRESET(SCIPbanditResetEpsgreedy)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   SCIP_Real* weights;
   int w;
   int nactions;
   SCIP_RANDNUMGEN* rng;

   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);

   weights = banditdata->weights;
   nactions = SCIPbanditGetNActions(bandit);
   assert(weights != NULL);
   assert(banditdata->priorities != NULL);
   assert(nactions > 0);

   rng = SCIPbanditGetRandnumgen(bandit);
   assert(rng != NULL);

   /* alter priorities slightly to make them unique */
   if( priorities != NULL )
   {
      for( w = 1; w < nactions; ++w )
      {
         assert(priorities[w] >= 0);
         banditdata->priorities[w] = priorities[w] + SCIPrandomGetReal(rng, -EPSGREEDY_SMALL, EPSGREEDY_SMALL);
      }
   }
   else
   {
      /* use random priorities */
      for( w = 0; w < nactions; ++w )
         banditdata->priorities[w] = SCIPrandomGetReal(rng, 0.0, 1.0);
   }

   /* reset weights and selection counters to 0 */
   BMSclearMemoryArray(weights, nactions);
   BMSclearMemoryArray(banditdata->sels, nactions);

   banditdata->nselections = 0;

   return SCIP_OKAY;
}

/*
 * interface methods of the Epsilon Greedy bandit algorithm
 */

/** internal method to create and reset epsilon greedy bandit algorithm */
SCIP_RETCODE SCIPbanditCreateEpsgreedy(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_BANDITVTABLE*    vtable,             /**< virtual function table with epsilon greedy callbacks */
   SCIP_BANDIT**         epsgreedy,          /**< pointer to store the epsilon greedy bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   SCIP_Real             eps,                /**< parameter to increase probability for exploration between all actions */
   SCIP_Bool             preferrecent,       /**< should the weights be updated in an exponentially decaying way? */
   SCIP_Real             decayfactor,        /**< the factor to reduce the weight of older observations if exponential decay is enabled */
   int                   avglim,             /**< nonnegative limit on observation number before the exponential decay starts,
                                              *   only relevant if exponential decay is enabled */
   int                   nactions,           /**< the positive number of possible actions */
   unsigned int          initseed            /**< initial random seed */
   )
{
   SCIP_BANDITDATA* banditdata;

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &banditdata) );
   assert(banditdata != NULL);
   assert(eps >= 0.0);

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &banditdata->weights, nactions) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &banditdata->priorities, nactions) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &banditdata->sels, nactions) );
   banditdata->eps = eps;
   banditdata->nselections = 0;
   banditdata->preferrecent = preferrecent;
   banditdata->decayfactor = decayfactor;
   banditdata->avglim = avglim;

   SCIP_CALL( SCIPbanditCreate(epsgreedy, vtable, blkmem, bufmem, priorities, nactions, initseed, banditdata) );

   return SCIP_OKAY;
}

/** create and resets an epsilon greedy bandit algorithm */
SCIP_RETCODE SCIPcreateBanditEpsgreedy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         epsgreedy,          /**< pointer to store the epsilon greedy bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   SCIP_Real             eps,                /**< parameter to increase probability for exploration between all actions */
   SCIP_Bool             preferrecent,       /**< should the weights be updated in an exponentially decaying way? */
   SCIP_Real             decayfactor,        /**< the factor to reduce the weight of older observations if exponential decay is enabled */
   int                   avglim,             /**< nonnegative limit on observation number before the exponential decay starts,
                                              *   only relevant if exponential decay is enabled */
   int                   nactions,           /**< the positive number of possible actions */
   unsigned int          initseed            /**< initial seed for random number generation */
   )
{
   SCIP_BANDITVTABLE* vtable;
   assert(scip != NULL);
   assert(epsgreedy != NULL);

   vtable = SCIPfindBanditvtable(scip, BANDIT_NAME);
   if( vtable == NULL )
   {
      SCIPerrorMessage("Could not find virtual function table for %s bandit algorithm\n", BANDIT_NAME);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbanditCreateEpsgreedy(SCIPblkmem(scip), SCIPbuffer(scip), vtable, epsgreedy,
         priorities, eps, preferrecent, decayfactor, avglim, nactions, SCIPinitializeRandomSeed(scip, initseed)) );

   return SCIP_OKAY;
}

/** get weights array of epsilon greedy bandit algorithm */
SCIP_Real* SCIPgetWeightsEpsgreedy(
   SCIP_BANDIT*          epsgreedy           /**< epsilon greedy bandit algorithm */
   )
{
   SCIP_BANDITDATA* banditdata;
   assert(epsgreedy != NULL);
   banditdata = SCIPbanditGetData(epsgreedy);
   assert(banditdata != NULL);

   return banditdata->weights;
}

/** set epsilon parameter of epsilon greedy bandit algorithm */
void SCIPsetEpsilonEpsgreedy(
   SCIP_BANDIT*          epsgreedy,          /**< epsilon greedy bandit algorithm */
   SCIP_Real             eps                 /**< parameter to increase probability for exploration between all actions */
   )
{
   SCIP_BANDITDATA* banditdata;
   assert(epsgreedy != NULL);
   assert(eps >= 0);

   banditdata = SCIPbanditGetData(epsgreedy);

   banditdata->eps = eps;
}


/** creates the epsilon greedy bandit algorithm includes it in SCIP */
SCIP_RETCODE SCIPincludeBanditvtableEpsgreedy(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BANDITVTABLE* banditvtable;

   SCIP_CALL( SCIPincludeBanditvtable(scip, &banditvtable, BANDIT_NAME,
         SCIPbanditFreeEpsgreedy, SCIPbanditSelectEpsgreedy, SCIPbanditUpdateEpsgreedy, SCIPbanditResetEpsgreedy) );

   return SCIP_OKAY;
}
