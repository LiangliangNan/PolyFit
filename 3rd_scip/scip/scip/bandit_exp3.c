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

/**@file   bandit_exp3.c
 * @brief  methods for Exp.3 bandit selection
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/bandit_exp3.h"

#define BANDIT_NAME "exp3"
#define NUMTOL 1e-6

/*
 * Data structures
 */

/** implementation specific data of Exp.3 bandit algorithm */
struct SCIP_BanditData
{
   SCIP_Real*            weights;            /**< exponential weight for each arm */
   SCIP_Real             weightsum;          /**< the sum of all weights */
   SCIP_Real             gamma;              /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution */
   SCIP_Real             beta;               /**< gain offset between 0 and 1 at every observation */
};

/*
 * Local methods
 */

/*
 * Callback methods of bandit algorithm
 */

/** callback to free bandit specific data structures */
SCIP_DECL_BANDITFREE(SCIPbanditFreeExp3)
{  /*lint --e{715}*/

   SCIP_BANDITDATA* banditdata;
   int nactions;
   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   BMSfreeBlockMemoryArray(blkmem, &banditdata->weights, nactions);

   BMSfreeBlockMemory(blkmem, &banditdata);

   SCIPbanditSetData(bandit, NULL);

   return SCIP_OKAY;
}

/** selection callback for bandit selector */
SCIP_DECL_BANDITSELECT(SCIPbanditSelectExp3)
{  /*lint --e{715}*/

   SCIP_BANDITDATA* banditdata;
   SCIP_RANDNUMGEN* rng;
   SCIP_Real randnr;
   SCIP_Real psum;
   SCIP_Real gammaoverk;
   SCIP_Real oneminusgamma;
   SCIP_Real* weights;
   SCIP_Real weightsum;
   int i;
   int nactions;

   assert(bandit != NULL);
   assert(selection != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   rng = SCIPbanditGetRandnumgen(bandit);
   assert(rng != NULL);
   nactions = SCIPbanditGetNActions(bandit);


   /* draw a random number between 0 and 1 */
   randnr = SCIPrandomGetReal(rng, 0.0, 1.0);

   /* initialize some local variables to speed up probability computations */
   oneminusgamma = 1 - banditdata->gamma;
   gammaoverk = banditdata->gamma / (SCIP_Real)nactions;
   weightsum = banditdata->weightsum;
   weights = banditdata->weights;
   psum = 0.0;

   /* loop over probability distribution until rand is reached
    * the loop terminates without looking at the last action,
    * which is then selected automatically if the target probability
    * is not reached earlier
    */
   for( i = 0; i < nactions - 1; ++i )
   {
      SCIP_Real prob;

      /* compute the probability for arm i as convex kombination of a uniform distribution and a weighted distribution */
      prob = oneminusgamma * weights[i] / weightsum + gammaoverk;
      psum += prob;

      /* break and select element if target probability is reached */
      if( randnr <= psum )
         break;
   }

   /* select element i, which is the last action in case that the break statement hasn't been reached */
   *selection = i;

   return SCIP_OKAY;
}

/** update callback for bandit algorithm */
SCIP_DECL_BANDITUPDATE(SCIPbanditUpdateExp3)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   SCIP_Real eta;
   SCIP_Real gainestim;
   SCIP_Real beta;
   SCIP_Real weightsum;
   SCIP_Real newweightsum;
   SCIP_Real* weights;
   SCIP_Real oneminusgamma;
   SCIP_Real gammaoverk;
   int nactions;

   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   assert(selection >= 0);
   assert(selection < nactions);

   /* the learning rate eta */
   eta = 1.0 / (SCIP_Real)nactions;

   beta = banditdata->beta;
   oneminusgamma = 1.0 - banditdata->gamma;
   gammaoverk = banditdata->gamma * eta;
   weights = banditdata->weights;
   weightsum = banditdata->weightsum;
   newweightsum = weightsum;

   /* if beta is zero, only the observation for the current arm needs an update */
   if( EPSZ(beta, NUMTOL) )
   {
      SCIP_Real probai;
      probai = oneminusgamma * weights[selection] / weightsum + gammaoverk;

      assert(probai > 0.0);

      gainestim = score / probai;
      newweightsum -= weights[selection];
      weights[selection] *= exp(eta * gainestim);
      newweightsum += weights[selection];
   }
   else
   {
      int j;
      newweightsum = 0.0;

      /* loop over all items and update their weights based on the influence of the beta parameter */
      for( j = 0; j < nactions; ++j )
      {
         SCIP_Real probaj;
         probaj = oneminusgamma * weights[j] / weightsum + gammaoverk;

         assert(probaj > 0.0);

         /* consider the score only for the chosen arm i, use constant beta offset otherwise */
         if( j == selection )
            gainestim = (score + beta) / probaj;
         else
            gainestim = beta / probaj;

         weights[j] *= exp(eta * gainestim);
         newweightsum += weights[j];
      }
   }

   banditdata->weightsum = newweightsum;

   return SCIP_OKAY;
}

/** reset callback for bandit algorithm */
SCIP_DECL_BANDITRESET(SCIPbanditResetExp3)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   SCIP_Real* weights;
   int nactions;
   int i;

   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);
   weights = banditdata->weights;

   assert(nactions > 0);

   banditdata->weightsum = (1.0 + NUMTOL) * (SCIP_Real)nactions;

   /* in case of priorities, weights are normalized to sum up to nactions */
   if( priorities != NULL )
   {
      SCIP_Real normalization;
      SCIP_Real priosum;
      priosum = 0.0;

      /* compute sum of priorities */
      for( i = 0; i < nactions; ++i )
      {
         assert(priorities[i] >= 0);
         priosum += priorities[i];
      }

      /* if there are positive priorities, normalize the weights */
      if( priosum > 0.0 )
      {
         normalization = nactions / priosum;
         for( i = 0; i < nactions; ++i )
            weights[i] = (priorities[i] * normalization) + NUMTOL;
      }
      else
      {
         /* use uniform distribution in case of all priorities being 0.0 */
         for( i = 0; i < nactions; ++i )
            weights[i] = 1.0 + NUMTOL;
      }
   }
   else
   {
      /* use uniform distribution in case of unspecified priorities */
      for( i = 0; i < nactions; ++i )
         weights[i] = 1.0 + NUMTOL;
   }

   return SCIP_OKAY;
}


/*
 * bandit algorithm specific interface methods
 */

/** direct bandit creation method for the core where no SCIP pointer is available */
SCIP_RETCODE SCIPbanditCreateExp3(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_BANDITVTABLE*    vtable,             /**< virtual function table for callback functions of Exp.3 */
   SCIP_BANDIT**         exp3,               /**< pointer to store bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   SCIP_Real             gammaparam,         /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution */
   SCIP_Real             beta,               /**< gain offset between 0 and 1 at every observation */
   int                   nactions,           /**< the positive number of actions for this bandit algorithm */
   unsigned int          initseed            /**< initial random seed */
   )
{
   SCIP_BANDITDATA* banditdata;

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &banditdata) );
   assert(banditdata != NULL);

   banditdata->gamma = gammaparam;
   banditdata->beta = beta;
   assert(gammaparam >= 0 && gammaparam <= 1);
   assert(beta >= 0 && beta <= 1);

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &banditdata->weights, nactions) );

   SCIP_CALL( SCIPbanditCreate(exp3, vtable, blkmem, bufmem, priorities, nactions, initseed, banditdata) );

   return SCIP_OKAY;
}

/** creates and resets an Exp.3 bandit algorithm using \p scip pointer */
SCIP_RETCODE SCIPcreateBanditExp3(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         exp3,               /**< pointer to store bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   SCIP_Real             gammaparam,         /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution */
   SCIP_Real             beta,               /**< gain offset between 0 and 1 at every observation */
   int                   nactions,           /**< the positive number of actions for this bandit algorithm */
   unsigned int          initseed            /**< initial seed for random number generation */
   )
{
   SCIP_BANDITVTABLE* vtable;

   vtable = SCIPfindBanditvtable(scip, BANDIT_NAME);
   if( vtable == NULL )
   {
      SCIPerrorMessage("Could not find virtual function table for %s bandit algorithm\n", BANDIT_NAME);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbanditCreateExp3(SCIPblkmem(scip), SCIPbuffer(scip), vtable, exp3,
         priorities, gammaparam, beta, nactions, SCIPinitializeRandomSeed(scip, (int)(initseed % INT_MAX))) );

   return SCIP_OKAY;
}

/** set gamma parameter of Exp.3 bandit algorithm to increase weight of uniform distribution */
void SCIPsetGammaExp3(
   SCIP_BANDIT*          exp3,               /**< bandit algorithm */
   SCIP_Real             gammaparam          /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution */
   )
{
   SCIP_BANDITDATA* banditdata = SCIPbanditGetData(exp3);

   assert(gammaparam >= 0 && gammaparam <= 1);

   banditdata->gamma = gammaparam;
}

/** set beta parameter of Exp.3 bandit algorithm to increase gain offset for actions that were not played */
void SCIPsetBetaExp3(
   SCIP_BANDIT*          exp3,               /**< bandit algorithm */
   SCIP_Real             beta                /**< gain offset between 0 and 1 at every observation */
   )
{
   SCIP_BANDITDATA* banditdata = SCIPbanditGetData(exp3);

   assert(beta >= 0 && beta <= 1);

   banditdata->beta = beta;
}

/** returns probability to play an action */
SCIP_Real SCIPgetProbabilityExp3(
   SCIP_BANDIT*          exp3,               /**< bandit algorithm */
   int                   action              /**< index of the requested action */
   )
{
   SCIP_BANDITDATA* banditdata = SCIPbanditGetData(exp3);

   assert(banditdata->weightsum > 0.0);
   assert(SCIPbanditGetNActions(exp3) > 0);

   return (1.0 - banditdata->gamma) * banditdata->weights[action] / banditdata->weightsum + banditdata->gamma / (SCIP_Real)SCIPbanditGetNActions(exp3);
}

/* include virtual function table for Exp.3 bandit algorithms */
SCIP_RETCODE SCIPincludeBanditvtableExp3(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BANDITVTABLE* vtable;

   SCIP_CALL( SCIPincludeBanditvtable(scip, &vtable, BANDIT_NAME,
         SCIPbanditFreeExp3, SCIPbanditSelectExp3, SCIPbanditUpdateExp3, SCIPbanditResetExp3) );
   assert(vtable != NULL);

   return SCIP_OKAY;
}
