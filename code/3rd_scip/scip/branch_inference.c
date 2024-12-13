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

/**@file   branch_inference.c
 * @ingroup DEFPLUGINS_BRANCH
 * @brief  inference history branching rule
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/branch_inference.h"
#include "scip/pub_branch.h"
#include "scip/pub_history.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_message.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_var.h"
#include <string.h>


/**@name Branching rule properties
 *
 * @{
 */

#define BRANCHRULE_NAME          "inference"
#define BRANCHRULE_DESC          "inference history branching"
#define BRANCHRULE_PRIORITY      1000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

/**@} */

/**@name Default parameter values
 *
 * @{
 */

#define DEFAULT_CONFLICTWEIGHT  1000.0  /**< weight in score calculations for conflict score */
#define DEFAULT_CUTOFFWEIGHT       1.0  /**< weight in score calculations for cutoff score */
#define DEFAULT_INFERENCEWEIGHT    1.0  /**< weight in score calculations for inference score */
#define DEFAULT_RELIABLESCORE    0.001  /**< score which is seen to be reliable for a branching decision */
#define DEFAULT_FRACTIONALS        TRUE /**< should branching on LP solution be restricted to the fractional variables? */
#define DEFAULT_USEWEIGHTEDSUM     TRUE /**< should a weighted sum of inference, conflict and cutoff weights be used? */

#define DEFAULT_CONFLICTPRIO         1  /**< priority value for using conflict weights in lex. order */
#define DEFAULT_CUTOFFPRIO           1  /**< priority value for using cutoff weights in lex. order */

/**@} */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Real             conflictweight;     /**< weight in score calculations for conflict score */
   SCIP_Real             cutoffweight;       /**< weight in score calculations for cutoff score */
   SCIP_Real             inferenceweight;    /**< weight in score calculations for inference score */
   SCIP_Real             reliablescore;      /**< score which is seen to be reliable for a branching decision */
   SCIP_Bool             fractionals;        /**< should branching on LP solution be restricted to the fractional variables? */
   SCIP_Bool             useweightedsum;     /**< should a weighted sum of inference, conflict and cutoff weights be used? */
   int                   conflictprio;       /**< priority value for using conflict weights in lexicographic order */
   int                   cutoffprio;         /**< priority value for using cutoff weights in lexicographic order */
};

/** evaluate the given candidate with the given score against the currently best know candidate, tiebreaking included */
static
void evaluateValueCand(
   SCIP_VAR*             cand,               /**< candidate to be checked */
   SCIP_Real             score,              /**< score of the candidate */
   SCIP_Real             branchpoint,        /**< potential branching point */
   SCIP_BRANCHDIR        branchdir,          /**< potential branching direction */
   SCIP_VAR**            bestcand,           /**< pointer to the currently best candidate */
   SCIP_Real*            bestscore,          /**< pointer to the score of the currently best candidate */
   SCIP_Real*            bestbranchpoint,    /**< pointer to store the (best) branching point */
   SCIP_BRANCHDIR*       bestbranchdir       /**< pointer to store the branching direction relative to the branching point */
   )
{
   /* evaluate the candidate against the currently best candidate */
   if( *bestscore < score )
   {
      /* the score of the candidate is better than the currently best know candidate */
      *bestscore = score;
      *bestcand = cand;
      *bestbranchpoint = branchpoint;
      *bestbranchdir = branchdir;
   }
   else if( (*bestscore) == score ) /*lint !e777*/
   {
      SCIP_Real bestobj;
      SCIP_Real candobj;

      bestobj = REALABS(SCIPvarGetObj(*bestcand));
      candobj = REALABS(SCIPvarGetObj(cand));

      /* the candidate has the same score as the best known candidate; therefore we use a second and third
       * criteria to detect a unique best candidate;
       *
       * - the second criteria prefers the candidate with a larger absolute value of its objective coefficient
       *   since branching on that variable might trigger further propagation w.r.t. objective function
       * - if the absolute values of the objective coefficient are equal the variable index is used to define a
       *   unique best candidate
       *
       * @note It is very important to select a unique best candidate. Otherwise the solver might vary w.r.t. the
       *       performance to much since the candidate array which is used here (SCIPgetPseudoBranchCands() or
       *       SCIPgetLPBranchCands()) gets dynamically changed during the solution process. In particular,
       *       starting a probing mode might already change the order of these arrays. To be independent of that
       *       the selection should be unique. Otherwise, to selection process can get influenced by starting a
       *       probing or not.
       */
      if( bestobj < candobj || (bestobj == candobj && SCIPvarGetIndex(*bestcand) < SCIPvarGetIndex(cand)) ) /*lint !e777*/
      {
         *bestcand = cand;
         *bestbranchpoint = branchpoint;
         *bestbranchdir = branchdir;
      }
   }
}

/** evaluate the given candidate with the given score against the currently best know candidate */
static
void evaluateAggrCand(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             cand,               /**< candidate to be checked */
   SCIP_Real             score,              /**< score of the candidate */
   SCIP_Real             val,                /**< solution value of the candidate */
   SCIP_VAR**            bestcand,           /**< pointer to the currently best candidate */
   SCIP_Real*            bestscore,          /**< pointer to the score of the currently best candidate */
   SCIP_Real*            bestval,            /**< pointer to the solution value of the currently best candidate */
   SCIP_VAR**            bestcands,          /**< buffer array to return selected candidates */
   int*                  nbestcands          /**< pointer to return number of selected candidates */
   )
{
   /* evaluate the candidate against the currently best candidate */
   /** TODO: consider a weaker comparison of some kind */
   if( *bestscore < score )
   {
      /* the score of the candidate is better than the currently best know candidate, so it should be the first candidate in bestcands and nbestcands should be set to 1*/
      *bestscore = score;
      *bestcand = cand;
      *bestval = val;
      *nbestcands = 1;
      bestcands[0] = cand;
   }
   /** TODO: consider a weaker comparison of some kind */
   else if( SCIPisEQ(scip, *bestscore, score) )
   {
      /* the score of the candidate is comparable to the currently known best, so we add it to bestcands and increase nbestcands by 1*/
      bestcands[*nbestcands] = cand;
      (*nbestcands)++;
   }
}

/** choose a singular best candidate from bestcands and move it to the beginning of the candidate array */
static
void tiebreakAggrCand(
   SCIP_VAR**            bestcands,          /**< buffer array to return selected candidates */
   int                   nbestcands          /**< number of selected candidates */
   )
{
   int c;

   for( c = 0; c < nbestcands; ++c )
   {
      SCIP_Real bestobj;
      SCIP_Real candobj;

      bestobj = REALABS(SCIPvarGetObj(bestcands[0]));
      candobj = REALABS(SCIPvarGetObj(bestcands[c]));

      /* the candidate has the same score as the best known candidate; therefore we use a second and third
       * criteria to detect a unique best candidate;
       *
       * - the second criteria prefers the candidate with a larger absolute value of its objective coefficient
       *   since branching on that variable might trigger further propagation w.r.t. objective function
       * - if the absolute values of the objective coefficient are equal the variable index is used to define a
       *   unique best candidate
       *
       * @note It is very important to select a unique best candidate. Otherwise the solver might vary w.r.t. the
       *       performance too much since the candidate array which is used here (SCIPgetPseudoBranchCands() or
       *       SCIPgetLPBranchCands()) gets dynamically changed during the solution process. In particular,
       *       starting a probing mode might already change the order of these arrays. To be independent of that
       *       the selection should be unique. Otherwise, to selection process can get influenced by starting a
       *       probing or not.
       */
      if( bestobj < candobj || (bestobj == candobj && SCIPvarGetIndex(bestcands[0]) < SCIPvarGetIndex(bestcands[c])) ) /*lint !e777*/
      {
         bestcands[0] = bestcands[c];
      }
   }
}

/** check if the score for the given domain value and variable domain value is better than the current best know one */
static
void checkValueScore(
   SCIP_Real             value,              /**< domain value */
   SCIP_HISTORY*         history,            /**< variable history for given donain value */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_Real             conflictweight,     /**< weight in score calculations for conflict score */
   SCIP_Real             cutoffweight,       /**< weight in score calculations for cutoff score */
   SCIP_Real             reliablescore,      /**< score which is seen to be reliable for a branching decision */
   SCIP_Real*            bestscore,          /**< pointer to store the best score */
   SCIP_Real*            branchpoint,        /**< pointer to store the (best) branching point */
   SCIP_BRANCHDIR*       branchdir           /**< pointer to store the branching direction relative to the branching point */
   )
{
   SCIP_Real conflictscore;
   SCIP_Real cutoffscore;
   SCIP_Real score;

   conflictscore = SCIPhistoryGetVSIDS(history, dir);
   cutoffscore = SCIPhistoryGetCutoffSum(history, dir);

   /* in case the conflict score is below the reliable score we set it to zero since it is seen to be
    * unreliable
    */
   if( conflictscore  < reliablescore )
      conflictscore = 0.0;

   /* in case the cutoff score is below the reliable score we set it to zero since it is seen to be unreliable */
   if( cutoffscore < reliablescore )
      cutoffscore = 0.0;

   /* compute weight score */
   score = conflictweight * conflictscore + cutoffweight * cutoffscore;

   if( score > *bestscore )
   {
      (*bestscore) = score;
      (*branchpoint) = value;
      (*branchdir) = dir;
   }
}

/** return an aggregated score for the given variable using the conflict score and cutoff score */
static
SCIP_Real getAggrScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             conflictweight,     /**< weight in score calculations for conflict score */
   SCIP_Real             inferenceweight,    /**< weight in score calculations for inference score */
   SCIP_Real             cutoffweight,       /**< weight in score calculations for cutoff score */
   SCIP_Real             reliablescore       /**< score which is seen to be reliable for a branching decision */
   )
{
   SCIP_Real conflictscore;
   SCIP_Real cutoffscore;

   conflictscore = SCIPgetVarConflictScore(scip, var);
   cutoffscore = SCIPgetVarAvgInferenceCutoffScore(scip, var, cutoffweight);

   /* in case the conflict score is below the reliable score we set it to zero since it is seen to be
    * unreliable
    */
   if( conflictscore  < reliablescore )
      conflictscore = 0.0;

   /* in case the cutoff score is below the reliable score we set it to zero since it is seen to be unreliable */
   if( cutoffscore < reliablescore )
      cutoffscore = 0.0;

   /* compute weighted score for the candidate */
   return (conflictweight * conflictscore + inferenceweight * cutoffscore);
}

/** return an aggregated score for the given variable using the conflict score and cutoff score */
static
SCIP_Real getValueScore(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             conflictweight,     /**< weight in score calculations for conflict score */
   SCIP_Real             cutoffweight,       /**< weight in score calculations for cutoff score */
   SCIP_Real             reliablescore,      /**< score which is seen to be reliable for a branching decision */
   SCIP_Real*            branchpoint,        /**< pointer to store the branching point */
   SCIP_BRANCHDIR*       branchdir           /**< pointer to store the branching direction relative to the branching point */
   )
{
   SCIP_VALUEHISTORY* valuehistory;
   SCIP_Real bestscore;

   (*branchpoint) = SCIP_UNKNOWN;
   (*branchdir) = SCIP_BRANCHDIR_UPWARDS;

   valuehistory = SCIPvarGetValuehistory(var);
   bestscore = 0.0;

   if( valuehistory != NULL )
   {
      SCIP_HISTORY** histories;
      SCIP_Real* values;
      int nvalues;
      int v;

      histories = SCIPvaluehistoryGetHistories(valuehistory);
      values = SCIPvaluehistoryGetValues(valuehistory);
      nvalues = SCIPvaluehistoryGetNValues(valuehistory);

      for( v = 0; v < nvalues; ++v )
      {
         SCIP_Real value;

         value = values[v];

         /* skip all domain values which are smaller or equal to the lower bound */
         if( value <= SCIPvarGetLbLocal(var) )
            continue;

         /* skip all domain values which are larger or equal to the upper bound */
         if( value >= SCIPvarGetUbLocal(var) )
            break;

         /* check var <= value */
         checkValueScore(value, histories[v], SCIP_BRANCHDIR_DOWNWARDS, conflictweight, cutoffweight, reliablescore, &bestscore, branchpoint, branchdir);

         /* check var >= value */
         checkValueScore(value, histories[v], SCIP_BRANCHDIR_UPWARDS, conflictweight, cutoffweight, reliablescore, &bestscore, branchpoint, branchdir);
      }
   }

   return bestscore;
}

static
void selectBestCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            cands,              /**< candidate array */
   SCIP_Real*            candsols,           /**< array of candidate solution values, or NULL */
   int                   ncands,             /**< number of candidates */
   SCIP_Real             conflictweight,     /**< weight in score calculations for conflict score */
   SCIP_Real             inferenceweight,    /**< weight in score calculations for inference score */
   SCIP_Real             cutoffweight,       /**< weight in score calculations for cutoff score */
   SCIP_Real             reliablescore,      /**< score which is seen to be reliable for a branching decision */
   SCIP_VAR**            bestcands,          /**< buffer array to return selected candidates */
   int*                  nbestcands          /**< pointer to return number of selected candidates */
   )
{
   SCIP_VAR* bestaggrcand;
   SCIP_Real bestval;
   SCIP_Real bestaggrscore;
   int c;

   bestaggrcand = cands[0];
   assert(cands[0] != NULL);

   bestval = candsols[0];
   bestcands[0] = cands[0];
   *nbestcands = 1;

   /* get aggregated score for the first candidate */
   bestaggrscore = getAggrScore(scip, cands[0], conflictweight, inferenceweight, cutoffweight, reliablescore);

   for( c = 1; c < ncands; ++c )
   {
      SCIP_VAR* cand;
      SCIP_Real val;
      SCIP_Real aggrscore;

      cand = cands[c];
      assert(cand != NULL);

      val = candsols[c];

      /* get score for the candidate */
      aggrscore = getAggrScore(scip, cand, conflictweight, inferenceweight, cutoffweight, reliablescore);

      /*lint -e777*/
      SCIPdebugMsg(scip, " -> cand <%s>: prio=%d, solval=%g, score=%g\n", SCIPvarGetName(cand), SCIPvarGetBranchPriority(cand),
         val == SCIP_UNKNOWN ? SCIPgetVarSol(scip, cand) : val, aggrscore);

      /* evaluate the candidate against the currently best candidate w.r.t. aggregated score */
      evaluateAggrCand(scip, cand, aggrscore, val, &bestaggrcand, &bestaggrscore, &bestval, bestcands, nbestcands);
   }
}  /*lint --e{438}*/


/** selects a variable out of the given candidate array and performs the branching */
static
SCIP_RETCODE performBranchingSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            cands,              /**< candidate array */
   SCIP_Real*            candsols,           /**< array of candidate solution values, or NULL */
   int                   ncands,             /**< number of candidates */
   SCIP_Real             conflictweight,     /**< weight in score calculations for conflict score */
   SCIP_Real             inferenceweight,    /**< weight in score calculations for inference score */
   SCIP_Real             cutoffweight,       /**< weight in score calculations for cutoff score */
   SCIP_Real             reliablescore,      /**< score which is seen to be reliable for a branching decision */
   SCIP_Bool             useweightedsum,     /**< should a weighted sum of inference, conflict and cutoff weights be used? */
   SCIP_RESULT*          result,             /**< buffer to store result (branched, reduced domain, ...) */
   int                   conflictprio,       /**< priority value for using conflict weights in lex. order */
   int                   cutoffprio          /**< priority value for using conflict weights in lex. order */
   )
{
   SCIP_VAR* bestaggrcand;
   SCIP_Real bestval;
   SCIP_NODE* downchild;
   SCIP_NODE* eqchild;
   SCIP_NODE* upchild;
   SCIP_VAR** bestcands;
   int nbestcands;
   int c;

   assert(ncands > 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* check if conflict score, inferences, and cutoff score should be used in combination; otherwise just use
    * inference */
   if( useweightedsum == FALSE )
   {
      conflictprio = 0;
      cutoffprio = 0;
      conflictweight = 0.0;
      inferenceweight = 1.0;
      cutoffweight = 0.0;
   }

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &bestcands, ncands) );
   nbestcands = 0;

   if( conflictprio > cutoffprio )
   {
      /* select the best candidates w.r.t. the first criterion */
      selectBestCands(scip, cands, candsols, ncands, conflictweight, 0.0, 0.0, reliablescore,
            bestcands, &nbestcands);

      /* select the best candidates w.r.t. the second criterion; we use bestcands and nbestcands as input and
       * output, so the method must make sure to overwrite the last argument only at the very end */
      if( nbestcands > 1 )
      {
         selectBestCands(scip, bestcands, candsols, nbestcands, 0.0, inferenceweight, cutoffweight, reliablescore,
               bestcands, &nbestcands);
      }
   }
   else if( conflictprio == cutoffprio )
   {
      /* select the best candidates w.r.t. weighted sum of both criteria */
      selectBestCands(scip, cands, candsols, ncands, conflictweight, inferenceweight, cutoffweight, reliablescore,
            bestcands, &nbestcands);
   }
   else
   {
      assert(conflictprio < cutoffprio);

      /* select the best candidates w.r.t. the first criterion */
      selectBestCands(scip, cands, candsols, ncands, 0.0, inferenceweight, cutoffweight, reliablescore,
            bestcands, &nbestcands);

      /* select the best candidates w.r.t. the second criterion; we use bestcands and nbestcands as input and
       * output, so the method must make sure to overwrite the last argument only at the very end */
      if( nbestcands > 1 )
      {
         /* select the best candidates w.r.t. the first criterion */
         selectBestCands(scip, bestcands, candsols, nbestcands, conflictweight, 0.0, 0.0, reliablescore,
               bestcands, &nbestcands);
      }
   }

   assert(nbestcands == 0 || bestcands[0] != NULL);

   /* final tie breaking */
   if( nbestcands > 1 )
   {
      tiebreakAggrCand(bestcands, nbestcands);
      nbestcands = 1;
   }

   assert(nbestcands == 1);

   bestaggrcand = bestcands[0];
   bestval = -SCIP_INVALID;

   /* loop over cands, find bestcands[0], and store corresponding candsols value in bestval */
   for( c = 0; c < ncands; ++c )
   {
      if( bestaggrcand == cands[c] )
      {
         bestval = candsols[c];
         break;
      }
   }

   assert(bestval != -SCIP_INVALID);

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &bestcands);

   assert(bestaggrcand != NULL);

   SCIPdebugMsg(scip, " -> %d candidates, selected variable <%s>[%g,%g] (prio=%d, solval=%.12f, conflict=%g cutoff=%g, inference=%g)\n",
      ncands, SCIPvarGetName(bestaggrcand), SCIPvarGetLbLocal (bestaggrcand), SCIPvarGetUbLocal(bestaggrcand), SCIPvarGetBranchPriority(bestaggrcand),
      bestval == SCIP_UNKNOWN ? SCIPgetVarSol(scip, bestaggrcand) : bestval, /*lint !e777*/
      SCIPgetVarConflictScore(scip, bestaggrcand),  SCIPgetVarAvgInferenceCutoffScore(scip, bestaggrcand, cutoffweight),
      SCIPgetVarAvgInferenceScore(scip, bestaggrcand));

   assert(candsols != NULL);
   /* perform the branching */
   SCIP_CALL( SCIPbranchVarVal(scip, bestaggrcand, SCIPgetBranchingPoint(scip, bestaggrcand, bestval), &downchild, &eqchild, &upchild) );

   if( downchild != NULL || eqchild != NULL || upchild != NULL )
   {
      *result = SCIP_BRANCHED;
   }
   else
   {
      /* if there are no children, then variable should have been fixed by SCIPbranchVar(Val) */
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(bestaggrcand), SCIPvarGetUbLocal(bestaggrcand)));
      *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}


/** selects a variable out of the given candidate array and performs the branching */
static
SCIP_RETCODE performBranchingNoSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            cands,              /**< candidate array */
   int                   ncands,             /**< number of candidates */
   SCIP_Real             conflictweight,     /**< weight in score calculations for conflict score */
   SCIP_Real             inferenceweight,    /**< weight in score calculations for inference score */
   SCIP_Real             cutoffweight,       /**< weight in score calculations for cutoff score */
   SCIP_Real             reliablescore,      /**< score which is seen to be reliable for a branching decision */
   SCIP_Bool             useweightedsum,     /**< should a weighted sum of inference, conflict and cutoff weights be used? */
   SCIP_RESULT*          result              /**< buffer to store result (branched, reduced domain, ...) */
   )
{
   SCIP_VAR* bestaggrcand;
   SCIP_VAR* bestvaluecand;
   SCIP_Real bestval;
   SCIP_Real bestaggrscore;
   SCIP_Real bestvaluescore;
   SCIP_Real bestbranchpoint;
   SCIP_BRANCHDIR bestbranchdir;
   SCIP_NODE* downchild;
   SCIP_NODE* eqchild;
   SCIP_NODE* upchild;
   SCIP_VAR** bestcands;
   int nbestcands;

   bestbranchpoint = SCIP_UNKNOWN;
   bestbranchdir = SCIP_BRANCHDIR_DOWNWARDS;
   bestvaluescore = 0.0;
   bestvaluecand = NULL;

   assert(ncands > 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;


   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &bestcands, ncands) );
   nbestcands = 0;

   /* check if the weighted sum between the average inferences and conflict score should be used */
   if( useweightedsum )
   {
      int c;

      bestaggrcand = cands[0];
      bestvaluecand = cands[0];
      assert(cands[0] != NULL);

      bestval = SCIP_UNKNOWN;

      /* get domain value score for the first candidate */
      bestvaluescore = getValueScore(cands[0], conflictweight, cutoffweight, reliablescore, &bestbranchpoint, &bestbranchdir);
      SCIPdebugMsg(scip, "current best value candidate <%s>[%g,%g] %s <%g> (value %g)\n",
         SCIPvarGetName(bestvaluecand), SCIPvarGetLbLocal(bestvaluecand), SCIPvarGetUbLocal(bestvaluecand),
         bestbranchdir == SCIP_BRANCHDIR_DOWNWARDS ? "<=" : ">=", bestbranchpoint, bestvaluescore);

      /* get aggregated score for the first candidate */
      bestaggrscore = getAggrScore(scip, cands[0], conflictweight, inferenceweight, cutoffweight, reliablescore);

      for( c = 1; c < ncands; ++c )
      {
         SCIP_VAR* cand;
         SCIP_Real val;
         SCIP_Real aggrscore;
         SCIP_Real branchpoint;
         SCIP_BRANCHDIR branchdir;
         SCIP_Real valuescore;

         cand = cands[c];
         assert(cand != NULL);

         val = SCIP_UNKNOWN;

         /* get domain value score for the candidate */
         valuescore = getValueScore(cand, conflictweight, cutoffweight, reliablescore, &branchpoint, &branchdir);

         /* evaluate the candidate against the currently best candidate w.r.t. domain value score */
         evaluateValueCand(cand, valuescore, branchpoint, branchdir, &bestvaluecand, &bestvaluescore, &bestbranchpoint, &bestbranchdir);

         SCIPdebugMsg(scip, "current best value candidate <%s>[%g,%g] %s <%g> (value %g)\n",
            SCIPvarGetName(bestvaluecand), SCIPvarGetLbLocal(bestvaluecand), SCIPvarGetUbLocal(bestvaluecand),
            bestbranchdir == SCIP_BRANCHDIR_DOWNWARDS ? "<=" : ">=", bestbranchpoint, bestvaluescore);

         /* get aggregated score for the candidate */
         aggrscore = getAggrScore(scip, cand, conflictweight, inferenceweight, cutoffweight, reliablescore);

         /*lint -e777*/
         SCIPdebugMsg(scip, " -> cand <%s>: prio=%d, solval=%g, score=%g\n", SCIPvarGetName(cand), SCIPvarGetBranchPriority(cand),
            val == SCIP_UNKNOWN ? SCIPgetVarSol(scip, cand) : val, aggrscore);

         /* evaluate the candidate against the currently best candidate w.r.t. aggregated score */
         evaluateAggrCand(scip, cand, aggrscore, val, &bestaggrcand, &bestaggrscore, &bestval, bestcands, &nbestcands);
      }
   }
   else
   {
      int c;

      bestaggrcand = cands[0];
      assert(cands[0] != NULL);

      bestval = SCIP_UNKNOWN;

      bestaggrscore = SCIPgetVarAvgInferenceScore(scip, cands[0]);

      /* search for variable with best score w.r.t. average inferences per branching */
      for( c = 1; c < ncands; ++c )
      {
         SCIP_VAR* cand;
         SCIP_Real val;
         SCIP_Real aggrscore;

         cand = cands[c];
         assert(cand != NULL);

         val = SCIP_UNKNOWN;

         aggrscore = SCIPgetVarAvgInferenceScore(scip, cand);

         /* in case the average inferences score is below the reliable score we set it to zero since it is seen to be
          * unreliable
          */
         if( aggrscore < reliablescore )
            aggrscore = 0.0;

         SCIPdebugMsg(scip, " -> cand <%s>: prio=%d, solval=%g, score=%g\n", SCIPvarGetName(cand), SCIPvarGetBranchPriority(cand),
            val == SCIP_UNKNOWN ? SCIPgetVarSol(scip, cand) : val, aggrscore); /*lint !e777*/

         /* evaluate the candidate against the currently best candidate */
         evaluateAggrCand(scip, cand, aggrscore, val, &bestaggrcand, &bestaggrscore, &bestval, bestcands, &nbestcands);
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &bestcands);

   assert(bestaggrcand != NULL);

   SCIPdebugMsg(scip, " -> %d candidates, selected variable <%s>[%g,%g] (prio=%d, solval=%.12f, score=%g, conflict=%g cutoff=%g, inference=%g)\n",
      ncands, SCIPvarGetName(bestaggrcand), SCIPvarGetLbLocal (bestaggrcand), SCIPvarGetUbLocal(bestaggrcand), SCIPvarGetBranchPriority(bestaggrcand),
      bestval == SCIP_UNKNOWN ? SCIPgetVarSol(scip, bestaggrcand) : bestval, bestaggrscore, /*lint !e777*/
      SCIPgetVarConflictScore(scip, bestaggrcand),  SCIPgetVarAvgInferenceCutoffScore(scip, bestaggrcand, cutoffweight),
      SCIPgetVarAvgInferenceScore(scip, bestaggrcand));

   if( bestbranchpoint == SCIP_UNKNOWN ) /*lint !e777*/
   {
      SCIP_CALL( SCIPbranchVar(scip, bestaggrcand, &downchild, &eqchild, &upchild) );
   }
   else
   {
      /* perform the branching */
      SCIP_Real estimate;
      SCIP_Real downprio;
      SCIP_Real upprio;
      SCIP_Real downub;
      SCIP_Real uplb;

      assert(bestvaluecand != NULL);

      downprio = 0.0;
      upprio = 0.0;

      if( bestbranchdir == SCIP_BRANCHDIR_DOWNWARDS )
      {
         downprio = 1.0;
         downub = bestbranchpoint;
         uplb = bestbranchpoint + 1.0;
      }
      else
      {
         upprio = 1.0;
         downub = bestbranchpoint - 1.0;
         uplb = bestbranchpoint;
      }

      /* calculate the child estimate */
      estimate = SCIPcalcChildEstimate(scip, bestvaluecand, downub);

      /* create down child */
      SCIP_CALL( SCIPcreateChild(scip, &downchild, downprio, estimate) );

      /* change upper bound in down child */
      SCIP_CALL( SCIPchgVarUbNode(scip, downchild, bestvaluecand, downub) );

      /* calculate the child estimate */
      estimate = SCIPcalcChildEstimate(scip, bestvaluecand, uplb);

      /* create up child */
      SCIP_CALL( SCIPcreateChild(scip, &upchild, upprio, estimate) );

      /* change lower bound in up child */
      SCIP_CALL( SCIPchgVarLbNode(scip, upchild, bestvaluecand, uplb) );

      SCIPdebugMsg(scip, "branch on variable <%s> and value <%g>\n", SCIPvarGetName(bestvaluecand), bestbranchpoint);

      eqchild = NULL;
   }
   if( downchild != NULL || eqchild != NULL || upchild != NULL )
   {
      *result = SCIP_BRANCHED;
   }
   else
   {
      /* if there are no children, then variable should have been fixed by SCIPbranchVar(Val) */
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(bestaggrcand), SCIPvarGetUbLocal(bestaggrcand)));
      *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyInference)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleInference(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeInference)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeBlockMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpInference)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** cands;
   int ncands;

   SCIPdebugMsg(scip, "Execlp method of inference branching\n");

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( branchruledata->fractionals )
   {
      /* get LP candidates (fractional integer variables) */
      SCIP_CALL( SCIPgetLPBranchCands(scip, &cands, NULL, NULL, NULL, &ncands, NULL) );
   }
   else
   {
      /* get pseudo candidates (non-fixed integer variables) */
      SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, NULL, &ncands) );
   }

   /* perform the branching */
   SCIP_CALL( performBranchingNoSol(scip, cands, ncands, branchruledata->conflictweight,
         branchruledata->inferenceweight, branchruledata->cutoffweight, branchruledata->reliablescore,
         branchruledata->useweightedsum, result) );

   return SCIP_OKAY;
}


/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextInference)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** cands;
   SCIP_Real* candsols;
   int ncands;

   SCIPdebugMsg(scip, "Execext method of inference branching\n");

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get branching candidates */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &cands, &candsols, NULL, &ncands, NULL, NULL, NULL, NULL) );
   assert(ncands > 0);

   /* perform the branching */
   SCIP_CALL( performBranchingSol(scip, cands, candsols, ncands, branchruledata->conflictweight,
         branchruledata->inferenceweight, branchruledata->cutoffweight, branchruledata->reliablescore,
         branchruledata->useweightedsum, result, branchruledata->conflictprio, branchruledata->cutoffprio) );

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsInference)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** cands;
   int ncands;

   SCIPdebugMsg(scip, "Execps method of inference branching\n");

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get pseudo candidates (non-fixed integer variables) */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, NULL, &ncands) );

   /* perform the branching */
   SCIP_CALL( performBranchingNoSol(scip, cands, ncands, branchruledata->conflictweight,
         branchruledata->inferenceweight, branchruledata->cutoffweight, branchruledata->reliablescore,
         branchruledata->useweightedsum, result) );

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the inference history branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleInference(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create inference branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non-fundamental callbacks via specific setter functions*/
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyInference) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeInference) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpInference) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextInference) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsInference) );

   /* inference branching rule parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/inference/conflictweight",
         "weight in score calculations for conflict score",
         &branchruledata->conflictweight, TRUE, DEFAULT_CONFLICTWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/inference/inferenceweight",
         "weight in score calculations for inference score",
         &branchruledata->inferenceweight, TRUE, DEFAULT_INFERENCEWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/inference/cutoffweight",
         "weight in score calculations for cutoff score",
         &branchruledata->cutoffweight, TRUE, DEFAULT_CUTOFFWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/inference/fractionals",
         "should branching on LP solution be restricted to the fractional variables?",
         &branchruledata->fractionals, TRUE, DEFAULT_FRACTIONALS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/inference/useweightedsum",
         "should a weighted sum of inference, conflict and cutoff weights be used?",
         &branchruledata->useweightedsum, FALSE, DEFAULT_USEWEIGHTEDSUM, NULL, NULL) );
   /* inference branching rule parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/inference/reliablescore",
         "weight in score calculations for conflict score",
         &branchruledata->reliablescore, TRUE, DEFAULT_RELIABLESCORE, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   /* parameters for lexicographical ordering */
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/inference/conflictprio",
         "priority value for using conflict weights in lex. order",
         &branchruledata->conflictprio, FALSE, DEFAULT_CONFLICTPRIO, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/inference/cutoffprio",
         "priority value for using cutoff weights in lex. order",
         &branchruledata->cutoffprio, FALSE, DEFAULT_CUTOFFPRIO, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
