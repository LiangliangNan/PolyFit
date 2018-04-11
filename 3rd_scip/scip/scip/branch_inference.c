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

/**@file   branch_inference.c
 * @brief  inference history branching rule
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_inference.h"


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
};

/** evaluate the given candidate with the given score against the currently best know candidate */
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
   if( (*bestscore) < score )
   {
      /* the score of the candidate is better than the currently best know candidate */
      (*bestscore) = score;
      (*bestcand) = cand;
      (*bestbranchpoint) = branchpoint;
      (*bestbranchdir) = branchdir;
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
         (*bestcand) = cand;
         (*bestbranchpoint) = branchpoint;
         (*bestbranchdir) = branchdir;
      }
   }
}

/** evaluate the given candidate with the given score against the currently best know candidate */
static
void evaluateAggrCand(
   SCIP_VAR*             cand,               /**< candidate to be checked */
   SCIP_Real             score,              /**< score of the candidate */
   SCIP_Real             val,                /**< solution value of the candidate */
   SCIP_VAR**            bestcand,           /**< pointer to the currently best candidate */
   SCIP_Real*            bestscore,          /**< pointer to the score of the currently best candidate */
   SCIP_Real*            bestval             /**< pointer to the solution value of the currently best candidate */
   )
{
   /* evaluate the candidate against the currently best candidate */
   if( (*bestscore) < score )
   {
      /* the score of the candidate is better than the currently best know candidate */
      (*bestscore) = score;
      (*bestcand) = cand;
      (*bestval) = val;
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
         (*bestcand) = cand;
         (*bestval) = val;
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
   SCIP*                 scip,               /**< SCIP data structure */
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


/** selects a variable out of the given candidate array and performs the branching */
static
SCIP_RETCODE performBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            cands,              /**< candidate array */
   SCIP_Real*            candsols,           /**< array of candidate solution values, or NULL */
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

   bestbranchpoint = SCIP_UNKNOWN;
   bestbranchdir = SCIP_BRANCHDIR_DOWNWARDS;
   bestvaluescore = 0.0;
   bestvaluecand = NULL;

   assert(ncands > 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* check if the weighted sum between the average inferences and conflict score should be used */
   if( useweightedsum )
   {
      int c;

      bestaggrcand = cands[0];
      bestvaluecand = cands[0];
      assert(cands[0] != NULL);

      if( candsols == NULL )
      {
         bestval = SCIP_UNKNOWN;

         /* get domain value score for the first candidate */
         bestvaluescore = getValueScore(scip, cands[0], conflictweight, cutoffweight, reliablescore, &bestbranchpoint, &bestbranchdir);
         SCIPdebugMsg(scip, "current best value candidate <%s>[%g,%g] %s <%g> (value %g)\n",
            SCIPvarGetName(bestvaluecand), SCIPvarGetLbLocal(bestvaluecand), SCIPvarGetUbLocal(bestvaluecand),
            bestbranchdir == SCIP_BRANCHDIR_DOWNWARDS ? "<=" : ">=", bestbranchpoint, bestvaluescore);
      }
      else
         bestval = candsols[0];

      /* get aggregated score for the first candidate */
      bestaggrscore = getAggrScore(scip, cands[0], conflictweight, inferenceweight, cutoffweight, reliablescore);

      for( c = 1; c < ncands; ++c )
      {
         SCIP_VAR* cand;
         SCIP_Real val;
         SCIP_Real aggrscore;
         SCIP_Real branchpoint;
         SCIP_BRANCHDIR branchdir;

         cand = cands[c];
         assert(cand != NULL);

         if( candsols == NULL )
         {
            SCIP_Real valuescore;

            val = SCIP_UNKNOWN;

            /* get domain value score for the candidate */
            valuescore = getValueScore(scip, cand, conflictweight, cutoffweight, reliablescore, &branchpoint, &branchdir);

            /* evaluate the candidate against the currently best candidate w.r.t. domain value score */
            evaluateValueCand(cand, valuescore, branchpoint, branchdir, &bestvaluecand, &bestvaluescore, &bestbranchpoint, &bestbranchdir);

            SCIPdebugMsg(scip, "current best value candidate <%s>[%g,%g] %s <%g> (value %g)\n",
               SCIPvarGetName(bestvaluecand), SCIPvarGetLbLocal(bestvaluecand), SCIPvarGetUbLocal(bestvaluecand),
               bestbranchdir == SCIP_BRANCHDIR_DOWNWARDS ? "<=" : ">=", bestbranchpoint, bestvaluescore);
         }
         else
            val = candsols[c];

         /* get aggregated score for the candidate */
         aggrscore = getAggrScore(scip, cand, conflictweight, inferenceweight, cutoffweight, reliablescore);

         SCIPdebugMsg(scip, " -> cand <%s>: prio=%d, solval=%g, score=%g\n", SCIPvarGetName(cand), SCIPvarGetBranchPriority(cand),
            val == SCIP_UNKNOWN ? SCIPgetVarSol(scip, cand) : val, aggrscore); /*lint !e777*/

         /* evaluate the candidate against the currently best candidate w.r.t. aggregated score */
         evaluateAggrCand(cand, aggrscore, val, &bestaggrcand, &bestaggrscore, &bestval);
      }
   }
   else
   {
      int c;

      bestaggrcand = cands[0];
      assert(cands[0] != NULL);

      if( candsols != NULL )
         bestval = candsols[0];
      else
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

         if( candsols != NULL )
            val = candsols[c];
         else
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
         evaluateAggrCand(cand, aggrscore, val, &bestaggrcand, &bestaggrscore, &bestval);
      }
   }

   assert(bestaggrcand != NULL);

   SCIPdebugMsg(scip, " -> %d candidates, selected variable <%s>[%g,%g] (prio=%d, solval=%.12f, score=%g, conflict=%g cutoff=%g, inference=%g)\n",
      ncands, SCIPvarGetName(bestaggrcand), SCIPvarGetLbLocal (bestaggrcand), SCIPvarGetUbLocal(bestaggrcand), SCIPvarGetBranchPriority(bestaggrcand),
      bestval == SCIP_UNKNOWN ? SCIPgetVarSol(scip, bestaggrcand) : bestval, bestaggrscore, /*lint !e777*/
      SCIPgetVarConflictScore(scip, bestaggrcand),  SCIPgetVarAvgInferenceCutoffScore(scip, bestaggrcand, cutoffweight),
      SCIPgetVarAvgInferenceScore(scip, bestaggrcand));

   /* perform the branching */
   if( candsols != NULL )
   {
      SCIP_CALL( SCIPbranchVarVal(scip, bestaggrcand, SCIPgetBranchingPoint(scip, bestaggrcand, bestval), &downchild, &eqchild, &upchild) );
   }
   else if( bestbranchpoint == SCIP_UNKNOWN ) /*lint !e777*/
   {
      SCIP_CALL( SCIPbranchVar(scip, bestaggrcand, &downchild, &eqchild, &upchild) );
   }
   else
   {
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
   SCIP_CALL( performBranching(scip, cands, NULL, ncands, branchruledata->conflictweight,
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
   SCIP_CALL( performBranching(scip, cands, candsols, ncands, branchruledata->conflictweight,
         branchruledata->inferenceweight, branchruledata->cutoffweight, branchruledata->reliablescore,
         branchruledata->useweightedsum, result) );

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
   SCIP_CALL( performBranching(scip, cands, NULL, ncands, branchruledata->conflictweight,
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

   return SCIP_OKAY;
}
