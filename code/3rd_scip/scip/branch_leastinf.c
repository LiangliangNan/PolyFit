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

/**@file   branch_leastinf.c
 * @brief  least infeasible LP branching rule
 * @author Tobias Achterberg
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_leastinf.h"

#define BRANCHRULE_NAME          "leastinf"
#define BRANCHRULE_DESC          "least infeasible branching"
#define BRANCHRULE_PRIORITY      50
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

/*
 * Local methods
 */

/** compares the so far best branching candidate with a new candidate and updates best candidate, if new candidate is better */
static
void updateBestCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            bestvar,            /**< best branching candidate */
   SCIP_Real*            bestscore,          /**< score of best branching candidate */
   SCIP_Real*            bestobj,            /**< absolute objective value of best branching candidate */
   SCIP_Real*            bestsol,            /**< proposed branching point of best branching candidate */
   SCIP_VAR*             cand,               /**< branching candidate to consider */
   SCIP_Real             candscore,          /**< scoring of branching candidate */
   SCIP_Real             candsol             /**< proposed branching point of branching candidate */
   )
{
   SCIP_Real obj;

   assert(scip != NULL);
   assert(bestvar != NULL);
   assert(bestscore != NULL);
   assert(bestobj != NULL);
   assert(*bestobj >= 0.0);
   assert(cand != NULL);

   /* a branching variable candidate should either be an active problem variable or a multi-aggregated variable */
   assert(SCIPvarIsActive(SCIPvarGetProbvar(cand)) ||
      SCIPvarGetStatus(SCIPvarGetProbvar(cand)) == SCIP_VARSTATUS_MULTAGGR);

   if( SCIPvarGetStatus(SCIPvarGetProbvar(cand)) == SCIP_VARSTATUS_MULTAGGR )
   {
      /* for a multi-aggregated variable, we call updateBestCandidate function recursively with all variables in the multi-aggregation */
      SCIP_VAR** multvars;
      int nmultvars;
      int i;
      SCIP_Bool success;
      SCIP_Real multvarlb;
      SCIP_Real multvarub;

      cand = SCIPvarGetProbvar(cand);
      multvars = SCIPvarGetMultaggrVars(cand);
      nmultvars = SCIPvarGetMultaggrNVars(cand);

      /* if we have a candidate branching point, then first register only aggregation variables
       * for which we can compute a corresponding branching point too (see also comments below)
       * if this fails, then register all (unfixed) aggregation variables, thereby forgetting about candsol
       */
      success = FALSE;
      if( candsol != SCIP_INVALID ) /*lint !e777*/
      {
         SCIP_Real* multscalars;
         SCIP_Real minact;
         SCIP_Real maxact;
         SCIP_Real aggrvarsol;
         SCIP_Real aggrvarsol1;
         SCIP_Real aggrvarsol2;

         multscalars = SCIPvarGetMultaggrScalars(cand);

         /* for computing the branching point, we need the current bounds of the multi-aggregated variable */
         minact = SCIPcomputeVarLbLocal(scip, cand);
         maxact = SCIPcomputeVarUbLocal(scip, cand);

         for( i = 0; i < nmultvars; ++i )
         {
            /* skip fixed variables */
            multvarlb = SCIPcomputeVarLbLocal(scip, multvars[i]);
            multvarub = SCIPcomputeVarUbLocal(scip, multvars[i]);
            if( SCIPisEQ(scip, multvarlb, multvarub) )
               continue;

            assert(multscalars != NULL);
            assert(multscalars[i] != 0.0);

            /* we cannot ensure that both the upper bound in the left node and the lower bound in the right node
             * will be candsol by a clever choice for the branching point of multvars[i],
             * but we can try to ensure that at least one of them will be at candsol
             */
            if( multscalars[i] > 0.0 )
            {
               /*    cand >= candsol
                * if multvars[i] >= (candsol - (maxact - multscalars[i] * ub(multvars[i]))) / multscalars[i]
                *                 = (candsol - maxact) / multscalars[i] + ub(multvars[i])
                */
               aggrvarsol1 = (candsol - maxact) / multscalars[i] + multvarub;

               /*     cand <= candsol
                * if multvars[i] <= (candsol - (minact - multscalar[i] * lb(multvars[i]))) / multscalars[i]
                *                 = (candsol - minact) / multscalars[i] + lb(multvars[i])
                */
               aggrvarsol2 = (candsol - minact) / multscalars[i] + multvarlb;
            }
            else
            {
               /*    cand >= candsol
                * if multvars[i] <= (candsol - (maxact - multscalars[i] * lb(multvars[i]))) / multscalars[i]
                *                 = (candsol - maxact) / multscalars[i] + lb(multvars[i])
                */
               aggrvarsol2 = (candsol - maxact) / multscalars[i] + multvarlb;

               /*    cand <= candsol
                * if multvars[i] >= (candsol - (minact - multscalar[i] * ub(multvars[i]))) / multscalars[i]
                *                 = (candsol - minact) / multscalars[i] + ub(multvars[i])
                */
               aggrvarsol1 = (candsol - minact) / multscalars[i] + multvarub;
            }

            /* by the above choice, aggrvarsol1 <= ub(multvars[i]) and aggrvarsol2 >= lb(multvars[i])
             * if aggrvarsol1 <= lb(multvars[i]) or aggrvarsol2 >= ub(multvars[i]), then choose the other one
             * if both are out of bounds, then give up
             * if both are inside bounds, then choose the one closer to 0.0 (someone has better idea???)
             */
            if( SCIPisFeasLE(scip, aggrvarsol1, multvarlb) )
            {
               if( SCIPisFeasGE(scip, aggrvarsol2, multvarub) )
                  continue;
               else
                  aggrvarsol = aggrvarsol2;
            }
            else
            {
               if( SCIPisFeasGE(scip, aggrvarsol2, multvarub) )
                  aggrvarsol = aggrvarsol1;
               else
                  aggrvarsol = REALABS(aggrvarsol1) < REALABS(aggrvarsol2) ? aggrvarsol1 : aggrvarsol2;
            }
            success = TRUE;

            updateBestCandidate(scip, bestvar, bestscore, bestobj, bestsol,
                  multvars[i], candscore, aggrvarsol);
         }
      }

      if( !success )
         for( i = 0; i < nmultvars; ++i )
         {
            /* skip fixed variables */
            multvarlb = SCIPcomputeVarLbLocal(scip, multvars[i]);
            multvarub = SCIPcomputeVarUbLocal(scip, multvars[i]);
            if( SCIPisEQ(scip, multvarlb, multvarub) )
               continue;

            updateBestCandidate(scip, bestvar, bestscore, bestobj, bestsol,
               multvars[i], candscore, SCIP_INVALID);
         }

      assert(*bestvar != NULL); /* if all variables were fixed, something is strange */

      return;
   }

   candscore *= SCIPvarGetBranchFactor(cand);
   obj = SCIPvarGetObj(cand);
   obj = REALABS(obj);
   if( SCIPisInfinity(scip, *bestscore)
      || (!SCIPisInfinity(scip, candscore) && 
          (SCIPisLT(scip, candscore, *bestscore) || (SCIPisLE(scip, candscore, *bestscore) && obj > *bestobj))) )
   {
      *bestvar = cand;
      *bestscore = candscore;
      *bestobj = obj;
      *bestsol = candsol;
   }
}

/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyLeastinf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleLeastinf(scip) );

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpLeastinf)
{  /*lint --e{715}*/
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   SCIP_Real infeasibility;
   SCIP_Real score;
   SCIP_Real obj;
   SCIP_Real bestscore;
   SCIP_Real bestobj;
   int bestcand;
   int i;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Execlp method of leastinf branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands, NULL) );
   assert(nlpcands > 0);

   /* search the least infeasible candidate */
   bestscore = SCIP_REAL_MIN;
   bestobj = 0.0;
   bestcand = -1;
   for( i = 0; i < nlpcands; ++i )
   {
      assert(lpcands[i] != NULL);

      infeasibility = lpcandsfrac[i];
      infeasibility = MIN(infeasibility, 1.0-infeasibility);
      score = 1.0 - infeasibility;
      score *= SCIPvarGetBranchFactor(lpcands[i]);
      obj = SCIPvarGetObj(lpcands[i]);
      obj = REALABS(obj);
      if( SCIPisGT(scip, score, bestscore)
         || (SCIPisGE(scip, score, bestscore) && obj > bestobj) )
      {
         bestscore = score;
         bestobj = obj;
         bestcand = i;
      }
   }
   assert(bestcand >= 0);

   SCIPdebugMsg(scip, " -> %d candidates, selected candidate %d: variable <%s> (frac=%g, obj=%g, factor=%g, score=%g)\n",
      nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandsfrac[bestcand], bestobj,
      SCIPvarGetBranchFactor(lpcands[bestcand]), bestscore);

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextLeastinf)
{  /*lint --e{715}*/
   SCIP_VAR** externcands;
   SCIP_Real* externcandssol;
   SCIP_Real* externcandsscore;
   int nexterncands;
   SCIP_VAR* bestcand;
   SCIP_Real bestscore;
   SCIP_Real bestobj;
   SCIP_Real bestsol;
   SCIP_Real brpoint;
   int i;
   SCIP_NODE* downchild;
   SCIP_NODE* eqchild;
   SCIP_NODE* upchild;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Execext method of leastinf branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &externcands, &externcandssol, &externcandsscore, NULL, &nexterncands, NULL, NULL, NULL) );
   assert(nexterncands > 0);

   /* search the least infeasible candidate */
   bestscore = SCIPinfinity(scip);
   bestobj  = 0.0;
   bestcand = NULL;
   bestsol = SCIP_INVALID;
   for( i = 0; i < nexterncands; ++i )
   {
      updateBestCandidate(scip, &bestcand, &bestscore, &bestobj, &bestsol, externcands[i], externcandsscore[i], externcandssol[i]);
   }

   if( bestcand == NULL )
   {
      SCIPerrorMessage("branchExecextLeastinf failed to select a branching variable from %d candidates\n", nexterncands);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   brpoint = SCIPgetBranchingPoint(scip, bestcand, bestsol);

   SCIPdebugMsg(scip, " -> %d candidates, selected variable <%s> (infeas=%g, obj=%g, factor=%g, score=%g), branching point=%g\n",
      nexterncands, SCIPvarGetName(bestcand), bestsol, bestobj,
      SCIPvarGetBranchFactor(bestcand), bestscore, brpoint);

   /* perform the branching */
   SCIP_CALL( SCIPbranchVarVal(scip, bestcand, brpoint, &downchild, &eqchild, &upchild) );

   if( downchild != NULL || eqchild != NULL || upchild != NULL )
   {
      *result = SCIP_BRANCHED;
   }
   else
   {
      /* if there are no children, then variable should have been fixed by SCIPbranchVarVal */
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(bestcand), SCIPvarGetUbLocal(bestcand)));
      *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the least infeasible LP branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleLeastinf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULE* branchrule;

   /* include branching rule */
   branchrule = NULL;
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, NULL) );
   assert(branchrule != NULL);

   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyLeastinf) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpLeastinf) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextLeastinf) );

   return SCIP_OKAY;
}
