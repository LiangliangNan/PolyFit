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

/**@file   branch_relpscost.c
 * @brief  reliable pseudo costs branching rule
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_relpscost.h"
#include "scip/cons_and.h"
#include "scip/pub_misc.h"

#define BRANCHRULE_NAME          "relpscost"
#define BRANCHRULE_DESC          "reliability branching on pseudo cost values"
#define BRANCHRULE_PRIORITY      10000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define DEFAULT_CONFLICTWEIGHT   0.01        /**< weight in score calculations for conflict score */
#define DEFAULT_CONFLENGTHWEIGHT 0.0         /**< weight in score calculations for conflict length score*/
#define DEFAULT_INFERENCEWEIGHT  0.0001      /**< weight in score calculations for inference score */
#define DEFAULT_CUTOFFWEIGHT     0.0001      /**< weight in score calculations for cutoff score */
#define DEFAULT_PSCOSTWEIGHT     1.0         /**< weight in score calculations for pseudo cost score */
#define DEFAULT_NLSCOREWEIGHT    0.1         /**< weight in score calculations for nlcount score */
#define DEFAULT_MINRELIABLE      1.0         /**< minimal value for minimum pseudo cost size to regard pseudo cost value as reliable */
#define DEFAULT_MAXRELIABLE      5.0         /**< maximal value for minimum pseudo cost size to regard pseudo cost value as reliable */
#define DEFAULT_SBITERQUOT       0.5         /**< maximal fraction of strong branching LP iterations compared to normal iterations */
#define DEFAULT_SBITEROFS        100000      /**< additional number of allowed strong branching LP iterations */
#define DEFAULT_MAXLOOKAHEAD     9           /**< maximal number of further variables evaluated without better score */
#define DEFAULT_INITCAND         100         /**< maximal number of candidates initialized with strong branching per node */
#define DEFAULT_INITITER         0           /**< iteration limit for strong branching initialization of pseudo cost entries (0: auto) */
#define DEFAULT_MAXBDCHGS        5           /**< maximal number of bound tightenings before the node is reevaluated (-1: unlimited) */
#define DEFAULT_MAXPROPROUNDS   -2           /**< maximum number of propagation rounds to be performed during strong branching
                                              *   before solving the LP (-1: no limit, -2: parameter settings) */
#define DEFAULT_PROBINGBOUNDS    TRUE        /**< should valid bounds be identified in a probing-like fashion during strong
                                              *   branching (only with propagation)? */
#define DEFAULT_USERELERRORFORRELIABILITY FALSE /**< should reliability be based on relative errors? */
#define DEFAULT_LOWERRORTOL      0.05        /**< lowest tolerance beneath which relative errors are reliable */
#define DEFAULT_HIGHERRORTOL     1.0         /**< highest tolerance beneath which relative errors are reliable */
#define DEFAULT_USEHYPTESTFORRELIABILITY FALSE /**< should the strong branching decision be based on a hypothesis test? */
#define DEFAULT_USEDYNAMICCONFIDENCE FALSE   /**< should the confidence level be adjusted dynamically? */
#define DEFAULT_STORESEMIINITCOSTS FALSE     /**< should strong branching result be considered for pseudo costs if the other direction was infeasible? */
#define DEFAULT_USESBLOCALINFO   FALSE       /**< should the scoring function use only local cutoff and inference information obtained for strong branching candidates? */
#define DEFAULT_CONFIDENCELEVEL  2           /**< The confidence level for statistical methods, between 0 (Min) and 4 (Max). */
#define DEFAULT_SKIPBADINITCANDS TRUE        /**< should branching rule skip candidates that have a low probability to be
                                              *  better than the best strong-branching or pseudo-candidate? */
#define DEFAULT_STARTRANDSEED    5           /**< start random seed for random number generation */
#define DEFAULT_RANDINITORDER    FALSE       /**< should slight perturbation of scores be used to break ties in the prior scores? */
#define DEFAULT_USESMALLWEIGHTSITLIM FALSE   /**< should smaller weights be used for pseudo cost updates after hitting the LP iteration limit? */
#define DEFAULT_DYNAMICWEIGHTS   TRUE        /**< should the weights of the branching rule be adjusted dynamically during solving based
                                              *   infeasible and objective leaf counters? */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Real             conflictweight;     /**< weight in score calculations for conflict score */
   SCIP_Real             conflengthweight;   /**< weight in score calculations for conflict length score */
   SCIP_Real             inferenceweight;    /**< weight in score calculations for inference score */
   SCIP_Real             cutoffweight;       /**< weight in score calculations for cutoff score */
   SCIP_Real             pscostweight;       /**< weight in score calculations for pseudo cost score */
   SCIP_Real             nlscoreweight;      /**< weight in score calculations for nlcount score */
   SCIP_Real             minreliable;        /**< minimal value for minimum pseudo cost size to regard pseudo cost value as reliable */
   SCIP_Real             maxreliable;        /**< maximal value for minimum pseudo cost size to regard pseudo cost value as reliable */
   SCIP_Real             sbiterquot;         /**< maximal fraction of strong branching LP iterations compared to normal iterations */
   int                   sbiterofs;          /**< additional number of allowed strong branching LP iterations */
   int                   maxlookahead;       /**< maximal number of further variables evaluated without better score */
   int                   initcand;           /**< maximal number of candidates initialized with strong branching per node */
   int                   inititer;           /**< iteration limit for strong branching initialization of pseudo cost entries (0: auto) */
   int                   maxbdchgs;          /**< maximal number of bound tightenings before the node is reevaluated (-1: unlimited) */
   int                   maxproprounds;      /**< maximum number of propagation rounds to be performed during strong branching
                                              *   before solving the LP (-1: no limit, -2: parameter settings) */
   SCIP_Bool             probingbounds;      /**< should valid bounds be identified in a probing-like fashion during strong
                                              *   branching (only with propagation)? */
   SCIP_Bool             userelerrorforreliability; /**< should reliability be based on relative errors? */
   SCIP_Real             lowerrortol;        /**< lowest tolerance beneath which relative errors are reliable */
   SCIP_Real             higherrortol;       /**< highest tolerance beneath which relative errors are reliable */
   SCIP_Bool             usehyptestforreliability; /**< should the strong branching decision be based on a hypothesis test? */
   SCIP_Bool             usedynamicconfidence; /**< should the confidence level be adjusted dynamically? */
   SCIP_Bool             storesemiinitcosts; /**< should strong branching result be considered for pseudo costs if the other direction was infeasible? */
   SCIP_Bool             usesblocalinfo;     /**< should the scoring function disregard cutoffs for variable if sb-lookahead was feasible ? */
   SCIP_Bool             skipbadinitcands;   /**< should branching rule skip candidates that have a low probability to be
                                               *  better than the best strong-branching or pseudo-candidate? */
   SCIP_Bool             dynamicweights;     /**< should the weights of the branching rule be adjusted dynamically during solving based on objective and infeasible leaf counters? */
   int                   confidencelevel;    /**< The confidence level for statistical methods, between 0 (Min) and 4 (Max). */
   int*                  nlcount;            /**< array to store nonlinear count values */
   int                   nlcountsize;        /**< length of nlcount array */
   int                   nlcountmax;         /**< maximum entry in nlcount array or 1 if NULL */
   SCIP_Bool             randinitorder;      /**< should slight perturbation of scores be used to break ties in the prior scores? */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   int                   startrandseed;      /**< start random seed for random number generation */
   SCIP_Bool             usesmallweightsitlim; /**< should smaller weights be used for pseudo cost updates after hitting the LP iteration limit? */
};

/*
 * local methods
 */

/** return probindex of variable or corresponding active variable (if negated or aggregated) or -1 (if multiaggregated) */
static
SCIP_RETCODE binvarGetActiveProbindex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable */
   int*                  probindex           /**< buffer to store probindex */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(SCIPvarIsBinary(var));
   assert(probindex != NULL);

   *probindex = SCIPvarGetProbindex(var);

   /* if variable is not active, try to find active representative */
   if( *probindex == -1 )
   {
      SCIP_VAR* repvar;
      SCIP_Bool negated;

      SCIP_CALL( SCIPgetBinvarRepresentative(scip, var, &repvar, &negated) );
      assert(repvar != NULL);
      assert(SCIPvarGetStatus(repvar) != SCIP_VARSTATUS_FIXED);

      if( SCIPvarIsActive(repvar) )
         *probindex = SCIPvarGetProbindex(repvar);
      else if( SCIPvarIsNegated(repvar) )
         *probindex = SCIPvarGetProbindex(SCIPvarGetNegationVar(repvar));
   }

   return SCIP_OKAY;
}

/** counts number of nonlinear constraints in which each variable appears */
static
SCIP_RETCODE countNonlinearities(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  nlcount,            /**< pointer to array for storing count values */
   int                   nlcountsize,        /**< buffer for storing length of nlcount array */
   int*                  nlcountmax          /**< buffer for storing maximum value in nlcount array */
   )
{
   SCIP_CONSHDLR* andconshdlr;
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(nlcount != NULL);
   assert(nlcountmax != NULL);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(nlcountsize >= nvars);

   /* get nonlinearity for constraints in NLP */
   if( SCIPisNLPConstructed(scip) )
   {
      assert(SCIPgetNNLPVars(scip) == nvars);
      SCIP_CALL( SCIPgetNLPVarsNonlinearity(scip, nlcount) );
   }
   else
   {
      BMSclearMemoryArray(nlcount, nvars);
   }

   /* increase counters for and constraints */
   andconshdlr = SCIPfindConshdlr(scip, "and");
   if( andconshdlr != NULL )
   {
      int c;

      for( c = 0; c < SCIPconshdlrGetNActiveConss(andconshdlr); c++ )
      {
         SCIP_CONS* andcons;
         SCIP_VAR** andvars;
         SCIP_VAR* andres;
         int probindex;
         int nandvars;
         int v;

         /* get constraint and variables */
         andcons = SCIPconshdlrGetConss(andconshdlr)[c];
         nandvars = SCIPgetNVarsAnd(scip, andcons);
         andvars = SCIPgetVarsAnd(scip, andcons);
         andres = SCIPgetResultantAnd(scip, andcons);

         probindex = -1;
         for( v = 0; v < nandvars; v++ )
         {
            /* don't rely on the and conshdlr removing fixed variables
             * @todo fix the and conshdlr in that respect
             */
            if( SCIPvarGetStatus(andvars[v]) != SCIP_VARSTATUS_FIXED )
            {
               SCIP_CALL( binvarGetActiveProbindex(scip, andvars[v], &probindex) );
               if( probindex >= 0 )
                  nlcount[probindex]++;
            }
         }

         SCIP_CALL( binvarGetActiveProbindex(scip, andres, &probindex) );
         if( probindex >= 0 )
            nlcount[probindex]++;
      }
   }

   /* compute maximum count value */
   *nlcountmax = 1;
   for( i = 0; i < nvars; i++ )
   {
      if( *nlcountmax < nlcount[i] )
         *nlcountmax = nlcount[i];
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE branchruledataEnsureNlcount(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   int nvars;

   assert(scip != NULL);
   assert(branchruledata != NULL);

   nvars = SCIPgetNVars(scip);

   /**@todo test whether we want to apply this as if problem has only and constraints */
   /**@todo update changes in and constraints */
   if( branchruledata->nlscoreweight > 0.0 ) /*  && SCIPisNLPConstructed(scip) */
   {
      if( branchruledata->nlcount == NULL )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->nlcount, nvars) );
         branchruledata->nlcountsize = nvars;

         SCIP_CALL( countNonlinearities(scip, branchruledata->nlcount, branchruledata->nlcountsize, &branchruledata->nlcountmax) );
      }
      else if( branchruledata->nlcountsize < nvars )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &branchruledata->nlcount, branchruledata->nlcountsize, nvars) );
         /**@todo should we update nlcounts for new variables? */
         BMSclearMemoryArray(&(branchruledata->nlcount[branchruledata->nlcountsize]), nvars - branchruledata->nlcountsize); /*lint !e866*/
         branchruledata->nlcountsize = nvars;
      }
      assert(branchruledata->nlcount != NULL);
      assert(branchruledata->nlcountsize == nvars);
      assert(branchruledata->nlcountmax >= 1);
   }
   else
   {
      SCIPfreeBlockMemoryArrayNull(scip, &branchruledata->nlcount, branchruledata->nlcountsize);
      branchruledata->nlcountsize = 0;
      branchruledata->nlcountmax = 1;
   }

   return SCIP_OKAY;
}


/** calculates nlscore value between 0 and 1 */
static
SCIP_Real calcNlscore(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  nlcount,            /**< array to store count values */
   int                   nlcountmax,         /**< maximum value in nlcount array */
   int                   probindex           /**< index of branching candidate */
   )
{
   if( nlcountmax >= 1 && nlcount != NULL )
   {
      SCIP_Real nlscore;

      assert(scip != NULL);
      assert(probindex >= 0);
      assert(probindex < SCIPgetNVars(scip));

      nlscore = nlcount[probindex] / (SCIP_Real)nlcountmax;

      assert(nlscore >= 0.0);
      assert(nlscore <= 1.0);
      return nlscore;
   }
   else
      return 0.0;
}

/** calculates an overall score value for the given individual score values */
static
SCIP_Real calcScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< branching rule data */
   SCIP_Real             conflictscore,      /**< conflict score of current variable */
   SCIP_Real             avgconflictscore,   /**< average conflict score */
   SCIP_Real             conflengthscore,    /**< conflict length score of current variable */
   SCIP_Real             avgconflengthscore, /**< average conflict length score */
   SCIP_Real             inferencescore,     /**< inference score of current variable */
   SCIP_Real             avginferencescore,  /**< average inference score */
   SCIP_Real             cutoffscore,        /**< cutoff score of current variable */
   SCIP_Real             avgcutoffscore,     /**< average cutoff score */
   SCIP_Real             pscostscore,        /**< pscost score of current variable */
   SCIP_Real             avgpscostscore,     /**< average pscost score */
   SCIP_Real             nlscore,            /**< nonlinear score of current variable between 0 and 1 */
   SCIP_Real             frac                /**< fractional value of variable in current solution */
   )
{
   SCIP_Real score;
   SCIP_Real dynamicfactor;

   assert(branchruledata != NULL);
   assert(0.0 < frac && frac < 1.0);

   if( branchruledata->dynamicweights )
   {
      dynamicfactor = (SCIPgetNInfeasibleLeaves(scip) + 1.0) / (SCIPgetNObjlimLeaves(scip) + 1.0);
   }
   else
      dynamicfactor = 1.0;

   score = dynamicfactor * (branchruledata->conflictweight * (1.0 - 1.0/(1.0+conflictscore/avgconflictscore))
            + branchruledata->conflengthweight * (1.0 - 1.0/(1.0+conflengthscore/avgconflengthscore))
            + branchruledata->inferenceweight * (1.0 - 1.0/(1.0+inferencescore/avginferencescore))
            + branchruledata->cutoffweight * (1.0 - 1.0/(1.0+cutoffscore/avgcutoffscore)))
         + branchruledata->pscostweight / dynamicfactor * (1.0 - 1.0/(1.0+pscostscore/avgpscostscore))
         + branchruledata->nlscoreweight * nlscore;

   /* avoid close to integral variables */
   if( MIN(frac, 1.0 - frac) < 10.0 * SCIPfeastol(scip) )
      score *= 1e-6;

   return score;
}

/** adds given index and direction to bound change arrays */
static
SCIP_RETCODE addBdchg(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 bdchginds,          /**< pointer to bound change index array */
   SCIP_BOUNDTYPE**      bdchgtypes,         /**< pointer to bound change types array */
   SCIP_Real**           bdchgbounds,        /**< pointer to bound change new bounds array */
   int*                  nbdchgs,            /**< pointer to number of bound changes */
   int                   ind,                /**< index to store in bound change index array */
   SCIP_BOUNDTYPE        type,               /**< type of the bound change to store in bound change type array */
   SCIP_Real             bound               /**< new bound to store in bound change new bounds array */
   )
{
   assert(bdchginds != NULL);
   assert(bdchgtypes != NULL);
   assert(bdchgbounds != NULL);
   assert(nbdchgs != NULL);

   SCIP_CALL( SCIPreallocBufferArray(scip, bdchginds, (*nbdchgs) + 1) );
   SCIP_CALL( SCIPreallocBufferArray(scip, bdchgtypes, (*nbdchgs) + 1) );
   SCIP_CALL( SCIPreallocBufferArray(scip, bdchgbounds, (*nbdchgs) + 1) );
   assert(*bdchginds != NULL);
   assert(*bdchgtypes != NULL);
   assert(*bdchgbounds != NULL);

   (*bdchginds)[*nbdchgs] = ind;
   (*bdchgtypes)[*nbdchgs] = type;
   (*bdchgbounds)[*nbdchgs] = bound;
   (*nbdchgs)++;

   return SCIP_OKAY;
}

/** frees bound change arrays */
static
void freeBdchgs(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 bdchginds,          /**< pointer to bound change index array */
   SCIP_BOUNDTYPE**      bdchgtypes,         /**< pointer to bound change types array */
   SCIP_Real**           bdchgbounds,        /**< pointer to bound change new bounds array */
   int*                  nbdchgs             /**< pointer to number of bound changes */
   )
{
   assert(bdchginds != NULL);
   assert(bdchgtypes != NULL);
   assert(bdchgbounds != NULL);
   assert(nbdchgs != NULL);

   SCIPfreeBufferArrayNull(scip, bdchgbounds);
   SCIPfreeBufferArrayNull(scip, bdchgtypes);
   SCIPfreeBufferArrayNull(scip, bdchginds);
   *nbdchgs = 0;
}

/** applies bound changes stored in bound change arrays */
static
SCIP_RETCODE applyBdchgs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   int*                  bdchginds,          /**< bound change index array */
   SCIP_BOUNDTYPE*       bdchgtypes,         /**< bound change types array */
   SCIP_Real*            bdchgbounds,        /**< bound change new bound array */
   int                   nbdchgs,            /**< number of bound changes */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
#ifndef NDEBUG
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
#endif
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   int i;

   assert(vars != NULL);

#ifndef NDEBUG
   /* find branching rule */
   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
#endif

   SCIPdebugMsg(scip, "applying %d bound changes\n", nbdchgs);

   for( i = 0; i < nbdchgs; ++i )
   {
      int v;

      v = bdchginds[i];

      SCIPdebugMsg(scip, " -> <%s> [%g,%g]\n",
         SCIPvarGetName(vars[v]), SCIPvarGetLbLocal(vars[v]), SCIPvarGetUbLocal(vars[v]));

      if( bdchgtypes[i] == SCIP_BOUNDTYPE_LOWER )
      {
         /* change lower bound of variable to given bound */
         SCIP_CALL( SCIPtightenVarLb(scip, vars[v], bdchgbounds[i], TRUE, &infeasible, &tightened) );
         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         /* if we did propagation, the bound change might already have been added */
         assert(tightened || (branchruledata->maxproprounds != 0));
      }
      else
      {
         assert(bdchgtypes[i] == SCIP_BOUNDTYPE_UPPER);

         /* change upper bound of variable to given bound */
         SCIP_CALL( SCIPtightenVarUb(scip, vars[v], bdchgbounds[i], TRUE, &infeasible, &tightened) );
         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         /* if we did propagation, the bound change might already have been added */
         assert(tightened || (branchruledata->maxproprounds != 0));
      }
      SCIPdebugMsg(scip, "  -> [%g,%g]\n", SCIPvarGetLbLocal(vars[v]), SCIPvarGetUbLocal(vars[v]));
   }

   return SCIP_OKAY;
}

/** execute reliability pseudo cost branching */
static
SCIP_RETCODE execRelpscost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_VAR**            branchcands,        /**< branching candidates */
   SCIP_Real*            branchcandssol,     /**< solution value for the branching candidates */
   SCIP_Real*            branchcandsfrac,    /**< fractional part of the branching candidates */
   int                   nbranchcands,       /**< number of branching candidates */
   SCIP_Bool             executebranch,      /**< execute a branching step or run probing only */
   SCIP_RESULT*          result              /**< pointer to the result of the execution */
   )
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_Real lpobjval;
   SCIP_Real bestsbdown;
   SCIP_Real bestsbup;
   SCIP_Real provedbound;
   SCIP_Bool bestsbdownvalid;
   SCIP_Bool bestsbupvalid;
   SCIP_Bool bestsbdowncutoff;
   SCIP_Bool bestsbupcutoff;
   SCIP_Bool bestisstrongbranch;
   SCIP_Bool allcolsinlp;
   SCIP_Bool exactsolve;
   int ninitcands;
   int bestcand;

   *result = SCIP_DIDNOTRUN;

   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get current LP objective bound of the local sub problem and global cutoff bound */
   lpobjval = SCIPgetLPObjval(scip);

   /* check, if we want to solve the problem exactly, meaning that strong branching information is not useful
    * for cutting off sub problems and improving lower bounds of children
    */
   exactsolve = SCIPisExactSolve(scip);

   /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
   allcolsinlp = SCIPallColsInLP(scip);

   bestcand = -1;
   bestisstrongbranch = FALSE;
   bestsbdown = SCIP_INVALID;
   bestsbup = SCIP_INVALID;
   bestsbdownvalid = FALSE;
   bestsbupvalid = FALSE;
   bestsbdowncutoff = FALSE;
   bestsbupcutoff = FALSE;
   provedbound = lpobjval;

   if( nbranchcands == 1 )
   {
      /* only one candidate: nothing has to be done */
      bestcand = 0;
      SCIPdebug(ninitcands = 0);
   }
   else
   {
      SCIP_VAR** vars;
      int* initcands;
      SCIP_Real* initcandscores;
      SCIP_Real* newlbs = NULL;
      SCIP_Real* newubs = NULL;
      int* bdchginds;
      SCIP_BOUNDTYPE* bdchgtypes;
      SCIP_Real* bdchgbounds;
      int maxninitcands;
      int nuninitcands;
      int nbdchgs;
      int nbdconflicts;
      SCIP_Real avgconflictscore;
      SCIP_Real avgconflengthscore;
      SCIP_Real avginferencescore;
      SCIP_Real avgcutoffscore;
      SCIP_Real avgpscostscore;
      SCIP_Real bestpsscore;
      SCIP_Real bestpsfracscore;
      SCIP_Real bestpsdomainscore;
      SCIP_Real bestsbscore;
      SCIP_Real bestuninitsbscore;
      SCIP_Real bestsbfracscore;
      SCIP_Real bestsbdomainscore;
      SCIP_Real prio;
      SCIP_Real reliable;
      SCIP_Real maxlookahead;
      SCIP_Real lookahead;
      SCIP_Real relerrorthreshold;
      SCIP_Bool initstrongbranching;
      SCIP_Bool propagate;
      SCIP_Bool probingbounds;
      SCIP_Longint nodenum;
      SCIP_Longint nlpiterationsquot;
      SCIP_Longint nsblpiterations;
      SCIP_Longint maxnsblpiterations;
      int bestsolidx;
      int maxbdchgs;
      int bestpscand;
      int bestsbcand;
      int bestuninitsbcand;
      int inititer;
      int nvars;
      int i;
      int c;
      SCIP_CONFIDENCELEVEL clevel;

      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      bestsolidx = SCIPgetBestSol(scip) == NULL ? -1 : SCIPsolGetIndex(SCIPgetBestSol(scip));

      /* get average conflict, inference, and pseudocost scores */
      avgconflictscore = SCIPgetAvgConflictScore(scip);
      avgconflictscore = MAX(avgconflictscore, 0.1);
      avgconflengthscore = SCIPgetAvgConflictlengthScore(scip);
      avgconflengthscore = MAX(avgconflengthscore, 0.1);
      avginferencescore = SCIPgetAvgInferenceScore(scip);
      avginferencescore = MAX(avginferencescore, 0.1);
      avgcutoffscore = SCIPgetAvgCutoffScore(scip);
      avgcutoffscore = MAX(avgcutoffscore, 0.1);
      avgpscostscore = SCIPgetAvgPseudocostScore(scip);
      avgpscostscore = MAX(avgpscostscore, 0.1);

      /* get nonlinear counts according to parameters */
      SCIP_CALL( branchruledataEnsureNlcount(scip, branchruledata) );

      initstrongbranching = FALSE;

      /* check whether propagation should be performed */
      propagate = (branchruledata->maxproprounds != 0);

      /* check whether valid bounds should be identified in probing-like fashion */
      probingbounds = propagate && branchruledata->probingbounds;

      /* get maximal number of candidates to initialize with strong branching; if the current solutions is not basic,
       * we cannot warmstart the simplex algorithm and therefore don't initialize any candidates
       */
      maxninitcands = MIN(nbranchcands, branchruledata->initcand);
      if( !SCIPisLPSolBasic(scip) )
         maxninitcands = 0;

      /* calculate maximal number of strong branching LP iterations; if we used too many, don't apply strong branching
       * any more
       */
      nlpiterationsquot = (SCIP_Longint)(branchruledata->sbiterquot * SCIPgetNNodeLPIterations(scip));
      maxnsblpiterations = nlpiterationsquot + branchruledata->sbiterofs + SCIPgetNRootStrongbranchLPIterations(scip);
      nsblpiterations = SCIPgetNStrongbranchLPIterations(scip);
      if( nsblpiterations > maxnsblpiterations )
         maxninitcands = 0;

      /* get buffer for storing the unreliable candidates */
      SCIP_CALL( SCIPallocBufferArray(scip, &initcands, maxninitcands+1) ); /* allocate one additional slot for convenience */
      SCIP_CALL( SCIPallocBufferArray(scip, &initcandscores, maxninitcands+1) );
      ninitcands = 0;

      /* get current node number */
      nodenum = SCIPgetNNodes(scip);

      /* initialize bound change arrays */
      bdchginds = NULL;
      bdchgtypes = NULL;
      bdchgbounds = NULL;
      nbdchgs = 0;
      nbdconflicts = 0;
      maxbdchgs = branchruledata->maxbdchgs;
      if( maxbdchgs == -1 )
         maxbdchgs = INT_MAX;

      /* calculate value used as reliability */
      prio = (maxnsblpiterations - nsblpiterations)/(nsblpiterations + 1.0);
      prio = MIN(prio, 1.0);
      prio = MAX(prio, (nlpiterationsquot - nsblpiterations)/(nsblpiterations + 1.0));
      reliable = (1.0-prio) * branchruledata->minreliable + prio * branchruledata->maxreliable;

      /* calculate the threshold for the relative error in the same way; low tolerance is more strict than higher tolerance */
      relerrorthreshold = (1.0 - prio) * branchruledata->higherrortol + prio * branchruledata->lowerrortol;

      clevel = (SCIP_CONFIDENCELEVEL)branchruledata->confidencelevel;
      /* determine the confidence level for hypothesis testing based on value of prio */
      if( branchruledata->usedynamicconfidence )
      {
         /* with decreasing priority, use a less strict confidence level */
         if( prio >= 0.9 )
            clevel = SCIP_CONFIDENCELEVEL_MAX;
         else if( prio >= 0.7 )
            clevel = SCIP_CONFIDENCELEVEL_HIGH;
         else if( prio >= 0.5 )
            clevel = SCIP_CONFIDENCELEVEL_MEDIUM;
         else if( prio >= 0.3 )
            clevel = SCIP_CONFIDENCELEVEL_LOW;
         else
            clevel = SCIP_CONFIDENCELEVEL_MIN;
      }

      /* search for the best pseudo cost candidate, while remembering unreliable candidates in a sorted buffer */
      nuninitcands = 0;
      bestpscand = -1;
      bestpsscore = -SCIPinfinity(scip);
      bestpsfracscore = -SCIPinfinity(scip);
      bestpsdomainscore = -SCIPinfinity(scip);

      /* search for the best candidate first */
      if( branchruledata->usehyptestforreliability )
      {
         for( c = 0; c < nbranchcands; ++c )
         {
            SCIP_Real conflictscore;
            SCIP_Real conflengthscore;
            SCIP_Real inferencescore;
            SCIP_Real cutoffscore;
            SCIP_Real pscostscore;
            SCIP_Real nlscore;
            SCIP_Real score;

            conflictscore = SCIPgetVarConflictScore(scip, branchcands[c]);
            conflengthscore = SCIPgetVarConflictlengthScore(scip, branchcands[c]);
            inferencescore = SCIPgetVarAvgInferenceScore(scip, branchcands[c]);
            cutoffscore = SCIPgetVarAvgCutoffScore(scip, branchcands[c]);
            nlscore = calcNlscore(scip, branchruledata->nlcount, branchruledata->nlcountmax, SCIPvarGetProbindex(branchcands[c]));
            pscostscore = SCIPgetVarPseudocostScore(scip, branchcands[c], branchcandssol[c]);

            /* replace the pseudo cost score with the already calculated one;
             * @todo: use old data for strong branching with propagation?
             */
            if( SCIPgetVarStrongbranchNode(scip, branchcands[c]) == nodenum )
            {
               SCIP_Real down;
               SCIP_Real up;
               SCIP_Real lastlpobjval;
               SCIP_Real downgain;
               SCIP_Real upgain;

               /* use the score of the strong branching call at the current node */
               SCIP_CALL( SCIPgetVarStrongbranchLast(scip, branchcands[c], &down, &up, NULL, NULL, NULL, &lastlpobjval) );
               downgain = MAX(down - lastlpobjval, 0.0);
               upgain = MAX(up - lastlpobjval, 0.0);
               pscostscore = SCIPgetBranchScore(scip, branchcands[c], downgain, upgain);

               SCIPdebugMsg(scip, " -> strong branching on variable <%s> already performed (down=%g (%+g), up=%g (%+g), pscostscore=%g)\n",
                  SCIPvarGetName(branchcands[c]), down, downgain, up, upgain, pscostscore);
            }

            score = calcScore(scip, branchruledata, conflictscore, avgconflictscore, conflengthscore, avgconflengthscore,
               inferencescore, avginferencescore, cutoffscore, avgcutoffscore, pscostscore, avgpscostscore, nlscore, branchcandsfrac[c]);

            /* check for better score of candidate */
            if( SCIPisSumGE(scip, score, bestpsscore) )
            {
               SCIP_Real fracscore;
               SCIP_Real domainscore;

               fracscore = MIN(branchcandsfrac[c], 1.0 - branchcandsfrac[c]);
               domainscore = -(SCIPvarGetUbLocal(branchcands[c]) - SCIPvarGetLbLocal(branchcands[c]));
               if( SCIPisSumGT(scip, score, bestpsscore)
                     || SCIPisSumGT(scip, fracscore, bestpsfracscore)
                     || (SCIPisSumGE(scip, fracscore, bestpsfracscore) && domainscore > bestpsdomainscore) )
               {
                  bestpscand = c;
                  bestpsscore = score;
                  bestpsfracscore = fracscore;
                  bestpsdomainscore = domainscore;
               }
            }
         }
      }

      for( c = 0; c < nbranchcands; ++c )
      {
         SCIP_Real conflictscore;
         SCIP_Real conflengthscore;
         SCIP_Real inferencescore;
         SCIP_Real cutoffscore;
         SCIP_Real pscostscore;
         SCIP_Real nlscore;
         SCIP_Real score;
         SCIP_Bool usesb;

         assert(branchcands[c] != NULL);
         assert(!SCIPisFeasIntegral(scip, branchcandssol[c]));

         /* get conflict, inference, cutoff, nonlinear, and pseudo cost scores for candidate */
         conflictscore = SCIPgetVarConflictScore(scip, branchcands[c]);
         conflengthscore = SCIPgetVarConflictlengthScore(scip, branchcands[c]);
         inferencescore = SCIPgetVarAvgInferenceScore(scip, branchcands[c]);
         cutoffscore = SCIPgetVarAvgCutoffScore(scip, branchcands[c]);
         nlscore = calcNlscore(scip, branchruledata->nlcount, branchruledata->nlcountmax, SCIPvarGetProbindex(branchcands[c]));
         pscostscore = SCIPgetVarPseudocostScore(scip, branchcands[c], branchcandssol[c]);
         usesb = FALSE;

         /* don't use strong branching on variables that have already been initialized at the current node;
          * instead replace the pseudo cost score with the already calculated one;
          * @todo: use old data for strong branching with propagation?
          */
         if( SCIPgetVarStrongbranchNode(scip, branchcands[c]) == nodenum )
         {
            SCIP_Real down;
            SCIP_Real up;
            SCIP_Real lastlpobjval;
            SCIP_Real downgain;
            SCIP_Real upgain;

            /* use the score of the strong branching call at the current node */
            SCIP_CALL( SCIPgetVarStrongbranchLast(scip, branchcands[c], &down, &up, NULL, NULL, NULL, &lastlpobjval) );
            downgain = MAX(down - lastlpobjval, 0.0);
            upgain = MAX(up - lastlpobjval, 0.0);
            pscostscore = SCIPgetBranchScore(scip, branchcands[c], downgain, upgain);

            SCIPdebugMsg(scip, " -> strong branching on variable <%s> already performed (down=%g (%+g), up=%g (%+g), pscostscore=%g)\n",
               SCIPvarGetName(branchcands[c]), down, downgain, up, upgain, pscostscore);
         }
         else if( maxninitcands > 0 )
         {
            SCIP_Real downsize;
            SCIP_Real upsize;
            SCIP_Real size;

            /* check, if the pseudo cost score of the variable is reliable */
            downsize = SCIPgetVarPseudocostCountCurrentRun(scip, branchcands[c], SCIP_BRANCHDIR_DOWNWARDS);
            upsize = SCIPgetVarPseudocostCountCurrentRun(scip, branchcands[c], SCIP_BRANCHDIR_UPWARDS);
            size = MIN(downsize, upsize);

            /* determine if variable is considered reliable based on the current reliability setting */
            /* check fixed number threshold (aka original) reliability first */
            assert(!branchruledata->usehyptestforreliability || bestpscand >= 0);
            usesb = FALSE;
            if( size < reliable )
               usesb = TRUE;
            else if( branchruledata->userelerrorforreliability && branchruledata->usehyptestforreliability )
            {
               if( !SCIPisVarPscostRelerrorReliable(scip, branchcands[c], relerrorthreshold, clevel) &&
                     !SCIPsignificantVarPscostDifference(scip, branchcands[bestpscand], branchcandsfrac[bestpscand],
                        branchcands[c], branchcandsfrac[c], SCIP_BRANCHDIR_DOWNWARDS, clevel, TRUE) &&
                     !SCIPsignificantVarPscostDifference(scip, branchcands[bestpscand], 1 - branchcandsfrac[bestpscand],
                        branchcands[c], 1 - branchcandsfrac[c], SCIP_BRANCHDIR_UPWARDS, clevel, TRUE) )
                  usesb = TRUE;
            }
            /* check if relative error is tolerable */
            else if( branchruledata->userelerrorforreliability &&
                  !SCIPisVarPscostRelerrorReliable(scip, branchcands[c], relerrorthreshold, clevel))
               usesb = TRUE;
            /* check if best pseudo-candidate is significantly better in both directions, use strong-branching otherwise */
            else if( branchruledata->usehyptestforreliability &&
                  !SCIPsignificantVarPscostDifference(scip, branchcands[bestpscand], branchcandsfrac[bestpscand],
                        branchcands[c], branchcandsfrac[c], SCIP_BRANCHDIR_DOWNWARDS, clevel, TRUE) &&
                  !SCIPsignificantVarPscostDifference(scip, branchcands[bestpscand], 1 - branchcandsfrac[bestpscand],
                        branchcands[c], 1 - branchcandsfrac[c], SCIP_BRANCHDIR_UPWARDS, clevel, TRUE))
               usesb = TRUE;

            /* count the number of variables that are completely uninitialized */
            if( size < 0.1 )
               nuninitcands++;
         }

         /* combine the five score values */
         score = calcScore(scip, branchruledata, conflictscore, avgconflictscore, conflengthscore, avgconflengthscore,
            inferencescore, avginferencescore, cutoffscore, avgcutoffscore, pscostscore, avgpscostscore, nlscore, branchcandsfrac[c]);

         if( usesb )
         {
            int j;

            /* assign a random score to this uninitialized candidate */
            if( branchruledata->randinitorder )
               score += SCIPrandomGetReal(branchruledata->randnumgen, 0.0, 1e-4);

            /* pseudo cost of variable is not reliable: insert candidate in initcands buffer */
            for( j = ninitcands; j > 0 && score > initcandscores[j-1]; --j )
            {
               initcands[j] = initcands[j-1];
               initcandscores[j] = initcandscores[j-1];
            }
            initcands[j] = c;
            initcandscores[j] = score;
            ninitcands++;
            ninitcands = MIN(ninitcands, maxninitcands);
         }
         /* in the case of hypothesis reliability, the best pseudo candidate has been determined already */
         else if( !branchruledata->usehyptestforreliability )
         {
            /* variable will keep it's pseudo cost value: check for better score of candidate */
            if( SCIPisSumGE(scip, score, bestpsscore) )
            {
               SCIP_Real fracscore;
               SCIP_Real domainscore;

               fracscore = MIN(branchcandsfrac[c], 1.0 - branchcandsfrac[c]);
               domainscore = -(SCIPvarGetUbLocal(branchcands[c]) - SCIPvarGetLbLocal(branchcands[c]));
               if( SCIPisSumGT(scip, score, bestpsscore)
                  || SCIPisSumGT(scip, fracscore, bestpsfracscore)
                  || (SCIPisSumGE(scip, fracscore, bestpsfracscore) && domainscore > bestpsdomainscore) )
               {
                  bestpscand = c;
                  bestpsscore = score;
                  bestpsfracscore = fracscore;
                  bestpsdomainscore = domainscore;
               }
            }
         }
      }

      /* in the special case that only the best pseudo candidate was selected for strong branching, skip the strong branching */
      if( branchruledata->usehyptestforreliability && ninitcands == 1 )
      {
         ninitcands = 0;
         SCIPdebugMsg(scip, "Only one single candidate for initialization-->Skipping strong branching\n");
      }

      /* initialize unreliable candidates with strong branching until maxlookahead is reached,
       * search best strong branching candidate
       */
      maxlookahead = (SCIP_Real)branchruledata->maxlookahead * (1.0 + (SCIP_Real)nuninitcands/(SCIP_Real)nbranchcands);
      inititer = branchruledata->inititer;
      if( inititer == 0 )
      {
         SCIP_Longint nlpiterations;
         SCIP_Longint nlps;

         /* iteration limit is set to twice the average number of iterations spent to resolve a dual feasible SCIP_LP;
          * at the first few nodes, this average is not very exact, so we better increase the iteration limit on
          * these very important nodes
          */
         nlpiterations = SCIPgetNDualResolveLPIterations(scip);
         nlps = SCIPgetNDualResolveLPs(scip);
         if( nlps == 0 )
         {
            nlpiterations = SCIPgetNNodeInitLPIterations(scip);
            nlps = SCIPgetNNodeInitLPs(scip);
            if( nlps == 0 )
            {
               nlpiterations = 1000;
               nlps = 1;
            }
         }
         assert(nlps >= 1);
         inititer = (int)(2*nlpiterations / nlps);
         inititer = (int)((SCIP_Real)inititer * (1.0 + 20.0/nodenum));
         inititer = MAX(inititer, 10);
         inititer = MIN(inititer, 500);
      }

      SCIPdebugMsg(scip, "strong branching (reliable=%g, %d/%d cands, %d uninit, maxcands=%d, maxlookahead=%g, maxbdchgs=%d, inititer=%d, iterations:%" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT ", basic:%u)\n",
         reliable, ninitcands, nbranchcands, nuninitcands, maxninitcands, maxlookahead, maxbdchgs, inititer,
         SCIPgetNStrongbranchLPIterations(scip), maxnsblpiterations, SCIPisLPSolBasic(scip));

      bestsbcand = -1;
      bestsbscore = -SCIPinfinity(scip);
      bestsbfracscore = -SCIPinfinity(scip);
      bestsbdomainscore = -SCIPinfinity(scip);
      bestuninitsbscore = -SCIPinfinity(scip);
      bestuninitsbcand = -1;
      lookahead = 0.0;
      for( i = 0; i < ninitcands && lookahead < maxlookahead && nbdchgs + nbdconflicts < maxbdchgs
              && (i < (int) maxlookahead || SCIPgetNStrongbranchLPIterations(scip) < maxnsblpiterations); ++i )
      {
         SCIP_Real down;
         SCIP_Real up;
         SCIP_Real downgain;
         SCIP_Real upgain;
         SCIP_Bool downvalid;
         SCIP_Bool upvalid;
         SCIP_Longint ndomredsdown;
         SCIP_Longint ndomredsup;
         SCIP_Bool lperror;
         SCIP_Bool downinf;
         SCIP_Bool upinf;
         SCIP_Bool downconflict;
         SCIP_Bool upconflict;

         /* get candidate number to initialize */
         c = initcands[i];
         assert(!SCIPisFeasIntegral(scip, branchcandssol[c]));

         if( branchruledata->skipbadinitcands )
         {
            SCIP_Bool skipsb = FALSE;
            /* if the current best candidate is a candidate found by strong branching, determine if candidate pseudo-costs are
             * significantly smaller in at least one direction, in which case we safe the execution of strong-branching for now
             */
            if( bestsbscore > bestpsscore && bestsbscore > bestuninitsbscore && bestsbupvalid && bestsbdownvalid )
            {
               assert(bestsbcand != -1);
               assert(bestsbup != SCIP_INVALID && bestsbdown != SCIP_INVALID); /*lint !e777 lint doesn't like comparing floats */

               /* test if the variable is unlikely to produce a better gain than the currently best one. Skip strong-branching
                * in such a case
                */
               if( SCIPpscostThresholdProbabilityTest(scip, branchcands[c], branchcandsfrac[c], bestsbdown,
                     SCIP_BRANCHDIR_DOWNWARDS, clevel)
                     || SCIPpscostThresholdProbabilityTest(scip, branchcands[c], 1.0 - branchcandsfrac[c], bestsbup,
                           SCIP_BRANCHDIR_UPWARDS, clevel) )
                  skipsb = TRUE;
            }
            /* the currently best candidate is also a pseudo-candidate; apply significance test and skip candidate if it
             * is significantly worse in at least one direction
             */
            else if( bestpscand != -1 && bestpsscore > bestuninitsbscore )
            {
               if( SCIPsignificantVarPscostDifference(scip, branchcands[bestpscand], branchcandsfrac[bestpscand],
                     branchcands[c], branchcandsfrac[c], SCIP_BRANCHDIR_DOWNWARDS, clevel, TRUE)
                     || SCIPsignificantVarPscostDifference(scip, branchcands[bestpscand], 1.0 - branchcandsfrac[bestpscand],
                           branchcands[c], 1.0 - branchcandsfrac[c], SCIP_BRANCHDIR_UPWARDS, clevel, TRUE) )
                  skipsb = TRUE;
            }
            /* compare against the best init cand that has been skipped already */
            else if( bestuninitsbcand != -1 )
            {
               if( SCIPsignificantVarPscostDifference(scip, branchcands[bestuninitsbcand], branchcandsfrac[bestuninitsbcand],
                     branchcands[c], branchcandsfrac[c], SCIP_BRANCHDIR_DOWNWARDS, clevel, TRUE)
                     || SCIPsignificantVarPscostDifference(scip, branchcands[bestuninitsbcand], 1.0 - branchcandsfrac[bestuninitsbcand],
                           branchcands[c], 1.0 - branchcandsfrac[c], SCIP_BRANCHDIR_UPWARDS, clevel, TRUE) )
                  skipsb = TRUE;
            }

            /* skip candidate, update the best score of an unitialized candidate */
            if( skipsb )
            {
               if( bestuninitsbcand == -1 )
               {
                  bestuninitsbcand = c;
                  bestuninitsbscore = initcandscores[i];
               }
               continue;
            }
         }
         SCIPdebugMsg(scip, "init pseudo cost (%g/%g) of <%s> at %g (score:%g) with strong branching (%d iterations) -- %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT " iterations\n",
            SCIPgetVarPseudocostCountCurrentRun(scip, branchcands[c], SCIP_BRANCHDIR_DOWNWARDS),
            SCIPgetVarPseudocostCountCurrentRun(scip, branchcands[c], SCIP_BRANCHDIR_UPWARDS),
            SCIPvarGetName(branchcands[c]), branchcandssol[c], initcandscores[i],
            inititer, SCIPgetNStrongbranchLPIterations(scip), maxnsblpiterations);

         /* use strong branching on candidate */
         if( !initstrongbranching )
         {
            initstrongbranching = TRUE;

            SCIP_CALL( SCIPstartStrongbranch(scip, propagate) );

            /* create arrays for probing-like bound tightening */
            if( probingbounds )
            {
               SCIP_CALL( SCIPallocBufferArray(scip, &newlbs, nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &newubs, nvars) );
            }
         }

         if( propagate )
         {
            /* apply strong branching */
            SCIP_CALL( SCIPgetVarStrongbranchWithPropagation(scip, branchcands[c], branchcandssol[c], lpobjval, inititer,
                  branchruledata->maxproprounds, &down, &up, &downvalid, &upvalid, &ndomredsdown, &ndomredsup, &downinf, &upinf,
                  &downconflict, &upconflict, &lperror, newlbs, newubs) );
         }
         else
         {
            /* apply strong branching */
            SCIP_CALL( SCIPgetVarStrongbranchFrac(scip, branchcands[c], inititer,
                  &down, &up, &downvalid, &upvalid, &downinf, &upinf, &downconflict, &upconflict, &lperror) );

            ndomredsdown = ndomredsup = 0;
         }

         /* check for an error in strong branching */
         if( lperror )
         {
            if( !SCIPisStopped(scip) )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                  "(node %" SCIP_LONGINT_FORMAT ") error in strong branching call for variable <%s> with solution %g\n",
                  SCIPgetNNodes(scip), SCIPvarGetName(branchcands[c]), branchcandssol[c]);
            }
            break;
         }

         /* Strong branching might have found a new primal solution which updated the cutoff bound. In this case, the
          * provedbound computed before can be higher than the cutoffbound and the current node can be cut off.
          * Additionally, also if the value for the current best candidate is valid and exceeds the new cutoff bound,
          * we want to change the domain of this variable rather than branching on it.
          */
         if( SCIPgetBestSol(scip) != NULL && SCIPsolGetIndex(SCIPgetBestSol(scip)) != bestsolidx )
         {
            bestsolidx = SCIPsolGetIndex(SCIPgetBestSol(scip));

            SCIPdebugMsg(scip, " -> strong branching on variable <%s> lead to a new incumbent\n",
               SCIPvarGetName(branchcands[c]));

            /* proved bound for current node is larger than new cutoff bound -> cut off current node */
            if( SCIPisGE(scip, provedbound, SCIPgetCutoffbound(scip)) )
            {
               SCIPdebugMsg(scip, " -> node can be cut off (provedbound=%g, cutoff=%g)\n", provedbound, SCIPgetCutoffbound(scip));

               *result = SCIP_CUTOFF;
               break; /* terminate initialization loop, because node is infeasible */
            }
            /* proved bound for down child of best candidate is larger than cutoff bound
             *  -> increase lower bound of best candidate
             * we must only do this if the LP is complete, i.e., we are not doing column generation
             */

            else if( bestsbcand != -1  && allcolsinlp )
            {
               if( !bestsbdowncutoff && bestsbdownvalid && SCIPisGE(scip, bestsbdown, SCIPgetCutoffbound(scip)) )
               {
                  bestsbdowncutoff = TRUE;

                  SCIPdebugMsg(scip, " -> valid dual bound for down child of best candidate <%s> is higher than new cutoff bound (valid=%u, bestsbdown=%g, cutoff=%g)\n",
                     SCIPvarGetName(branchcands[bestsbcand]), bestsbdownvalid, bestsbdown, SCIPgetCutoffbound(scip));

                  SCIPdebugMsg(scip, " -> increase lower bound of best candidate <%s> to %g\n",
                     SCIPvarGetName(branchcands[bestsbcand]), SCIPfeasCeil(scip, branchcandssol[bestsbcand]));

                  SCIP_CALL( addBdchg(scip, &bdchginds, &bdchgtypes, &bdchgbounds, &nbdchgs, SCIPvarGetProbindex(branchcands[bestsbcand]),
                        SCIP_BOUNDTYPE_LOWER, SCIPfeasCeil(scip, branchcandssol[bestsbcand])) );
               }
               /* proved bound for up child of best candidate is larger than cutoff bound -> decrease upper bound of best candidate */
               else if( !bestsbupcutoff && bestsbupvalid && SCIPisGE(scip, bestsbup, SCIPgetCutoffbound(scip)) )
               {
                  bestsbupcutoff = TRUE;

                  SCIPdebugMsg(scip, " -> valid dual bound for up child of best candidate <%s> is higher than new cutoff bound (valid=%u, bestsbup=%g, cutoff=%g)\n",
                     SCIPvarGetName(branchcands[bestsbcand]), bestsbupvalid, bestsbup, SCIPgetCutoffbound(scip));

                  SCIPdebugMsg(scip, " -> decrease upper bound of best candidate <%s> to %g\n",
                     SCIPvarGetName(branchcands[bestsbcand]), SCIPfeasFloor(scip, branchcandssol[bestsbcand]));

                  SCIP_CALL( addBdchg(scip, &bdchginds, &bdchgtypes, &bdchgbounds, &nbdchgs, SCIPvarGetProbindex(branchcands[bestsbcand]),
                        SCIP_BOUNDTYPE_UPPER, SCIPfeasFloor(scip, branchcandssol[bestsbcand])) );
               }
            }
         }

         /* evaluate strong branching */
         down = MAX(down, lpobjval);
         up = MAX(up, lpobjval);
         downgain = down - lpobjval;
         upgain = up - lpobjval;
         assert(!allcolsinlp || exactsolve || !downvalid || downinf == SCIPisGE(scip, down, SCIPgetCutoffbound(scip)));
         assert(!allcolsinlp || exactsolve || !upvalid || upinf == SCIPisGE(scip, up, SCIPgetCutoffbound(scip)));
         assert(downinf || !downconflict);
         assert(upinf || !upconflict);

         /* @todo: store pseudo cost only for valid bounds?
          * depending on the user parameter choice of storesemiinitcosts, pseudo costs are also updated in single directions,
          * if the node in the other direction was infeasible or cut off
          */
         if( !downinf
#ifdef WITH_LPSOLSTAT
               && SCIPgetLastStrongbranchLPSolStat(scip, SCIP_BRANCHDIR_DOWNWARDS) != SCIP_LPSOLSTAT_ITERLIMIT
#endif
               && (!upinf || branchruledata->storesemiinitcosts) )
         {
            SCIP_Real weight;

            /* smaller weights are given if the strong branching hit the time limit in the corresponding direction */
            if( branchruledata->usesmallweightsitlim )
               weight = SCIPgetLastStrongbranchLPSolStat(scip, SCIP_BRANCHDIR_DOWNWARDS) != SCIP_LPSOLSTAT_ITERLIMIT ? 1.0 : 0.5;
            else
               weight = 1.0;

            /* update pseudo cost values */
            SCIP_CALL( SCIPupdateVarPseudocost(scip, branchcands[c], 0.0 - branchcandsfrac[c], downgain, weight) );
         }
         if( !upinf
#ifdef WITH_LPSOLSTAT
               && SCIPgetLastStrongbranchLPSolStat(scip, SCIP_BRANCHDIR_UPWARDS) != SCIP_LPSOLSTAT_ITERLIMIT
#endif
               && (!downinf || branchruledata->storesemiinitcosts)  )
         {
            SCIP_Real weight;

            /* smaller weights are given if the strong branching hit the time limit in the corresponding direction */
            if( branchruledata->usesmallweightsitlim )
               weight = SCIPgetLastStrongbranchLPSolStat(scip, SCIP_BRANCHDIR_UPWARDS) != SCIP_LPSOLSTAT_ITERLIMIT ? 1.0 : 0.5;
            else
               weight = 1.0;

            SCIP_CALL( SCIPupdateVarPseudocost(scip, branchcands[c], 1.0 - branchcandsfrac[c], upgain, weight) );
         }

         /* the minimal lower bound of both children is a proved lower bound of the current subtree */
         if( allcolsinlp && !exactsolve && downvalid && upvalid )
         {
            SCIP_Real minbound;

            minbound = MIN(down, up);
            provedbound = MAX(provedbound, minbound);
            assert((downinf && upinf) || SCIPisLT(scip, provedbound, SCIPgetCutoffbound(scip)));

            /* save probing-like bounds detected during strong branching */
            if( probingbounds )
            {
               int v;

               assert(newlbs != NULL);
               assert(newubs != NULL);

               for( v = 0; v < nvars; ++v )
               {
                  if( SCIPisGT(scip, newlbs[v], SCIPvarGetLbLocal(vars[v])) )
                  {
                     SCIPdebugMsg(scip, "better lower bound for variable <%s>: %.9g -> %.9g (by strongbranching on <%s>)\n",
                        SCIPvarGetName(vars[v]), SCIPvarGetLbLocal(vars[v]), newlbs[v], SCIPvarGetName(branchcands[c]));

                     SCIP_CALL( addBdchg(scip, &bdchginds, &bdchgtypes, &bdchgbounds, &nbdchgs, v,
                           SCIP_BOUNDTYPE_LOWER, newlbs[v]) );
                  }
                  if( SCIPisLT(scip, newubs[v], SCIPvarGetUbLocal(vars[v])) )
                  {
                     SCIPdebugMsg(scip, "better upper bound for variable <%s>: %.9g -> %.9g (by strongbranching on <%s>)\n",
                        SCIPvarGetName(vars[v]), SCIPvarGetUbLocal(vars[v]), newubs[v], SCIPvarGetName(branchcands[c]));

                     SCIP_CALL( addBdchg(scip, &bdchginds, &bdchgtypes, &bdchgbounds, &nbdchgs, v,
                           SCIP_BOUNDTYPE_UPPER, newubs[v]) );
                  }
               }
            }
         }

         /* check if there are infeasible roundings */
         if( downinf || upinf )
         {
            assert(allcolsinlp || propagate);
            assert(!exactsolve);

            if( downinf && upinf )
            {
               /* both roundings are infeasible -> node is infeasible */
               SCIPdebugMsg(scip, " -> variable <%s> is infeasible in both directions (conflict: %u/%u)\n",
                  SCIPvarGetName(branchcands[c]), downconflict, upconflict);
               *result = SCIP_CUTOFF;
               break; /* terminate initialization loop, because node is infeasible */
            }
            else
            {
               /* rounding is infeasible in one direction -> round variable in other direction */
               SCIPdebugMsg(scip, " -> variable <%s> is infeasible in %s branch (conflict: %u/%u)\n",
                  SCIPvarGetName(branchcands[c]), downinf ? "downward" : "upward", downconflict, upconflict);
               SCIP_CALL( addBdchg(scip, &bdchginds, &bdchgtypes, &bdchgbounds, &nbdchgs, SCIPvarGetProbindex(branchcands[c]),
                     (downinf ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER),
                     (downinf ? SCIPfeasCeil(scip, branchcandssol[c]) : SCIPfeasFloor(scip, branchcandssol[c]))) );
               if( nbdchgs + nbdconflicts >= maxbdchgs )
                  break; /* terminate initialization loop, because enough roundings are performed */
            }
         }
         else
         {
            SCIP_Real conflictscore;
            SCIP_Real conflengthscore;
            SCIP_Real inferencescore;
            SCIP_Real cutoffscore;
            SCIP_Real pscostscore;
            SCIP_Real nlscore;
            SCIP_Real score;

            /* check for a better score */
            conflictscore = SCIPgetVarConflictScore(scip, branchcands[c]);
            conflengthscore = SCIPgetVarConflictlengthScore(scip, branchcands[c]);
            nlscore = calcNlscore(scip, branchruledata->nlcount, branchruledata->nlcountmax, SCIPvarGetProbindex(branchcands[c]));

            /* optionally, use only local information obtained via strong branching for this candidate, i.e., local
             * domain reductions and no cutoff score
             */
            inferencescore = branchruledata->usesblocalinfo ? SCIPgetBranchScore(scip, branchcands[c], (SCIP_Real)ndomredsdown, (SCIP_Real)ndomredsup)
                  : SCIPgetVarAvgInferenceScore(scip, branchcands[c]);
            cutoffscore = branchruledata->usesblocalinfo ? 0.0 : SCIPgetVarAvgCutoffScore(scip, branchcands[c]);
            pscostscore = SCIPgetBranchScore(scip, branchcands[c], downgain, upgain);

            score = calcScore(scip, branchruledata, conflictscore, avgconflictscore, conflengthscore, avgconflengthscore,
               inferencescore, avginferencescore, cutoffscore, avgcutoffscore, pscostscore, avgpscostscore, nlscore, branchcandsfrac[c]);

            if( SCIPisSumGE(scip, score, bestsbscore) )
            {
               SCIP_Real fracscore;
               SCIP_Real domainscore;

               fracscore = MIN(branchcandsfrac[c], 1.0 - branchcandsfrac[c]);
               domainscore = -(SCIPvarGetUbLocal(branchcands[c]) - SCIPvarGetLbLocal(branchcands[c]));
               if( SCIPisSumGT(scip, score, bestsbscore)
                  || SCIPisSumGT(scip, fracscore, bestsbfracscore)
                  || (SCIPisSumGE(scip, fracscore, bestsbfracscore) && domainscore > bestsbdomainscore) )
               {
                  bestsbcand = c;
                  bestsbscore = score;
                  bestsbdown = down;
                  bestsbup = up;
                  bestsbdownvalid = downvalid;
                  bestsbupvalid = upvalid;
                  bestsbdowncutoff = FALSE;
                  bestsbupcutoff = FALSE;
                  bestsbfracscore = fracscore;
                  bestsbdomainscore = domainscore;
                  lookahead = 0.0;
               }
               else
                  lookahead += 0.5;
            }
            else
               lookahead += 1.0;

            SCIPdebugMsg(scip, " -> variable <%s> (solval=%g, down=%g (%+g,valid=%u), up=%g (%+g,valid=%u), score=%g/ %g/%g %g/%g -> %g)\n",
               SCIPvarGetName(branchcands[c]), branchcandssol[c], down, downgain, downvalid, up, upgain, upvalid,
               pscostscore, conflictscore, conflengthscore, inferencescore, cutoffscore,  score);
         }
      }
#ifdef SCIP_DEBUG
      if( bestsbcand >= 0 )
      {
         SCIPdebugMsg(scip, " -> best: <%s> (%g / %g / %g), lookahead=%g/%g\n",
            SCIPvarGetName(branchcands[bestsbcand]), bestsbscore, bestsbfracscore, bestsbdomainscore,
            lookahead, maxlookahead);
      }
#endif

      if( initstrongbranching )
      {
         if( probingbounds )
         {
            assert(newlbs != NULL);
            assert(newubs != NULL);

            SCIPfreeBufferArray(scip, &newubs);
            SCIPfreeBufferArray(scip, &newlbs);
         }

         SCIP_CALL( SCIPendStrongbranch(scip) );

         if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OBJLIMIT || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE )
         {
            assert(SCIPhasCurrentNodeLP(scip));
            *result = SCIP_CUTOFF;
         }
      }

      if( *result != SCIP_CUTOFF )
      {
         /* get the score of the best uninitialized strong branching candidate */
         if( i < ninitcands && bestuninitsbcand == -1 )
            bestuninitsbscore = initcandscores[i];

         /* if the best pseudo cost candidate is better than the best uninitialized strong branching candidate,
          * compare it to the best initialized strong branching candidate
          */
         if( bestpsscore > bestuninitsbscore && SCIPisSumGT(scip, bestpsscore, bestsbscore) )
         {
            bestcand = bestpscand;
            bestisstrongbranch = FALSE;
         }
         else if( bestsbcand >= 0 )
         {
            bestcand = bestsbcand;
            bestisstrongbranch = TRUE;
         }
         else
         {
            /* no candidate was initialized, and the best score is the one of the first candidate in the initialization
             * queue
             */
            assert(ninitcands >= 1);
            bestcand = initcands[0];
            bestisstrongbranch = FALSE;
         }

         /* update lower bound of current node */
         if( allcolsinlp && !exactsolve )
         {
            assert(SCIPisLT(scip, provedbound, SCIPgetCutoffbound(scip)));
            SCIP_CALL( SCIPupdateNodeLowerbound(scip, SCIPgetCurrentNode(scip), provedbound) );
         }
      }

      /* apply domain reductions */
      if( nbdchgs > 0 )
      {
         if( *result != SCIP_CUTOFF )
         {
            SCIP_CALL( applyBdchgs(scip, vars, bdchginds, bdchgtypes, bdchgbounds, nbdchgs, result) );
            if( *result != SCIP_CUTOFF )
               *result = SCIP_REDUCEDDOM;
         }
         freeBdchgs(scip, &bdchginds, &bdchgtypes, &bdchgbounds, &nbdchgs);
      }

      /* free buffer for the unreliable candidates */
      SCIPfreeBufferArray(scip, &initcandscores);
      SCIPfreeBufferArray(scip, &initcands);
   }

   /* if no domain could be reduced, create the branching */
   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED && executebranch )
   {
      SCIP_NODE* downchild;
      SCIP_NODE* upchild;
      SCIP_VAR* var;
      SCIP_Real val;

      assert(*result == SCIP_DIDNOTRUN);
      assert(0 <= bestcand && bestcand < nbranchcands);
      assert(!SCIPisFeasIntegral(scip, branchcandssol[bestcand]));
      assert(SCIPisLT(scip, provedbound, SCIPgetCutoffbound(scip)));
      assert(!bestsbdowncutoff && !bestsbupcutoff);

      var = branchcands[bestcand];
      val = branchcandssol[bestcand];

      /* perform the branching */
      SCIPdebugMsg(scip, " -> %d (%d) cands, sel cand %d: var <%s> (sol=%g, down=%g (%+g), up=%g (%+g), sb=%u, psc=%g/%g [%g])\n",
         nbranchcands, ninitcands, bestcand, SCIPvarGetName(var), branchcandssol[bestcand],
         bestsbdown, bestsbdown - lpobjval, bestsbup, bestsbup - lpobjval, bestisstrongbranch,
         SCIPgetVarPseudocostCurrentRun(scip, var, SCIP_BRANCHDIR_DOWNWARDS),
         SCIPgetVarPseudocostCurrentRun(scip, var, SCIP_BRANCHDIR_UPWARDS),
         SCIPgetVarPseudocostScoreCurrentRun(scip, var, branchcandssol[bestcand]));
      SCIP_CALL( SCIPbranchVarVal(scip, var, val, &downchild, NULL, &upchild) );
      assert(downchild != NULL);
      assert(upchild != NULL);

      /* update the lower bounds in the children */
      if( bestisstrongbranch && allcolsinlp && !exactsolve )
      {
         if( bestsbdownvalid )
         {
            assert(SCIPisLT(scip, bestsbdown, SCIPgetCutoffbound(scip)));
            SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, bestsbdown) );
            assert(SCIPisGE(scip, SCIPgetNodeLowerbound(scip, downchild), provedbound));
         }
         if( bestsbupvalid )
         {
            assert(SCIPisLT(scip, bestsbup, SCIPgetCutoffbound(scip)));
            SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, bestsbup) );
            assert(SCIPisGE(scip, SCIPgetNodeLowerbound(scip, upchild), provedbound));
         }
      }

      SCIPdebugMsg(scip, " -> down child's lowerbound: %g\n", SCIPnodeGetLowerbound(downchild));
      SCIPdebugMsg(scip, " -> up child's lowerbound  : %g\n", SCIPnodeGetLowerbound(upchild));

      assert(SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_INFEASIBLE && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OBJLIMIT);

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyRelpscost)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleRelpscost(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeRelpscost)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeBlockMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
static
SCIP_DECL_BRANCHINITSOL(branchInitsolRelpscost)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* initialize branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   branchruledata->nlcount = NULL;
   branchruledata->nlcountsize = 0;
   branchruledata->nlcountmax = 1;
   assert(branchruledata->startrandseed >= 0);

   /* create a random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &branchruledata->randnumgen,
         (unsigned int)branchruledata->startrandseed) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolRelpscost)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free memory in branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeBlockMemoryArrayNull(scip, &branchruledata->nlcount, branchruledata->nlcountsize);

   /* free random number generator */
   SCIPfreeRandom(scip, &branchruledata->randnumgen);

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpRelpscost)
{  /*lint --e{715}*/
   SCIP_VAR** tmplpcands;
   SCIP_VAR** lpcands;
   SCIP_Real* tmplpcandssol;
   SCIP_Real* lpcandssol;
   SCIP_Real* tmplpcandsfrac;
   SCIP_Real* lpcandsfrac;
   int nlpcands;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Execlp method of relpscost branching in node %llu\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      *result = SCIP_DIDNOTRUN;
      SCIPdebugMsg(scip, "Could not apply relpscost branching, as the current LP was not solved to optimality.\n");

      return SCIP_OKAY;
   }

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &tmplpcands, &tmplpcandssol, &tmplpcandsfrac, NULL, &nlpcands, NULL) );
   assert(nlpcands > 0);

   /* copy LP banching candidates and solution values, because they will be updated w.r.t. the strong branching LP
    * solution
    */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcands, tmplpcands, nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandssol, tmplpcandssol, nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandsfrac, tmplpcandsfrac, nlpcands) );

   /* execute branching rule */
   SCIP_CALL( execRelpscost(scip, branchrule, lpcands, lpcandssol, lpcandsfrac, nlpcands, TRUE, result) );

   SCIPfreeBufferArray(scip, &lpcandsfrac);
   SCIPfreeBufferArray(scip, &lpcandssol);
   SCIPfreeBufferArray(scip, &lpcands);

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the reliable pseudo cost branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleRelpscost(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create relpscost branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non-fundamental callbacks via specific setter functions*/
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyRelpscost) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeRelpscost) );
   SCIP_CALL( SCIPsetBranchruleInitsol(scip, branchrule, branchInitsolRelpscost) );
   SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolRelpscost) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpRelpscost) );

   /* relpscost branching rule parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/relpscost/conflictweight",
         "weight in score calculations for conflict score",
         &branchruledata->conflictweight, TRUE, DEFAULT_CONFLICTWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/relpscost/conflictlengthweight",
         "weight in score calculations for conflict length score",
         &branchruledata->conflengthweight, TRUE, DEFAULT_CONFLENGTHWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/relpscost/inferenceweight",
         "weight in score calculations for inference score",
         &branchruledata->inferenceweight, TRUE, DEFAULT_INFERENCEWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/relpscost/cutoffweight",
         "weight in score calculations for cutoff score",
         &branchruledata->cutoffweight, TRUE, DEFAULT_CUTOFFWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/relpscost/pscostweight",
         "weight in score calculations for pseudo cost score",
         &branchruledata->pscostweight, TRUE, DEFAULT_PSCOSTWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/relpscost/nlscoreweight",
         "weight in score calculations for nlcount score",
         &branchruledata->nlscoreweight, TRUE, DEFAULT_NLSCOREWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/relpscost/minreliable",
         "minimal value for minimum pseudo cost size to regard pseudo cost value as reliable",
         &branchruledata->minreliable, TRUE, DEFAULT_MINRELIABLE, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/relpscost/maxreliable",
         "maximal value for minimum pseudo cost size to regard pseudo cost value as reliable",
         &branchruledata->maxreliable, TRUE, DEFAULT_MAXRELIABLE, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/relpscost/sbiterquot",
         "maximal fraction of strong branching LP iterations compared to node relaxation LP iterations",
         &branchruledata->sbiterquot, FALSE, DEFAULT_SBITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/relpscost/sbiterofs",
         "additional number of allowed strong branching LP iterations",
         &branchruledata->sbiterofs, FALSE, DEFAULT_SBITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/relpscost/maxlookahead",
         "maximal number of further variables evaluated without better score",
         &branchruledata->maxlookahead, TRUE, DEFAULT_MAXLOOKAHEAD, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/relpscost/initcand",
         "maximal number of candidates initialized with strong branching per node",
         &branchruledata->initcand, FALSE, DEFAULT_INITCAND, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/relpscost/inititer",
         "iteration limit for strong branching initializations of pseudo cost entries (0: auto)",
         &branchruledata->inititer, FALSE, DEFAULT_INITITER, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/relpscost/maxbdchgs",
         "maximal number of bound tightenings before the node is reevaluated (-1: unlimited)",
         &branchruledata->maxbdchgs, TRUE, DEFAULT_MAXBDCHGS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/relpscost/maxproprounds",
         "maximum number of propagation rounds to be performed during strong branching before solving the LP (-1: no limit, -2: parameter settings)",
         &branchruledata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -2, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/relpscost/probingbounds",
         "should valid bounds be identified in a probing-like fashion during strong branching (only with propagation)?",
         &branchruledata->probingbounds, TRUE, DEFAULT_PROBINGBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/relpscost/userelerrorreliability",
         "should reliability be based on relative errors?", &branchruledata->userelerrorforreliability, TRUE, DEFAULT_USERELERRORFORRELIABILITY,
         NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/relpscost/lowerrortol", "low relative error tolerance for reliability",
         &branchruledata->lowerrortol, TRUE, DEFAULT_LOWERRORTOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/relpscost/higherrortol", "high relative error tolerance for reliability",
         &branchruledata->higherrortol, TRUE, DEFAULT_HIGHERRORTOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/relpscost/storesemiinitcosts",
         "should strong branching result be considered for pseudo costs if the other direction was infeasible?",
         &branchruledata->storesemiinitcosts, TRUE, DEFAULT_STORESEMIINITCOSTS,
         NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/relpscost/usesblocalinfo",
         "should the scoring function use only local cutoff and inference information obtained for strong branching candidates?",
         &branchruledata->usesblocalinfo, TRUE, DEFAULT_USESBLOCALINFO,
         NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/relpscost/usehyptestforreliability",
         "should the strong branching decision be based on a hypothesis test?",
         &branchruledata->usehyptestforreliability, TRUE, DEFAULT_USEHYPTESTFORRELIABILITY,
         NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/relpscost/usedynamicconfidence",
         "should the confidence level be adjusted dynamically?",
         &branchruledata->usedynamicconfidence, TRUE, DEFAULT_USEDYNAMICCONFIDENCE,
         NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "branching/relpscost/skipbadinitcands",
         "should branching rule skip candidates that have a low probability to "
         "be better than the best strong-branching or pseudo-candidate?",
         &branchruledata->skipbadinitcands, TRUE, DEFAULT_SKIPBADINITCANDS,
         NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/relpscost/confidencelevel",
         "the confidence level for statistical methods, between 0 (Min) and 4 (Max).",
         &branchruledata->confidencelevel, TRUE, DEFAULT_CONFIDENCELEVEL, 0, 4, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/relpscost/randinitorder",
         "should candidates be initialized in randomized order?",
         &branchruledata->randinitorder, TRUE, DEFAULT_RANDINITORDER,
         NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/relpscost/usesmallweightsitlim",
         "should smaller weights be used for pseudo cost updates after hitting the LP iteration limit?",
         &branchruledata->usesmallweightsitlim, TRUE, DEFAULT_USESMALLWEIGHTSITLIM,
         NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/relpscost/dynamicweights",
         "should the weights of the branching rule be adjusted dynamically during solving based on objective and infeasible leaf counters?",
         &branchruledata->dynamicweights, TRUE, DEFAULT_DYNAMICWEIGHTS,
         NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/relpscost/startrandseed", "start seed for random number generation",
         &branchruledata->startrandseed, TRUE, DEFAULT_STARTRANDSEED, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** execution reliability pseudo cost branching with the given branching candidates */
SCIP_RETCODE SCIPexecRelpscostBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            branchcands,        /**< branching candidates */
   SCIP_Real*            branchcandssol,     /**< solution value for the branching candidates */
   SCIP_Real*            branchcandsfrac,    /**< fractional part of the branching candidates */
   int                   nbranchcands,       /**< number of branching candidates */
   SCIP_Bool             executebranching,   /**< perform a branching step after probing */
   SCIP_RESULT*          result              /**< pointer to the result of the execution */
   )
{
   SCIP_BRANCHRULE* branchrule;

   assert(scip != NULL);
   assert(result != NULL);

   /* find branching rule */
   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   /* execute branching rule */
   SCIP_CALL( execRelpscost(scip, branchrule, branchcands, branchcandssol, branchcandsfrac, nbranchcands, executebranching, result) );

   return SCIP_OKAY;
}
