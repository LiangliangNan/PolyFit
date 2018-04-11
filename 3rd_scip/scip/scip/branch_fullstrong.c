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

/**@file   branch_fullstrong.c
 * @brief  full strong LP branching rule
 * @author Tobias Achterberg
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_fullstrong.h"


#define BRANCHRULE_NAME          "fullstrong"
#define BRANCHRULE_DESC          "full strong branching"
#define BRANCHRULE_PRIORITY      0
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define DEFAULT_REEVALAGE        10LL        /**< number of intermediate LPs solved to trigger reevaluation of strong branching
                                              *   value for a variable that was already evaluated at the current node */
#define DEFAULT_MAXPROPROUNDS      -2        /**< maximum number of propagation rounds to be performed during strong branching
                                              *   before solving the LP (-1: no limit, -2: parameter settings) */
#define DEFAULT_PROBINGBOUNDS    TRUE        /**< should valid bounds be identified in a probing-like fashion during strong
                                              *   branching (only with propagation)? */
#define DEFAULT_FORCESTRONGBRANCH FALSE      /**< should strong branching be applied even if there is just a single candidate? */


/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Longint          reevalage;          /**< number of intermediate LPs solved to trigger reevaluation of strong branching
                                              *   value for a variable that was already evaluated at the current node */
   int                   maxproprounds;      /**< maximum number of propagation rounds to be performed during strong branching
                                              *   before solving the LP (-1: no limit, -2: parameter settings) */
   SCIP_Bool             probingbounds;      /**< should valid bounds be identified in a probing-like fashion during strong
                                              *   branching (only with propagation)? */
   SCIP_Bool             forcestrongbranch;  /**< should strong branching be applied even if there is just a single candidate? */
   int                   lastcand;           /**< last evaluated candidate of last branching rule execution */
   int                   skipsize;           /**< size of skipdown and skipup array */
   SCIP_Bool*            skipdown;           /**< should be branching on down child be skipped? */
   SCIP_Bool*            skipup;             /**< should be branching on up child be skipped? */
};


/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyFullstrong)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleFullstrong(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeFullstrong)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &branchruledata->skipdown, branchruledata->skipsize);
   SCIPfreeBlockMemoryArrayNull(scip, &branchruledata->skipup, branchruledata->skipsize);

   SCIPfreeBlockMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitFullstrong)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* initialize branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   branchruledata->lastcand = 0;

   return SCIP_OKAY;
}

/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitFullstrong)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* initialize branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert((branchruledata->skipdown != NULL) == (branchruledata->skipup != NULL));

   if( branchruledata->skipdown != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &branchruledata->skipup, branchruledata->skipsize);
      SCIPfreeBlockMemoryArray(scip, &branchruledata->skipdown, branchruledata->skipsize);
      branchruledata->skipdown = NULL;
      branchruledata->skipup = NULL;
      branchruledata->skipsize = 0;
   }

   return SCIP_OKAY;
}

/**
 * Selects a variable from a set of candidates by strong branching
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 * @note The variables in the lpcands array must have a fractional value in the current LP solution
 */
SCIP_RETCODE SCIPselectVarStrongBranching(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP_VAR**            lpcands,            /**< branching candidates                                */
   SCIP_Real*            lpcandssol,         /**< solution values of the branching candidates         */
   SCIP_Real*            lpcandsfrac,        /**< fractional values of the branching candidates       */
   SCIP_Bool*            skipdown,           /**< should down branchings be skipped? */
   SCIP_Bool*            skipup,             /**< should up branchings be skipped? */
   int                   nlpcands,           /**< number of branching candidates                      */
   int                   npriolpcands,       /**< number of priority branching candidates             */
   int                   ncomplete,          /**< number of branching candidates without skip         */
   int*                  start,              /**< starting index in lpcands                           */
   int                   maxproprounds,      /**< maximum number of propagation rounds to be performed during strong
                                              *   branching before solving the LP (-1: no limit, -2: parameter settings) */
   SCIP_Bool             probingbounds,      /**< should valid bounds be identified in a probing-like fashion during
                                              *   strong branching (only with propagation)? */
   SCIP_Bool             forcestrongbranch,  /**< should strong branching be applied even if there is just a single candidate? */
   int*                  bestcand,           /**< best candidate for branching                        */
   SCIP_Real*            bestdown,           /**< objective value of the down branch for bestcand     */
   SCIP_Real*            bestup,             /**< objective value of the up branch for bestcand       */
   SCIP_Real*            bestscore,          /**< score for bestcand                                  */
   SCIP_Bool*            bestdownvalid,      /**< is bestdown a valid dual bound for the down branch? */
   SCIP_Bool*            bestupvalid,        /**< is bestup a valid dual bound for the up branch?     */
   SCIP_Real*            provedbound,        /**< proved dual bound for current subtree               */
   SCIP_RESULT*          result              /**< result pointer                                      */
   )
{  /*lint --e{715}*/
   SCIP_VAR** vars = NULL;
   SCIP_Real* newlbs = NULL;
   SCIP_Real* newubs = NULL;
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_Longint reevalage;
   SCIP_Longint nodenum;
   SCIP_Real down;
   SCIP_Real up;
   SCIP_Real downgain;
   SCIP_Real upgain;
   SCIP_Real score;
   SCIP_Real lpobjval;
   SCIP_Bool exactsolve;
   SCIP_Bool lperror;
   SCIP_Bool allcolsinlp;
   SCIP_Bool downvalid;
   SCIP_Bool upvalid;
   SCIP_Bool downinf;
   SCIP_Bool upinf;
   SCIP_Bool downconflict;
   SCIP_Bool upconflict;
   SCIP_Bool bothgains;
   SCIP_Bool propagate;
   int nvars = 0;
   int nsbcalls;
   int i;
   int c;

   assert(scip != NULL);
   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(lpcandsfrac != NULL);
   assert(bestcand != NULL);
   assert(skipdown != NULL);
   assert(skipup != NULL);
   assert(bestdown != NULL);
   assert(bestup != NULL);
   assert(bestscore != NULL);
   assert(bestdownvalid != NULL);
   assert(bestupvalid != NULL);
   assert(provedbound != NULL);
   assert(result != NULL);
   assert(nlpcands > 0);

   /* check, if we want to solve the problem exactly, meaning that strong branching information is not useful
    * for cutting off sub problems and improving lower bounds of children
    */
   exactsolve = SCIPisExactSolve(scip);

   /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
   allcolsinlp = SCIPallColsInLP(scip);

   /* get current node number */
   nodenum = SCIPgetNNodes(scip);

   /* get current LP objective bound of the local sub problem and global cutoff bound */
   lpobjval = SCIPgetLPObjval(scip);
   *provedbound = lpobjval;

   *bestcand = 0;
   *bestdown = lpobjval;
   *bestup = lpobjval;
   *bestdownvalid = TRUE;
   *bestupvalid = TRUE;
   *bestscore = -SCIPinfinity(scip);

   /* if only one candidate exists, choose this one without applying strong branching; also, when SCIP is about to be
    * stopped, all strongbranching evaluations will be aborted anyway, thus we can return immediately
    */
   if( (!forcestrongbranch && nlpcands == 1) || SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* this assert may not hold if SCIP is stopped, thus we only check it here */
   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

   /* get branching rule */
   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* auto-setting for reevalage */
   reevalage = branchruledata->reevalage;

   /* check whether propagation should be performed */
   propagate = (maxproprounds != 0);

   /* if we don't do propagation, we cannot identify valid bounds in a probing-like fashion */
   if( !propagate )
      probingbounds = FALSE;

   /* create arrays for probing-like bound tightening */
   if( probingbounds )
   {
      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      SCIP_CALL( SCIPallocBufferArray(scip, &newlbs, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newubs, nvars) );
   }

    /* initialize strong branching */
   SCIP_CALL( SCIPstartStrongbranch(scip, propagate) );

   /* search the full strong candidate
    * cycle through the candidates, starting with the position evaluated in the last run
    */
   nsbcalls = 0;
   bothgains = TRUE;
   for( i = 0, c = *start; i < nlpcands && (!bothgains || i < ncomplete); ++i, ++c )
   {
      c = c % nlpcands;
      assert(lpcands[c] != NULL);

      /* don't use strong branching on variables that have already been initialized at the current node,
       * and that were evaluated not too long ago
       */
      if( SCIPgetVarStrongbranchNode(scip, lpcands[c]) == nodenum
         && SCIPgetVarStrongbranchLPAge(scip, lpcands[c]) < reevalage )
      {
         SCIP_Real lastlpobjval;

         /* use the score of the strong branching call at the current node */
         SCIP_CALL( SCIPgetVarStrongbranchLast(scip, lpcands[c], &down, &up, NULL, NULL, NULL, &lastlpobjval) );
         downgain = MAX(down - lastlpobjval, 0.0);
         upgain = MAX(up - lastlpobjval, 0.0);
         downvalid = FALSE;
         upvalid = FALSE;
         downinf = FALSE;
         upinf = FALSE;
         downconflict = FALSE;
         upconflict = FALSE;
         lperror = FALSE;
         SCIPdebugMsg(scip, "strong branching on variable <%s> already performed (lpage=%" SCIP_LONGINT_FORMAT ", down=%g (%+g), up=%g (%+g))\n",
            SCIPvarGetName(lpcands[c]), SCIPgetVarStrongbranchLPAge(scip, lpcands[c]), down, downgain, up, upgain);
      }
      else
      {
         SCIPdebugMsg(scip, "applying strong branching%s on variable <%s> with solution %g\n",
            propagate ? "with propagation" : "", SCIPvarGetName(lpcands[c]), lpcandssol[c]);
         assert(i >= ncomplete || (!skipdown[i] && !skipup[i]));

         /* apply strong branching */
         up = -SCIPinfinity(scip);
         down = -SCIPinfinity(scip);

         if( propagate )
         {
            SCIP_CALL( SCIPgetVarStrongbranchWithPropagation(scip, lpcands[c], lpcandssol[c], lpobjval, INT_MAX,
                  maxproprounds, skipdown[i] ? NULL : &down, skipup[i] ? NULL : &up, &downvalid,
                  &upvalid, NULL, NULL, &downinf, &upinf, &downconflict, &upconflict, &lperror, newlbs, newubs) );

            SCIPdebugMsg(scip, "-> down=%.9g (gain=%.9g, valid=%u, inf=%u, conflict=%u), up=%.9g (gain=%.9g, valid=%u, inf=%u, conflict=%u)\n",
               down, down - lpobjval, downvalid, downinf, downconflict, up, up - lpobjval, upvalid, upinf, upconflict);
         }
         else
         {
            SCIP_CALL( SCIPgetVarStrongbranchFrac(scip, lpcands[c], INT_MAX,
                  skipdown[i] ? NULL : &down, skipup[i] ? NULL : &up, &downvalid, &upvalid, &downinf, &upinf,
                  &downconflict, &upconflict, &lperror) );
         }
         nsbcalls++;

         /* display node information line */
         if( SCIPgetDepth(scip) == 0 && nsbcalls % 100 == 0 )
         {
            SCIP_CALL( SCIPprintDisplayLine(scip, NULL, SCIP_VERBLEVEL_HIGH, TRUE) );
         }

         /* check for an error in strong branching */
         if( lperror )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "(node %" SCIP_LONGINT_FORMAT ") error in strong branching call%s for variable <%s> with solution %g\n",
               SCIPgetNNodes(scip), propagate ? " with propagation" : "", SCIPvarGetName(lpcands[c]), lpcandssol[c]);
            break;
         }

         /* evaluate strong branching */
         down = MAX(down, lpobjval);
         up = MAX(up, lpobjval);
         downgain = down - lpobjval;
         upgain = up - lpobjval;
         if( !SCIPisFeasZero(scip,downgain) && !SCIPisFeasZero(scip,upgain) )
            bothgains = TRUE;

         assert(!allcolsinlp || exactsolve || !downvalid || downinf == SCIPisGE(scip, down, SCIPgetCutoffbound(scip)));
         assert(!allcolsinlp || exactsolve || !upvalid || upinf == SCIPisGE(scip, up, SCIPgetCutoffbound(scip)));
         assert(downinf || !downconflict);
         assert(upinf || !upconflict);

         /* check if there are infeasible roundings */
         if( downinf || upinf )
         {
            /* if we didn't do propagation, we can only detect infeasibility if the LP is a valid relaxation */
            assert(allcolsinlp || propagate);
            assert(!exactsolve);

            if( downinf && upinf )
            {
               /* both roundings are infeasible -> node is infeasible */
               *result = SCIP_CUTOFF;
               SCIPdebugMsg(scip, " -> variable <%s> is infeasible in both directions\n", SCIPvarGetName(lpcands[c]));
               break; /* terminate initialization loop, because node is infeasible */
            }
            else if( downinf )
            {
               SCIP_Bool infeasible;
               SCIP_Bool tightened;

               /* downwards rounding is infeasible -> change lower bound of variable to upward rounding */
               SCIP_CALL( SCIPtightenVarLb(scip, lpcands[c], SCIPfeasCeil(scip, lpcandssol[c]), TRUE, &infeasible, &tightened) );
               assert(!infeasible);

               /* if we did propagation, the bound change might already have been added */
               assert(tightened || propagate);

               *result = SCIP_REDUCEDDOM;
               SCIPdebugMsg(scip, " -> variable <%s> is infeasible in downward branch\n", SCIPvarGetName(lpcands[c]));
               break; /* terminate initialization loop, because LP was changed */
            }
            else
            {
               SCIP_Bool infeasible;
               SCIP_Bool tightened;

               /* upwards rounding is infeasible -> change upper bound of variable to downward rounding */
               assert(upinf);
               SCIP_CALL( SCIPtightenVarUb(scip, lpcands[c], SCIPfeasFloor(scip, lpcandssol[c]), TRUE, &infeasible, &tightened) );
               assert(!infeasible);

               /* if we did propagation, the bound change might already have been added */
               assert(tightened || propagate);

               *result = SCIP_REDUCEDDOM;
               SCIPdebugMsg(scip, " -> variable <%s> is infeasible in upward branch\n", SCIPvarGetName(lpcands[c]));
               break; /* terminate initialization loop, because LP was changed */
            }
         }
         else if( allcolsinlp && !exactsolve && downvalid && upvalid )
         {
            SCIP_Real minbound;

            /* the minimal lower bound of both children is a proved lower bound of the current subtree */
            minbound = MIN(down, up);
            *provedbound = MAX(*provedbound, minbound);

            /* apply probing-like bounds detected during strong branching */
            if( probingbounds )
            {
               int nboundchgs;
               int v;

               assert(vars != NULL);
               assert(nvars > 0);
               assert(newlbs != NULL);
               assert(newubs != NULL);

               nboundchgs = 0;

               for( v = 0; v < nvars; ++v )
               {
                  if( SCIPisGT(scip, newlbs[v], SCIPvarGetLbLocal(vars[v])) )
                  {
                     SCIPdebugMsg(scip, "better lower bound for variable <%s>: %.9g -> %.9g (strongbranching on var <%s>\n",
                        SCIPvarGetName(vars[v]), SCIPvarGetLbLocal(vars[v]), newlbs[v], SCIPvarGetName(lpcands[c]));

                     SCIP_CALL( SCIPchgVarLb(scip, vars[v], newlbs[v]) );
                     ++nboundchgs;
                  }
                  if( SCIPisLT(scip, newubs[v], SCIPvarGetUbLocal(vars[v])) )
                  {
                     SCIPdebugMsg(scip, "better upper bound for variable <%s>: %.9g -> %.9g (strongbranching on var <%s>\n",
                        SCIPvarGetName(vars[v]), SCIPvarGetUbLocal(vars[v]), newubs[v], SCIPvarGetName(lpcands[c]));

                     SCIP_CALL( SCIPchgVarUb(scip, vars[v], newubs[v]) );
                     ++nboundchgs;
                  }
               }

               if( nboundchgs > 0 )
               {
                  *result = SCIP_REDUCEDDOM;
                  SCIPdebugMsg(scip, " -> strong branching with propagation on variable <%s> led to %d bound changes\n",
                     SCIPvarGetName(lpcands[c]), nboundchgs);
                  break; /* terminate initialization loop, because LP was changed */
               }
            }
         }

         /* update pseudo cost values */
         assert(!downinf); /* otherwise, we would have terminated the initialization loop */
         assert(!upinf);
         SCIP_CALL( SCIPupdateVarPseudocost(scip, lpcands[c], 0.0-lpcandsfrac[c], downgain, 1.0) );
         SCIP_CALL( SCIPupdateVarPseudocost(scip, lpcands[c], 1.0-lpcandsfrac[c], upgain, 1.0) );
      }

      /* check for a better score, if we are within the maximum priority candidates */
      if( c < npriolpcands )
      {
         score = SCIPgetBranchScore(scip, lpcands[c], downgain, upgain);
         if( score > *bestscore )
         {
            *bestcand = c;
            *bestdown = down;
            *bestup = up;
            *bestdownvalid = downvalid;
            *bestupvalid = upvalid;
            *bestscore = score;
         }
      }
      else
      {
         SCIPdebug(score = 0.0;)
      }

      SCIPdebugMsg(scip, " -> cand %d/%d (prio:%d) var <%s> (solval=%g, downgain=%g, upgain=%g, score=%g) -- best: <%s> (%g)\n",
         c, nlpcands, npriolpcands, SCIPvarGetName(lpcands[c]), lpcandssol[c], downgain, upgain, score,
         SCIPvarGetName(lpcands[*bestcand]), *bestscore);
   }

   /* end strong branching */
   SCIP_CALL( SCIPendStrongbranch(scip) );

   *start = c;

   if( probingbounds )
   {
      assert(newlbs != NULL);
      assert(newubs != NULL);

      SCIPfreeBufferArray(scip, &newlbs);
      SCIPfreeBufferArray(scip, &newubs);
   }

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpFullstrong)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** tmplpcands;
   SCIP_VAR** lpcands;
   SCIP_Real* tmplpcandssol;
   SCIP_Real* lpcandssol;
   SCIP_Real* tmplpcandsfrac;
   SCIP_Real* lpcandsfrac;
   SCIP_Real bestdown;
   SCIP_Real bestup;
   SCIP_Real bestscore;
   SCIP_Real provedbound;
   SCIP_Bool bestdownvalid;
   SCIP_Bool bestupvalid;
   int nlpcands;
   int npriolpcands;
   int bestcand;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Execlp method of fullstrong branching\n");

   *result = SCIP_DIDNOTRUN;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &tmplpcands, &tmplpcandssol, &tmplpcandsfrac, &nlpcands, &npriolpcands, NULL) );
   assert(nlpcands > 0);
   assert(npriolpcands > 0);

   /* copy LP banching candidates and solution values, because they will be updated w.r.t. the strong branching LP
    * solution
    */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcands, tmplpcands, nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandssol, tmplpcandssol, nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandsfrac, tmplpcandsfrac, nlpcands) );

   if( branchruledata->skipdown == NULL )
   {
      assert(branchruledata->skipup == NULL);

      branchruledata->skipsize = SCIPgetNVars(scip);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->skipdown, branchruledata->skipsize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->skipup, branchruledata->skipsize) );
      BMSclearMemoryArray(branchruledata->skipdown, branchruledata->skipsize);
      BMSclearMemoryArray(branchruledata->skipup, branchruledata->skipsize);
   }

   SCIP_CALL( SCIPselectVarStrongBranching(scip, lpcands, lpcandssol, lpcandsfrac, branchruledata->skipdown,
         branchruledata->skipup, nlpcands, npriolpcands, nlpcands, &branchruledata->lastcand,
         branchruledata->maxproprounds, branchruledata->probingbounds, branchruledata->forcestrongbranch, &bestcand,
         &bestdown, &bestup, &bestscore, &bestdownvalid, &bestupvalid, &provedbound, result) );

   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED )
   {
      SCIP_NODE* downchild;
      SCIP_NODE* upchild;
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_Bool allcolsinlp;
      SCIP_Bool exactsolve;

      assert(*result == SCIP_DIDNOTRUN);
      assert(0 <= bestcand && bestcand < nlpcands);
      assert(SCIPisLT(scip, provedbound, SCIPgetCutoffbound(scip)));

      var = lpcands[bestcand];
      val = lpcandssol[bestcand];

      /* perform the branching */
      SCIPdebugMsg(scip, " -> %d candidates, selected candidate %d: variable <%s> (solval=%g, down=%g, up=%g, score=%g)\n",
         nlpcands, bestcand, SCIPvarGetName(var), lpcandssol[bestcand], bestdown, bestup, bestscore);
      SCIP_CALL( SCIPbranchVarVal(scip, var, val, &downchild, NULL, &upchild) );
      assert(downchild != NULL);
      assert(upchild != NULL);

      /* check, if we want to solve the problem exactly, meaning that strong branching information is not useful
       * for cutting off sub problems and improving lower bounds of children
       */
      exactsolve = SCIPisExactSolve(scip);

      /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
      allcolsinlp = SCIPallColsInLP(scip);

      /* update the lower bounds in the children */
      if( allcolsinlp && !exactsolve )
      {
         SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, bestdownvalid ? MAX(bestdown, provedbound) : provedbound) );
         SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, bestupvalid ? MAX(bestup, provedbound) : provedbound) );
      }
      SCIPdebugMsg(scip, " -> down child's lowerbound: %g\n", SCIPnodeGetLowerbound(downchild));
      SCIPdebugMsg(scip, " -> up child's lowerbound: %g\n", SCIPnodeGetLowerbound(upchild));

      *result = SCIP_BRANCHED;
   }

   SCIPfreeBufferArray(scip, &lpcandsfrac);
   SCIPfreeBufferArray(scip, &lpcandssol);
   SCIPfreeBufferArray(scip, &lpcands);

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the full strong LP branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleFullstrong(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create fullstrong branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );
   branchruledata->lastcand = 0;
   branchruledata->skipsize = 0;
   branchruledata->skipup = NULL;
   branchruledata->skipdown = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non-fundamental callbacks via specific setter functions*/
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyFullstrong) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeFullstrong) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitFullstrong) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitFullstrong) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpFullstrong) );

   /* fullstrong branching rule parameters */
   SCIP_CALL( SCIPaddLongintParam(scip,
         "branching/fullstrong/reevalage",
         "number of intermediate LPs solved to trigger reevaluation of strong branching value for a variable that was already evaluated at the current node",
         &branchruledata->reevalage, TRUE, DEFAULT_REEVALAGE, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/fullstrong/maxproprounds",
         "maximum number of propagation rounds to be performed during strong branching before solving the LP (-1: no limit, -2: parameter settings)",
         &branchruledata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -2, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/fullstrong/probingbounds",
         "should valid bounds be identified in a probing-like fashion during strong branching (only with propagation)?",
         &branchruledata->probingbounds, TRUE, DEFAULT_PROBINGBOUNDS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/fullstrong/forcestrongbranch",
         "should strong branching be applied even if there is just a single candidate?",
         &branchruledata->forcestrongbranch, TRUE, DEFAULT_FORCESTRONGBRANCH, NULL, NULL) );

   return SCIP_OKAY;
}
