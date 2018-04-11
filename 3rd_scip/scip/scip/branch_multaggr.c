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
/**@file   branch_multaggr.c
 * @brief  fullstrong branching on fractional and multi-aggregated variables
 * @author Anna Melchiori
 * @author Gerald Gamrath
 *
 * This branching rule uses all fractional binary and integer variables as candidates,
 * as well as fractional multiaggregated binary and integer variables. Although not
 * directly contained in the presolved problem anymore, the multi-aggregation provides
 * an affine linear sum of integer variables, on which branching can be performed.
 *
 * For more details, see
 * G.Gamrath, A.Melchiori, T.Berthold, A.M.Gleixner, D.Salvagnin: Branching on Multi-aggregated Variables
 * (http://dx.doi.org/10.1007/978-3-319-18008-3_10)
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_multaggr.h"
#include "scip/branch_fullstrong.h"
#include "scip/cons_linear.h"
#include "scip/var.h"
#include "scip/set.h"
#include "scip/pub_tree.h"
#include "scip/struct_scip.h"
#include "scip/clock.h"

#define BRANCHRULE_NAME            "multaggr"
#define BRANCHRULE_DESC            "fullstrong branching on fractional and multi-aggregated variables"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0


#define DEFAULT_REEVALAGE         0LL        /**< number of intermediate LPs solved to trigger reevaluation of strong branching
                                              *   value for a variable that was already evaluated at the current node */
#define DEFAULT_MAXPROPROUNDS       0        /**< maximum number of propagation rounds to be performed during multaggr branching
                                              *   before solving the LP (-1: no limit, -2: parameter settings) */
#define DEFAULT_PROBINGBOUNDS    TRUE        /**< should valid bounds be identified in a probing-like fashion during multi-aggr
                                              *   branching (only with propagation)? */

/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Longint          reevalage;          /**< number of intermediate LPs solved to trigger reevaluation of strong branching
                                              *   value for a variable that was already evaluated at the current node */
   SCIP_Bool             probingbounds;      /**< should valid bounds be identified in a probing-like fashion during strong
                                              *   branching (only with propagation)? */
   int                   lastcand;           /**< last evaluated candidate of last branching rule execution */
   int                   maxproprounds;      /**< maximum number of propagation rounds to be performed during strong branching
                                              *   before solving the LP (-1: no limit, -2: parameter settings) */
   int                   skipsize;           /**< size of skipdown and skipup array */
   SCIP_Bool*            skipdown;           /**< should be branching on down child be skipped? */
   SCIP_Bool*            skipup;             /**< should be branching on up child be skipped? */
#ifdef SCIP_STATISTIC
   SCIP_CLOCK*           clckstrongbr;       /**< clock to store the time spent inside the strong branching function on fractional variables */
   SCIP_CLOCK*           clckmultaggrbr;     /**< clock to store the time spent inside the strong branching function on multi-aggragated variables */
   SCIP_Real*            ratioggain;         /**< for each occurence of a branching on a multi-aggregated variable we store the ratio of the gain that
                                              *   we would have obtained branching on the best fractional variable over the gain obtained
                                              *   branching on the current multi-aggregated variable */
   SCIP_Real             ameanratio;         /**< arithmetic mean of the ratioggain array */
   SCIP_Bool             noupdate;           /**< pointer to store if the probing LP has not been solved so we do not want to
                                              *   update statistics */
   int                   firstmultaggrdepth; /**< depth of the first branching on a multi-aggregated variable */
   int                   rundepth;           /**< the run of the first multi-aggregated branching */
   int                   nmultaggrbranch;    /**< number of branchings on multi-aggregated variables */
   int                   nfracbranch;        /**< number of branchings on fractional variables */
   int                   nEscore;            /**< number of times that the bestscore over all multi-aggregated variables is equal to the best
                                              *   fractional variables score and thus we do not branch on the multi-aggregate variable */
   int                   nmultaggrcutoff;    /**< number of cutoffs detected during the probing mode on multi-aggregated variables */
   int                   nmultaggrconsadd;   /**< number of times that a probing constraint of a multi-aggregated variable has been
                                              *   added to the original problem */
   int                   nfractcutoff;       /**< number of cutoffs detected during strong branching on fractional variables */
   int                   nfractconsadd;      /**< number of times that during strong branching on fractional variables a constraint has been
                                              *   added to the original problem or a variables domain has been reduced */
   int                   nmultaggrvars;      /**< number of multi-aggregated variables in the problem of the last run */
   int                   nrun;               /**< number of restarts */
   int                   size;               /**< size of the provided array to store the ratio gain */
   int                   nstrongbrcall;      /**< number of times that the selectVarstrongBranching function has been called */
   int                   nmultaggrbrcall;    /**< number of times that the selectVarMultAggrBranching function has been called */
   int                   totallpcands;       /**< total number of observed lpcands over all selectVarstrongBranching function calls */
   int                   totalmultaggrcands; /**< total number of observed multi-aggregregated candidates over all selectVarMultAggrBranching
                                              *   function calls */
#endif
};


/*
 * Local methods
 */

/* this function ensures that the allocated memory is enough to store statistics data */
#ifdef SCIP_STATISTIC
static
SCIP_RETCODE ensureArraySize(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(branchruledata->ratioggain != NULL);
   assert(branchruledata->nmultaggrbranch >= 0);
   assert(branchruledata->size >= 0);

   /* check whether the size of the array is big enough; reallocate memory if needed */
   if( branchruledata->nmultaggrbranch + 1 > branchruledata->size )
   {
      int newsize = SCIPcalcMemGrowSize(scip, branchruledata->nmultaggrbranch + 1);
      assert(newsize >= branchruledata->nmultaggrbranch + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &branchruledata->ratioggain, branchruledata->size, newsize) );
      branchruledata->size = newsize;
   }
   return SCIP_OKAY;
}
#endif

/* this function gives us the best candidate for branching among the multi-aggregated variables of the problem
 * and the best fractional integer variable already selected by strong branching
*/
static
SCIP_RETCODE selectVarMultAggrBranching(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_VAR**            bestcand,           /**< the best candidate variable selected by strong branching */
   SCIP_Real*            bestscore,          /**< score of the best branching candidate */
   SCIP_Real*            bestsol,            /**< solution value of the best branching candidate */
   SCIP_Real*            bestdown,           /**< objective value of the down node when branching on bestcand */
   SCIP_Real*            bestup,             /**< objective value of the up node when branching on bestcand */
   SCIP_Bool*            bestdownvalid,      /**< is bestdown a valid dual bound for the down branch? */
   SCIP_Bool*            bestupvalid,        /**< is bestup a valid dual bound for the up branch? */
   SCIP_Real*            provedbound,        /**< proved dual bound for the current subtree */
   SCIP_Real*            estimatedown,       /**< pointer to store the down child nodes estimate */
   SCIP_Real*            estimateup,         /**< pointer to store the up child nodes estimate */
#ifdef SCIP_STATISTIC
   SCIP_Real*            bestmultaggrscore,  /**< pointer to store the multi aggregated score */
#endif
   SCIP_RESULT*          result              /**< pointer to store results of branching */
   )
{
   SCIP_VAR** fixvars;
   SCIP_CONS* probingconsdown;
   SCIP_CONS* probingconsup;
   SCIP_NODE* node;
   SCIP_Real* fixvarssols;
   SCIP_Real fixvarssol;
   SCIP_Real lpobjval;
   SCIP_Bool exactsolve;
   SCIP_Bool allcolsinlp;
   SCIP_Bool downnodeinf = FALSE;
   SCIP_Bool startprobing = TRUE;
   SCIP_Bool endprobing = FALSE;
   int nfixvars;
   int i;
   int j;
   int k;

   /* import branchrule data for statistics */
#ifdef SCIP_STATISTIC
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
#endif

   assert(scip != NULL);
   assert(bestcand != NULL);
   assert(bestscore != NULL);

   /* check, if we want to solve the problem exactly, meaning that strong branching information is not useful
    * for cutting off sub problems and improving lower bounds of children
    */
   exactsolve = SCIPisExactSolve(scip);

   /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
   allcolsinlp = SCIPallColsInLP(scip);

   /* get fixed variables */
   fixvars = SCIPgetFixedVars(scip);
   nfixvars = SCIPgetNFixedVars(scip);
   SCIPdebugMsg(scip, " fractional variable: <%s> with value: %f is selected by strong branching\n", SCIPvarGetName(*bestcand), *bestsol);

   /* check if we would exceed the depth limit */
   if( SCIP_MAXTREEDEPTH <= SCIPgetDepth(scip) )
   {
      SCIPdebugMsg(scip, "cannot perform probing in selectVarMultAggrBranching, depth limit reached.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( nfixvars != 0 )
   {
      assert(fixvars != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &fixvarssols, nfixvars) );
      lpobjval = SCIPgetLPObjval(scip);

      /* store the values of the fixed variables at the current optimal solution */
      for( i = 0; i < nfixvars; i++ )
      {
         assert(fixvars[i] != NULL);
         fixvarssols[i] = SCIPvarGetLPSol(fixvars[i]);
      }

      for( i = 0; i < nfixvars; i++ )
      {
         assert(fixvars[i] != NULL);

         /* only integer and binary multi-aggregated variables are potential branching candidates */
         if( SCIPvarGetStatus(fixvars[i]) == SCIP_VARSTATUS_MULTAGGR && (SCIPvarGetType(fixvars[i]) == SCIP_VARTYPE_INTEGER ||
               SCIPvarGetType(fixvars[i]) == SCIP_VARTYPE_BINARY) )
         {
            fixvarssol = fixvarssols[i];

            /* start probing mode for the fractional multi-aggregated variable */
            if( !SCIPisFeasIntegral(scip, fixvarssol) )
            {
               SCIP_VAR** downvars = NULL;
               SCIP_VAR** upvars = NULL;
               SCIP_Real* downvarssols = NULL;
               SCIP_Real* upvarssols = NULL;
               SCIP_LPSOLSTAT solstatdown;
               SCIP_LPSOLSTAT solstatup;
               SCIP_Real downobjval;
               SCIP_Real upobjval;
               SCIP_Real estimateprobdown = 0.0;
               SCIP_Real estimateprobup = 0.0;
               SCIP_Bool downinf;
               SCIP_Bool upinf;
               SCIP_Bool lperror;
               int ndownvars;
               int nupvars;

               /* start the probing mode if this is the first entrance */
               if( startprobing )
               {
                  SCIP_CALL( SCIPstartProbing(scip) );
                  startprobing = FALSE;
                  endprobing = TRUE;

                  SCIPdebugMsg(scip, "PROBING MODE:\n");
               }

               SCIPdebugMsg(scip, " multi-aggregated variable: <%s> with value: %f\n", SCIPvarGetName(fixvars[i]), fixvarssol);

               SCIPstatistic(branchruledata->totalmultaggrcands += 1);

               /* create the multi-aggregated rounded down constraint */
               SCIP_CALL( SCIPcreateConsLinear(scip, &probingconsdown, "probingconsdown", SCIPvarGetMultaggrNVars(fixvars[i]),
                     SCIPvarGetMultaggrVars(fixvars[i]), SCIPvarGetMultaggrScalars(fixvars[i]), -SCIPinfinity(scip),
                     SCIPfeasFloor(scip, fixvarssol) - SCIPvarGetMultaggrConstant(fixvars[i]), TRUE, TRUE, FALSE, FALSE,
                     TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
               assert(probingconsdown != NULL);

               /* create the down child probing node */
               SCIP_CALL( SCIPnewProbingNode(scip) );
               node = SCIPgetCurrentNode(scip);
               assert(node != NULL);

               SCIP_CALL( SCIPaddConsNode(scip, node, probingconsdown, NULL) );
               SCIP_CALL( SCIPreleaseCons(scip, &probingconsdown) );

#ifdef PRINTNODECONS
               SCIPdebugMsg(scip, " created down probing node with constraint:\n");
               SCIP_CALL( SCIPprintCons(scip, probingconsdown, NULL) );
               SCIPinfoMessage(scip, NULL, "\n");
#endif

               /* solve the down child probing node */
               SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &downinf) );
               solstatdown = SCIPgetLPSolstat(scip);
               lperror = lperror || (solstatdown == SCIP_LPSOLSTAT_NOTSOLVED && downinf == 0) || (solstatdown == SCIP_LPSOLSTAT_ITERLIMIT) ||
                  (solstatdown == SCIP_LPSOLSTAT_TIMELIMIT);
               assert(solstatdown != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

               /* break the branching rule if an error occurred, problem was not solved, iteration or time limit was reached */
               if( lperror )
               {
                  SCIPdebugMsg(scip, "error solving down node probing LP: status=%d\n", solstatdown);
                  SCIPstatistic(branchruledata->noupdate = TRUE);
                  break;
               }

               downobjval = SCIPgetLPObjval(scip);
               downinf = downinf || SCIPisGE(scip, downobjval, SCIPgetCutoffbound(scip));
               assert(((solstatdown != SCIP_LPSOLSTAT_INFEASIBLE) && (solstatdown != SCIP_LPSOLSTAT_OBJLIMIT)) || downinf);

               if( !downinf )
               {
                  /* when an optimal solution has been found calculate down child's estimate based on pseudo costs */
                  /* estimate = lowerbound + sum(min{f_j * pscdown_j, (1-f_j) * pscup_j}) */
                  estimateprobdown = SCIPnodeGetLowerbound(node);
                  SCIP_CALL( SCIPgetLPBranchCands(scip, &downvars, &downvarssols, NULL, &ndownvars, NULL, NULL) );

                  for( j = 0 ; j < ndownvars; j++ )
                  {
                     SCIP_Real estimateincr;
                     SCIP_Real pscdown;
                     SCIP_Real pscup;

                     assert(downvars != NULL);
                     assert(downvars[j] != NULL);

                     pscdown = SCIPvarGetPseudocost(downvars[j], scip->stat, SCIPsetFeasFloor(scip->set, downvarssols[j]) - downvarssols[j]);
                     pscup = SCIPvarGetPseudocost(downvars[j], scip->stat, SCIPsetFeasCeil(scip->set, downvarssols[j]) - downvarssols[j]);
                     estimateincr = MIN(pscdown, pscup);

                     estimateprobdown += estimateincr;
                  }
               }
               SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

               /* create the multi-aggregated rounded up constraint */
               SCIP_CALL( SCIPcreateConsLinear(scip, &probingconsup, "probingconsup", SCIPvarGetMultaggrNVars(fixvars[i]), SCIPvarGetMultaggrVars(fixvars[i]),
                     SCIPvarGetMultaggrScalars(fixvars[i]),  SCIPfeasCeil(scip, fixvarssol) -  SCIPvarGetMultaggrConstant(fixvars[i]), SCIPinfinity(scip),
                     TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
               assert(probingconsup != NULL);

               /* create the up child probing node */
               SCIP_CALL( SCIPnewProbingNode(scip) );
               node = SCIPgetCurrentNode(scip);

               SCIP_CALL( SCIPaddConsNode(scip, node, probingconsup, NULL) );
               SCIP_CALL( SCIPreleaseCons(scip, &probingconsup) );

#ifdef PRINTNODECONS
               SCIPdebugMsg(scip, " created up probing node with constraint:\n");
               SCIP_CALL( SCIPprintCons(scip, probingconsup, NULL) );
               SCIPinfoMessage(scip, NULL, "\n");
#endif
               /* solve the up child probing node */
               SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &upinf) );
               solstatup = SCIPgetLPSolstat(scip);
               lperror = lperror ||  (solstatup == SCIP_LPSOLSTAT_NOTSOLVED && upinf == 0) || (solstatup == SCIP_LPSOLSTAT_ITERLIMIT) ||
                  (solstatup == SCIP_LPSOLSTAT_TIMELIMIT);
               assert(solstatup != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

               /* break the branching rule if an error occurred, problem was not solved, iteration or time limit was reached */
               if( lperror )
               {
                  SCIPdebugMsg(scip, "error solving up node probing LP: status=%d\n", solstatup);
                  SCIPstatistic(branchruledata->noupdate = TRUE);
                  break;
               }

               upobjval = SCIPgetLPObjval(scip);
               upinf = upinf || SCIPisGE(scip, upobjval, SCIPgetCutoffbound(scip));
               assert(((solstatup != SCIP_LPSOLSTAT_INFEASIBLE) && (solstatup != SCIP_LPSOLSTAT_OBJLIMIT)) || upinf);

               SCIPdebugMsg(scip, " down node objval: %g up node objval: %g\n", downobjval, upobjval);

               if( !upinf )
               {
                  /* when an optimal solution has been found  calculate up child's estimate based on pseudo costs */
                  /* estimate = lowerbound + sum(min{f_j * pscdown_j, (1-f_j) * pscup_j}) */
                  estimateprobup = SCIPnodeGetLowerbound(node);
                  SCIP_CALL( SCIPgetLPBranchCands(scip, &upvars, &upvarssols, NULL, &nupvars, NULL, NULL) );

                  for( k = 0 ; k < nupvars; k++ )
                  {
                     SCIP_Real estimateincr;
                     SCIP_Real pscdown;
                     SCIP_Real pscup;

                     assert(upvars != NULL);
                     assert(upvars[k] != NULL);

                     pscdown = SCIPvarGetPseudocost(upvars[k], scip->stat, SCIPsetFeasFloor(scip->set, upvarssols[k]) - upvarssols[k]);
                     pscup = SCIPvarGetPseudocost(upvars[k], scip->stat, SCIPsetFeasCeil(scip->set, upvarssols[k]) - upvarssols[k]);
                     estimateincr = MIN(pscdown, pscup);
                     estimateprobup += estimateincr;
                  }
               }
               SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

               /* check whether the children nodes are solved to optimality and give a valid new lower bound or not */
               if( downinf || upinf )
               {
                  /* check if the LP is a valid relaxation and we can then collect new information */
                  if( allcolsinlp )
                  {
                     /* cut off the node either when both children are infeasible or the objective limit was reached;
                      * if only one child is feasible with LP value smaller than objective limit, add the corresponding
                      * constraint to the problem and break the branching rule in order to solve the updated LP
                      */
                     if( downinf && upinf )
                     {
                        SCIPdebugMsg(scip, "node can be cut off due to strong branching on multi-aggregated variable <%s>\n",
                           SCIPvarGetName(fixvars[i]));
                        SCIPstatistic(branchruledata->nmultaggrcutoff += 1);

                        *result = SCIP_CUTOFF;
                        break;
                     }
                     else
                     {
                        assert(!lperror);

                        if( downinf )
                           downnodeinf = TRUE;

                        SCIPdebugMsg(scip, "%s child of multi-aggregated variable <%s> is infeasible\n",
                           downinf ? "down" : "up", SCIPvarGetName(fixvars[i]) );
                        SCIPstatistic(branchruledata->nmultaggrconsadd += 1);

                        *result = SCIP_CONSADDED;
                        break;
                     }
                  }
               }
               else
               {
                  /* if both children are solved to optimality and they both give a new valid bound, calculate the score of the
                   * multi-aggregated variable
                   */
                  SCIP_Real downgain;
                  SCIP_Real upgain;
                  SCIP_Real down;
                  SCIP_Real up;
                  SCIP_Real score;
                  SCIP_Real minbound;

                  assert(!downinf);
                  assert(!upinf);
                  assert(!lperror);

                  SCIPdebugMsg(scip, " both probing nodes are valid while branching on multi-aggregated variable: <%s>\n ", SCIPvarGetName(fixvars[i]));

                  down = MAX(downobjval, lpobjval);
                  up = MAX(upobjval, lpobjval);
                  downgain = down - lpobjval;
                  upgain = up - lpobjval;
                  score = SCIPgetBranchScore(scip, NULL, downgain, upgain);

                  if( allcolsinlp && !exactsolve )
                  {
                     /* the minimal lower bound of both children is a proved lower bound of the current subtree */
                     minbound = MIN(downobjval, upobjval);
                     *provedbound = MAX(*provedbound, minbound);
                  }

                  SCIPstatistic(
                     if( score > *bestmultaggrscore )
                        *bestmultaggrscore = score;
                     );

                  /* update the best branching candidate and all its values if a strictly greater score has been found */
                  if( score > *bestscore )
                  {
                     SCIPstatistic(
                        if( branchruledata->nmultaggrbranch == 0 )
                        {
                           branchruledata->rundepth = SCIPgetNRuns(scip);
                           branchruledata->firstmultaggrdepth = SCIPgetFocusDepth(scip);
                        }
                        )

                     SCIPdebugMsg(scip, " <%s> is a better candidate for branching\n", SCIPvarGetName(fixvars[i]));

                     *bestscore = MAX(score, *bestscore);
                     *bestcand = fixvars[i];
                     *bestsol = fixvarssol;
                     *bestdown = downobjval;
                     *bestup = upobjval;
                     *bestdownvalid = TRUE;
                     *bestupvalid = TRUE;
                     *estimatedown = estimateprobdown;
                     *estimateup = estimateprobup;
                  }
                  assert(bestscore != NULL);
                  assert(bestcand != NULL);
                  assert(bestup != NULL);
                  assert(bestdown != NULL);
               }
            }
         }
      }

      /* end probing mode */
      if( endprobing )
      {
         SCIP_CALL( SCIPendProbing(scip) );
      }

      SCIPdebugMsg(scip, "\n");

      /* one of the child nodes was infeasible, add the other constraint to the current node */
      if( *result == SCIP_CONSADDED )
      {
         node = SCIPgetCurrentNode(scip);
         if( downnodeinf )
         {
            SCIP_CALL( SCIPcreateConsLinear(scip, &probingconsup, "infconsup", SCIPvarGetMultaggrNVars(fixvars[i]),
                  SCIPvarGetMultaggrVars(fixvars[i]), SCIPvarGetMultaggrScalars(fixvars[i]),
                  SCIPfeasCeil(scip, fixvarssols[i]) - SCIPvarGetMultaggrConstant(fixvars[i]), SCIPinfinity(scip),
                  TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
            assert(probingconsup != NULL);
            SCIP_CALL( SCIPaddConsNode(scip, node, probingconsup, NULL) );
            SCIPdebugMsg(scip, " <%s> new valid constraint has been added to the original problem\n", SCIPconsGetName(probingconsup));
            SCIP_CALL( SCIPreleaseCons(scip, &probingconsup) );
         }
         else
         {
            SCIP_CALL( SCIPcreateConsLinear(scip, &probingconsdown, "infconsdown", SCIPvarGetMultaggrNVars(fixvars[i]),
                  SCIPvarGetMultaggrVars(fixvars[i]), SCIPvarGetMultaggrScalars(fixvars[i]), - SCIPinfinity(scip),
                  SCIPfeasFloor(scip, fixvarssols[i]) - SCIPvarGetMultaggrConstant(fixvars[i]), TRUE, TRUE, FALSE, FALSE,
                  TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
            assert(probingconsdown != NULL);
            SCIP_CALL( SCIPaddConsNode(scip, node, probingconsdown, NULL) );
            SCIPdebugMsg(scip, " <%s> new valid constraint has been added to the original problem\n", SCIPconsGetName(probingconsdown));
            SCIP_CALL( SCIPreleaseCons(scip, &probingconsdown) );
         }
      }
      SCIPfreeBufferArray(scip, &fixvarssols);
   }
   return SCIP_OKAY;
}


/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyMultAggr)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleMultAggr(scip) ) ;

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeMultAggr)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPstatistic(SCIPfreeBlockMemoryArrayNull(scip , &branchruledata->ratioggain, branchruledata->size));
   SCIPfreeBlockMemoryArrayNull(scip, &branchruledata->skipdown, branchruledata->skipsize);
   SCIPfreeBlockMemoryArrayNull(scip, &branchruledata->skipup, branchruledata->skipsize);

   SCIPfreeBlockMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitMultAggr)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   branchruledata->lastcand = 0;
   SCIPstatistic(
      branchruledata->firstmultaggrdepth = 0;
      branchruledata->nmultaggrbranch = 0;
      branchruledata->nfracbranch = 0;
      branchruledata->nEscore = 0;
      branchruledata->nmultaggrcutoff = 0;
      branchruledata->nmultaggrconsadd = 0;
      branchruledata->nfractcutoff = 0;
      branchruledata->nfractconsadd = 0;
      branchruledata->nrun = 0;
      branchruledata->size = 100;
      branchruledata->ameanratio = 0.0;
      branchruledata->noupdate = FALSE;
      branchruledata->clckstrongbr = NULL;
      branchruledata->clckmultaggrbr = NULL;
      branchruledata->nstrongbrcall = 0;
      branchruledata->nmultaggrbrcall = 0;
      branchruledata->totalmultaggrcands = 0;
      branchruledata->totallpcands = 0;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->ratioggain, branchruledata->size) );
      BMSclearMemoryArray(branchruledata->ratioggain, branchruledata->size);
      SCIP_CALL( SCIPcreateClock(scip, &branchruledata->clckstrongbr) );
      SCIP_CALL( SCIPcreateClock(scip, &branchruledata->clckmultaggrbr) );
   )
   return SCIP_OKAY;
}

/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitMultAggr)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIPstatistic(int j = 0);

   /* initialize branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert((branchruledata->skipdown != NULL) == (branchruledata->skipup != NULL));

   /* print statistics */
   SCIPstatistic(
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Multi-aggregated branching stats          : \n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  nmultaggrvars                           :  %d    (last run)\n",
         branchruledata->nmultaggrvars);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  firstmultaggrbranchdepth                :  %d    (in run %d)\n",
         branchruledata->firstmultaggrdepth,
         branchruledata->rundepth);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  nmultaggrbranch                         :  %d    (tot %d)\n",
         branchruledata->nmultaggrbranch, branchruledata->nmultaggrbranch +  branchruledata->nfracbranch);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  nmultaggrcutoff                         :  %d\n", branchruledata->nmultaggrcutoff);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  nmultaggrconsadd                        :  %d\n", branchruledata->nmultaggrconsadd);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  nfractcutoff                            :  %d\n", branchruledata->nfractcutoff);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  nfractconsadd                           :  %d\n", branchruledata->nfractconsadd);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  nEscore                                 :  %d\n", branchruledata->nEscore);

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Branching Time                          :    \n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  nstrongbrcall                           :  %d\n", branchruledata->nstrongbrcall);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  totalstrongbrtime                       :  %g\n",
         SCIPgetClockTime(scip, branchruledata->clckstrongbr));
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  totallpcands                            :  %d\n", branchruledata->totallpcands);

      if( branchruledata->totallpcands != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  averagetimestrongbr                     :  %g\n",
         SCIPgetClockTime(scip, branchruledata->clckstrongbr) / branchruledata->totallpcands);
      }
      else
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  averagetimestrongbr                     :  %s\n", "--");
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  nmultaggrbrcall                         :  %d\n", branchruledata->nmultaggrbrcall);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  totalmultaggrbrtime                     :  %g\n",
         SCIPgetClockTime(scip, branchruledata->clckmultaggrbr));
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  totalmultaggrcands                      :  %d\n", branchruledata->totalmultaggrcands);

     if( branchruledata->totalmultaggrcands != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  averagetimemultaggrbr                   :  %g\n",
         SCIPgetClockTime(scip, branchruledata->clckmultaggrbr) / branchruledata->totalmultaggrcands);
      }
      else
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  averagetimemultaggrbr                   :  %s\n", "--");
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Ratioggain                              :\n");
      if( branchruledata->nmultaggrbranch != 0 )
      {
         for( j = 0; j < branchruledata->nmultaggrbranch; j++ )
         {
            branchruledata->ameanratio += branchruledata->ratioggain[j];
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  %g", branchruledata->ratioggain[j]);
         }

         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
         branchruledata->ameanratio = branchruledata->ameanratio / branchruledata->nmultaggrbranch;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  ameanratio                              :  %4.2f\n", branchruledata->ameanratio);
      }
      else
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  ameanratio                              :  %s\n", "--");
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");

      /* free arrays */
      if( branchruledata->ratioggain != NULL )
      {
         SCIPfreeMemoryArray(scip, &branchruledata->ratioggain);
         branchruledata->ratioggain = NULL;
      }
      SCIP_CALL( SCIPfreeClock(scip, &branchruledata->clckstrongbr) );
      SCIP_CALL( SCIPfreeClock(scip, &branchruledata->clckmultaggrbr) );
   )
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

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpMultAggr)
{  /*lint --e{715}*/

   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** lpcands;
   SCIP_VAR** tmplpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   SCIP_Real* tmplpcandssol;
   SCIP_Real* tmplpcandsfrac;
   SCIP_NODE* downchild;
   SCIP_NODE* upchild;
   SCIP_Real bestup;
   SCIP_Real bestdown;
   SCIP_Real bestscore;
   SCIP_Real provedbound;
   SCIP_Real estimatedown = 0.0;
   SCIP_Real estimateup = 0.0;
   SCIP_Bool bestdownvalid;
   SCIP_Bool bestupvalid;
   SCIP_Longint oldreevalage;
   int bestcandpos;
   int nlpcands;
   int npriolpcands;
   SCIPstatistic(
      SCIP_Real lpobjval;
      SCIP_Bool reoptimize;
   )

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Execlp method of mult-aggreg branching\n ");
   *result = SCIP_DIDNOTRUN;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIP_CALL( SCIPgetLongintParam(scip, "branching/fullstrong/reevalage", &oldreevalage) );
   SCIP_CALL( SCIPsetLongintParam(scip, "branching/fullstrong/reevalage", branchruledata->reevalage) );

   /* get the lpobjval and the number of multi aggregated variables of the problem as a statistic counter */
   SCIPstatistic(
      reoptimize = FALSE;
      lpobjval = SCIPgetLPObjval(scip);

      if( SCIPgetNRuns(scip) != branchruledata->nrun )
      {
         SCIP_VAR** fixvars;
         int nfixvars;
         int i;

         branchruledata->nmultaggrvars = 0;
         fixvars = SCIPgetFixedVars(scip);
         nfixvars = SCIPgetNFixedVars(scip);

         if( nfixvars != 0 )
         {
            for( i = 0; i < nfixvars; i++ )
            {
               if( SCIPvarGetStatus(fixvars[i]) == SCIP_VARSTATUS_MULTAGGR && (SCIPvarGetType(fixvars[i]) == SCIP_VARTYPE_INTEGER ||
                     SCIPvarGetType(fixvars[i]) == SCIP_VARTYPE_BINARY) )
               {
                  branchruledata->nmultaggrvars += 1;
               }
            }
         }
         branchruledata->nrun = SCIPgetNRuns(scip);
      }
   )

   /* get all branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &tmplpcands, &tmplpcandssol, &tmplpcandsfrac, &nlpcands, &npriolpcands, NULL) );
   assert(nlpcands > 0);
   assert(npriolpcands > 0);

   /* copy LP branching candidates and solution values, because they will be updated w.r.t. the strong branching LP
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

   /* start the clock to get the time spent inside the function */
   SCIPstatistic(
      SCIP_CALL( SCIPstartClock(scip, branchruledata->clckstrongbr) );
      );

   /* compute strong branching among the array of fractional variables in order to get the best one */
   SCIP_CALL( SCIPselectVarStrongBranching(scip, lpcands, lpcandssol, lpcandsfrac, branchruledata->skipdown,
         branchruledata->skipup, nlpcands, npriolpcands, nlpcands, &branchruledata->lastcand,
         branchruledata->maxproprounds, branchruledata->probingbounds, TRUE,
         &bestcandpos, &bestdown, &bestup, &bestscore, &bestdownvalid, &bestupvalid, &provedbound, result) );

   SCIPstatistic(
      SCIP_CALL( SCIPstopClock(scip, branchruledata->clckstrongbr) );
      branchruledata->totallpcands += SCIPgetNLPBranchCands(scip);
      branchruledata->nstrongbrcall += 1;
      )

   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED )
   {
      SCIP_VAR* bestcand = lpcands[bestcandpos];
      SCIP_Real bestsol = lpcandssol[bestcandpos];
      SCIPstatistic( SCIP_Real bestmultaggrscore = -SCIPinfinity(scip); )

      SCIPstatistic(
         SCIP_Real fdowngain = 0.0;
         SCIP_Real fupgain = 0.0;

	 /* reoptimize is set to true if strong branching on fractional variables did not explicitly evaluate the objective
          * values of the probing child nodes and thus we do not have updated information
          */
         if( SCIPisLT(scip, SCIPgetVarStrongbranchLPAge(scip, bestcand), branchruledata->reevalage)
            || branchruledata->maxproprounds != 0 )
            reoptimize = TRUE;

         /* store values needed for the ratioggain statistic */
         if( !reoptimize )
         {
            SCIP_Real fdown;
            SCIP_Real fup;

            fdown = MAX(bestdown, lpobjval);
            fup = MAX(bestup, lpobjval);
            fdowngain = fdown - lpobjval;
            fupgain = fup - lpobjval;
         }

         /* start and then stop the clock to get the time spent inside the function */
         SCIP_CALL( SCIPstartClock(scip, branchruledata->clckmultaggrbr) );
      )

      /* compute strong branching among the multi-aggregated variables and the best fractional variable */
#ifdef SCIP_STATISTIC
      SCIP_CALL( selectVarMultAggrBranching(scip, &bestcand, &bestscore, &bestsol, &bestdown, &bestup, &bestdownvalid, &bestupvalid, &provedbound,
            &estimatedown, &estimateup, &bestmultaggrscore, result) );
#else
      SCIP_CALL( selectVarMultAggrBranching(scip, &bestcand, &bestscore, &bestsol, &bestdown, &bestup, &bestdownvalid, &bestupvalid, &provedbound,
            &estimatedown, &estimateup, result) );
#endif
      SCIPstatistic(
         SCIP_CALL( SCIPstopClock(scip, branchruledata->clckmultaggrbr) );
         branchruledata->nmultaggrbrcall += 1;
      )

      if( *result != SCIP_CUTOFF && *result != SCIP_CONSADDED )
      {
         SCIPstatistic(
            if( !(branchruledata->noupdate) )
            {
               if( SCIPisEQ(scip, bestmultaggrscore, bestscore) )
                  branchruledata->nEscore += 1;
            }
            )

         assert(bestcand != NULL);
         SCIPdebugMsg(scip, "BRANCHING MODE:\n");

         /* perform branching on the best found candidate */
         if( SCIPvarGetStatus(bestcand) == SCIP_VARSTATUS_MULTAGGR )
         {
            SCIP_CONS* multaggrconsdown;
            SCIP_CONS* multaggrconsup;

            SCIPstatistic(
               if( !(branchruledata->noupdate) )
               {
                  branchruledata->nmultaggrbranch += 1;

                  if( !reoptimize )
                  {
                     SCIP_Real gfractbranch;
                     SCIP_Real gmultaggrbranch;
                     SCIP_Real downgain;
                     SCIP_Real upgain;
                     SCIP_Real down;
                     SCIP_Real up;
                     int nmultaggrbranch;

                     down = MAX(bestdown, lpobjval);
                     up = MAX(bestup, lpobjval);
                     downgain = down - lpobjval;
                     upgain = up - lpobjval;

                     SCIP_CALL( ensureArraySize(scip, branchruledata) );

                     gfractbranch= SQRT(MAX(fdowngain,1e-06) * MAX(fupgain,1e-06));
                     gmultaggrbranch = SQRT(MAX(downgain,1e-06) * MAX(upgain,1e-06));

                     nmultaggrbranch = branchruledata->nmultaggrbranch;

                     if( gmultaggrbranch == 0.0 )
                     {
                        branchruledata->ratioggain[nmultaggrbranch - 1] = 1;
                     }
                     else
                     {
                        branchruledata->ratioggain[nmultaggrbranch - 1] = gfractbranch / gmultaggrbranch;
                     }
                  }
               }
               )

            /* create the multi-aggregated constraints rounded up and down */
            SCIP_CALL( SCIPcreateConsLinear(scip, &multaggrconsdown, "consdown", SCIPvarGetMultaggrNVars(bestcand),
                  SCIPvarGetMultaggrVars(bestcand), SCIPvarGetMultaggrScalars(bestcand), - SCIPinfinity(scip),
                  SCIPfeasFloor(scip, bestsol) - SCIPvarGetMultaggrConstant(bestcand),
                  TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

            SCIP_CALL( SCIPcreateConsLinear(scip, &multaggrconsup, "consup", SCIPvarGetMultaggrNVars(bestcand),
                  SCIPvarGetMultaggrVars(bestcand), SCIPvarGetMultaggrScalars(bestcand),
                  SCIPfeasCeil(scip, bestsol) -  SCIPvarGetMultaggrConstant(bestcand), SCIPinfinity(scip),
                  TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

            /* create the child nodes */
            SCIP_CALL( SCIPcreateChild(scip, &downchild, 1.0, estimatedown) );
            SCIPdebugMsg(scip, " down node: lowerbound %f estimate %f\n", SCIPnodeGetLowerbound(downchild), SCIPnodeGetEstimate(downchild));

            SCIP_CALL( SCIPcreateChild(scip, &upchild, 1.0, estimateup) );
            SCIPdebugMsg(scip, " up node: lowerbound %f estimate %f\n", SCIPnodeGetLowerbound(upchild), SCIPnodeGetEstimate(upchild));

            assert(downchild != NULL);
            assert(upchild != NULL);

            SCIP_CALL( SCIPaddConsNode(scip, downchild, multaggrconsdown, NULL) );
            SCIP_CALL( SCIPaddConsNode(scip, upchild, multaggrconsup, NULL) );

#ifdef PRINTNODECONS
            SCIPdebugMsg(scip, "branching at node %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

            SCIPdebugMsg(scip, "created child node %lld with constraint:\n", SCIPnodeGetNumber(downchild));
            SCIP_CALL( SCIPprintCons(scip, multaggrconsdown, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");

            SCIPdebugMsg(scip, "created child node %lld with constraint:\n", SCIPnodeGetNumber(upchild));
            SCIP_CALL( SCIPprintCons(scip, multaggrconsup, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif

            /* relase constraints */
            SCIP_CALL( SCIPreleaseCons(scip, &multaggrconsdown) );
            SCIP_CALL( SCIPreleaseCons(scip, &multaggrconsup) );

            SCIPdebugMsg(scip, "BRANCHED on multi-aggregated variable <%s>\n", SCIPvarGetName(bestcand));

            *result = SCIP_BRANCHED;
         }
         else
         {
            SCIPstatistic(
               if( !(branchruledata->noupdate) )
                  branchruledata->nfracbranch += 1
             );

            assert(*result == SCIP_DIDNOTRUN);
            assert(SCIPisLT(scip, provedbound, SCIPgetCutoffbound(scip)));

            SCIP_CALL( SCIPbranchVarVal(scip, bestcand, bestsol, &downchild, NULL, &upchild) );

            assert(downchild != NULL);
            assert(upchild != NULL);

            SCIPdebugMsg(scip, "BRANCHED on fractional variable <%s>\n", SCIPvarGetName(bestcand));

            *result = SCIP_BRANCHED;
         }

         /* update the lower bounds in the children; we must not do this if columns are missing in the LP
          * (e.g., because we are doing branch-and-price) or the problem should be solved exactly
          */
         if( SCIPallColsInLP(scip) && !SCIPisExactSolve(scip) )
         {
            SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, bestdownvalid ? MAX(bestdown, provedbound) : provedbound) );
            SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, bestupvalid ? MAX(bestup, provedbound) : provedbound) );
         }
         SCIPdebugMsg(scip, " -> down child's lowerbound: %g\n", SCIPnodeGetLowerbound(downchild));
         SCIPdebugMsg(scip, " -> up child's lowerbound: %g\n", SCIPnodeGetLowerbound(upchild));
      }
   }
   else
   {
      SCIPdebugMsg(scip, "strong branching breaks\n" );

      SCIPstatistic(
         if( *result == SCIP_CUTOFF )
         {
            branchruledata->nfractcutoff += 1;
         }
         else
         {
            branchruledata->nfractconsadd += 1;
         }
      )
   }

   SCIPfreeBufferArray(scip, &lpcandsfrac);
   SCIPfreeBufferArray(scip, &lpcandssol);
   SCIPfreeBufferArray(scip, &lpcands);

   SCIP_CALL( SCIPsetLongintParam(scip, "branching/fullstrong/reevalage", oldreevalage) );

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the multi-aggregated branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleMultAggr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create multaggr branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );
   branchruledata->lastcand = 0;
   branchruledata->skipsize = 0;
   branchruledata->skipup = NULL;
   branchruledata->skipdown = NULL;
   SCIPstatistic(branchruledata->ratioggain = NULL);

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyMultAggr) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeMultAggr) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitMultAggr) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitMultAggr) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpMultAggr) );

   /* multi-aggregated branching rule parameters */
   SCIP_CALL( SCIPaddLongintParam(scip,
         "branching/multaggr/reevalage",
         "number of intermediate LPs solved to trigger reevaluation of strong branching value for a variable that was already evaluated at the current node",
         &branchruledata->reevalage, TRUE, DEFAULT_REEVALAGE, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/multaggr/maxproprounds",
         "maximum number of propagation rounds to be performed during multaggr branching before solving the LP (-1: no limit, -2: parameter settings)",
         &branchruledata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -2, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/multaggr/probingbounds",
         "should valid bounds be identified in a probing-like fashion during multaggr branching (only with propagation)?",
         &branchruledata->probingbounds, TRUE, DEFAULT_PROBINGBOUNDS, NULL, NULL) );

   return SCIP_OKAY;
}
