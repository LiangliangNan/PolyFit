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

/**@file   scip_branch.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for branching rule plugins and branching
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/branch.h"
#include "scip/debug.h"
#include "scip/lp.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/var.h"
#include "scip/scip_branch.h"
#include "scip/scip_numerics.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_primal.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_var.h"
#include "scip/tree.h"


/** creates a branching rule and includes it in SCIP
 *
 *  @note method has all branching rule callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeBranchruleBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of branching rule */
   const char*           desc,               /**< description of branching rule */
   int                   priority,           /**< priority of the branching rule */
   int                   maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying branching rule
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
   SCIP_DECL_BRANCHCOPY  ((*branchcopy)),    /**< copy method of branching rule or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   SCIP_DECL_BRANCHINIT  ((*branchinit)),    /**< initialize branching rule */
   SCIP_DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialize branching rule */
   SCIP_DECL_BRANCHINITSOL((*branchinitsol)),/**< solving process initialization method of branching rule */
   SCIP_DECL_BRANCHEXITSOL((*branchexitsol)),/**< solving process deinitialization method of branching rule */
   SCIP_DECL_BRANCHEXECLP((*branchexeclp)),  /**< branching execution method for fractional LP solutions */
   SCIP_DECL_BRANCHEXECEXT((*branchexecext)),/**< branching execution method for external candidates */
   SCIP_DECL_BRANCHEXECPS((*branchexecps)),  /**< branching execution method for not completely fixed pseudo solutions */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   SCIP_BRANCHRULE* branchrule;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeBranchrule", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether branching rule is already present */
   if( SCIPfindBranchrule(scip, name) != NULL )
   {
      SCIPerrorMessage("branching rule <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbranchruleCreate(&branchrule, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, priority, maxdepth,
         maxbounddist, branchcopy, branchfree, branchinit, branchexit, branchinitsol, branchexitsol,
         branchexeclp, branchexecext, branchexecps, branchruledata) );
   SCIP_CALL( SCIPsetIncludeBranchrule(scip->set, branchrule) );

   return SCIP_OKAY;
}

/** creates a branching rule and includes it in SCIP. All non-fundamental (or optional) callbacks will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetBranchruleInit(), SCIPsetBranchruleExit(),
 *  SCIPsetBranchruleCopy(), SCIPsetBranchruleFree(), SCIPsetBranchruleInitsol(), SCIPsetBranchruleExitsol(),
 *  SCIPsetBranchruleExecLp(), SCIPsetBranchruleExecExt(), and SCIPsetBranchruleExecPs().
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeBranchrule() instead
 */
SCIP_RETCODE SCIPincludeBranchruleBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE**     branchruleptr,      /**< pointer to branching rule, or NULL */
   const char*           name,               /**< name of branching rule */
   const char*           desc,               /**< description of branching rule */
   int                   priority,           /**< priority of the branching rule */
   int                   maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying branching rule
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   SCIP_BRANCHRULE* branchrule;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeBranchruleBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether branching rule is already present */
   if( SCIPfindBranchrule(scip, name) != NULL )
   {
      SCIPerrorMessage("branching rule <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbranchruleCreate(&branchrule, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority, maxdepth,
         maxbounddist, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, NULL, NULL, branchruledata) );

   SCIP_CALL( SCIPsetIncludeBranchrule(scip->set, branchrule) );

   if( branchruleptr != NULL )
      *branchruleptr = branchrule;

   return SCIP_OKAY;
}

/** sets copy method of branching rule */
SCIP_RETCODE SCIPsetBranchruleCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHCOPY  ((*branchcopy))     /**< copy method of branching rule or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBranchruleCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(branchrule != NULL);

   SCIPbranchruleSetCopy(branchrule, branchcopy);

   return SCIP_OKAY;
}

/** sets destructor method of branching rule */
SCIP_RETCODE SCIPsetBranchruleFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHFREE  ((*branchfree))     /**< destructor of branching rule */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBranchruleFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(branchrule != NULL);

   SCIPbranchruleSetFree(branchrule, branchfree);

   return SCIP_OKAY;
}

/** sets initialization method of branching rule */
SCIP_RETCODE SCIPsetBranchruleInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHINIT  ((*branchinit))     /**< initialize branching rule */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBranchruleInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(branchrule != NULL);

   SCIPbranchruleSetInit(branchrule, branchinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of branching rule */
SCIP_RETCODE SCIPsetBranchruleExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXIT  ((*branchexit))     /**< deinitialize branching rule */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBranchruleExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(branchrule != NULL);

   SCIPbranchruleSetExit(branchrule, branchexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of branching rule */
SCIP_RETCODE SCIPsetBranchruleInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHINITSOL((*branchinitsol)) /**< solving process initialization method of branching rule */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBranchruleInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(branchrule != NULL);

   SCIPbranchruleSetInitsol(branchrule, branchinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of branching rule */
SCIP_RETCODE SCIPsetBranchruleExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXITSOL((*branchexitsol)) /**< solving process deinitialization method of branching rule */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBranchruleExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(branchrule != NULL);

   SCIPbranchruleSetExitsol(branchrule, branchexitsol);

   return SCIP_OKAY;
}

/** sets branching execution method for fractional LP solutions */
SCIP_RETCODE SCIPsetBranchruleExecLp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXECLP((*branchexeclp))   /**< branching execution method for fractional LP solutions */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBranchruleExecLp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(branchrule != NULL);

   SCIPbranchruleSetExecLp(branchrule, branchexeclp);

   return SCIP_OKAY;
}

/** sets branching execution method for external candidates  */
SCIP_RETCODE SCIPsetBranchruleExecExt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXECEXT((*branchexecext)) /**< branching execution method for external candidates */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBranchruleExecExt", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(branchrule != NULL);

   SCIPbranchruleSetExecExt(branchrule, branchexecext);

   return SCIP_OKAY;
}

/** sets branching execution method for not completely fixed pseudo solutions */
SCIP_RETCODE SCIPsetBranchruleExecPs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXECPS((*branchexecps))   /**< branching execution method for not completely fixed pseudo solutions */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBranchruleExecPs", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(branchrule != NULL);

   SCIPbranchruleSetExecPs(branchrule, branchexecps);

   return SCIP_OKAY;
}

/** returns the branching rule of the given name, or NULL if not existing */
SCIP_BRANCHRULE* SCIPfindBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of branching rule */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   SCIPsetSortBranchrules(scip->set);

   return SCIPsetFindBranchrule(scip->set, name);
}

/** returns the array of currently available branching rules */
SCIP_BRANCHRULE** SCIPgetBranchrules(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->branchrules;
}

/** returns the number of currently available branching rules */
int SCIPgetNBranchrules(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nbranchrules;
}

/** sets the priority of a branching rule */
SCIP_RETCODE SCIPsetBranchrulePriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   int                   priority            /**< new priority of the branching rule */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPbranchruleSetPriority(branchrule, scip->set, priority);

   return SCIP_OKAY;
}

/** sets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
SCIP_RETCODE SCIPsetBranchruleMaxdepth(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   int                   maxdepth            /**< new maxdepth of the branching rule */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPbranchruleSetMaxdepth(branchrule, maxdepth);

   return SCIP_OKAY;
}

/** sets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
SCIP_RETCODE SCIPsetBranchruleMaxbounddist(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_Real             maxbounddist        /**< new maxbounddist of the branching rule */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPbranchruleSetMaxbounddist(branchrule, maxbounddist);

   return SCIP_OKAY;
}

/** gets branching candidates for LP solution branching (fractional variables) along with solution values,
 *  fractionalities, and number of branching candidates; The number of branching candidates does NOT
 *  account for fractional implicit integer variables which should not be used for branching decisions.
 *
 *  Fractional implicit integer variables are stored at the positions *nlpcands to *nlpcands + *nfracimplvars - 1
 *
 *  branching rules should always select the branching candidate among the first npriolpcands of the candidate
 *  list
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetLPBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   SCIP_Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   SCIP_Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*                  nlpcands,           /**< pointer to store the number of LP branching candidates, or NULL */
   int*                  npriolpcands,       /**< pointer to store the number of candidates with maximal priority, or NULL */
   int*                  nfracimplvars       /**< pointer to store the number of fractional implicit integer variables, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_OPTIMAL && SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      SCIPerrorMessage("LP not solved to optimality - solstat=%d\n", SCIPlpGetSolstat(scip->lp));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbranchcandGetLPCands(scip->branchcand, scip->set, scip->stat, scip->lp,
         lpcands, lpcandssol, lpcandsfrac, nlpcands, npriolpcands, nfracimplvars) );

   return SCIP_OKAY;
}

/** gets number of branching candidates for LP solution branching (number of fractional variables)
 *
 *  @return the number of branching candidates for LP solution branching (number of fractional variables).
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNLPBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RETCODE retcode;
   int nlpcands;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_OPTIMAL && SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      SCIPerrorMessage("LP not solved to optimality\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }

   retcode = SCIPbranchcandGetLPCands(scip->branchcand, scip->set, scip->stat, scip->lp,
      NULL, NULL, NULL, &nlpcands, NULL, NULL);

   if( retcode != SCIP_OKAY )
   {
      SCIPerrorMessage("Error <%d> during computation of the number of LP branching candidates\n", retcode);
      SCIPABORT();
      return 0; /*lint !e527*/
   }

   return nlpcands;
}

/** gets number of branching candidates with maximal priority for LP solution branching
 *
 *  @return the number of branching candidates with maximal priority for LP solution branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNPrioLPBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RETCODE retcode;
   int npriolpcands;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrioLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_OPTIMAL && SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      SCIPerrorMessage("LP not solved to optimality\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }

   retcode = SCIPbranchcandGetLPCands(scip->branchcand, scip->set, scip->stat, scip->lp,
      NULL, NULL, NULL, NULL, &npriolpcands, NULL);

   if( retcode != SCIP_OKAY )
   {
      SCIPerrorMessage("Error <%d> during computation of the number of LP branching candidates with maximal priority\n", retcode);
      SCIPABORT();
      return 0; /*lint !e527*/
   }

   return npriolpcands;
}

/** gets external branching candidates along with solution values, scores, and number of branching candidates;
 *  these branching candidates can be used by relaxations or nonlinear constraint handlers;
 *  branching rules should always select the branching candidate among the first nprioexterncands of the candidate
 *  list
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note Candidate variables with maximal priority are ordered: binaries first, then integers, implicit integers and
 *        continuous last.
 */
SCIP_RETCODE SCIPgetExternBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           externcands,        /**< pointer to store the array of extern branching candidates, or NULL */
   SCIP_Real**           externcandssol,     /**< pointer to store the array of extern candidate solution values, or NULL */
   SCIP_Real**           externcandsscore,   /**< pointer to store the array of extern candidate scores, or NULL */
   int*                  nexterncands,       /**< pointer to store the number of extern branching candidates, or NULL */
   int*                  nprioexterncands,   /**< pointer to store the number of candidates with maximal priority, or NULL */
   int*                  nprioexternbins,    /**< pointer to store the number of binary candidates with maximal priority, or NULL */
   int*                  nprioexternints,    /**< pointer to store the number of integer candidates with maximal priority, or NULL */
   int*                  nprioexternimpls    /**< pointer to store the number of implicit integer candidates with maximal priority,
                                              *   or NULL */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetExternBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchcandGetExternCands(scip->branchcand, externcands, externcandssol, externcandsscore, nexterncands,
         nprioexterncands, nprioexternbins, nprioexternints, nprioexternimpls) );

   return SCIP_OKAY;
}

/** gets number of external branching candidates
 *
 *  @return the number of external branching candidates.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNExternBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNExternBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNExternCands(scip->branchcand);
}

/** gets number of external branching candidates with maximal branch priority
 *
 *  @return the number of external branching candidates with maximal branch priority.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNPrioExternBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrioExternBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioExternCands(scip->branchcand);
}

/** gets number of binary external branching candidates with maximal branch priority
 *
 *  @return the number of binary external branching candidates with maximal branch priority.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNPrioExternBranchBins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrioExternBranchBins", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioExternBins(scip->branchcand);
}

/** gets number of integer external branching candidates with maximal branch priority
 *
 *  @return the number of integer external branching candidates with maximal branch priority.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNPrioExternBranchInts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrioExternBranchInts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioExternInts(scip->branchcand);
}

/** gets number of implicit integer external branching candidates with maximal branch priority
 *
 *  @return the number of implicit integer external branching candidates with maximal branch priority.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNPrioExternBranchImpls(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrioExternBranchImpls", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioExternImpls(scip->branchcand);
}

/** gets number of continuous external branching candidates with maximal branch priority
 *
 *  @return the number of continuous external branching candidates with maximal branch priority.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNPrioExternBranchConts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrioExternBranchConts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioExternConts(scip->branchcand);
}

/** insert variable, its score and its solution value into the external branching candidate storage
 * the relative difference of the current lower and upper bounds of a continuous variable must be at least epsilon
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPaddExternBranchCand(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to insert */
   SCIP_Real             score,              /**< score of external candidate, e.g. infeasibility */
   SCIP_Real             solval              /**< value of the variable in the current solution */
   )
{
   assert(scip != NULL);
   assert(var->scip == scip);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddExternBranchCand", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchcandAddExternCand(scip->branchcand, scip->set, var, score, solval) );

   return SCIP_OKAY;
}

/** removes all external candidates from the storage for external branching
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
void SCIPclearExternBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPclearExternBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPbranchcandClearExternCands(scip->branchcand);
}

/** checks whether the given variable is contained in the candidate storage for external branching
 *
 *  @return whether the given variable is contained in the candidate storage for external branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPcontainsExternBranchCand(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to look for */
   )
{
   assert(scip != NULL);
   assert(var->scip == scip);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPcontainsExternBranchCand", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandContainsExternCand(scip->branchcand, var);
}

/** gets branching candidates for pseudo solution branching (non-fixed variables) along with the number of candidates
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetPseudoBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*                  npseudocands,       /**< pointer to store the number of pseudo branching candidates, or NULL */
   int*                  npriopseudocands    /**< pointer to store the number of candidates with maximal priority, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchcandGetPseudoCands(scip->branchcand, scip->set, scip->transprob,
         pseudocands, npseudocands, npriopseudocands) );

   return SCIP_OKAY;
}

/** gets number of branching candidates for pseudo solution branching (non-fixed variables)
 *
 *  @return the number branching candidates for pseudo solution branching (non-fixed variables).
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNPseudoBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPseudoCands(scip->branchcand);
}

/** gets number of branching candidates with maximal branch priority for pseudo solution branching
 *
 *  @return the number of branching candidates with maximal branch priority for pseudo solution branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNPrioPseudoBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrioPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoCands(scip->branchcand);
}

/** gets number of binary branching candidates with maximal branch priority for pseudo solution branching
 *
 *  @return the number of binary branching candidates with maximal branch priority for pseudo solution branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNPrioPseudoBranchBins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrioPseudoBranchBins", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoBins(scip->branchcand);
}

/** gets number of integer branching candidates with maximal branch priority for pseudo solution branching
 *
 *  @return the number of integer branching candidates with maximal branch priority for pseudo solution branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNPrioPseudoBranchInts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrioPseudoBranchInts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoInts(scip->branchcand);
}

/** gets number of implicit integer branching candidates with maximal branch priority for pseudo solution branching
 *
 *  @return the number of implicit integer branching candidates with maximal branch priority for pseudo solution branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNPrioPseudoBranchImpls(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrioPseudoBranchImpls", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoImpls(scip->branchcand);
}

/** calculates the branching score out of the gain predictions for a binary branching
 *
 *  @return the branching score out of the gain predictions for a binary branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetBranchScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   SCIP_Real             downgain,           /**< prediction of objective gain for rounding downwards */
   SCIP_Real             upgain              /**< prediction of objective gain for rounding upwards */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetBranchScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var == NULL || var->scip == scip );

   return SCIPbranchGetScore(scip->set, var, downgain, upgain);
}

/** calculates the branching score out of the gain predictions for a branching with arbitrary many children
 *
 *  @return the branching score out of the gain predictions for a branching with arbitrary many children.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetBranchScoreMultiple(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   int                   nchildren,          /**< number of children that the branching will create */
   SCIP_Real*            gains               /**< prediction of objective gain for each child */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetBranchScoreMultiple", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPbranchGetScoreMultiple(scip->set, var, nchildren, gains);
}

/** computes a branching point for a continuous or discrete variable
 *
 *  @see SCIPbranchGetBranchingPoint
 *
 *  @return the branching point for a continuous or discrete variable.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetBranchingPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable, of which the branching point should be computed */
   SCIP_Real             suggestion          /**< suggestion for branching point, or SCIP_INVALID if no suggestion */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetBranchingPoint", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPbranchGetBranchingPoint(scip->set, scip->tree, var, suggestion);
}

/** calculates the node selection priority for moving the given variable's LP value to the given target value;
 *  this node selection priority can be given to the SCIPcreateChild() call
 *
 *  @return the node selection priority for moving the given variable's LP value to the given target value.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPcalcNodeselPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable on which the branching is applied */
   SCIP_BRANCHDIR        branchdir,          /**< type of branching that was performed: upwards, downwards, or fixed;
                                              *   fixed should only be used, when both bounds changed
                                              */
   SCIP_Real             targetvalue         /**< new value of the variable in the child node */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPcalcNodeselPriority", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPtreeCalcNodeselPriority(scip->tree, scip->set, scip->stat, var, branchdir, targetvalue);
}

/** calculates an estimate for the objective of the best feasible solution contained in the subtree after applying the given
 *  branching; this estimate can be given to the SCIPcreateChild() call
 *
 *  @return the estimate for the objective of the best feasible solution contained in the subtree after applying the given
 *  branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPcalcChildEstimate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable on which the branching is applied */
   SCIP_Real             targetvalue         /**< new value of the variable in the child node */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPcalcChildEstimate", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPtreeCalcChildEstimate(scip->tree, scip->set, scip->stat, var, targetvalue);
}

/** calculates the increase of the estimate for the objective of the best feasible solution contained in the subtree
 *  after applying the given branching
 *
 *  @return the increase of the estimate for the objective of the best feasible solution contained in the subtree after
 *          applying the given branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPcalcChildEstimateIncrease(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable on which the branching is applied */
   SCIP_Real             varsol,             /**< solution value of variable */
   SCIP_Real             targetvalue         /**< new value of the variable in the child node */
   )
{
   SCIP_Real estimateinc;

   assert(scip != NULL);
   assert(var != NULL);

   /* compute increase above parent node's (i.e., focus node's) estimate value */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      estimateinc = SCIPvarGetPseudocost(var, scip->stat, targetvalue - varsol);
   else
   {
      SCIP_Real pscdown;
      SCIP_Real pscup;

      /* calculate estimate based on pseudo costs:
       *   estimate = lowerbound + sum(min{f_j * pscdown_j, (1-f_j) * pscup_j})
       *            = parentestimate - min{f_b * pscdown_b, (1-f_b) * pscup_b} + (targetvalue-oldvalue)*{pscdown_b or pscup_b}
       */
      pscdown = SCIPvarGetPseudocost(var, scip->stat, SCIPsetFeasFloor(scip->set, varsol) - varsol);
      pscup = SCIPvarGetPseudocost(var, scip->stat, SCIPsetFeasCeil(scip->set, varsol) - varsol);
      estimateinc = SCIPvarGetPseudocost(var, scip->stat, targetvalue - varsol) - MIN(pscdown, pscup);
   }

   /* due to rounding errors estimateinc might be slightly negative */
   if( estimateinc > 0.0 )
      estimateinc = 0.0;

   return estimateinc;
}

/** creates a child node of the focus node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcreateChild(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE**           node,               /**< pointer to node data structure */
   SCIP_Real             nodeselprio,        /**< node selection priority of new node */
   SCIP_Real             estimate            /**< estimate for(transformed) objective value of best feasible solution in subtree */
   )
{
   assert(node != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateChild", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnodeCreateChild(node, scip->mem->probmem, scip->set, scip->stat, scip->tree, nodeselprio, estimate) );

   return SCIP_OKAY;
}

/** branches on a non-continuous variable v using the current LP or pseudo solution;
 *  if solution value x' is fractional, two child nodes will be created
 *  (x <= floor(x'), x >= ceil(x')),
 *  if solution value is integral, the x' is equal to lower or upper bound of the branching
 *  variable and the bounds of v are finite, then two child nodes will be created
 *  (x <= x'', x >= x''+1 with x'' = floor((lb + ub)/2)),
 *  otherwise (up to) three child nodes will be created
 *  (x <= x'-1, x == x', x >= x'+1)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPbranchVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_NODE**           downchild,          /**< pointer to return the left child with variable rounded down, or NULL */
   SCIP_NODE**           eqchild,            /**< pointer to return the middle child with variable fixed, or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child with variable rounded up, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPbranchVar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      SCIPerrorMessage("cannot branch on continuous variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   if( SCIPsetIsEQ(scip->set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      SCIPerrorMessage("cannot branch on variable <%s> with fixed domain [%.15g,%.15g]\n",
         SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPtreeBranchVar(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
         scip->lp, scip->branchcand, scip->eventqueue, var, SCIP_INVALID, downchild, eqchild, upchild) );

   return SCIP_OKAY;
}

/** branches a variable x using a given domain hole; two child nodes (x <= left, x >= right) are created
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPbranchVarHole(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             left,               /**< left side of the domain hole */
   SCIP_Real             right,              /**< right side of the domain hole */
   SCIP_NODE**           downchild,          /**< pointer to return the left child (x <= left), or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child (x >= right), or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPbranchVarHole", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   SCIP_CALL( SCIPtreeBranchVarHole(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
         scip->origprob, scip->lp, scip->branchcand, scip->eventqueue, var, left, right, downchild, upchild) );

   return SCIP_OKAY;
}

/** branches on a variable x using a given value x';
 *  for continuous variables with relative domain width larger epsilon, x' must not be one of the bounds;
 *  two child nodes (x <= x', x >= x') are created;
 *  for integer variables, if solution value x' is fractional, two child nodes are created
 *  (x <= floor(x'), x >= ceil(x')),
 *  if x' is integral, three child nodes are created
 *  (x <= x'-1, x == x', x >= x'+1)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPbranchVarVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             val,                /**< value to branch on */
   SCIP_NODE**           downchild,          /**< pointer to return the left child with variable rounded down, or NULL */
   SCIP_NODE**           eqchild,            /**< pointer to return the middle child with variable fixed, or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child with variable rounded up, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPbranchVarVal", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   /* A continuous variable will be fixed if SCIPisRelEQ(lb,ub) is true. Otherwise, the given branching value should be
    * such that its value is not equal to one of the bounds. We assert this by requiring that it is at least eps/2 away
    * from each bound. The 2.1 is there, because ub-lb may be in (eps, 2*eps], in which case there is no way to choose a
    * branching value that is at least eps away from both bounds. However, if the variable bounds are below/above
    * -/+infinity * 2.1, then SCIPisLT will give an assert, so we omit the check in this case.
    */
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS ||
      SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) ||
      SCIPisInfinity(scip, -2.1*SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, 2.1*SCIPvarGetUbLocal(var)) ||
      (SCIPisLT(scip, 2.1*SCIPvarGetLbLocal(var), 2.1*val) && SCIPisLT(scip, 2.1*val, 2.1*SCIPvarGetUbLocal(var)) ) );

   if( SCIPsetIsEQ(scip->set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      SCIPerrorMessage("cannot branch on variable <%s> with fixed domain [%.15g,%.15g]\n",
         SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPtreeBranchVar(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
         scip->lp, scip->branchcand, scip->eventqueue, var, val, downchild, eqchild, upchild) );

   return SCIP_OKAY;
}

/** n-ary branching on a variable x using a given value
 *
 *  Branches on variable x such that up to n/2 children are created on each side of the usual branching value.
 *  The branching value is selected as in SCIPbranchVarVal().
 *  The parameters minwidth and widthfactor determine the domain width of the branching variable in the child nodes.
 *  If n is odd, one child with domain width 'width' and having the branching value in the middle is created.
 *  Otherwise, two children with domain width 'width' and being left and right of the branching value are created.
 *  Next further nodes to the left and right are created, where width is multiplied by widthfactor with increasing distance
 *  from the first nodes.
 *  The initial width is calculated such that n/2 nodes are created to the left and to the right of the branching value.
 *  If this value is below minwidth, the initial width is set to minwidth, which may result in creating less than n nodes.
 *
 *  Giving a large value for widthfactor results in creating children with small domain when close to the branching value
 *  and large domain when closer to the current variable bounds. That is, setting widthfactor to a very large value and n to 3
 *  results in a ternary branching where the branching variable is mostly fixed in the middle child.
 *  Setting widthfactor to 1.0 results in children where the branching variable always has the same domain width
 *  (except for one child if the branching value is not in the middle).
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPbranchVarValNary(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             val,                /**< value to branch on */
   int                   n,                  /**< attempted number of children to be created, must be >= 2 */
   SCIP_Real             minwidth,           /**< minimal domain width in children */
   SCIP_Real             widthfactor,        /**< multiplier for children domain width with increasing distance from val, must be >= 1.0 */
   int*                  nchildren           /**< pointer to store number of created children, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPbranchVarValNary", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   /* see comment in SCIPbranchVarVal */
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS ||
      SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) ||
      SCIPisInfinity(scip, -2.1*SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, 2.1*SCIPvarGetUbLocal(var)) ||
      (SCIPisLT(scip, 2.1*SCIPvarGetLbLocal(var), 2.1*val) && SCIPisLT(scip, 2.1*val, 2.1*SCIPvarGetUbLocal(var)) ) );

   if( SCIPsetIsEQ(scip->set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      SCIPerrorMessage("cannot branch on variable <%s> with fixed domain [%.15g,%.15g]\n",
         SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPtreeBranchVarNary(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
         scip->origprob, scip->lp, scip->branchcand, scip->eventqueue, var, val, n, minwidth, widthfactor, nchildren) );

   return SCIP_OKAY;
}

/** calls branching rules to branch on an LP solution; if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *  if the branch priority of an unfixed variable is larger than the maximal branch priority of the fractional
 *  variables, pseudo solution branching is applied on the unfixed variables with maximal branch priority
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPbranchLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPbranchLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchExecLP(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
         scip->tree, scip->reopt, scip->lp, scip->sepastore, scip->branchcand, scip->eventqueue, scip->primal->cutoffbound,
         TRUE, result) );

   return SCIP_OKAY;
}

/** calls branching rules to branch on a external candidates; if no such candidates exist, the result is SCIP_DIDNOTRUN
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPbranchExtern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPbranchExtern", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchExecExtern(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
         scip->tree, scip->reopt, scip->lp, scip->sepastore, scip->branchcand, scip->eventqueue, scip->primal->cutoffbound,
         TRUE, result) );

   return SCIP_OKAY;
}

/** calls branching rules to branch on a pseudo solution; if no unfixed variables exist, the result is SCIP_DIDNOTRUN
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPbranchPseudo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPbranchPseudo", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchExecPseudo(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
         scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->primal->cutoffbound, TRUE, result) );

   return SCIP_OKAY;
}
