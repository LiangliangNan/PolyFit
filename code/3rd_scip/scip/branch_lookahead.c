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

/**@file   branch_lookahead.c
 * @ingroup DEFPLUGINS_BRANCH
 * @ingroup BRANCHINGRULES
 * @brief  lookahead LP branching rule
 * @author Christoph Schubert
 * @author Gerald Gamrath
 *
 * The (multi-level) lookahead branching rule applies strong branching to every fractional value of the LP solution
 * at the current node of the branch-and-bound tree, as well as recursivly to every temporary child problem created by this
 * strong branching. The rule selects the candidate with the best proven dual bound.
 *
 * The branching rule was motivated by the following technical report:
 *
 * @par
 * Wasu Glankwamdee and Jeff Linderoth@n
 * Lookahead Branching for Mixed Integer Programming@n
 * Technical Report 06T-004, Department of Industrial and Systems Engineering, Lehigh University.
 *
 * For a more mathematical description and a comparison between lookahead branching and other branching rules
 * in SCIP, we refer to
 *
 * @par
 * Christoph Schubert@n
 * Multi-Level Lookahead Branching@n
 * Master Thesis, Technische Universität Berlin, 2017@n
 */

/* Supported defines:
 * PRINTNODECONS: prints the binary constraints added
 * SCIP_DEBUG: prints detailed execution information
 * SCIP_STATISTIC: prints some statistics after the branching rule is freed */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "lpi/lpi.h"
#include "scip/branch_lookahead.h"
#include "scip/cons_logicor.h"
#include "scip/pub_branch.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cons.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include <string.h>

#define BRANCHRULE_NAME            "lookahead"
#define BRANCHRULE_DESC            "full strong branching over multiple levels"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define DEFAULT_USEBINARYCONSTRAINTS   FALSE /**< should binary constraints be collected and applied? */
#define DEFAULT_ADDCLIQUE              FALSE /**< add binary constraints with two variables found at the root node also as a clique? */
#define DEFAULT_ADDBINCONSROW          0     /**< should binary constraints be added as rows to the base LP?
                                              *   (0: no, 1: separate, 2: as initial rows) */
#define DEFAULT_USEDOMAINREDUCTION     TRUE  /**< Should domain reductions be collected and applied? */
#define DEFAULT_MERGEDOMAINREDUCTIONS  FALSE /**< should domain reductions of feasible siblings should be merged? */
#define DEFAULT_PREFERSIMPLEBOUNDS     FALSE /**< should domain reductions only be applied if there are simple bound changes? */
#define DEFAULT_ONLYVIOLDOMREDS        FALSE /**< Should only domain reductions that violate the LP solution be applied? */
#define DEFAULT_MAXNVIOLATEDCONS       1     /**< How many constraints that are violated by the base lp solution
                                              *   should be gathered until the rule is stopped and they are added? */
#define DEFAULT_MAXNVIOLATEDBINCONS    0     /**< How many binary constraints that are violated by the base lp
                                              *   solution should be gathered until the rule is stopped and they are
                                              *   added? */
#define DEFAULT_MAXNVIOLATEDDOMREDS    1     /**< How many domain reductions that are violated by the base lp solution
                                              *   should be gathered until the rule is stopped and they are added? */
#define DEFAULT_STOREUNVIOLATEDSOL     TRUE  /**< If only non violating constraints are added, should the branching
                                              *   decision be stored till the next call? */
#define DEFAULT_REEVALAGE              10LL  /**< Max number of LPs solved after which a previous prob branching
                                              *   result is recalculated. */
#define DEFAULT_REEVALAGEFSB           10LL  /**< Max number of LPs solved after which a previous FSB scoring
                                              *   result is recalculated. */
#define DEFAULT_RECURSIONDEPTH         2     /**< The max depth of LAB. */
#define DEFAULT_ADDNONVIOCONS          FALSE /**< Should binary constraints, that are not violated by the base LP, be
                                              *   collected and added? */
#define DEFAULT_PROPAGATE              TRUE  /**< Should domain propagation be executed before each temporary node is
                                              *   solved? */
#define DEFAULT_USELEVEL2DATA          TRUE  /**< should branching data generated at depth level 2 be stored for re-using it? */
#define DEFAULT_APPLYCHILDBOUNDS       FALSE /**< should bounds known for child nodes be applied? */
#define DEFAULT_ENFORCEMAXDOMREDS      FALSE /**< should the maximum number of domain reductions maxnviolateddomreds be enforced? */
#define DEFAULT_UPDATEBRANCHINGRESULTS FALSE /**< should branching results (and scores) be updated w.r.t. proven dual bounds? */
#define DEFAULT_MAXPROPROUNDS          0     /**< maximum number of propagation rounds to perform at temporary
                                              *   nodes (-1: unlimited, 0: SCIP default) */
#define DEFAULT_ABBREVIATED            TRUE  /**< Toggles the abbreviated LAB. */
#define DEFAULT_MAXNCANDS              4     /**< If abbreviated: The max number of candidates to consider at the base node */
#define DEFAULT_MAXNDEEPERCANDS        2     /**< If abbreviated: The max number of candidates to consider per deeper node
                                              *   (0: same as base node) */
#define DEFAULT_REUSEBASIS             TRUE  /**< If abbreviated: Should the information gathered to obtain the best
                                              *   candidates be reused? */
#define DEFAULT_ABBREVPSEUDO           FALSE /**< If abbreviated: Use pseudo costs to estimate the score of a
                                              *   candidate. */
#define DEFAULT_LEVEL2AVGSCORE         FALSE /**< should the average score be used for uninitialized scores in level 2? */
#define DEFAULT_LEVEL2ZEROSCORE        FALSE /**< should uninitialized scores be set to 0? */
#define DEFAULT_SCORINGFUNCTION        'a'   /**< scoring function to be used at the base level */
#define DEFAULT_DEEPERSCORINGFUNCTION  'x'   /**< scoring function to be used at deeper levels */
#define DEFAULT_SCORINGSCORINGFUNCTION 'd'   /**< scoring function to be used for FSB scoring */
#define DEFAULT_MINWEIGHT              0.8   /**< default value for the weight of the minimum in the convex combination of two
                                              *   child gains (taken from the paper) */
#define DEFAULT_WORSEFACTOR           -1.0   /**< if the FSB score is of a candidate is worse than the best by this factor, skip this candidate (-1: disable) */
#define DEFAULT_FILTERBYMAXGAIN       FALSE  /**< should lookahead branching only be applied if the max gain in level 1 is not uniquely that of the best candidate? */

#ifdef SCIP_DEBUG
/* Adjusted debug message that also prints the current probing depth. */
#define LABdebugMessage(scip,lvl,...)        do                                                                            \
                                             {                                                                             \
                                                SCIP_STAGE stage;                                                          \
                                                SCIPverbMessage(scip, lvl, NULL, "[%s:%-4d] ", __FILE__, __LINE__);        \
                                                stage = SCIPgetStage(scip);                                                \
                                                if( stage == SCIP_STAGE_INIT )                                             \
                                                {                                                                          \
                                                   SCIPverbMessage(scip, lvl, NULL, "Init   : ");                          \
                                                }                                                                          \
                                                else if( stage == SCIP_STAGE_FREE )                                        \
                                                {                                                                          \
                                                   SCIPverbMessage(scip, lvl, NULL, "Free   : ");                          \
                                                }                                                                          \
                                                else if( SCIPinProbing(scip) )                                             \
                                                {                                                                          \
                                                   SCIPverbMessage(scip, lvl, NULL, "%*sDepth %i: ",                       \
                                                      2 * SCIPgetProbingDepth(scip), "", SCIPgetProbingDepth(scip));       \
                                                }                                                                          \
                                                else                                                                       \
                                                {                                                                          \
                                                   SCIPverbMessage(scip, lvl, NULL, "Base   : ");                          \
                                                }                                                                          \
                                                SCIPverbMessage(scip, lvl, NULL, __VA_ARGS__);                             \
                                             }                                                                             \
                                             while( FALSE )

/* Writes a debug message without the leading information. Can be used to append something to an output of LABdebugMessage*/
#define LABdebugMessagePrint(scip,lvl,...)   do                                                                            \
                                             {                                                                             \
                                                SCIPverbMessage(scip, lvl, NULL, __VA_ARGS__);                             \
                                             }                                                                             \
                                             while( FALSE )
#else
#define LABdebugMessage(scip,lvl,...)        /**/
/*#define LABdebugMessagePrint(scip,lvl,...)   only used with SCIP_DEBUG defined */
#endif

/*
 * Data structures
 */

/** A struct holding information to speed up the solving time for solving a problem again. This is filled by the FSB
 *  scoring routine that is run to get the best candidates. It is then read by the actual ALAB routine. */
typedef struct
{
   SCIP_LPISTATE*        lpistate;           /**< the basis information that may be set before another solve lp call */
   SCIP_LPINORMS*        lpinorms;           /**< the norms that may be set before another solve lp call */
   SCIP_Bool             primalfeas;         /**< indicates whether the solution was primal feasible */
   SCIP_Bool             dualfeas;           /**< indicates whether the solution was dual feasible */
} WARMSTARTINFO;

/** Allocates the warm start information on the buffer and initializes it with default values. */
static
SCIP_RETCODE warmStartInfoCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   WARMSTARTINFO**       warmstartinfo       /**< the warmstartinfo to allocate and initialize */
   )
{
   assert(scip != NULL);
   assert(warmstartinfo != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, warmstartinfo) );

   (*warmstartinfo)->lpistate = NULL;
   (*warmstartinfo)->lpinorms = NULL;
   (*warmstartinfo)->primalfeas = FALSE;
   (*warmstartinfo)->dualfeas = FALSE;

   return SCIP_OKAY;
}

/** checks that the warm start info can be read into the lp solver. */
static
SCIP_Bool warmStartInfoIsAvailable(
   WARMSTARTINFO*        warmstartinfo       /**< the warm start info to check (may be NULL) */
   )
{
   return warmstartinfo != NULL && warmstartinfo->lpistate != NULL;
}

/** Frees the allocated buffer memory of the warm start info. */
static
SCIP_RETCODE warmStartInfoFree(
   SCIP*                 scip,               /**< SCIP data structure */
   WARMSTARTINFO**       warmstartinfo       /**< the warm start info to free */
   )
{
   SCIP_LPI* lpi;
   BMS_BLKMEM* blkmem;

   assert(scip != NULL);
   assert(warmstartinfo != NULL);

   SCIP_CALL( SCIPgetLPI(scip, &lpi) );
   blkmem = SCIPblkmem(scip);

   if( (*warmstartinfo)->lpistate != NULL )
   {
      SCIP_CALL( SCIPlpiFreeState(lpi, blkmem, &(*warmstartinfo)->lpistate) );
   }

   if( (*warmstartinfo)->lpinorms != NULL )
   {
      SCIP_CALL( SCIPlpiFreeNorms(lpi, blkmem, &(*warmstartinfo)->lpinorms) );
   }

   SCIPfreeBlockMemory(scip, warmstartinfo);

   return SCIP_OKAY;
}

/** A struct containing all information needed to branch on a variable. */
typedef struct
{
   SCIP_VAR*             branchvar;          /**< the variable to branch on */
   SCIP_Real             branchval;          /**< the fractional value to branch on */
   SCIP_Real             fracval;            /**< the fractional part of the value to branch on (val - floor(val)) */
   WARMSTARTINFO*        downwarmstartinfo;  /**< the warm start info containing the lp data from a previous down branch */
   WARMSTARTINFO*        upwarmstartinfo;    /**< the warm start info containing the lp data from a previous up branch */
} CANDIDATE;

/** Allocates the candidate on the buffer and initializes it with default values. */
static
SCIP_RETCODE candidateCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATE**           candidate           /**< the candidate to allocate and initialize */
   )
{
   assert(scip != NULL);
   assert(candidate != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, candidate) );

   (*candidate)->downwarmstartinfo = NULL;
   (*candidate)->upwarmstartinfo = NULL;
   (*candidate)->branchvar = NULL;

   return SCIP_OKAY;
}

/** free the warm starting information for the given candidate */
static
SCIP_RETCODE candidateFreeWarmStartInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATE*            candidate           /**< the candidate to free the warm starting information for */
   )
{
   assert(scip != NULL);
   assert(candidate != NULL);

   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "freeing warmstart info of candidate <%s>(%u/%u)...\n",
      SCIPvarGetName(candidate->branchvar),
      candidate->upwarmstartinfo != NULL, candidate->downwarmstartinfo != NULL);

   if( candidate->upwarmstartinfo != NULL )
   {
      SCIP_CALL( warmStartInfoFree(scip, &candidate->upwarmstartinfo) );
   }
   if( candidate->downwarmstartinfo != NULL )
   {
      SCIP_CALL( warmStartInfoFree(scip, &candidate->downwarmstartinfo) );
   }

   return SCIP_OKAY;
}


/** Frees the allocated buffer memory of the candidate and clears the contained lpi memories. */
static
SCIP_RETCODE candidateFree(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATE**           candidate           /**< the candidate to free */
   )
{
   assert(scip != NULL);
   assert(candidate != NULL);

   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "freeing candidate <%s>(%u/%u)...\n",
      (*candidate) != NULL ? SCIPvarGetName((*candidate)->branchvar) : "none",
      (*candidate)->upwarmstartinfo != NULL, (*candidate)->downwarmstartinfo != NULL);

   /* if a candidate is freed, we no longer need the content of the warm start info */
   SCIP_CALL( candidateFreeWarmStartInfo(scip, *candidate) );

   SCIPfreeBlockMemory(scip, candidate);
   return SCIP_OKAY;
}

/** Store the current lp solution in the warm start info for further usage. */
static
SCIP_RETCODE candidateStoreWarmStartInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATE*            candidate,          /**< the branching candidate */
   SCIP_Bool             down                /**< is the info for down branching? */
   )
{
   SCIP_LPI* lpi;
   BMS_BLKMEM* blkmem;
   WARMSTARTINFO* warmstartinfo;

   assert(scip != NULL);
   assert(candidate != NULL);

   SCIP_CALL( SCIPgetLPI(scip, &lpi) );
   blkmem = SCIPblkmem(scip);

   if( down )
   {
      if( candidate->downwarmstartinfo == NULL )
      {
         SCIP_CALL( warmStartInfoCreate(scip, &candidate->downwarmstartinfo) );
      }
      warmstartinfo = candidate->downwarmstartinfo;
   }
   else
   {
      if( candidate->upwarmstartinfo == NULL )
      {
         SCIP_CALL( warmStartInfoCreate(scip, &candidate->upwarmstartinfo) );
      }
      warmstartinfo = candidate->upwarmstartinfo;
   }

   SCIP_CALL( SCIPlpiGetState(lpi, blkmem, &warmstartinfo->lpistate) );

   SCIP_CALL( SCIPlpiGetNorms(lpi, blkmem, &warmstartinfo->lpinorms) );

   warmstartinfo->primalfeas = SCIPlpiIsPrimalFeasible(lpi);
   warmstartinfo->dualfeas = SCIPlpiIsDualFeasible(lpi);

   assert(warmstartinfo->lpistate != NULL);
   /* warmstartinfo->lpinorms may be NULL */

   return SCIP_OKAY;
}

/** returns whether the candidate has stored warm starting information for the given direction */
static
SCIP_Bool candidateHasWarmStartInfo(
   CANDIDATE*            candidate,          /**< the branching candidate */
   SCIP_Bool             down                /**< is the info for down branching? */
   )
{
   assert(candidate != NULL);

   return warmStartInfoIsAvailable(down ? candidate->downwarmstartinfo : candidate->upwarmstartinfo);
}


/** loads the warm starting information of the candidate for the given direction */
static
SCIP_RETCODE candidateLoadWarmStartInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATE*            candidate,          /**< the branching candidate */
   SCIP_Bool             down                /**< is the info for down branching? */
   )
{
   WARMSTARTINFO* warmstartinfo;

   assert(scip != NULL);
   assert(candidate != NULL);

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "loading basis...\n");

   if( down )
      warmstartinfo = candidate->downwarmstartinfo;
   else
      warmstartinfo = candidate->upwarmstartinfo;

   /* As we solved the very same LP some time earlier and stored the state (the basis) and the norms, we can now set those in
    * the LP solver, such that the solution does not (in best case) need any further calculation.
    * Some iterations may occur, as the conflict analysis may have added some constraints in the meantime. */
   SCIP_CALL( SCIPsetProbingLPState(scip, &(warmstartinfo->lpistate), &(warmstartinfo->lpinorms), warmstartinfo->primalfeas,
         warmstartinfo->dualfeas) );

   /* The state and norms will be freed later by the SCIP framework. Therefore they are set to NULL to enforce that we won't
    * free them on our own. */
   assert(warmstartinfo->lpistate == NULL);
   assert(warmstartinfo->lpinorms == NULL);

   return SCIP_OKAY;
}


/** Holds the information needed for branching on a variable. */
typedef struct
{
   SCIP_VAR*             branchvar;          /**< the variable to branch on, may be NULL */
   SCIP_Real             branchval;          /**< the fractional value to branch on */
   SCIP_Real*            downlowerbounds;    /**< variable lower bounds for down child */
   SCIP_Real*            downupperbounds;    /**< variable upper bounds for down child */
   SCIP_Real*            uplowerbounds;      /**< variable lower bounds for up child */
   SCIP_Real*            upupperbounds;      /**< variable upper bounds for up child */
   SCIP_Real             downdb;             /**< dual bound for down branch */
   SCIP_Real             updb;               /**< dual bound for the up branch */
   SCIP_Real             proveddb;           /**< proven dual bound for the current node */
   SCIP_Real             score;              /**< score of the branching decision */
   SCIP_Bool             downdbvalid;        /**< Indicator for the validity of the downdb value. Is FALSE, if no actual
                                              *   branching occurred or the value was determined by an LP not solved to
                                              *   optimality. */
   SCIP_Bool             updbvalid;          /**< Indicator for the validity of the updb value. Is FALSE, if no actual
                                              *   branching occurred or the value was determined by an LP not solved to
                                              *   optimality. */
   SCIP_Bool             boundsvalid;        /**< are variable bounds for down and up child valid? */
   int                   boundssize;         /**< size of bounds arrays */
} BRANCHINGDECISION;

/** initialize a branching decsion with default values */
static
void branchingDecisionInit(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGDECISION*    decision            /**< the decision to initialize */
   )
{
   assert(scip != NULL);
   assert(decision != NULL);

   decision->branchvar = NULL;
   decision->branchval = SCIP_INVALID;
   decision->downlowerbounds = NULL;
   decision->downupperbounds = NULL;
   decision->uplowerbounds = NULL;
   decision->upupperbounds = NULL;
   decision->downdb = -SCIPinfinity(scip);
   decision->downdbvalid = FALSE;
   decision->updb = -SCIPinfinity(scip);
   decision->updbvalid = FALSE;
   decision->boundsvalid = FALSE;
   decision->proveddb = -SCIPinfinity(scip);
   decision->score = -SCIPinfinity(scip);
   decision->boundssize = 0;
}


/** allocates a branching decision in the buffer and initializes it with default values. */
static
SCIP_RETCODE branchingDecisionCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGDECISION**   decision            /**< pointer to the decision to allocate and initialize */
   )
{
   assert(scip != NULL);
   assert(decision != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, decision) );

   branchingDecisionInit(scip, *decision);

   return SCIP_OKAY;
}

/** copies the data from the source branching decision storage to the target storage;
 *  this is used to store the most important information (i.e., the dual bounds obtained) so that it can be used in a
 *  subsequent call in case the LP solution did not change because we only added bound changes that did not forbid the
 *  current LP solution;
 *  however, we do not want to store all the domain changes for the two potential child nodes for this rare case, they
 *  will be identified when processing the child nodes anyway
 */
static
void branchingDecisionCopy(
   BRANCHINGDECISION*    sourcedecision,     /**< the source branching decision */
   BRANCHINGDECISION*    targetdecision      /**< the target branching decision */
   )
{
   assert(sourcedecision != NULL);
   assert(targetdecision != NULL);

   targetdecision->branchvar = sourcedecision->branchvar;
   targetdecision->branchval = sourcedecision->branchval;
   targetdecision->downdb = sourcedecision->downdb;
   targetdecision->downdbvalid = sourcedecision->downdbvalid;
   targetdecision->updb = sourcedecision->updb;
   targetdecision->updbvalid = sourcedecision->updbvalid;
   targetdecision->proveddb = sourcedecision->proveddb;
   targetdecision->score = sourcedecision->score;

   assert(targetdecision->downlowerbounds == NULL);
   assert(targetdecision->downupperbounds == NULL);
   assert(targetdecision->uplowerbounds == NULL);
   assert(targetdecision->upupperbounds == NULL);
   assert(targetdecision->boundsvalid == FALSE);
   assert(targetdecision->boundssize == 0);
}

/** Checks whether the given branching decision can be used to branch on. */
static
SCIP_Bool branchingDecisionIsValid(
   BRANCHINGDECISION*    decision            /**< the branching decision to check */
   )
{
   assert(decision != NULL);

   /* a branching decision is deemed valid, if the var pointer is not on the default NULL value (see the allocate method) */
   return decision->branchvar != NULL;
}

/* ensure that the array that stores the bounds for both child nodes is large enough */
static
SCIP_RETCODE branchingDecisionEnsureBoundArraysSize(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGDECISION*    decision,           /**< branching decision */
   int                   nvars               /**< number of problem variables */
   )
{
   assert(decision != NULL);

   if( decision->boundssize == 0 )
   {
      decision->boundssize = nvars;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &decision->downlowerbounds, decision->boundssize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &decision->downupperbounds, decision->boundssize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &decision->uplowerbounds, decision->boundssize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &decision->upupperbounds, decision->boundssize) );
   }
   assert(decision->boundssize == nvars);

   return SCIP_OKAY;
}

/** Frees the allocated memory of the branching decision. */
static
void branchingDecisionFree(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGDECISION**   decision            /**< pointer to the decision to be freed */
   )
{
   assert(scip != NULL);
   assert(decision != NULL);

   if( (*decision)->boundssize != 0 )
   {
      assert((*decision)->downlowerbounds != NULL);
      assert((*decision)->downupperbounds != NULL);
      assert((*decision)->uplowerbounds != NULL);
      assert((*decision)->upupperbounds != NULL);

      SCIPfreeBlockMemoryArray(scip, &(*decision)->downlowerbounds, (*decision)->boundssize);
      SCIPfreeBlockMemoryArray(scip, &(*decision)->downupperbounds, (*decision)->boundssize);
      SCIPfreeBlockMemoryArray(scip, &(*decision)->uplowerbounds, (*decision)->boundssize);
      SCIPfreeBlockMemoryArray(scip, &(*decision)->upupperbounds, (*decision)->boundssize);
   }

   SCIPfreeBuffer(scip, decision);
}

/** A container to hold the result of a branching. */
typedef struct
{
   SCIP_Real             objval;             /**< The objective value of the solved lp. Only contains meaningful data, if
                                              *   cutoff == FALSE. */
   SCIP_Real             dualbound;          /**< The best dual bound for this branching, may be changed by deeper level
                                              *   branchings. */
   SCIP_Longint          niterations;        /**< The number of probing iterations needed in sub branch. */
   SCIP_Bool             cutoff;             /**< Indicates whether the node was infeasible and was cutoff. */
   SCIP_Bool             dualboundvalid;     /**< Is the value of the dual bound valid? That means, was the according LP
                                              *   or the sub problems solved to optimality? */
   int                   ndeepestcutoffs;    /**< number of cutoffs on the lowest level below this child */
   SCIP_Real             deeperscore;        /**< best score computed for the deeper lookahead level */
   SCIP_Real             bestgain;           /**< best gain (w.r.t. to the base lp) on the lowest level below this child */
   SCIP_Real             totalgains;         /**< sum over all gains that are valid in both children */
   int                   ntotalgains;        /**< number of gains summed in totalgains */
   int                   ndeepestnodes;      /**< number of nodes processed in the deepest level */
} BRANCHINGRESULTDATA;

/** Allocates a branching result in the buffer. */
static
SCIP_RETCODE branchingResultDataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGRESULTDATA** resultdata          /**< pointer to the result to be allocated */
   )
{
   assert(scip != NULL);
   assert(resultdata != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, resultdata) );

   return SCIP_OKAY;
}

/** Initiates the branching result with default values. */
static
void branchingResultDataInit(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGRESULTDATA*  resultdata          /**< pointer to the result to be initialized */
   )
{
   assert(scip != NULL);
   assert(resultdata != NULL);

   resultdata->objval = -SCIPinfinity(scip);
   resultdata->dualbound = -SCIPinfinity(scip);
   resultdata->cutoff = FALSE;
   resultdata->dualboundvalid = FALSE;
   resultdata->niterations = 0;
   resultdata->ndeepestcutoffs = 0;
   resultdata->deeperscore = -SCIPinfinity(scip);
   resultdata->bestgain = 0.;
   resultdata->totalgains = 0.;
   resultdata->ntotalgains = 0;
   resultdata->ndeepestnodes = 0;
}

/** Copies the data from the source to the target. */
static
void branchingResultDataCopy(
   BRANCHINGRESULTDATA*  sourcedata,         /**< the source branching result */
   BRANCHINGRESULTDATA*  targetdata          /**< the target branching result */
   )
{
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   targetdata->cutoff = sourcedata->cutoff;
   targetdata->objval = sourcedata->objval;
   targetdata->dualbound = sourcedata->dualbound;
   targetdata->dualboundvalid = sourcedata->dualboundvalid;
   targetdata->niterations = sourcedata->niterations;
   targetdata->ndeepestcutoffs = sourcedata->ndeepestcutoffs;
   targetdata->deeperscore = sourcedata->deeperscore;
   targetdata->bestgain = sourcedata->bestgain;
   targetdata->totalgains = sourcedata->totalgains;
   targetdata->ntotalgains = sourcedata->ntotalgains;
   targetdata->ndeepestnodes = sourcedata->ndeepestnodes;
}

/** Frees the allocated buffer memory of the branching result. */
static
void branchingResultDataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGRESULTDATA** resultdata          /**< pointer to the result to be freed */
   )
{
   assert(scip != NULL);
   assert(resultdata != NULL);

   SCIPfreeBuffer(scip, resultdata);
}

/** a container to hold the result of a second-level LP */
typedef struct
{
   SCIP_Real             lpobjval;           /**< the objective value of the solved lp; only contains meaningful data, if
                                              *   cutoff == FALSE. */
   SCIP_Real             branchval1;         /**< new bound for first branching variable */
   SCIP_Real             branchval2;         /**< new bound for second branching variable */
   unsigned int          branchvar1:30;      /**< problem index of first branching variable */
   unsigned int          branchvar2:30;      /**< problem index of second branching variable */
   unsigned int          branchdir1:1;       /**< branching direction for first branching variable (0:down, 1:up) */
   unsigned int          branchdir2:1;       /**< branching direction for second branching variable (0:down, 1:up) */
   unsigned int          cutoff:1;           /**< indicates whether the node was infeasible and was cut off. */
   unsigned int          valid:1;            /**< is the lpobjval a valid dual bound? */
} LEVEL2RESULT;

/** a container to hold the results of all second-level LPs */
typedef struct
{
   LEVEL2RESULT**        level2results;      /**< array with all level2 results */
   SCIP_Real             branchval1;         /**< new bound for first branching variable */
   SCIP_Real             branchval2;         /**< new bound for second branching variable */
   int                   nlevel2results;     /**< number of level2 results stored */
   int                   level2resultssize;  /**< size of level2results array */
   unsigned int          branchvar1:30;      /**< problem index of first branching variable */
   unsigned int          branchvar2:30;      /**< problem index of second branching variable */
   unsigned int          branchdir1:1;       /**< branching direction for first branching variable (0:down, 1:up) */
   unsigned int          branchdir2:1;       /**< branching direction for second branching variable (0:down, 1:up) */
} LEVEL2DATA;

/** allocates a double branching result in the memory and fills it with the information stored in the level 2 data */
static
SCIP_RETCODE level2resultCreateFromData(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVEL2DATA*           data,               /**< level2 data */
   LEVEL2RESULT**        result              /**< pointer to the result to be allocated */
   )
{
   assert(scip != NULL);
   assert(data != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, result) );

   if( data->branchvar1 < data->branchvar2 )
   {
      (*result)->branchval1 = data->branchval1;
      (*result)->branchval2 = data->branchval2;
      (*result)->branchvar1 = data->branchvar1; /*lint !e732*/
      (*result)->branchvar2 = data->branchvar2; /*lint !e732*/
      (*result)->branchdir1 = data->branchdir1;
      (*result)->branchdir2 = data->branchdir2;
   }
   else
   {
      (*result)->branchval1 = data->branchval2;
      (*result)->branchval2 = data->branchval1;
      (*result)->branchvar1 = data->branchvar2; /*lint !e732*/
      (*result)->branchvar2 = data->branchvar1; /*lint !e732*/
      (*result)->branchdir1 = data->branchdir2;
      (*result)->branchdir2 = data->branchdir1;
   }

   return SCIP_OKAY;
}


#ifdef SCIP_DEBUG
/** prints the double branching result */
static
void level2resultPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVEL2RESULT*         result              /**< pointer to the result to be initialized */
   )
{
   SCIP_VAR** vars;

   assert(result != NULL);

   vars = SCIPgetVars(scip);

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH,
      "level 2 result: <%s> %s %g + <%s> %s %g: lpval: %.9g, inf: %d, valid: %d\n",
      SCIPvarGetName(vars[result->branchvar1]), result->branchdir1 ? ">=" : "<=", result->branchval1,
      SCIPvarGetName(vars[result->branchvar2]), result->branchdir2 ? ">=" : "<=", result->branchval2,
      result->lpobjval, result->cutoff, result->valid);
}
#else
#define level2resultPrint(scip,result) /**/
#endif

/** frees the allocated memory of the double branching result */
static
void level2resultFree(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVEL2RESULT**        result              /**< pointer to the result to be freed */
   )
{
   assert(scip != NULL);
   assert(result != NULL);

   SCIPfreeBlockMemory(scip, result);
}

/** returns TRUE iff both level 2 results are equal; two branchings are equal if they branched on the same variables
 *  with the same values
 */
static
SCIP_Bool level2resultEqual(
   LEVEL2RESULT*         result1,            /**< first level 2 result */
   LEVEL2RESULT*         result2             /**< second level 2 result */
   )
{
   assert(result1->branchvar1 < result1->branchvar2);
   assert(result2->branchvar1 < result2->branchvar2);

   /* check all cases */
   if( result1->branchvar1 != result2->branchvar1
      || result1->branchvar2 != result2->branchvar2
      || result1->branchdir1 != result2->branchdir1
      || result1->branchdir2 != result2->branchdir2
      || result1->branchval1 > result2->branchval1 + 0.5
      || result1->branchval1 < result2->branchval1 - 0.5
      || result1->branchval2 > result2->branchval2 + 0.5
      || result1->branchval2 < result2->branchval2 - 0.5)
      return FALSE;

   return TRUE;
}

/** allocates the level2 data */
static
SCIP_RETCODE level2dataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVEL2DATA**          data                /**< pointer to the data to be allocated */
   )
{
   assert(scip != NULL);
   assert(data != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, data) );

   (*data)->level2results = NULL;
   (*data)->branchval1 = -SCIPinfinity(scip);
   (*data)->branchval2 = -SCIPinfinity(scip);
   (*data)->nlevel2results = 0;
   (*data)->level2resultssize = 0;
   (*data)->branchvar1 = 0;
   (*data)->branchvar2 = 0;
   (*data)->branchdir1 = 0;
   (*data)->branchdir2 = 0;

   return SCIP_OKAY;
}

/** frees the allocated memory of the level2 data */
static
void level2dataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVEL2DATA**          data                /**< pointer to the data to be freed */
   )
{
   assert(scip != NULL);
   assert(data != NULL);

   while( (*data)->nlevel2results > 0 )
   {
      --(*data)->nlevel2results;
      level2resultFree(scip, &((*data)->level2results[(*data)->nlevel2results]));
   }
   assert((*data)->nlevel2results == 0);

   if( (*data)->level2results != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*data)->level2results, (*data)->level2resultssize);
   }

   SCIPfreeBlockMemory(scip, data);
}

/** ensures that level2results can store at least one more element */
static
SCIP_RETCODE level2dataEnsureSize(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVEL2DATA*           data                /**< level2 data */
   )
{
   assert(scip != NULL);
   assert(data != NULL);

   if( data->nlevel2results >= data->level2resultssize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, data->level2resultssize + 1);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &data->level2results, data->level2resultssize, newsize) );
      data->level2resultssize = newsize;
   }

   return SCIP_OKAY;
}

/** get a result from the level 2 data */
static
SCIP_RETCODE level2dataGetResult(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVEL2DATA*           data,               /**< level2 data */
   LEVEL2RESULT**        result              /**< pointer to store result */
   )
{
   LEVEL2RESULT* tmpresult;
   int i;

   assert(data != NULL);
   assert(result != NULL);

   *result = NULL;

   /* we branched twice on the same variable; the result cannot be stored already */
   if( data->branchvar1 == data->branchvar2 )
   {
      assert(SCIPvarGetType(SCIPgetVars(scip)[data->branchvar1]) != SCIP_VARTYPE_BINARY);
      return SCIP_OKAY;
   }

   SCIP_CALL( level2resultCreateFromData(scip, data, &tmpresult) );

   /* search for a level 2 result with the same branching decisions */
   for( i = 0; i < data->nlevel2results; ++i )
   {
      if( level2resultEqual(data->level2results[i], tmpresult) )
      {
         *result = data->level2results[i];
      }
   }

   level2resultFree(scip, &tmpresult);

   return SCIP_OKAY;
}


/** store a new result in the level 2 data */
static
SCIP_RETCODE level2dataStoreResult(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVEL2DATA*           data,               /**< level2 data */
   SCIP_Real             lpobjval,           /**< LP objective value */
   SCIP_Bool             cutoff,             /**< was the LP infeasible? */
   SCIP_Bool             valid,              /**< is the LP value a valid dual bound? */
   SCIP_Bool*            duplicate           /**< pointer to store whether information for the same branching decisions was already stored */
   )
{
   LEVEL2RESULT* result;
   int i;

   assert(scip != NULL);
   assert(data != NULL);
   assert(duplicate != NULL);

   *duplicate = FALSE;

   /* we branched twice on the same variable; the result cannot be re-used lated */
   if( data->branchvar1 == data->branchvar2 )
   {
      assert(SCIPvarGetType(SCIPgetVars(scip)[data->branchvar1]) != SCIP_VARTYPE_BINARY);
      return SCIP_OKAY;
   }

   SCIP_CALL( level2dataEnsureSize(scip, data) );

   SCIP_CALL( level2resultCreateFromData(scip, data, &result) );

   result->lpobjval = lpobjval;
   result->cutoff = cutoff;
   result->valid = valid;

   /* search for a level 2 result with the same branching decisions*/
   for( i = 0; i < data->nlevel2results; ++i )
   {
      if( level2resultEqual( data->level2results[i], result) )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "same level2 node already processed:\n");
         level2resultPrint(scip, data->level2results[i]);
         level2resultPrint(scip, result);
         *duplicate = TRUE;
      }
   }

   data->level2results[data->nlevel2results] = result;
   ++data->nlevel2results;
   assert(data->nlevel2results <= data->level2resultssize);

   return SCIP_OKAY;
}


/** The data that is preserved over multiple runs of the branching rule. */
typedef struct
{
   BRANCHINGDECISION*    olddecision;        /**< The previous decision that gets used for the case that in the previous run
                                              *   only non-violating implied binary constraints were added.*/
   SCIP_Longint          oldnnodelpiterations; /**< node LP iterations when previous branching decision was stored */
   SCIP_Longint          oldnnodelps;        /**< node LPs when previous branching decision was stored */
   SCIP_Longint          oldntotalnodes;     /**< node at which previous branching decision was stored */
   SCIP_Longint*         lastbranchid;       /**< The node id at which the var was last branched on (for a given branching
                                              *   var). */
   SCIP_Longint*         lastbranchnlps;     /**< The number of (non-probing) LPs that where solved when the var was last
                                              *   branched on. */
   SCIP_Real*            lastbranchlpobjval; /**< The lp objval at which var was last branched on. */
   BRANCHINGRESULTDATA** lastbranchupres;    /**< The result of the last up branching for a given var. */
   BRANCHINGRESULTDATA** lastbranchdownres;  /**< The result of the last down branching for a given var. */
   int                   restartindex;       /**< The index at which the iteration over the number of candidates starts. */
   int                   nvars;              /**< The number of variables that can be stored in the arrays. */
} PERSISTENTDATA;

/** The parameter that can be changed by the user/caller and alter the behaviour of the lookahead branching. */
typedef struct
{
   SCIP_Longint          reevalage;          /**< The number of "normal" (not probing) lps that may have been solved before
                                              *   we stop using old data and start recalculating new first level data. */
   SCIP_Longint          reevalagefsb;       /**< The number of "normal" (not probing) lps that may have been solved before
                                              *   we stop using old FSB data and start recalculating new first level data. */
   int                   maxnviolatedcons;   /**< The number of constraints (domain reductions and binary constraints) we
                                              *   want to gather before restarting the run. Set to -1 for an unbounded
                                              *   number of constraints. */
   int                   maxnviolatedbincons;/**< The number of binary constraints we want to gather before restarting the
                                              *   run. Set to -1 for an undbounded number of binary constraints. */
   int                   maxnviolateddomreds;/**< The number of domain reductions we want to gather before restarting the
                                              *   run. Set to -1 for an undbounded number of domain reductions. */
   int                   recursiondepth;     /**< How deep should the recursion go? Default for Lookahead: 2 */
   int                   maxncands;          /**< If abbreviated == TRUE, at most how many candidates should be handled at the base node? */
   int                   maxndeepercands;    /**< If abbreviated == TRUE, at most how many candidates should be handled in deeper nodes? */
   SCIP_Bool             usedomainreduction; /**< indicates whether the data for domain reductions should be gathered and
                                              *   used. */
   SCIP_Bool             mergedomainreductions; /**< should domain reductions of feasible siblings should be merged? */
   SCIP_Bool             prefersimplebounds; /**<    should domain reductions only be applied if there are simple bound changes? */
   SCIP_Bool             onlyvioldomreds;    /**< Should only domain reductions that violate the LP solution be applied? */
   SCIP_Bool             usebincons;         /**< indicates whether the data for the implied binary constraints should
                                              *   be gathered and used */
   int                   addbinconsrow;      /**< should binary constraints be added as rows to the base LP?
                                              *   (0: no, 1: separate, 2: as initial rows) */
   SCIP_Bool             addnonviocons;      /**< Should constraints be added, that are not violated by the base LP? */
   SCIP_Bool             abbreviated;        /**< Should the abbreviated version be used? */
   SCIP_Bool             reusebasis;         /**< If abbreviated == TRUE, should the solution lp-basis of the FSB run be
                                              *   used in the first abbreviated level?  */
   SCIP_Bool             storeunviolatedsol; /**< Should a solution/decision be stored, to speed up the next iteration
                                              *   after adding the constraints/domreds? */
   SCIP_Bool             abbrevpseudo;       /**< If abbreviated == TRUE, should pseudocost values be used, to approximate
                                              *   the scoring? */
   SCIP_Bool             level2avgscore;     /**< should the average score be used for uninitialized scores in level 2? */
   SCIP_Bool             level2zeroscore;    /**< should uninitialized scores in level 2 be set to zero? */
   SCIP_Bool             addclique;          /**< add binary constraints with two variables found at the root node also as a clique? */
   SCIP_Bool             propagate;          /**< Should the problem be propagated before solving each inner node? */
   SCIP_Bool             uselevel2data;      /**< should branching data generated at depth level 2 be stored for re-using it? */
   SCIP_Bool             applychildbounds;   /**< should bounds known for child nodes be applied? */
   SCIP_Bool             enforcemaxdomreds;  /**< should the maximum number of domain reductions maxnviolateddomreds be enforced? */
   SCIP_Bool             updatebranchingresults; /**< should branching results (and scores) be updated w.r.t. proven dual bounds? */
   SCIP_Bool             inscoring;          /**< are we currently in FSB-scoring (only used internally) */
   int                   maxproprounds;      /**< maximum number of propagation rounds to perform at temporary nodes
                                              *   (-1: unlimited, 0: SCIP default) */
   char                  scoringfunction;    /**< scoring function at base level */
   char                  deeperscoringfunction; /**< scoring function at deeper levels */
   char                  scoringscoringfunction;/**< scoring function for FSB scoring */
   SCIP_Real             minweight;          /**< weight of the min gain of two child problems */
   SCIP_Real             worsefactor;        /**< if the FSB score is of a candidate is worse than the best by this factor, skip this candidate (-1: disable) */
   SCIP_Bool             filterbymaxgain;    /**< should lookahead branching only be applied if the max gain in level 1 is not uniquely that of the best candidate? */
} CONFIGURATION;


#if defined(SCIP_DEBUG) || defined(SCIP_STATISTIC)
#define MAXRESULT SCIP_DELAYNODE

/** returns a human readable name for the given result enum value */
static
const char* getStatusString(
   SCIP_RESULT           result              /**< enum value to get the string representation for */
   )
{
   assert(result >= 1);
   assert(result <= 18);

   switch( result )
   {
   case SCIP_DIDNOTRUN:
      return "SCIP_DIDNOTRUN";
   case SCIP_DELAYED:
      return "SCIP_DELAYED";
   case SCIP_DIDNOTFIND:
      return "SCIP_DIDNOTFIND";
   case SCIP_FEASIBLE:
      return "SCIP_FEASIBLE";
   case SCIP_INFEASIBLE:
      return "SCIP_INFEASIBLE";
   case SCIP_UNBOUNDED:
      return "SCIP_UNBOUNDED";
   case SCIP_CUTOFF:
      return "SCIP_CUTOFF";
   case SCIP_SEPARATED:
      return "SCIP_SEPARATED";
   case SCIP_NEWROUND:
      return "SCIP_NEWROUND";
   case SCIP_REDUCEDDOM:
      return "SCIP_REDUCEDDOM";
   case SCIP_CONSADDED:
      return "SCIP_CONSADDED";
   case SCIP_CONSCHANGED:
      return "SCIP_CONSCHANGED";
   case SCIP_BRANCHED:
      return "SCIP_BRANCHED";
   case SCIP_SOLVELP:
      return "SCIP_SOLVELP";
   case SCIP_FOUNDSOL:
      return "SCIP_FOUNDSOL";
   case SCIP_SUSPENDED:
      return "SCIP_SUSPENDED";
   case SCIP_SUCCESS:
      return "SCIP_SUCCESS";
   case SCIP_DELAYNODE:
      return "SCIP_DELAYNODE";
   default:
      SCIPerrorMessage("result code %d not treated in lookahead branching rule\n", result);
      SCIPABORT();
      return "UNKNOWN";
   }
}
#endif

#ifdef SCIP_STATISTIC
/** The data used for some statistical analysis. */
typedef struct
{
   int*                  nresults;           /**< Array of counters for each result state the lookahead branching finished.
                                              *   The first (0) entry is unused, as the result states are indexed 1-based
                                              *   and we use this index as our array index. */
   int*                  nsinglecutoffs;     /**< The number of single cutoffs on a (probing) node per probingdepth. */
   int*                  nfullcutoffs;       /**< The number of double cutoffs on a (probing) node per probingdepth. */
   int*                  nlpssolved;         /**< The number of all lps solved for a given probingdepth (incl. FSB). */
   int*                  nlpssolvedfsb;      /**< The number of lps solved by the initial FSB to get the FSB scores. */
   int*                  nduplicatelps;      /**< The number of lps solved for duplicate grand-child nodes. */
   SCIP_Longint*         nlpiterations;      /**< The number of all lp iterations needed for a given probingdepth
                                              *   (incl. FSB). */
   SCIP_Longint*         nlpiterationsfsb;   /**< The number of lp iterations needed to get the FSB scores. */
   int*                  npropdomred;        /**< The number of domain reductions based on domain propagation per
                                              *   progingdepth. */
   int*                  noldbranchused;     /**< The number of times old branching data is used (see the reevalage
                                              *   parameter in the CONFIGURATION struct) */
   int*                  noldbranchusedfsb;  /**< The number of times old FSB scoring data is used (see the reevalagefsb
                                              *   parameter in the CONFIGURATION struct) */
   int*                  chosenfsbcand;      /**< If abbreviated, this is the number of times each candidate was finally
                                              *   chosen by the following LAB */
   int*                  stopafterfsb;       /**< If abbreviated, this is the number of times the rule was stopped after
                                              *   scoring candidates by FSB, e.g., by adding constraints or domreds. */
   int*                  cutoffafterfsb;     /**< If abbreviated, this is the number of times the rule was stopped after
                                              *   scoring candidates by FSB because of a found cutoff. */
   int*                  domredafterfsb;     /**< If abbreviated, this is the number of times the rule was stopped after
                                              *   scoring candidates by FSB because of a found domain reduction. */
   int                   nsinglecandidate;   /**< number of times a single candidate was given to the recursion routine */
   int                   nsingleafterfilter; /**< number of times a single candidate remained after filtering */
   int                   noldcandidate;      /**< number of times the old candidate from last call with nonviolating
                                              *   reductions was branched on */
   int                   nlperrorcalls;      /**< number of times an LP error occured and LAB branched without completely
                                              *   evaluating all candidates */
   int                   nlimitcalls;        /**< number of times a time limit was reached and LAB branched without
                                              *   completely evaluating all candidates */
   int                   ntotalresults;      /**< The total sum of the entries in nresults. */
   int                   nbinconst;          /**< The number of binary constraints added to the base node. */
   int                   nbinconstvio;       /**< The number of binary constraints added to the base node, that are violated
                                              *   by the LP at that node. */
   int                   ndomred;            /**< The number of domain reductions added to the base node. */
   int                   ndomredvio;         /**< The number of domain reductions added to the base node, that are violated
                                              *   by the LP at that node. */
   int                   ndepthreached;      /**< The number of times the branching was aborted due to a too small depth. */
   int                   ndomredcons;        /**< The number of binary constraints ignored, as they would be dom reds. */
   int                   ncutoffproofnodes;  /**< The number of nodes needed to prove all found cutoffs. */
   int                   ndomredproofnodes;  /**< The number of nodes needed to prove all found domreds. */
   int                   ncliquesadded;      /**< The number of cliques added in the root node. */
   int                   maxnbestcands;      /**< if abbreviated, this is the maximum number of candidates to investigate */
   int                   recursiondepth;     /**< The recursiondepth of the LAB. Can be used to access the depth-dependent
                                              *   arrays contained in the statistics. */
} STATISTICS;

/** Initializes the statistics with the start values. */
static
void statisticsInit(
   STATISTICS*           statistics          /**< the statistics to be initialized */
   )
{
   int i;

   assert(statistics != NULL);
   assert(statistics->recursiondepth > 0);

   statistics->nsinglecandidate = 0;
   statistics->nsingleafterfilter = 0;
   statistics->noldcandidate = 0;
   statistics->nlperrorcalls = 0;
   statistics->nlimitcalls = 0;
   statistics->ntotalresults = 0;
   statistics->nbinconst = 0;
   statistics->nbinconstvio = 0;
   statistics->ndomredvio = 0;
   statistics->ndepthreached = 0;
   statistics->ndomred = 0;
   statistics->ndomredcons = 0;
   statistics->ncutoffproofnodes = 0;
   statistics->ndomredproofnodes = 0;
   statistics->ncliquesadded = 0;

   for( i = 0; i <= MAXRESULT; i++)
   {
      statistics->nresults[i] = 0;
   }

   for( i = 0; i < statistics->recursiondepth; i++ )
   {
      statistics->noldbranchused[i] = 0;
      statistics->noldbranchusedfsb[i] = 0;
      statistics->npropdomred[i] = 0;
      statistics->nfullcutoffs[i] = 0;
      statistics->nlpssolved[i] = 0;
      statistics->nlpssolvedfsb[i] = 0;
      statistics->nduplicatelps[i] = 0;
      statistics->nlpiterations[i] = 0;
      statistics->nlpiterationsfsb[i] = 0;
      statistics->nsinglecutoffs[i] = 0;
      statistics->stopafterfsb[i] = 0;
      statistics->cutoffafterfsb[i] = 0;
      statistics->domredafterfsb[i] = 0;
   }

   for( i = 0; i < statistics->maxnbestcands; i++ )
   {
      statistics->chosenfsbcand[i] = 0;
   }
}

/** Prints the content of the statistics to stdout. */
static
void statisticsPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   STATISTICS*           statistics          /**< the statistics to print */
   )
{
   assert(scip != NULL);
   assert(statistics != NULL);
   assert(statistics->recursiondepth > 0);

   /* only print something, if we have any statistics */
   if( statistics->ntotalresults > 0 )
   {
      int i;

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Lookahead Branching was called <%i> times.\n", statistics->ntotalresults);

      for( i = 1; i <= MAXRESULT; i++ )
      {
         SCIP_RESULT currentresult = (SCIP_RESULT)i;
         /* see type_result.h for the id <-> enum mapping */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Result <%s> was chosen <%i> times\n", getStatusString(currentresult),
            statistics->nresults[i]);
      }

      for( i = 0; i < statistics->maxnbestcands; i++ )
      {
         if( statistics->chosenfsbcand[i] > 0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "The %i. variable (w.r.t. the FSB score) was chosen as the final result %i times.\n",
               i+1, statistics->chosenfsbcand[i]);
         }
      }

      for( i = 0; i < statistics->recursiondepth; i++ )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "In depth <%i>, branching was stopped after the scoring FSB %i times, %i times because of a cutoff and %i times because of a domain reduction\n",
            i, statistics->stopafterfsb[i], statistics->cutoffafterfsb[i], statistics->domredafterfsb[i]);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "In depth <%i>, <%i> fullcutoffs and <%i> single cutoffs were found.\n",
            i, statistics->nfullcutoffs[i], statistics->nsinglecutoffs[i]);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "In depth <%i>, <%i> LPs were solved, <%i> of them to calculate the FSB score, <%i> were saved for duplicate grandchildren.\n",
            i, statistics->nlpssolved[i], statistics->nlpssolvedfsb[i], statistics->nduplicatelps[i]);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "In depth <%i>, <%" SCIP_LONGINT_FORMAT "> iterations were needed to solve the LPs, <%"
            SCIP_LONGINT_FORMAT "> of them to calculate the FSB score.\n", i, statistics->nlpiterations[i],
            statistics->nlpiterationsfsb[i]);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "In depth <%i>, a decision was discarded <%i> times due to domain reduction because of"
            " propagation.\n", i, statistics->npropdomred[i]);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "In depth <%i>, old LAB branching results were used in <%i> cases, old FSB scores in <%d> cases.\n",
            i, statistics->noldbranchused[i], statistics->noldbranchusedfsb[i]);
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "One single branching candidate was given <%i> times, after filtering, a single candidate remained <%i> times.\n",
         statistics->nsinglecandidate, statistics->nsingleafterfilter);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "The old branching candidate was used <%i> times.\n",
         statistics->noldcandidate);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "An LP error led to branching before all candidates were evaluated <%i> times.\n",
         statistics->nlperrorcalls);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "A reached (time) limit led to branching before all candidates were evaluated <%i> times.\n",
         statistics->nlimitcalls);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Depth limit was reached <%i> times.\n", statistics->ndepthreached);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Ignored <%i> binary constraints, that would be domain reductions.\n",
         statistics->ndomredcons);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Added <%i> binary constraints, of which <%i> where violated by the base LP.\n",
         statistics->nbinconst, statistics->nbinconstvio);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Reduced the domain of <%i> vars, <%i> of them where violated by the base LP.\n",
         statistics->ndomred, statistics->ndomredvio);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Added <%i> cliques found as binary constraint in the root node\n",
         statistics->ncliquesadded);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Needed <%i> additional nodes to prove the cutoffs of base nodes\n",
         statistics->ncutoffproofnodes);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Needed <%i> additional nodes to prove the domain reductions\n",
         statistics->ndomredproofnodes);
   }
}

/** Helper struct to store the statistical data needed in a single run. */
typedef struct
{
   int                   ncutoffproofnodes;  /**< The number of nodes needed to prove the current cutoff. */
} LOCALSTATISTICS;

/** Allocates the local statistics in buffer memory and initializes it with default values. */
static
SCIP_RETCODE localStatisticsAllocate(
   SCIP*                 scip,               /**< SCIP data structure */
   LOCALSTATISTICS**     localstats          /**< pointer to the local statistics to allocate and initialize */
   )
{
   assert(scip != NULL);
   assert(localstats != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, localstats) );

   (*localstats)->ncutoffproofnodes = 0;

   return SCIP_OKAY;
}

/** Frees the allocated buffer memory of the local statistics. */
static
void localStatisticsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   LOCALSTATISTICS**     localstats          /**< pointer to the local statistics to be freed */
   )
{
   assert(scip != NULL);
   assert(localstats != NULL);

   SCIPfreeBuffer(scip, localstats);
}
#endif

/** branching rule data */
struct SCIP_BranchruleData
{
   CONFIGURATION*        config;             /**< the parameter that influence the behaviour of the lookahead branching */
   PERSISTENTDATA*       persistent;         /**< the data that persists over multiple branching decisions */
   SCIP_Bool             isinitialized;      /**< indicates whether the fields in this struct are initialized */
#ifdef SCIP_STATISTIC
   STATISTICS*           statistics;         /**< statistical data container */
#endif
};

/** all constraints that were created and may be added to the base node */
typedef struct
{
   SCIP_VAR***           consvars;           /**< array containing the variables for each constraint to be created */
   int*                  nconsvars;          /**< number of vars in each element of 'consvars' */
   SCIP_Bool*            violated;           /**< indicating whether a constraint is violated by the base solution */
   int                   nelements;          /**< number of elements in 'consvars' and 'nconsvars' */
   int                   memorysize;         /**< number of entries that the array 'consvars' may hold before the
                                              *   array is reallocated. */
   int                   nviolatedcons;      /**< number of constraints that are violated by the base LP solution. */
} CONSTRAINTLIST;

/** Allocate and initialize the list holding the constraints. */
static
SCIP_RETCODE constraintListCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTRAINTLIST**      conslist,           /**< Pointer to the list to be allocated and initialized. */
   int                   startsize           /**< The number of entries the list initially can hold. */
   )
{
   assert(scip != NULL);
   assert(conslist != NULL);
   assert(startsize > 0);

   SCIP_CALL( SCIPallocBuffer(scip, conslist) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*conslist)->consvars, startsize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*conslist)->nconsvars, startsize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*conslist)->violated, startsize) );

   /* We start without any constraints */
   (*conslist)->nelements = 0;
   (*conslist)->memorysize = startsize;
   (*conslist)->nviolatedcons = 0;

   return SCIP_OKAY;
}

/** Append an element to the end of the list of constraints. */
static
SCIP_RETCODE constraintListAppend(
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTRAINTLIST*       list,               /**< list to add the consvars to */
   SCIP_VAR**            consvars,           /**< array of variables for the constraint to be created later */
   int                   nconsvars,          /**< number of elements in 'consvars' */
   SCIP_Bool             violated            /**< indicates whether the constraint is violated by the base lp */
   )
{
   assert(scip != NULL);
   assert(list != NULL);
   assert(consvars != NULL);
   assert(nconsvars > 0);

   /* In case the list tries to hold more elements than it has space, reallocate  */
   if( list->memorysize == list->nelements )
   {
      /* resize the array, such that it can hold the new element */
      int newmemsize = SCIPcalcMemGrowSize(scip, list->memorysize + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &list->consvars, list->memorysize, newmemsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &list->nconsvars, list->memorysize, newmemsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &list->violated, list->memorysize, newmemsize) );
      list->memorysize = newmemsize;
   }

   /* Set the new vars at the first unused place, which is the length used as index */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &list->consvars[list->nelements], consvars, nconsvars) ); /*lint !e866*/
   list->nconsvars[list->nelements] = nconsvars;
   list->violated[list->nelements] = violated;
   list->nelements++;

   return SCIP_OKAY;
}

/** Free all resources of a constraint list in opposite order to the allocation. */
static
void constraintListFree(
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTRAINTLIST**      conslist            /**< Pointer to the list to be freed. */
   )
{
   int i;

   assert(scip != NULL);
   assert(conslist != NULL);

   for( i = 0; i < (*conslist)->nelements; i++ )
   {
      SCIPfreeBlockMemoryArray(scip, &(*conslist)->consvars[i], (*conslist)->nconsvars[i]);
   }

   SCIPfreeBlockMemoryArray(scip, &(*conslist)->violated, (*conslist)->memorysize);
   SCIPfreeBlockMemoryArray(scip, &(*conslist)->nconsvars, (*conslist)->memorysize);
   SCIPfreeBlockMemoryArray(scip, &(*conslist)->consvars, (*conslist)->memorysize);
   SCIPfreeBuffer(scip, conslist);
}

/**
 * list of binary variables currently branched on
 * a down branching (x <= 0) is saved as the negated variable (1-x)
 * an up branching (x >= 1) is saved as the original variable (x)
 * these variables are used to build the binary constraint in case that a ('binary') branch is cut off
 */
typedef struct
{
   SCIP_VAR**            binaryvars;         /**< The binary variables currently branched on. */
   int                   nbinaryvars;        /**< The number of entries in 'nbinaryvars'. */
   int                   memorysize;         /**< The number of entries that the array 'binaryvars' may hold before the
                                              *   array is reallocated. */
} BINARYVARLIST;

/** Allocates and initializes the BINARYVARLIST struct. */
static
SCIP_RETCODE binaryVarListCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYVARLIST**       list,               /**< Pointer to the list to be allocated and initialized. */
   int                   startsize           /**< The number of entries the list initially can hold. */
   )
{
   assert(scip != NULL);
   assert(list != NULL);
   assert(startsize > 0);

   SCIP_CALL( SCIPallocBuffer(scip, list) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*list)->binaryvars, startsize) );

   /* We start with no entries and the (current) max length */
   (*list)->nbinaryvars = 0;
   (*list)->memorysize = startsize;

   return SCIP_OKAY;
}

/** Appends a binary variable to the list, reallocating the list if necessary. */
static
void binaryVarListAppend(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYVARLIST*        list,               /**< The list to add the var to. */
   SCIP_VAR*             vartoadd            /**< The binary var to add to the list. */
   )
{
   assert(scip != NULL);
   assert(list != NULL);
   assert(vartoadd != NULL);
   assert(SCIPvarIsBinary(vartoadd));
   assert(list->nbinaryvars < list->memorysize);

   /* Set the new var at the first unused place, which is the length used as index */
   list->binaryvars[list->nbinaryvars] = vartoadd;
   list->nbinaryvars++;
}

/** Remove the last element from the list. */
static
void binaryVarListDrop(
   BINARYVARLIST*        list                /**< The list to remove the last element from. */
   )
{
   assert(list != NULL);
   assert(list->nbinaryvars > 0);
   assert(list->binaryvars[list->nbinaryvars-1] != NULL);

   /* decrement the number of entries in the actual list */
   list->nbinaryvars--;
}

/** Frees all resources allocated by a BINARYVARLIST in opposite order of allocation. */
static
void binaryVarListFree(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYVARLIST**       list                /**< Pointer to the list to free */
   )
{
   assert(scip != NULL);
   assert(list != NULL);

   SCIPfreeBufferArray(scip, &(*list)->binaryvars);
   SCIPfreeBuffer(scip, list);
}

/** struct holding the relevant data for handling binary constraints */
typedef struct
{
   BINARYVARLIST*        binaryvars;         /**< current binary vars, used to fill the conslist */
   CONSTRAINTLIST*       conslist;           /**< list of constraints to be created */
} BINCONSDATA;

/** Allocate and initialize the BINCONSDATA struct. */
static
SCIP_RETCODE binConsDataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   BINCONSDATA**         consdata,           /**< Pointer to the struct to be allocated and initialized. */
   int                   maxdepth,           /**< The depth of the recursion as an upper bound of branch vars to hold. */
   int                   nstartcons          /**< The start size of the array containing the constraints. */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(maxdepth > 0);
   assert(nstartcons > 0);

   SCIP_CALL( SCIPallocBuffer(scip, consdata) );
   SCIP_CALL( binaryVarListCreate(scip, &(*consdata)->binaryvars, maxdepth) );
   SCIP_CALL( constraintListCreate(scip, &(*consdata)->conslist, nstartcons) );

   return SCIP_OKAY;
}

/** Free all resources in a BINCONSDATA in opposite order of allocation. */
static
void binConsDataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   BINCONSDATA**         consdata            /**< Pointer to the struct to be freed. */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   constraintListFree(scip, &(*consdata)->conslist);
   binaryVarListFree(scip, &(*consdata)->binaryvars);
   SCIPfreeBuffer(scip, consdata);
}

/** A struct acting as a fixed list of candidates */
typedef struct
{
   CANDIDATE**           candidates;         /**< the array of candidates */
   int                   ncandidates;        /**< the number of actual entries in candidates (without trailing NULLs); this
                                              *   is NOT the length of the candidates array, but the number of candidates in
                                              *   it */
} CANDIDATELIST;

/** allocates the candidate list on the buffer WITHOUT initializing the contained array of candidates. */
static
SCIP_RETCODE candidateListCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATELIST**       candidatelist,      /**< the candidate list to allocate */
   int                   ncandidates         /**< the number of candidates the list must hold */
   )
{
   assert(scip != NULL);
   assert(candidatelist != NULL);
   assert(ncandidates >= 0);

   SCIP_CALL( SCIPallocBuffer(scip, candidatelist) );

   if( ncandidates > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(*candidatelist)->candidates, ncandidates) );
   }
   else
      (*candidatelist)->candidates = NULL;

   (*candidatelist)->ncandidates = ncandidates;

   return SCIP_OKAY;
}

/** allocates the given list and fills it with all fractional candidates of the current LP solution. */
static
SCIP_RETCODE candidateListGetAllFractionalCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATELIST**       candidatelist       /**< the list to allocate and fill */
   )
{
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   int i;

   assert(scip != NULL);
   assert(candidatelist != NULL);

   /* get all fractional candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL, NULL) );

   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(lpcandsfrac != NULL);

   SCIP_CALL( candidateListCreate(scip, candidatelist, nlpcands) );

   for( i = 0; i < nlpcands; i++ )
   {
      CANDIDATE* candidate;

      SCIP_CALL( candidateCreate(scip, &candidate) );
      assert(candidate != NULL);

      candidate->branchvar = lpcands[i];
      candidate->branchval = lpcandssol[i];
      candidate->fracval = lpcandsfrac[i];

      (*candidatelist)->candidates[i] = candidate;

      LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "created candidate <%s>...\n",
         (candidate) != NULL ? SCIPvarGetName((candidate)->branchvar) : "none");
   }

   return SCIP_OKAY;
}

/** frees the allocated buffer memory of the candidate list and frees the contained candidates. */
static
SCIP_RETCODE candidateListFree(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATELIST**       candidatelist       /**< the list to be freed */
   )
{
   int i;

   assert(scip != NULL);
   assert(candidatelist != NULL);
   assert((*candidatelist)->ncandidates > 0 || (*candidatelist)->candidates == NULL);

   if( (*candidatelist)->candidates != NULL )
   {
      for( i = (*candidatelist)->ncandidates - 1; i >= 0; i-- )
      {
         CANDIDATE* cand = (*candidatelist)->candidates[i];
         if( cand != NULL )
         {
            SCIP_CALL( candidateFree(scip, &cand) );
         }
      }

      SCIPfreeBufferArray(scip, &(*candidatelist)->candidates);
   }
   SCIPfreeBuffer(scip, candidatelist);

   return SCIP_OKAY;
}

/** keeps only the first candidates and frees the remaining ones */
static
SCIP_RETCODE candidateListKeep(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATELIST*        candidatelist,      /**< the list to allocate and fill */
   int                   nindices            /**< the number of candidates to keep (starting from 0) */
   )
{
   int i;

   assert(scip != NULL);
   assert(candidatelist != NULL);
   assert(0 < nindices);
   assert(nindices <= candidatelist->ncandidates);

   /* only keep the first nindices candidates and free the remaining ones */
   for( i = nindices; i < candidatelist->ncandidates; i++ )
   {
      CANDIDATE* cand = candidatelist->candidates[i];
      if( cand != NULL )
      {
         SCIP_CALL( candidateFree(scip, &cand) );
         candidatelist->candidates[i] = NULL;
      }
   }
   candidatelist->ncandidates = nindices;

   return SCIP_OKAY;
}

/** all domain reductions found through cutoff of branches */
typedef struct
{
   SCIP_Real*            lowerbounds;        /**< The new lower bounds found for each variable in the problem. */
   SCIP_Real*            upperbounds;        /**< The new upper bounds found for each variable in the problem. */
   SCIP_Shortbool*       baselpviolated;     /**< Indicates whether the base lp solution violates the new bounds of a var.*/
   int                   nviolatedvars;      /**< Tracks the number of vars that have a violated (by the base lp) new lower
                                              *   or upper bound. */
   int                   nchangedvars;       /**< Tracks the number of vars, that have a changed domain. (a change on both,
                                              *   upper and lower bound, counts as one.) */
   int                   nsimplebounds;      /**< number of changed bounds resulting from infeasible child nodes */
#ifdef SCIP_STATISTIC
   int*                  lowerboundnproofs;  /**< The number of nodes needed to prove the lower bound for each variable. */
   int*                  upperboundnproofs;  /**< The number of nodes needed to prove the upper bound for each variable. */
#endif
} DOMAINREDUCTIONS;

/** allocate the struct on the buffer and initialize it with the default values */
static
SCIP_RETCODE domainReductionsCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   DOMAINREDUCTIONS**    domreds             /**< The struct that has to be allocated and initialized. */
   )
{
   SCIP_VAR** vars;
   int ntotalvars;
   int v;

   assert(scip != NULL);
   assert(domreds != NULL);

   /* The arrays saves the data for all variables in the problem via the ProbIndex. See SCIPvarGetProbindex() */
   vars = SCIPgetVars(scip);
   ntotalvars = SCIPgetNVars(scip);

   /* Allocate the struct and the contained arrays; initialize flags to FALSE */
   SCIP_CALL( SCIPallocBuffer(scip, domreds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*domreds)->lowerbounds, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*domreds)->upperbounds, ntotalvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &(*domreds)->baselpviolated, ntotalvars) );
#ifdef SCIP_STATISTIC
   SCIP_CALL( SCIPallocClearBufferArray(scip, &(*domreds)->lowerboundnproofs, ntotalvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &(*domreds)->upperboundnproofs, ntotalvars) );
#endif

   for( v = 0; v < ntotalvars; ++v )
   {
      (*domreds)->lowerbounds[v] = SCIPvarGetLbLocal(vars[v]);
      (*domreds)->upperbounds[v] = SCIPvarGetUbLocal(vars[v]);
   }

   /* At the start we have no domain reductions for any variable. */
   (*domreds)->nviolatedvars = 0;
   (*domreds)->nchangedvars = 0;
   (*domreds)->nsimplebounds = 0;

   return SCIP_OKAY;
}

/** frees the given DOMAINREDUCTIONS and all contained Arrays in the opposite order of allocation */
static
void domainReductionsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DOMAINREDUCTIONS**    domreds             /**< Pointer to the struct to be freed. */
   )
{
   assert(scip != NULL);
   assert(domreds != NULL);

#ifdef SCIP_STATISTIC
   SCIPfreeBufferArray(scip, &(*domreds)->upperboundnproofs);
   SCIPfreeBufferArray(scip, &(*domreds)->lowerboundnproofs);
#endif
   SCIPfreeBufferArray(scip, &(*domreds)->baselpviolated);
   SCIPfreeBufferArray(scip, &(*domreds)->upperbounds);
   SCIPfreeBufferArray(scip, &(*domreds)->lowerbounds);
   SCIPfreeBuffer(scip, domreds);
}

/** information about the current status of the branching */
typedef struct
{
   SCIP_Bool             addedbinconss;      /**< were binary constraints added? */
   SCIP_Bool             depthtoosmall;      /**< was the remaining depth too small to branch on? */
   SCIP_Bool             lperror;            /**< did an error occur while solving an LP */
   SCIP_Bool             cutoff;             /**< was the current node cut off? */
   SCIP_Bool             domredcutoff;       /**< was the current node cut off due to domain reductions? */
   SCIP_Bool             domred;             /**< were domain reductions added due to information obtained through
                                              *   branching? */
   SCIP_Bool             limitreached;       /**< was a limit (time, node, user, ...) reached? */
   SCIP_Bool             maxnconsreached;    /**< was the max number of constraints (bin conss and dom red) reached? */
} STATUS;

/** Allocates the status on the buffer memory and initializes it with default values. */
static
SCIP_RETCODE statusCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS**              status              /**< the status to be allocated */
   )
{
   assert(scip != NULL);
   assert(status != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, status) );

   (*status)->addedbinconss = FALSE;
   (*status)->depthtoosmall = FALSE;
   (*status)->lperror = FALSE;
   (*status)->cutoff = FALSE;
   (*status)->domred = FALSE;
   (*status)->domredcutoff = FALSE;
   (*status)->limitreached = FALSE;
   (*status)->maxnconsreached = FALSE;

   return SCIP_OKAY;
}

/** frees the allocated buffer memory of the status */
static
void statusFree(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS**              status              /**< the status to be freed */
   )
{
   assert(scip != NULL);
   assert(status != NULL);
   SCIPfreeBuffer(scip, status);
}

/** container struct to keep the calculated score for each variable */
typedef struct
{
   SCIP_Real*            scores;             /**< the scores for each problem variable */
   SCIP_Real*            downgains;          /**< the downgains for each problem variable */
   SCIP_Real*            upgains;            /**< the upgains for each problem variable */
   CANDIDATE**           bestsortedcands;    /**< array containing the best sorted variable indices w.r.t. their score */
   int                   nbestsortedcands;   /**< number of elements in bestsortedcands */
   SCIP_Real             scoresum;           /**< sum of set scores */
   int                   nsetscores;         /**< number of set scores */
} SCORECONTAINER;

/** resets the array containing the sorted indices w.r.t. their score. */
static
void scoreContainterResetBestSortedCands(
   SCORECONTAINER*       scorecontainer      /**< the score container to reset */
   )
{
   assert(scorecontainer != NULL);

   BMSclearMemoryArray(scorecontainer->bestsortedcands, scorecontainer->nbestsortedcands);
}

/** allocates the score container and inits it with default values */
static
SCIP_RETCODE scoreContainerCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCORECONTAINER**      scorecontainer,     /**< pointer to the score container to init */
   CONFIGURATION*        config              /**< config struct with the user configuration */
   )
{
   int ntotalvars;
   int ncands = config->maxncands;
   int i;

   assert(scip != NULL);
   assert(scorecontainer != NULL);
   assert(config != NULL);

   /* the container saves the score for all variables in the problem via the ProbIndex, see SCIPvarGetProbindex() */
   ntotalvars = SCIPgetNVars(scip);

   if( SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) < ncands )
      ncands = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);

   SCIP_CALL( SCIPallocBuffer(scip, scorecontainer) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*scorecontainer)->scores, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*scorecontainer)->downgains, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*scorecontainer)->upgains, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*scorecontainer)->bestsortedcands, ncands) );

   (*scorecontainer)->nbestsortedcands = ncands;
   (*scorecontainer)->scoresum = 0.0;
   (*scorecontainer)->nsetscores = 0;

   scoreContainterResetBestSortedCands(*scorecontainer);

   /* init the scores to something negative, as scores are always non negative */
   for( i = 0; i < ntotalvars; i++ )
   {
      (*scorecontainer)->scores[i] = -1.0;
      (*scorecontainer)->downgains[i] = -1.0;
      (*scorecontainer)->upgains[i] = -1.0;
   }

   return SCIP_OKAY;
}

/** Finds the insertion index for the given score in the candidate list. The score of each candidate is taken from the
 *  scorecontainer. The first elements of the candidate list have to be sorted, as this method uses binary search to find
 *  the correct insertion point
 */
static
int findInsertionPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCORECONTAINER*       scorecontainer,     /**< container with all the scores for each candidate */
   SCIP_Real             scoretoinsert,      /**< score to find the insertion index for */
   CANDIDATE**           candidates,         /**< candidate list where the first nsorted elements are sorted (w.r.t. their
                                              *   score) */
   int                   ncandidates         /**< number of elements in candidates to consider, starting from 0 */
   )
{
   int left = 0;
   int right = ncandidates - 1;

   assert(scip != NULL);
   assert(scorecontainer != NULL);
   assert(candidates != NULL);
   assert(ncandidates >= 0);

   while( left <= right )
   {
      int mid = left + ((right - left) / 2);
      SCIP_Real midscore = -SCIPinfinity(scip);
      CANDIDATE *midcand = candidates[mid];

      if( midcand != NULL)
      {
         SCIP_VAR* midvar;
         int midindex;

         midvar = midcand->branchvar;
         midindex = SCIPvarGetProbindex(midvar);
         midscore = scorecontainer->scores[midindex];
      }

      if( SCIPisGT(scip, scoretoinsert, midscore) )
         right = mid - 1;
      else
         left = mid + 1;
   }

   return right + 1;
}

/** Inserts the given probindex into the sorted array in the container, moving all indices after it to the right. Then
 *  returns the element that does not fit into the array any longer. */
static
CANDIDATE* scoreContainerUpdateSortOrder(
   SCORECONTAINER*       scorecontainer,     /**< container to insert the index into */
   CANDIDATE*            candidate,          /**< the probindex of a variable to store */
   int                   insertpoint         /**< point to store the index at */
   )
{
   int i;
   CANDIDATE* movecand = candidate;

   assert(scorecontainer != NULL);
   assert(candidate != NULL);
   assert(insertpoint >= 0);

   for( i = insertpoint; i < scorecontainer->nbestsortedcands; i++ )
   {
      CANDIDATE* oldcand = scorecontainer->bestsortedcands[i];
      scorecontainer->bestsortedcands[i] = movecand;
      movecand = oldcand;
   }

   return movecand;
}

/** sets the score for the variable in the score container */
static
SCIP_RETCODE scoreContainerSetScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCORECONTAINER*       scorecontainer,     /**< the container to write into */
   CANDIDATE*            cand,               /**< candidate to add the score for */
   SCIP_Real             score,              /**< score to add */
   SCIP_Real             downgain,           /**< LP gain in down child */
   SCIP_Real             upgain              /**< LP gain in up child */
   )
{
   CANDIDATE* droppedcandidate;
   int probindex;
   int insertpoint;

   assert(scip != NULL);
   assert(scorecontainer != NULL);
   assert(cand != NULL);
   assert(SCIPisGE(scip, score, -0.2));

   probindex = SCIPvarGetProbindex(cand->branchvar);
   assert(probindex >= 0);

   if( scorecontainer->scores[probindex] < -0.5 )
   {
      ++scorecontainer->nsetscores;
      scorecontainer->scoresum += score;
   }
   else
   {
      scorecontainer->scoresum += (score -  scorecontainer->scores[probindex]);
   }

   scorecontainer->scores[probindex] = score;
   scorecontainer->downgains[probindex] = downgain;
   scorecontainer->upgains[probindex] = upgain;

   /* find the point in the sorted array where the new score should be inserted */
   insertpoint =  findInsertionPoint(scip, scorecontainer, score, scorecontainer->bestsortedcands,
      scorecontainer->nbestsortedcands);

   /* insert the current variable (cand) at the position calculated above, returning the candidate that
    * was removed at the end of the list; this candidate can be the given candidate for the case that the score does not
    * belong to the best ones */
   droppedcandidate = scoreContainerUpdateSortOrder(scorecontainer, cand, insertpoint);

   /* remove the warm start info from the dropped candidate */
   if( droppedcandidate != NULL )
   {
      SCIP_CALL( candidateFreeWarmStartInfo(scip, droppedcandidate) );
   }

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Stored score <%.9g> for var <%s>.\n", score, SCIPvarGetName(cand->branchvar));

   return SCIP_OKAY;
}

/** Frees the score container and all of its contained arrays. */
static
SCIP_RETCODE scoreContainerFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCORECONTAINER**      scorecontainer      /**< score container to free */
   )
{
   assert(scip != NULL);
   assert(scorecontainer != NULL);

   /* don't free the candidates inside the cands array, as those are handled by the candidate list */
   SCIPfreeBufferArray(scip, &(*scorecontainer)->bestsortedcands);
   SCIPfreeBufferArray(scip, &(*scorecontainer)->upgains);
   SCIPfreeBufferArray(scip, &(*scorecontainer)->downgains);
   SCIPfreeBufferArray(scip, &(*scorecontainer)->scores);
   SCIPfreeBuffer(scip, scorecontainer);

   return SCIP_OKAY;
}

/*
 * Local methods for the logic
 */

/** branches recursively on all given candidates */
static
SCIP_RETCODE selectVarRecursive(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS*               status,             /**< current status */
   PERSISTENTDATA*       persistent,         /**< container to store data over multiple calls to the branching rule; or NULL */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   SCIP_SOL*             baselpsol,          /**< base lp solution */
   DOMAINREDUCTIONS*     domainreductions,   /**< container collecting all domain reductions found */
   BINCONSDATA*          binconsdata,        /**< container collecting all binary constraints */
   CANDIDATELIST*        candidatelist,      /**< list of candidates to branch on */
   BRANCHINGDECISION*    decision,           /**< struct to store the final decision */
   SCORECONTAINER*       scorecontainer,     /**< container to retrieve already calculated scores */
   LEVEL2DATA*           level2data,         /**< level 2 LP results data */
   int                   recursiondepth,     /**< remaining recursion depth */
   SCIP_Real             lpobjval,           /**< LP objective value of current probing node*/
   SCIP_Real             baselpobjval,       /**< LP objective value of focus node (not probing) */
   SCIP_Longint*         niterations,        /**< pointer to store the total number of iterations for this variable */
   int*                  ndeepestcutoffs,    /**< pointer to store the total number of cutoffs on the deepest level */
   SCIP_Real*            bestgain,           /**< pointer to store the best gain found with these candidates */
   SCIP_Real*            totalgains,         /**< pointer to store the sum over all gains that are valid in both children */
   int*                  ntotalgains,        /**< pointer to store the number of gains summed in totalgains */
   int*                  ndeepestnodes       /**< pointer to store the number of nodes processed in the deepest level */
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics          /**< general statistical data */
   ,LOCALSTATISTICS*     localstats          /**< local statistics, may be disregarded */
   ,SCIP_Real*           firstscoreptr       /**< pointer to store score of first candidate, or NULL */
   ,SCIP_Real*           bestscoreptr        /**< pointer to store best score, or NULL */
#endif
   );

/** Adds the given lower bound to the DOMAINREDUCTIONS struct. */
static
void addLowerBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< The variable the bound should be added for. */
   SCIP_Real             lowerbound,         /**< The new lower bound for the variable. */
   SCIP_SOL*             baselpsol,          /**< The LP solution of the base problem. Used to check whether the domain
                                              *   reduction is violated by it. */
   SCIP_Bool             simplechange,       /**< does the change result from an infeasible child node? */
   DOMAINREDUCTIONS*     domainreductions    /**< The struct the domain reduction should be added to. */
#ifdef SCIP_STATISTIC
   ,int                  nproofnodes         /**< The number of nodes needed to prove the new lower bound. */
   ,int                  force               /**< should the number of proof nodes be added even if the bound is known already? */
#endif
   )
{
   int varindex;
   SCIP_Real basesolutionval;

   assert(scip != NULL);
   assert(var != NULL);
   assert(baselpsol != NULL);
   assert(domainreductions != NULL);
#ifdef SCIP_STATISTIC
   assert(nproofnodes >= 0);
#endif

   /* The arrays inside DOMAINREDUCTIONS are indexed via the problem index. */
   varindex = SCIPvarGetProbindex(var);

   lowerbound = SCIPadjustedVarLb(scip, var, lowerbound);

   if( SCIPisLT(scip, domainreductions->lowerbounds[varindex], lowerbound) )
   {
      /* the new lower bound is stronger (greater) than the old one,
       * so we update the bound and number of proof nodes */
      domainreductions->lowerbounds[varindex] = lowerbound;
      domainreductions->nchangedvars++;
      if( simplechange )
         domainreductions->nsimplebounds++;
#ifdef SCIP_STATISTIC
      domainreductions->lowerboundnproofs[varindex] = nproofnodes;
   }
   else
   {
      /* if the given lower bound is equal to the old one we take the smaller number of proof nodes */
      if( SCIPisEQ(scip, domainreductions->lowerbounds[varindex], lowerbound) &&
         (force || domainreductions->lowerboundnproofs[varindex] > nproofnodes) )
         domainreductions->lowerboundnproofs[varindex] = nproofnodes;
#endif
   }

   /* we get the solution value to check whether the domain reduction is violated in the base LP */
   basesolutionval = SCIPgetSolVal(scip, baselpsol, var);

   /* in case the new lower bound is greater than the base solution val and the base solution val is not violated by a
    * previously found bound, we increment the nviolatedvars counter and set the baselpviolated flag */
   if( SCIPisFeasGT(scip, domainreductions->lowerbounds[varindex], basesolutionval)
       && !domainreductions->baselpviolated[varindex] )
   {
      domainreductions->baselpviolated[varindex] = TRUE;
      domainreductions->nviolatedvars++;
   }
}

/** Adds the given upper bound to the DOMAINREDUCTIONS struct. */
static
void addUpperBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< The variable the bound should be added for. */
   SCIP_Real             upperbound,         /**< The new upper bound for the variable. */
   SCIP_SOL*             baselpsol,          /**< The LP solution of the base problem. Used to check whether the domain
                                              *   reduction is violated by it. */
   SCIP_Bool             simplechange,       /**< does the change result from an infeasible child node? */
   DOMAINREDUCTIONS*     domainreductions    /**< The struct the domain reduction should be added to. */
#ifdef SCIP_STATISTIC
   ,int                  nproofnodes         /**< The number of nodes needed to prove the new lower bound. */
   ,int                  force               /**< should the number of proof nodes be added even if the bound is known already? */
#endif
   )
{
   int varindex;
   SCIP_Real basesolutionval;

   assert(scip != NULL);
   assert(var != NULL);
   assert(baselpsol != NULL);
   assert(domainreductions != NULL);
#ifdef SCIP_STATISTIC
   assert(nproofnodes >= 0);
#endif

   /* The arrays inside DOMAINREDUCTIONS are indexed via the problem index. */
   varindex = SCIPvarGetProbindex(var);

   upperbound = SCIPadjustedVarUb(scip, var, upperbound);

   if( SCIPisLE(scip, domainreductions->upperbounds[varindex], upperbound) )
   {
#ifdef SCIP_STATISTIC
      /* if the given upper bound is equal to the old one we take the smaller number of proof nodes */
      if( SCIPisEQ(scip, domainreductions->upperbounds[varindex], upperbound) &&
         (force || domainreductions->upperboundnproofs[varindex] > nproofnodes) )
         domainreductions->upperboundnproofs[varindex] = nproofnodes;
#endif
   }
   else
   {
      /* the new upper bound is stronger (smaller) than the old one,
       * so we update the bound and number of proof nodes */
      domainreductions->upperbounds[varindex] = upperbound;
      domainreductions->nchangedvars++;
      if( simplechange )
         domainreductions->nsimplebounds++;
#ifdef SCIP_STATISTIC
      domainreductions->upperboundnproofs[varindex] = nproofnodes;
#endif
   }

   /* We get the solution value to check whether the domain reduction is violated in the base LP */
   basesolutionval = SCIPgetSolVal(scip, baselpsol, var);

   /* In case the new upper bound is smaller than the base solution val and the base solution val is not violated by a
    * previously found bound, we increment the nviolatedvars counter and set the baselpviolated flag.  */
   if( SCIPisFeasLT(scip, domainreductions->upperbounds[varindex], basesolutionval)
       && !domainreductions->baselpviolated[varindex] )
   {
      domainreductions->baselpviolated[varindex] = TRUE;
      domainreductions->nviolatedvars++;
   }
}

/** apply the domain reductions from a single struct to another one; this may be used in case one of the two child
 *  problems of a variable is infeasible
 */
static
void applySingleDeeperDomainReductions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             baselpsol,          /**< The LP solution of the base problem. Used to check whether the domain
                                              *   reduction is violated by it. */
   int                   maxstoredomreds,    /**< maximum number of domain reductions to store */
   DOMAINREDUCTIONS*     targetdomreds,      /**< The target that should be filled with the merged data. */
   DOMAINREDUCTIONS*     domreds             /**< source domain reductions */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(baselpsol != NULL);
   assert(targetdomreds != NULL);
   assert(domreds != NULL);

   /* as the bounds are tracked for all vars we have to iterate over all vars */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   assert(vars != NULL);
   assert(nvars > 0);

   for( i = 0; i < nvars; i++ )
   {
      if( targetdomreds->nviolatedvars >= maxstoredomreds )
         return;

#ifdef SCIP_STATISTIC
      /* adjust the proof nodes */
      addLowerBound(scip, vars[i], domreds->lowerbounds[i], baselpsol, TRUE, targetdomreds,
         domreds->lowerboundnproofs[i], FALSE);
#else
      addLowerBound(scip, vars[i], domreds->lowerbounds[i], baselpsol, TRUE, targetdomreds);
#endif

      if( targetdomreds->nviolatedvars >= maxstoredomreds )
         return;

#ifdef SCIP_STATISTIC
      addUpperBound(scip, vars[i], domreds->upperbounds[i], baselpsol, TRUE, targetdomreds,
         domreds->upperboundnproofs[i], FALSE);
#else
      addUpperBound(scip, vars[i], domreds->upperbounds[i], baselpsol, TRUE, targetdomreds);
#endif
   }
}

/**
 * merges the domain reduction data from the two given branching children data into the target parent data
 */
static
void applyDeeperDomainReductions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             baselpsol,          /**< The LP solution of the base problem. Used to check whether the domain
                                              *   reduction is violated by it. */
   int                   maxstoredomreds,    /**< maximum number of domain reductions to store */
   DOMAINREDUCTIONS*     targetdomreds,      /**< The target that should be filled with the merged data. */
   DOMAINREDUCTIONS*     downdomreds,        /**< One of the source DOMAINREDUCTIONS. */
   DOMAINREDUCTIONS*     updomreds           /**< The other source DOMAINREDUCTIONS. */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(baselpsol != NULL);
   assert(targetdomreds != NULL);
   assert(downdomreds != NULL);
   assert(updomreds != NULL);

   /* as the bounds are tracked for all vars we have to iterate over all vars */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   assert(vars != NULL);
   assert(nvars > 0);

   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Combining domain reductions from up and down child.\n");
   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Previous number of changed variable domains: %d\n",
      targetdomreds->nchangedvars);

   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Number of changed variable domains in up child: %d\n",
      updomreds->nchangedvars);
   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Number of changed variable domains in down child: %d\n",
      downdomreds->nchangedvars);

   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real newlowerbound;
      SCIP_Real newupperbound;

      assert(vars[i] != NULL);

      if( targetdomreds->nviolatedvars >= maxstoredomreds )
         return;

      /* the MIN of both lower bounds represents a valid lower bound at the parent node */
      newlowerbound = MIN(downdomreds->lowerbounds[i], updomreds->lowerbounds[i]);

      /* This MIN can now be added via the default add method */
#ifdef SCIP_STATISTIC
      addLowerBound(scip, vars[i], newlowerbound, baselpsol, FALSE, targetdomreds,
         MIN(4, downdomreds->lowerboundnproofs[i] + updomreds->lowerboundnproofs[i] + 2), FALSE);
#else
      addLowerBound(scip, vars[i], newlowerbound, baselpsol, FALSE, targetdomreds);
#endif

      if( targetdomreds->nviolatedvars >= maxstoredomreds )
         return;

      /* the MAX of both upper bounds represents a valid upper bound at the parent node */
      newupperbound = MAX(downdomreds->upperbounds[i], updomreds->upperbounds[i]);

      /* This MAX can now be added via the default add method */
#ifdef SCIP_STATISTIC
      addUpperBound(scip, vars[i], newupperbound, baselpsol, FALSE, targetdomreds,
         MIN(4, downdomreds->upperboundnproofs[i] + updomreds->upperboundnproofs[i] + 2), FALSE);
#else
      addUpperBound(scip, vars[i], newupperbound, baselpsol, FALSE, targetdomreds);
#endif
   }

   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Subsequent number of changed variable domains: %d\n",
      targetdomreds->nchangedvars);
}

/** Applies the domain reductions to the current node. */
static
SCIP_RETCODE applyDomainReductions(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   SCIP_SOL*             baselpsol,          /**< The LP solution of the base problem. Used to check whether the domain
                                              *   reduction is violated by it. */
   DOMAINREDUCTIONS*     domreds,            /**< The domain reductions that should be applied to the current node. */
   SCIP_Bool*            domredcutoff,       /**< pointer to store whether a cutoff was found due to domain reductions */
   SCIP_Bool*            domred              /**< pointer to store whether a domain change was added */
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics          /**< The statistics container. */
#endif
   )
{
   int i;
   SCIP_VAR** probvars;
   int nprobvars;
#if defined(SCIP_DEBUG) || defined(SCIP_STATISTIC)
   int nboundsadded = 0;
   int nboundsaddedvio = 0;
#endif

   assert(scip != NULL);
   assert(baselpsol != NULL);
   assert(domreds != NULL);
   assert(domredcutoff != NULL);
   assert(domred != NULL);
#ifdef SCIP_STATISTIC
   assert(statistics != NULL);
#endif

   /* initially we have no cutoff */
   *domredcutoff = FALSE;

   /* as the bounds are tracked for all vars we have to iterate over all vars */
   probvars = SCIPgetVars(scip);
   nprobvars = SCIPgetNVars(scip);

   assert(probvars != NULL);
   assert(nprobvars > 0);

   if( config->prefersimplebounds && domreds->nsimplebounds == 0 )
      return SCIP_OKAY;

   for( i = 0; i < nprobvars && !(*domredcutoff); i++ )
   {
      SCIP_VAR* var;
      SCIP_Real proposedbound;
      SCIP_Real baselpval;
#ifdef SCIP_DEBUG
      SCIP_Real oldbound;
#endif
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      var = probvars[i];

      assert(var != NULL);

      baselpval = SCIPgetSolVal(scip, baselpsol, var);

      if( SCIPisGT(scip, domreds->lowerbounds[i], SCIPvarGetLbLocal(var)) )
      {
         /* apply lower bound */
#ifdef SCIP_DEBUG
         oldbound = SCIPvarGetLbLocal(var);
#endif
         proposedbound = domreds->lowerbounds[i];

         if( config->onlyvioldomreds && SCIPisGE(scip, baselpval, proposedbound) )
            continue;

         /* add the new bound */
         SCIP_CALL( SCIPtightenVarLb(scip, var, proposedbound, TRUE, &infeasible, &tightened) );

#if defined(SCIP_DEBUG) || defined(SCIP_STATISTIC)
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Variable <%s>, old lower bound <%g>, proposed lower bound <%g>, new "
            "lower bound <%g>\n", SCIPvarGetName(var), oldbound, proposedbound, SCIPvarGetLbLocal(var));
#endif

         if( infeasible )
         {
            /* the domain reduction may result in an empty model (ub < lb) */
            *domredcutoff = TRUE;
            LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The domain reduction of variable <%s> resulted in an empty "
               "model.\n", SCIPvarGetName(var));
         }
         else if( tightened )
         {
            /* the lb is now strictly greater than before */
            *domred = TRUE;
#if defined(SCIP_DEBUG) || defined(SCIP_STATISTIC)
            nboundsadded++;
#endif
#ifdef SCIP_STATISTIC
            LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The lower bound of variable <%s> was successfully tightened (%d).\n",
               SCIPvarGetName(var), domreds->lowerboundnproofs[i]);
            statistics->ndomredproofnodes += domreds->lowerboundnproofs[i];
#endif

#if defined(SCIP_DEBUG) || defined(SCIP_STATISTIC)
            if( SCIPisLT(scip, baselpval, SCIPvarGetLbLocal(var)) )
            {
               LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The lower bound of variable <%s> is violated by the base lp "
                  "value <%g>.\n", SCIPvarGetName(var), baselpval);

               nboundsaddedvio++;
            }
#endif
         }
      }

      if( SCIPisLT(scip, domreds->upperbounds[i], SCIPvarGetUbLocal(var)) )
      {
         /* apply upper bound */
#ifdef SCIP_DEBUG
         oldbound = SCIPvarGetUbLocal(var);
#endif
         proposedbound = domreds->upperbounds[i];

         if( config->onlyvioldomreds && SCIPisLE(scip, baselpval, proposedbound) )
            continue;

         /* add the new bound */
         SCIP_CALL( SCIPtightenVarUb(scip, var, proposedbound, TRUE, &infeasible, &tightened) );

#if defined(SCIP_DEBUG) || defined(SCIP_STATISTIC)
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Variable <%s>, old upper bound <%g>, proposed upper bound <%g>, new "
            "upper bound <%g>\n", SCIPvarGetName(var), oldbound, proposedbound, SCIPvarGetUbLocal(var));
#endif

         if( infeasible )
         {
            /* the domain reduction may result in an empty model (ub < lb) */
            *domredcutoff = TRUE;
            LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The domain reduction of variable <%s> resulted in an empty "
               "model.\n", SCIPvarGetName(var));
         }
         else if( tightened )
         {
            /* the ub is now strictly smaller than before */
            *domred = TRUE;
#if defined(SCIP_DEBUG) || defined(SCIP_STATISTIC)
            nboundsadded++;
#endif
#ifdef SCIP_STATISTIC
            LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The upper bound of variable <%s> was successfully tightened (%d).\n",
               SCIPvarGetName(var), domreds->upperboundnproofs[i]);
            statistics->ndomredproofnodes += domreds->upperboundnproofs[i];
#endif

#if defined(SCIP_DEBUG) || defined(SCIP_STATISTIC)
            if( SCIPisGT(scip, baselpval, SCIPvarGetUbLocal(var)) )
            {
               LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The upper bound of variable <%s> is violated by the base lp "
                  "value <%g>.\n", SCIPvarGetName(var), baselpval);

               nboundsaddedvio++;
            }
#endif
         }
      }
   }

#ifdef SCIP_STATISTIC
   statistics->ndomred += nboundsadded;
   statistics->ndomredvio += nboundsaddedvio;
#endif

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Truly changed <%d> domains of the problem, <%d> of them are violated by the "
      "base lp.\n", nboundsadded, nboundsaddedvio);
   return SCIP_OKAY;
}

/** Copies the current LP solution into the given pointer. Needs to be freed after usage! */
static
SCIP_RETCODE copyCurrentSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            lpsol               /**< pointer to store the solution into */
   )
{
   assert(scip != NULL);
   assert(lpsol != NULL);

   /* create temporary solution */
   SCIP_CALL( SCIPcreateLPSol(scip, lpsol, NULL) );

   /* unlink the solution, so that newly solved lps don't have any influence on our copy */
   SCIP_CALL( SCIPunlinkSol(scip, *lpsol) );

   return SCIP_OKAY;
}

/** Executes the branching on a given variable with a given value. */
static
SCIP_RETCODE branchOnVar(
   SCIP*                 scip                /**< SCIP data structure */,
   CONFIGURATION*        config,             /**< config struct with the user configuration */
   BRANCHINGDECISION*    decision            /**< the decision with all the needed data */
   )
{
   SCIP_VAR* bestvar;
   SCIP_Real bestval;
   SCIP_NODE* downchild = NULL;
   SCIP_NODE* upchild = NULL;

   assert(scip != NULL);
   assert(decision != NULL);
   assert(config != NULL);

   bestvar = decision->branchvar;
   bestval = decision->branchval;
   assert(bestvar != NULL);

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Effective branching on var <%s> with value <%g(%g)>. Old domain: [%g..%g].\n",
      SCIPvarGetName(bestvar), bestval, SCIPgetSolVal(scip, NULL, bestvar), SCIPvarGetLbLocal(bestvar), SCIPvarGetUbLocal(bestvar));

   assert(!SCIPisIntegral(scip, bestval));

   /* branch on the given variable */
   assert(SCIPisLT(scip, SCIPvarGetLbLocal(bestvar), bestval));
   assert(SCIPisLT(scip, bestval, SCIPvarGetUbLocal(bestvar)));
   SCIP_CALL( SCIPbranchVarVal(scip, bestvar, bestval, &downchild, NULL, &upchild) );

   SCIPdebugMsg(scip, "down child (node %" SCIP_LONGINT_FORMAT "): branching bound change <%s> <= %g\n",
      SCIPnodeGetNumber(downchild), SCIPvarGetName(bestvar), SCIPfeasFloor(scip, bestval));
   SCIPdebugMsg(scip, "up child (node %" SCIP_LONGINT_FORMAT "): branching bound change <%s> >= %g\n",
      SCIPnodeGetNumber(upchild), SCIPvarGetName(bestvar), SCIPfeasCeil(scip, bestval));

   assert(downchild != NULL);
   assert(upchild != NULL);

   /* update the lower bounds in the children; we must not do this if columns are missing in the LP
    * (e.g., because we are doing branch-and-price) or the problem should be solved exactly */
   if( SCIPallColsInLP(scip) && !SCIPisExactSolve(scip) )
   {
      SCIP_Real bestdown = decision->downdb;
      SCIP_Bool bestdownvalid = decision->downdbvalid;
      SCIP_Real bestup = decision->updb;
      SCIP_Bool bestupvalid = decision->updbvalid;
      SCIP_Real provedbound = decision->proveddb;

      /* update the lower bound for the LPs for further children of both created nodes */
      SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, bestdownvalid ? MAX(bestdown, provedbound) : provedbound) );
      SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, bestupvalid ? MAX(bestup, provedbound) : provedbound) );

      if( decision->boundsvalid && config->applychildbounds )
      {
         SCIP_VAR** vars;
         int nvars;
         int i;

         assert(decision->downlowerbounds != NULL);
         assert(decision->downupperbounds != NULL);
         assert(decision->uplowerbounds != NULL);
         assert(decision->upupperbounds != NULL);

         nvars = SCIPgetNVars(scip);
         vars = SCIPgetVars(scip);

         assert(nvars == decision->boundssize);

         for( i = 0; i < nvars; i++ )
         {
            SCIP_VAR* var = vars[i];
            SCIP_Real currentlb;
            SCIP_Real currentub;
            SCIP_Real newlb = decision->downlowerbounds[i];
            SCIP_Real newub = decision->downupperbounds[i];
            assert(var != NULL);

            currentlb = SCIPvarGetLbLocal(var);
            currentub = SCIPvarGetUbLocal(var);

            /* update the lower bound of the lower child in case it is better than the current one */
            if( SCIPisGT(scip, newlb, currentlb) )
            {
               SCIP_CALL( SCIPchgVarLbNode(scip, downchild, var, newlb) );

               SCIPdebugMsg(scip, "down child (node %" SCIP_LONGINT_FORMAT "): add bound change <%s> >= %g\n",
                  SCIPnodeGetNumber(downchild), SCIPvarGetName(var), newlb);
            }

            /* update the upper bound of the lower child in case it is better than the current one AND it is not the
             * branching variable, as its upper bound is already updated
             */
            if( SCIPisLT(scip, newub, currentub) && var != bestvar )
            {
               SCIP_CALL( SCIPchgVarUbNode(scip, downchild, var, newub) );

               SCIPdebugMsg(scip, "down child (node %" SCIP_LONGINT_FORMAT "): add bound change <%s> <= %g\n",
                  SCIPnodeGetNumber(downchild), SCIPvarGetName(var), newub);
            }

            newlb = decision->uplowerbounds[i];
            newub = decision->upupperbounds[i];

            /* update the lower bound of the upper child in case it is better than the current one AND it is not the
             * branching variable, as its lower bound is already updated
             */
            if( SCIPisGT(scip, newlb, currentlb) && var != bestvar)
            {
               SCIP_CALL( SCIPchgVarLbNode(scip, upchild, var, newlb) );

               SCIPdebugMsg(scip, "up child (node %" SCIP_LONGINT_FORMAT "): add bound change <%s> >= %g\n",
                  SCIPnodeGetNumber(upchild), SCIPvarGetName(var), newlb);
            }

            /* update the upper bound of the upper child in case it is better than the current one */
            if( SCIPisLT(scip, newub, currentub) )
            {
               SCIP_CALL( SCIPchgVarUbNode(scip, upchild, var, newub) );

               SCIPdebugMsg(scip, "up child (node %" SCIP_LONGINT_FORMAT "): add bound change <%s> <= %g\n",
                  SCIPnodeGetNumber(upchild), SCIPvarGetName(var), newub);
            }
         }
      }
   }
   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, " -> down child's lowerbound: %.9g, estimate: %.9g\n",
      SCIPnodeGetLowerbound(downchild), SCIPnodeGetEstimate(downchild));
   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, " -> up child's lowerbound: %.9g, estimate: %.9g\n",
      SCIPnodeGetLowerbound(upchild), SCIPnodeGetEstimate(upchild));

   return SCIP_OKAY;
}

/** Get the number of iterations the last LP needed */
static
SCIP_RETCODE getNIterationsLastLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint*         iterations          /**< pointer to store the number of iterations */
   )
{
   SCIP_LPI* lpi;
   int tmpiter;

   assert(scip != NULL);
   assert(iterations != NULL);

   /* get the LP interface of the last solved LP */
   SCIP_CALL( SCIPgetLPI(scip, &lpi) );

   /* get the number of iterations from the interface */
   SCIP_CALL( SCIPlpiGetIterations(lpi, &tmpiter) );

   *iterations = (SCIP_Longint)tmpiter;

   return SCIP_OKAY;
}

/** Creates a new probing node with a new bound for the given candidate and solves the corresponding LP. */
static
SCIP_RETCODE executeBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,             /**< configuration to control the behavior */
   SCIP_Bool             downbranching,      /**< the branching direction */
   CANDIDATE*            candidate,          /**< the candidate to branch on */
   BRANCHINGRESULTDATA*  resultdata,         /**< pointer to the result data which gets filled with the status */
   SCIP_SOL*             baselpsol,          /**< the base lp solution */
   DOMAINREDUCTIONS*     domreds,            /**< struct to store the domain reduction found during propagation */
   STATUS*               status              /**< status will contain updated lperror and limit fields */
   )
{
   SCIP_Real oldupperbound;
   SCIP_Real oldlowerbound;
   SCIP_Real newbound;
   SCIP_LPSOLSTAT solstat;
   SCIP_VAR* branchvar;
   SCIP_Real branchval;

   assert(scip != NULL);
   assert(candidate != NULL);
   assert(resultdata != NULL);
   assert(status != NULL);
   assert(config != NULL);
   assert(status != NULL);

   branchvar = candidate->branchvar;
   branchval = candidate->branchval;

   assert(branchvar != NULL);
   assert(!SCIPisFeasIntegral(scip, branchval));

   if( downbranching )
   {
      /* round the given value down, so that it can be used as the new upper bound */
      newbound = SCIPfeasFloor(scip, branchval);
   }
   else
   {
      /* round the given value up, so that it can be used as the new lower bound */
      newbound = SCIPfeasCeil(scip, branchval);
   }

   oldupperbound = SCIPvarGetUbLocal(branchvar);
   oldlowerbound = SCIPvarGetLbLocal(branchvar);

#ifdef SCIP_DEBUG
   if( downbranching )
   {
      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "DownBranching: Var=<%s>, Proposed upper bound=<%g>, "
         "old bounds=[<%g>..<%g>], new bounds=[<%g>..<%g>]\n", SCIPvarGetName(branchvar), newbound, oldlowerbound,
         oldupperbound, oldlowerbound, newbound);
   }
   else
   {
      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "UpBranching: Var=<%s>, Proposed lower bound=<%g>, "
         "old bounds=[<%g>..<%g>], new bounds=[<%g>..<%g>]\n", SCIPvarGetName(branchvar), newbound, oldlowerbound,
         oldupperbound, newbound, oldupperbound);
   }
#endif

   if( (downbranching && newbound < oldlowerbound - 0.5)
      || (!downbranching && newbound > oldupperbound + 0.5) )
   {
      /* if lb > ub we can cutoff this node */
      resultdata->cutoff = TRUE;

      return SCIP_OKAY;
   }

   assert(!resultdata->cutoff);

   SCIP_CALL( SCIPnewProbingNode(scip) );

   if( downbranching )
   {
      /* down branching preparations */
      if( SCIPisFeasLT(scip, newbound, oldupperbound) )
      {
         /* If the new upper bound is smaller than the old upper bound and also
          * greater than (or equal to) the old lower bound, we set the new upper bound.
          * oldLowerBound <= newUpperBound < oldUpperBound */
         SCIP_CALL( SCIPchgVarUbProbing(scip, branchvar, newbound) );
      }
   }
   else
   {
      /* up branching preparations */
      if( SCIPisFeasGT(scip, newbound, oldlowerbound) )
      {
         /* If the new lower bound is greater than the old lower bound and also
          * smaller than (or equal to) the old upper bound, we set the new lower bound.
          * oldLowerBound < newLowerBound <= oldUpperBound
          */
         SCIP_CALL( SCIPchgVarLbProbing(scip, branchvar, newbound) );
      }
   }

   /* restore the stored LP data (e.g., the basis) from a filtering run */
   if( candidateHasWarmStartInfo(candidate, downbranching) )
   {
      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Restoring lp information for %s branch of variable <%s>\n",
         downbranching ? "down" : "up", SCIPvarGetName(branchvar));
      SCIP_CALL( candidateLoadWarmStartInfo(scip, candidate, downbranching) );
   }

   /* apply domain propagation */
   if( config->propagate )
   {
      SCIP_Longint ndomredsfound = 0;

      SCIP_CALL( SCIPpropagateProbing(scip, config->maxproprounds, &resultdata->cutoff, &ndomredsfound) );

      if( ndomredsfound > 0 )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Found %" SCIP_LONGINT_FORMAT " domain reductions via propagation.\n", ndomredsfound);

         /* domreds != NULL iff config->usedomainreduction */
         if( domreds != NULL )
         {
            int i;
            SCIP_VAR** problemvars = SCIPgetVars(scip);
            int nproblemvars = SCIPgetNVars(scip);

            assert(problemvars != NULL);

            assert(config->usedomainreduction);

            for( i = 0; i < nproblemvars; i++ )
            {
               SCIP_Real lowerbound;
               SCIP_Real upperbound;
               SCIP_VAR* var = problemvars[i];
               assert(var != NULL);

               lowerbound = SCIPvarGetLbLocal(var);
               upperbound = SCIPvarGetUbLocal(var);
#ifdef SCIP_STATISTIC
               addLowerBound(scip, var, lowerbound, baselpsol, FALSE, domreds, 0, FALSE);
               addUpperBound(scip, var, upperbound, baselpsol, FALSE, domreds, 0, FALSE);
#else
               addLowerBound(scip, var, lowerbound, baselpsol, FALSE, domreds);
               addUpperBound(scip, var, upperbound, baselpsol, FALSE, domreds);
#endif
            }
         }
      }
   }

   if( !resultdata->cutoff )
   {
      /* solve the prepared probing LP */
      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &status->lperror, &resultdata->cutoff) );

      /* store the number of iterations needed */
      SCIP_CALL( getNIterationsLastLP(scip, &resultdata->niterations) );

      solstat = SCIPgetLPSolstat(scip);
      assert(solstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

      /* for us an error occurred, if an error during the solving occurred or the lp could not be solved but was not
       * cutoff */
      status->lperror = status->lperror || (solstat == SCIP_LPSOLSTAT_NOTSOLVED && resultdata->cutoff == FALSE);

      /* if we seem to have reached a {time, iteration}-limit or the user cancelled the execution, we want to stop
       * further calculations and instead return the current calculation state */
      status->limitreached = (solstat == SCIP_LPSOLSTAT_ITERLIMIT) || (solstat == SCIP_LPSOLSTAT_TIMELIMIT);

      if( resultdata->cutoff )
      {
         resultdata->objval = SCIPinfinity(scip);
         resultdata->dualbound = SCIPinfinity(scip);
         resultdata->dualboundvalid = TRUE;
      }
      else if( !status->limitreached && !status->lperror )
      {
         SCIP_Bool foundsol = FALSE;

         SCIP_CALL( SCIPtryStrongbranchLPSol(scip, &foundsol, &resultdata->cutoff) );

         /* if we have no error, we save the new objective value and the cutoff decision in the resultdata */
         resultdata->objval = SCIPgetLPObjval(scip);
         resultdata->dualbound = SCIPgetLPObjval(scip);
         resultdata->dualboundvalid = TRUE;
         resultdata->cutoff = resultdata->cutoff || SCIPisGE(scip, resultdata->objval, SCIPgetCutoffbound(scip));

         assert(solstat != SCIP_LPSOLSTAT_INFEASIBLE || resultdata->cutoff);
      }
   }

   return SCIP_OKAY;
}

/** Creates a logic or constraint based on the given 'consvars'. This array has to consist of the negated
 * versions of the variables present on a cutoff "path" (path means all variables from the root directly
 * to the cutoff node).
 * Let x_1, ..., x_n be the variables on a path to a cutoff with the branchings x_i <= 1 for all i.
 * Summed up the constraints would look like x_1 + ... x_n <= n-1.
 * Let y_i = 1 - x_i. Then we have y_1 + ... + y_n >= 1 which is a logic or constraint.
 */
static
SCIP_RETCODE createBinaryConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,             /**< configuration containing flags changing the behavior */
   SCIP_CONS**           constraint,         /**< pointer to store the created constraint in */
   char*                 constraintname,     /**< name of the new constraint */
   SCIP_VAR**            consvars,           /**< array containing the negated binary vars */
   int                   nconsvars           /**< the number of elements in 'consvars' */
   )
{
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool removable;
   SCIP_Bool enforce = FALSE;
   SCIP_Bool check = FALSE;
   SCIP_Bool propagate = TRUE;
   SCIP_Bool local = TRUE;
   SCIP_Bool modifiable = FALSE;
   SCIP_Bool dynamic = FALSE;
   SCIP_Bool stickingatnode = FALSE;

   assert(scip != NULL);
   assert(config != NULL);
   assert(constraint != NULL);
   assert(constraintname != NULL);
   assert(consvars != NULL);
   assert(nconsvars > 0);

   initial = (config->addbinconsrow == 2);
   separate = (config->addbinconsrow == 1);
   removable = (config->addbinconsrow == 1);

   /* creating a logic or constraint based on the list of vars in 'consvars'.
    * A logic or constraints looks like that: y_1 + ... + y_n >= 1.
    */
   SCIP_CALL( SCIPcreateConsLogicor(scip, constraint, constraintname, nconsvars, consvars, initial, separate, enforce,
         check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   return SCIP_OKAY;
}

/**
 * Create a name for the binary constraint.
 */
static
void createBinaryConstraintName(
   SCIP_VAR**            binaryvars,         /**< the variables contained in the constraint */
   int                   nbinaryvars,        /**< the number of elements in 'binaryvars' */
   char*                 constraintname      /**< the char pointer to store the name in */
   )
{
   int i;

   assert(binaryvars != NULL);
   assert(nbinaryvars > 0);
   assert(constraintname != NULL);
   assert(binaryvars[0] != NULL);

   (void) SCIPsnprintf(constraintname, SCIP_MAXSTRLEN, "lookahead_bin_%s", SCIPvarGetName(binaryvars[0]));

   for( i = 1; i < nbinaryvars; i++ )
   {
      size_t oldlen;
      SCIP_VAR* var = binaryvars[i];
      assert(var != NULL);

      oldlen = strlen(constraintname);
      (void) strncat(constraintname, "_", SCIP_MAXSTRLEN-oldlen);
      (void) strncat(constraintname, SCIPvarGetName(var), SCIP_MAXSTRLEN-oldlen-1);
   }
}

/**
 * Add the constraints found during the lookahead branching.
 * The implied binary bounds were found when two or more consecutive branchings of binary variables were cutoff. Then these
 * branching constraints can be combined into a single 'binary constraint'.
 */
static
SCIP_RETCODE addBinaryConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   BINCONSDATA*          binconsdata,        /**< collected binary constraints */
   SCIP_SOL*             baselpsol           /**< the original lp solution, used to check the violation of the constraint */
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics          /**< statistics data */
#endif
   )
{
   assert(scip != NULL);
   assert(config != NULL);
   assert(binconsdata != NULL);
   assert(baselpsol != NULL);
   assert(binconsdata->binaryvars != NULL);
   assert(binconsdata->binaryvars->nbinaryvars > 0);

   /* if we only have one var for the constraint, we can ignore it as it is already added as a domain reduction. */
   if( binconsdata->binaryvars->nbinaryvars > 1 )
   {
      int i;
      SCIP_VAR** negatedvars;
      SCIP_Real lhssum = 0.0;
      SCIP_Bool violated;

      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Adding binary constraint for <%i> vars.\n",
         binconsdata->binaryvars->nbinaryvars);

      SCIP_CALL( SCIPallocBufferArray(scip, &negatedvars, binconsdata->binaryvars->nbinaryvars) );

      for( i = 0; i < binconsdata->binaryvars->nbinaryvars; i++ )
      {
         SCIP_VAR* var = binconsdata->binaryvars->binaryvars[i];
         assert(var != NULL);
         assert(SCIPvarIsBinary(var));

         SCIP_CALL( SCIPgetNegatedVar(scip, var, &negatedvars[i]) );
         lhssum += SCIPgetSolVal(scip, baselpsol, negatedvars[i]);
      }

      violated = (lhssum < 1);

      if( config->addnonviocons || violated )
      {
         SCIP_CALL( constraintListAppend(scip, binconsdata->conslist, negatedvars,
            binconsdata->binaryvars->nbinaryvars, violated) );

         /* the constraint we will be building is a logic or: we have a list of binary variables that were
          * cutoff while we branched on with >= 1. So we have the constraint: x_1 + ... + x_n <= n-1.
          * Let y = (1-x), then we have an equivalent formulation: y_1 + ... + y_n >= 1. If the base lp
          * is violating this constraint we count this for our number of violated constraints and bounds. */
         if( violated )
            binconsdata->conslist->nviolatedcons++;
      }

      SCIPfreeBufferArray(scip, &negatedvars);
   }
#ifdef SCIP_STATISTIC
   else
   {
      assert(statistics != NULL);
      statistics->ndomredcons++;
   }
#endif

   return SCIP_OKAY;
}

/** applies the binary constraints to the original problem. */
static
SCIP_RETCODE applyBinaryConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            basenode,           /**< original branching node */
   CONSTRAINTLIST*       conslist,           /**< list of constraints to be added */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   SCIP_Bool*            consadded,          /**< pointer to store whether at least one constraint was added */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the original problem was made infeasible */
   SCIP_Bool*            boundchange         /**< pointer to store whether a bound change has been applied by adding the
                                              *   constraint as a clique */
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics          /**< statistics data */
#endif
   )
{
   int nconsadded = 0;
   int i;
#ifdef SCIP_STATISTIC
   int nvioconsadded = 0;

   assert(statistics != NULL);
#endif
   assert(basenode != NULL);
   assert(conslist != NULL);
   assert(config != NULL);
   assert(consadded != NULL);
   assert(cutoff != NULL);
   assert(boundchange != NULL);

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "processing %d binary constraints.\n", conslist->nelements);

   if( conslist->nelements == 0 )
      return SCIP_OKAY;

   for( i = 0; i < conslist->nelements; i++ )
   {
      SCIP_VAR** vars = conslist->consvars[i];
      int nvars = conslist->nconsvars[i];
      int v;
#ifdef SCIP_STATISTIC
      SCIP_Bool violated = conslist->violated[i];
#endif

      assert(vars != NULL);

      for( v = 0; v < nvars; ++v )
      {
         assert(vars[v] != NULL);
         assert(SCIPvarIsBinary(vars[v]));

         if( SCIPvarGetLbLocal(vars[v]) > 0.5 )
            break;
      }

      /* no variable is fixed to 1 yet, so constraint is not redundant */
      if( v == nvars )
      {
         SCIP_CONS* constraint;
         char constraintname[SCIP_MAXSTRLEN];

         /* create a name for the new constraint */
         createBinaryConstraintName(vars, nvars, constraintname);
         /* create the constraint with the freshly created name */
         SCIP_CALL( createBinaryConstraint(scip, config, &constraint, constraintname, vars, nvars) );

#ifdef PRINTNODECONS
         SCIPinfoMessage(scip, NULL, "Created constraint:\n");
         SCIP_CALL( SCIPprintCons(scip, constraint, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
#endif
         /* add the constraint to the given node */
         SCIP_CALL( SCIPaddConsNode(scip, basenode, constraint, NULL) );

         nconsadded++;

#ifdef SCIP_STATISTIC
         if( violated )
            nvioconsadded++;
#endif

         /* release the constraint, as it is no longer needed */
         SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

         /* a 2-variable logicor constraint can be expressend as a clique on the negated variables;
          * add it to the clique table if we are at the root node */
         if( nvars == 2 && config->addclique && SCIPgetNNodes(scip) == 1 )
         {
            SCIP_Bool* values;
            SCIP_Bool infeasible;
            int nbdchgs;

            SCIP_CALL( SCIPallocClearBufferArray(scip, &values, nvars) );

            /* a two-variable logicor constraint x + y >= 1 yields the implication x == 0 -> y == 1, and is represented
             * by the clique inequality ~x + ~y <= 1
             */
            SCIP_CALL( SCIPaddClique(scip, vars, values, nvars, FALSE, &infeasible, &nbdchgs) );

#ifdef SCIP_STATISTIC
            statistics->ncliquesadded++;
#endif

            if( infeasible )
               *cutoff = TRUE;

            if( nbdchgs > 0 )
               *boundchange = TRUE;

            SCIPfreeBufferArray(scip, &values);
         }
      }
   }

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "added %d/%d binary constraints.\n", nconsadded, conslist->nelements);

   if( nconsadded > 0 )
   {
      *consadded = TRUE;

#ifdef SCIP_STATISTIC
      statistics->nbinconst += nconsadded;
      statistics->nbinconstvio += nvioconsadded;
#endif
   }

   return SCIP_OKAY;
}

/** checks whether the given bounds are still the bounds of the given variable */
static
SCIP_Bool areBoundsChanged(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check the bounds of */
   SCIP_Real             lowerbound,         /**< reference lower bound */
   SCIP_Real             upperbound          /**< reference upper bound */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(SCIPisFeasIntegral(scip, lowerbound));
   assert(SCIPisFeasIntegral(scip, upperbound));
   assert(!SCIPisEQ(scip, lowerbound, upperbound));
   assert(SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS);

   /* due to roundings the value might have changed slightly without an actual influence on the integral value */
   return SCIPvarGetLbLocal(var) > lowerbound + 0.5 || SCIPvarGetUbLocal(var) < upperbound - 0.5;
}

/** Checks whether the branching rule should continue or terminate with the currently gathered data */
static
SCIP_Bool isBranchFurther(
   STATUS*               status,             /**< current status */
   SCIP_Bool             checkdomreds        /**< should domain reductions be checked? */
   )
{
   assert(status != NULL);

   return !status->lperror && !status->cutoff && !status->limitreached
      && !status->maxnconsreached && (!checkdomreds || !status->domred);
}

/** Checks whether the branching rule should continue or terminate with the currently gathered data. Additionally decrements
 * the given loopcounter. This is needed to better emulate the behavior of FSB by LAB with a depth of 1. */
static
SCIP_Bool isBranchFurtherLoopDecrement(
   STATUS*               status,             /**< current status */
   int*                  loopcounter         /**< the counter to decrement */
   )
{
   SCIP_Bool branchfurther;

   assert(status != NULL);
   assert(loopcounter != NULL);

   branchfurther = isBranchFurther(status, FALSE);

   if( !branchfurther )
      (*loopcounter)--;

   return branchfurther;
}

/** determines whether the previous LAB result of a variable should be reused */
static
SCIP_Bool isUseOldBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   PERSISTENTDATA*       persistent,         /**< data storage over multiple calls to the rule */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   SCIP_VAR*             branchvar           /**< variable to check */
   )
{
   assert(scip != NULL);
   assert(config != NULL);
   assert(branchvar != NULL);

   /* an old branching can be reused, if we are still at the same node and just a few LPs were solved in between */
   if( config->inscoring )
   {
       return SCIPgetVarStrongbranchNode(scip, branchvar) == SCIPgetNNodes(scip)
          && SCIPgetVarStrongbranchLPAge(scip, branchvar) < config->reevalagefsb;
   }
   else
   {
      return persistent->lastbranchid[SCIPvarGetProbindex(branchvar)] == SCIPgetNNodes(scip)
         && SCIPgetNLPs(scip) - persistent->lastbranchnlps[SCIPvarGetProbindex(branchvar)] < config->reevalage;
   }
}

/** retrieves previous LAB result for the given variable */
static
SCIP_RETCODE getOldBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   PERSISTENTDATA*       persistent,         /**< data storage over multiple calls to the rule */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   SCIP_VAR*             branchvar,          /**< variable to get previous results for */
   BRANCHINGRESULTDATA*  downbranchingresult,/**< pointer to store the previous down result in */
   BRANCHINGRESULTDATA*  upbranchingresult,  /**< pointer to store the previous up result in */
   SCIP_Real*            oldlpobjval         /**< pointer to store the previous base lp objval in */
   )
{
   assert(scip != NULL);
   assert(persistent != NULL);
   assert(branchvar != NULL);
   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);
   assert(oldlpobjval != NULL);

   if( config->inscoring )
   {
      SCIP_CALL( SCIPgetVarStrongbranchLast(scip, branchvar, &downbranchingresult->dualbound, &upbranchingresult->dualbound,
            &downbranchingresult->dualboundvalid, &upbranchingresult->dualboundvalid, NULL, oldlpobjval) );
      downbranchingresult->objval = downbranchingresult->dualbound;
      upbranchingresult->objval = upbranchingresult->dualbound;
   }
   else
   {
      int varindex = SCIPvarGetProbindex(branchvar);

      branchingResultDataCopy(persistent->lastbranchdownres[varindex], downbranchingresult);
      branchingResultDataCopy(persistent->lastbranchupres[varindex], upbranchingresult);
      *oldlpobjval = persistent->lastbranchlpobjval[varindex];
   }

#ifdef SCIP_DEBUG
   {
      SCIP_Real downgain;
      SCIP_Real upgain;

      downgain = MAX(downbranchingresult->dualbound - *oldlpobjval, 0);
      upgain = MAX(upbranchingresult->dualbound - *oldlpobjval, 0);

      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Lookahead branching on variable <%s> already performed (lpage=%"
         SCIP_LONGINT_FORMAT ", down=%.9g (%+g), up=%.9g (%+g))\n", SCIPvarGetName(branchvar),
         SCIPgetNLPs(scip) - persistent->lastbranchnlps[SCIPvarGetProbindex(branchvar)],
         downbranchingresult->dualbound, downgain, upbranchingresult->dualbound, upgain);
   }
#endif

   return SCIP_OKAY;
}

/** stores the LAB result for use in a later call to the branching rule */
static
SCIP_RETCODE updateOldBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   PERSISTENTDATA*       persistent,         /**< data storage over multiple calls to the rule */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   SCIP_VAR*             branchvar,          /**< variable to store previous results for */
   SCIP_Real             branchval,          /**< the value of branchvar */
   BRANCHINGRESULTDATA*  downbranchingresult,/**< down branching result to store */
   BRANCHINGRESULTDATA*  upbranchingresult,  /**< up branching result to store */
   SCIP_Real             lpobjval            /**< base lp obj val */
   )
{
   assert(scip != NULL);
   assert(persistent != NULL);
   assert(branchvar != NULL);
   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   if( config->inscoring )
   {
      SCIP_Longint niterations = downbranchingresult->niterations + upbranchingresult->niterations;

      SCIP_CALL( SCIPsetVarStrongbranchData(scip, branchvar, lpobjval, branchval, downbranchingresult->dualbound,
            upbranchingresult->dualbound, downbranchingresult->dualboundvalid, upbranchingresult->dualboundvalid, niterations,
            INT_MAX) );
   }
   else
   {
      int varindex = SCIPvarGetProbindex(branchvar);

      branchingResultDataCopy(downbranchingresult, persistent->lastbranchdownres[varindex]);
      branchingResultDataCopy(upbranchingresult, persistent->lastbranchupres[varindex]);
      persistent->lastbranchlpobjval[varindex] = lpobjval;
      persistent->lastbranchid[varindex] = SCIPgetNNodes(scip);
      persistent->lastbranchnlps[varindex] = SCIPgetNLPs(scip);
   }

   return SCIP_OKAY;
}

/** calculates the FSB scores for the given candidates */
static
SCIP_RETCODE getFSBResult(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS*               status,             /**< current status */
   PERSISTENTDATA*       persistent,         /**< container to store data over multiple calls to the branching rule; or NULL */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   SCIP_SOL*             baselpsol,          /**< base lp solution */
   DOMAINREDUCTIONS*     domainreductions,   /**< container collecting all domain reductions found; or NULL */
   BINCONSDATA*          binconsdata,        /**< container collecting all binary constraints; or NULL */
   CANDIDATELIST*        candidatelist,      /**< list containing all candidates to consider */
   BRANCHINGDECISION*    decision,           /**< struct to store the final decision */
   SCORECONTAINER*       scorecontainer,     /**< container to retrieve already calculated scores; or NULL */
   LEVEL2DATA*           level2data,         /**< level 2 LP results data */
   SCIP_Real             lpobjval            /**< base LP objective value */
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics          /**< general statistical data */
   ,LOCALSTATISTICS*     localstats          /**< local statistics, may be disregarded */
#endif
   )
{
   assert(scip != NULL);
   assert(config != NULL);
   assert(candidatelist != NULL);
   assert(status != NULL);
   assert(scorecontainer != NULL);
   assert(SCIPinProbing(scip));
   assert(SCIPgetProbingDepth(scip) >= 0 && SCIPgetProbingDepth(scip) < config->recursiondepth);

   /* inform configuration that we are in scoring mode now */
   config->inscoring = TRUE;

#ifdef SCIP_STATISTIC
   SCIP_CALL( selectVarRecursive(scip, status, persistent, config, baselpsol, domainreductions,
         binconsdata, candidatelist, decision, scorecontainer, level2data, 1,
         lpobjval, lpobjval, NULL, NULL, NULL, NULL, NULL, NULL,
         statistics, localstats, NULL, NULL) );
#else
   SCIP_CALL( selectVarRecursive(scip, status, persistent, config, baselpsol, domainreductions,
         binconsdata, candidatelist, decision, scorecontainer, level2data, 1,
         lpobjval, lpobjval, NULL, NULL, NULL, NULL, NULL, NULL) );
#endif

   /* inform configuration that we leave scoring mode now */
   config->inscoring = FALSE;

   return SCIP_OKAY;
}

#ifdef SCIP_DEBUG
/** prints the given candidate list */
static
void printCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VERBLEVEL        lvl,                /**< verbosity level to print the list in */
   CANDIDATELIST*        candidatelist       /**< the list to be printed */
   )
{
   int ncands;
   int i;

   assert(scip != NULL);
   assert(candidatelist != NULL);

   ncands = candidatelist->ncandidates;

   LABdebugMessagePrint(scip, lvl, "[");

   for( i = 0; i < ncands; i++ )
   {
      CANDIDATE* cand = candidatelist->candidates[i];

      assert(cand != NULL);
      assert(cand->branchvar != NULL);

      LABdebugMessagePrint(scip, lvl, "%s", SCIPvarGetName(cand->branchvar));
      if(i != ncands-1)
      {
         LABdebugMessagePrint(scip, lvl, ", ");
      }
   }
   LABdebugMessagePrint(scip, lvl, "]\n");
}
#endif

/** calculates the score based on the down and up branching result */
static
SCIP_Real calculateScoreFromResult(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< variable to get the score for */
   BRANCHINGRESULTDATA*  downbranchingresult,/**< branching result of the down branch */
   BRANCHINGRESULTDATA*  upbranchingresult,  /**< branching result of the up branch */
   SCIP_Real             lpobjval            /**< objective value to get difference to as gain */
   )
{
   SCIP_Real score;
   SCIP_Real downgain = SCIPsumepsilon(scip);
   SCIP_Real upgain = SCIPsumepsilon(scip);

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   /* the gain is the difference of the dualbound of a child and the reference objective value;
    * by bounding it by zero we are safe from numerical troubles
    */
   if( !downbranchingresult->cutoff )
      downgain = MAX(downgain, downbranchingresult->dualbound - lpobjval);
   if( !upbranchingresult->cutoff )
      upgain = MAX(upgain, upbranchingresult->dualbound - lpobjval);

   downgain = 100.0 * downgain;
   upgain = 100.0 * upgain;

   /* in case a child is infeasible and therefore cutoff we take the gain of the other child to receive a somewhat
    * realistic gain for the infeasible child;
    * if both children are infeasible we just reset the initial zero values again
    */
   if( downbranchingresult->cutoff )
      downgain = 2 * upgain;
   if( upbranchingresult->cutoff )
      upgain = 2 * downgain;

   score = SCIPgetBranchScore(scip, branchvar, downgain, upgain);

   return score;
}

/** calculates the score based on the down and up branching result */
static
SCIP_Real calculateScoreFromResult2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< variable to get the score for */
   BRANCHINGRESULTDATA*  downbranchingresult,/**< branching result of the down branch */
   BRANCHINGRESULTDATA*  upbranchingresult,  /**< branching result of the up branch */
   SCIP_Real             lpobjval            /**< objective value to get difference to as gain */
   )
{
   SCIP_Real lpscore;
   SCIP_Real dbscore;
   SCIP_Real score;
   SCIP_Real downgain = SCIPsumepsilon(scip);
   SCIP_Real upgain = SCIPsumepsilon(scip);

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   /* the gain is the difference of the dualbound of a child and the reference objective value;
    * by bounding it by zero we are safe from numerical troubles
    */
   if( !downbranchingresult->cutoff )
      downgain = MAX(downgain, downbranchingresult->objval - lpobjval);
   if( !upbranchingresult->cutoff )
      upgain = MAX(upgain, upbranchingresult->objval - lpobjval);

   downgain = 100.0 * downgain;
   upgain = 100.0 * upgain;

   /* in case a child is infeasible and therefore cutoff we take the gain of the other child to receive a somewhat
    * realistic gain for the infeasible child;
    * if both children are infeasible we just reset the initial zero values again
    */
   if( downbranchingresult->cutoff )
      downgain = 2 * upgain;
   if( upbranchingresult->cutoff )
      upgain = 2 * downgain;

   lpscore = SCIPgetBranchScore(scip, branchvar, downgain, upgain);

   /* the gain is the difference of the dualbound of a child and the reference objective value;
    * by bounding it by zero we are safe from numerical troubles
    */
   if( !downbranchingresult->cutoff )
      downgain = MAX(SCIPsumepsilon(scip), downbranchingresult->dualbound - lpobjval); /*lint !e666*/
   if( !upbranchingresult->cutoff )
      upgain = MAX(SCIPsumepsilon(scip), upbranchingresult->dualbound - lpobjval); /*lint !e666*/

   downgain = 100.0 * downgain;
   upgain = 100.0 * upgain;

   /* in case a child is infeasible and therefore cutoff we take the gain of the other child to receive a somewhat
    * realistic gain for the infeasible child;
    * if both children are infeasible we just reset the initial zero values again
    */
   if( downbranchingresult->cutoff )
      downgain = 2 * upgain;
   if( upbranchingresult->cutoff )
      upgain = 2 * downgain;

   dbscore = SCIPgetBranchScore(scip, branchvar, downgain, upgain);

   score = SCIPgetBranchScore(scip, branchvar, lpscore, dbscore);

   return score;
}

/** calculates the score based on the down and up branching scores */
static
SCIP_Real calculateScoreFromDeeperscore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< variable to get the score for */
   BRANCHINGRESULTDATA*  downbranchingresult,/**< branching result of the down branch */
   BRANCHINGRESULTDATA*  upbranchingresult   /**< branching result of the up branch */
   )
{
   SCIP_Real score;
   SCIP_Real downscore;
   SCIP_Real upscore;

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   assert(downbranchingresult->deeperscore >= -0.2 || downbranchingresult->cutoff || SCIPisStopped(scip));
   assert(upbranchingresult->deeperscore >= -0.2 || upbranchingresult->cutoff || SCIPisStopped(scip));

   downscore = sqrt(downbranchingresult->deeperscore);
   upscore = sqrt(upbranchingresult->deeperscore);

   downscore = MAX(downscore, SCIPsumepsilon(scip)); /*lint !e666*/
   upscore = MAX(upscore, SCIPsumepsilon(scip)); /*lint !e666*/

   if( downbranchingresult->cutoff )
      downscore = 2 * upscore;
   if( upbranchingresult->cutoff )
      upscore = 2 * downscore;

   score = SCIPgetBranchScore(scip, branchvar, downscore, upscore);

   return score;
}

/** calculates the score based on the down and up branching scores */
static
SCIP_Real calculateScoreFromDeeperscoreAndCutoffs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< variable to get the score for */
   BRANCHINGRESULTDATA*  downbranchingresult,/**< branching result of the down branch */
   BRANCHINGRESULTDATA*  upbranchingresult   /**< branching result of the up branch */
   )
{
   SCIP_Real score;
   SCIP_Real downscore;
   SCIP_Real upscore;
   SCIP_Real totaldowngains;
   SCIP_Real totalupgains;
   SCIP_Real nlowestlevelcutoffs;
   int ntotaldowngains;
   int ntotalupgains;

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   assert(downbranchingresult->deeperscore >= -0.2 || downbranchingresult->cutoff || SCIPisStopped(scip));
   assert(upbranchingresult->deeperscore >= -0.2 || upbranchingresult->cutoff || SCIPisStopped(scip));

   nlowestlevelcutoffs = (1.0 * downbranchingresult->ndeepestcutoffs + upbranchingresult->ndeepestcutoffs)/(MAX(1,downbranchingresult->ndeepestnodes + upbranchingresult->ndeepestnodes));
   totaldowngains = downbranchingresult->totalgains;
   totalupgains = upbranchingresult->totalgains;
   ntotaldowngains = MAX(1, downbranchingresult->ntotalgains);
   ntotalupgains = MAX(1, upbranchingresult->ntotalgains);

   downscore = sqrt(downbranchingresult->deeperscore);
   upscore = sqrt(upbranchingresult->deeperscore);

   downscore = MAX(downscore, SCIPsumepsilon(scip)); /*lint !e666*/
   upscore = MAX(upscore, SCIPsumepsilon(scip)); /*lint !e666*/

   if( downbranchingresult->cutoff )
      downscore = 2 * upscore;
   if( upbranchingresult->cutoff )
      upscore = 2 * downscore;

   score = SCIPgetBranchScore(scip, branchvar, downscore, upscore);

   downscore = sqrt(totaldowngains/ntotaldowngains);
   upscore = sqrt(totalupgains/ntotalupgains);

   downscore = MAX(downscore, SCIPsumepsilon(scip)); /*lint !e666*/
   upscore = MAX(upscore, SCIPsumepsilon(scip)); /*lint !e666*/

   score += SCIPgetBranchScore(scip, branchvar, downscore, upscore)*nlowestlevelcutoffs;

   return score;
}

/** calculates the combined gain, weighted with parameters given by the user configuration */
static
SCIP_Real calculateWeightedGain(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,             /**< LAB configuration */
   BRANCHINGRESULTDATA*  downbranchingresult,/**< branching result of the down branch */
   BRANCHINGRESULTDATA*  upbranchingresult,  /**< branching result of the up branch */
   SCIP_Real             lpobjval            /**< objective value to get difference to as gain */
   )
{
   SCIP_Real downgain = 0.0;
   SCIP_Real upgain = 0.0;

   assert(config != NULL);
   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   /* the gain is the difference of the dualbound of a child and the reference objective value;
    * by bounding it by zero we are safe from numerical troubles
    */
   if( !downbranchingresult->cutoff )
      downgain = MAX(0, downbranchingresult->dualbound - lpobjval);
   if( !upbranchingresult->cutoff )
      upgain = MAX(0, upbranchingresult->dualbound - lpobjval);

   if( config->scoringfunction == 's' )
   {
      if( downbranchingresult->cutoff )
         downgain = SCIPinfinity(scip);
      if( upbranchingresult->cutoff )
         upgain = SCIPinfinity(scip);
   }
   else
   {
      /* in case a child is infeasible and therefore cutoff we take the gain of the other child to receive a somewhat
       * realistic gain for the infeasible child;
       * if both children are infeasible we just reset the initial zero values again
       */
      if( downbranchingresult->cutoff )
         downgain = upgain;
      if( upbranchingresult->cutoff )
         upgain = downgain;
   }

   return config->minweight * MIN(downgain, upgain) + (1.0 - config->minweight) * MAX(downgain, upgain);
}

/** calculates the score as mentioned in the lookahead branching paper by Glankwamdee and Linderoth;
 *  their score scales the number of cutoffs on the last layer of a 2-level temporary branching tree with the average gain of
 *  every last level problem; together with the best gain for each branch of a variable we get the final score
 */
static
SCIP_Real calculateScaledCutoffScore(
   BRANCHINGRESULTDATA*  downbranchingresult,/**< branching result of the down branch */
   BRANCHINGRESULTDATA*  upbranchingresult   /**< branching result of the up branch */
   )
{
   SCIP_Real bestdowngain;
   SCIP_Real bestupgain;
   SCIP_Real totaldowngains;
   SCIP_Real totalupgains;
   int nlowestlevelcutoffs;
   int ntotaldowngains;
   int ntotalupgains;

   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   nlowestlevelcutoffs = downbranchingresult->ndeepestcutoffs + upbranchingresult->ndeepestcutoffs;
   bestdowngain = downbranchingresult->bestgain;
   bestupgain = upbranchingresult->bestgain;
   totaldowngains = downbranchingresult->totalgains;
   totalupgains = upbranchingresult->totalgains;
   ntotaldowngains = MAX(1, downbranchingresult->ntotalgains);
   ntotalupgains = MAX(1, upbranchingresult->ntotalgains);

   return bestdowngain + bestupgain + (totaldowngains/ntotaldowngains + totalupgains/ntotalupgains)*nlowestlevelcutoffs;
}

/** calculates the score as mentioned in the lookahead branching paper by Glankwamdee and Linderoth;
 *  their score scales the number of cutoffs on the last layer of a 2-level temporary branching tree with the average gain of
 *  every last level problem; together with the best gain for each branch of a variable we get the final score
 */
static
SCIP_Real calculateWeightedCutoffScore(
   CONFIGURATION*        config,             /**< LAB configuration */
   BRANCHINGRESULTDATA*  downbranchingresult,/**< branching result of the down branch */
   BRANCHINGRESULTDATA*  upbranchingresult   /**< branching result of the up branch */
   )
{
   SCIP_Real bestdowngain;
   SCIP_Real bestupgain;
   SCIP_Real totaldowngains;
   SCIP_Real totalupgains;
   SCIP_Real nlowestlevelcutoffs;
   int ntotaldowngains;
   int ntotalupgains;

   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   nlowestlevelcutoffs = (1.0 * downbranchingresult->ndeepestcutoffs + upbranchingresult->ndeepestcutoffs)/(downbranchingresult->ndeepestnodes + upbranchingresult->ndeepestnodes);
   bestdowngain = downbranchingresult->bestgain;
   bestupgain = upbranchingresult->bestgain;
   totaldowngains = downbranchingresult->totalgains;
   totalupgains = upbranchingresult->totalgains;
   ntotaldowngains = MAX(1, downbranchingresult->ntotalgains);
   ntotalupgains = MAX(1, upbranchingresult->ntotalgains);

   return config->minweight*MIN(bestdowngain, bestupgain) + (1.0 - config->minweight)*MAX(bestdowngain, bestupgain) + (totaldowngains/ntotaldowngains + totalupgains/ntotalupgains)*nlowestlevelcutoffs;
}

static
SCIP_Real calculateCutoffScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< variable to get the score for */
   BRANCHINGRESULTDATA*  downbranchingresult,/**< branching result of the down branch */
   BRANCHINGRESULTDATA*  upbranchingresult,  /**< branching result of the up branch */
   SCIP_Real             lpobjval            /**< objective value to get difference to as gain */
   )
{
   SCIP_Real score;
   SCIP_Real downgain = SCIPsumepsilon(scip);
   SCIP_Real upgain = SCIPsumepsilon(scip);
   SCIP_Real gap;
   int nlowestlevelcutoffs;

   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   nlowestlevelcutoffs = 0;

   /* the gain is the difference of the dualbound of a child and the reference objective value;
    * by bounding it by zero we are safe from numerical troubles
    */
   if( !downbranchingresult->cutoff )
   {
      nlowestlevelcutoffs += downbranchingresult->ndeepestcutoffs;
      downgain = MAX(downgain, downbranchingresult->dualbound - lpobjval);
   }
   if( !upbranchingresult->cutoff )
   {
      nlowestlevelcutoffs += upbranchingresult->ndeepestcutoffs;
      upgain = MAX(upgain, upbranchingresult->dualbound - lpobjval);
   }

   /* in case a child is infeasible and therefore cutoff we take the gain of the other child to receive a somewhat
    * realistic gain for the infeasible child;
    * if both children are infeasible we just reset the initial zero values again
    */
   if( downbranchingresult->cutoff )
   {
      nlowestlevelcutoffs += 2 * SCIPgetNPseudoBranchCands(scip);
      downgain = 2 * upgain;
   }
   if( upbranchingresult->cutoff )
   {
      nlowestlevelcutoffs += 2 * SCIPgetNPseudoBranchCands(scip);
      upgain = 2 * downgain;
   }

   gap = SCIPgetCutoffbound(scip) - lpobjval;

   downgain = downgain/gap;
   upgain = upgain/gap;

   score = 1.0 * nlowestlevelcutoffs + SCIPgetBranchScore(scip, branchvar, downgain, upgain);

   return score;
}

static
SCIP_Real calculateRelCutoffScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< variable to get the score for */
   BRANCHINGRESULTDATA*  downbranchingresult,/**< branching result of the down branch */
   BRANCHINGRESULTDATA*  upbranchingresult,  /**< branching result of the up branch */
   SCIP_Real             lpobjval            /**< objective value to get difference to as gain */
   )
{
   SCIP_Real score;
   SCIP_Real downgain = SCIPsumepsilon(scip);
   SCIP_Real upgain = SCIPsumepsilon(scip);
   SCIP_Real gap;
   int factor;
   SCIP_Real nlowestlevelcutoffs;

   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   assert(downbranchingresult->ndeepestnodes + upbranchingresult->ndeepestnodes > 0 || (downbranchingresult->cutoff && upbranchingresult->cutoff));

   nlowestlevelcutoffs = (1.0 * downbranchingresult->ndeepestcutoffs + upbranchingresult->ndeepestcutoffs)/(1 + downbranchingresult->ndeepestnodes + upbranchingresult->ndeepestnodes);

   factor = SCIPgetNPseudoBranchCands(scip);
   if( factor > SCIPgetNLPRows(scip) )
      factor = SCIPgetNLPRows(scip);
   factor = factor * factor;

   /* the gain is the difference of the dualbound of a child and the reference objective value;
    * by bounding it by zero we are safe from numerical troubles
    */
   if( !downbranchingresult->cutoff )
   {
      downgain = MAX(downgain, downbranchingresult->dualbound - lpobjval);
   }
   if( !upbranchingresult->cutoff )
   {
      upgain = MAX(upgain, upbranchingresult->dualbound - lpobjval);
   }

   /* in case a child is infeasible and therefore cutoff we take the gain of the other child to receive a somewhat
    * realistic gain for the infeasible child;
    * if both children are infeasible we just reset the initial zero values again
    */
   if( downbranchingresult->cutoff )
   {
      downgain = 2 * upgain;
   }
   if( upbranchingresult->cutoff )
   {
      upgain = 2 * downgain;
   }

   gap = SCIPgetCutoffbound(scip) - lpobjval;

   downgain = downgain/gap;
   upgain = upgain/gap;

   score = factor * nlowestlevelcutoffs + SCIPgetBranchScore(scip, branchvar, downgain, upgain);

   return score;
}

/** scoring method that selects an actual scoring method based on the user configuration */
static
SCIP_Real calculateScore(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,             /**< LAB configuration */
   SCIP_VAR*             branchvar,          /**< variable to get the score for */
   BRANCHINGRESULTDATA*  downbranchingresult,/**< branching result of the down branch */
   BRANCHINGRESULTDATA*  upbranchingresult,  /**< branching result of the up branch */
   SCIP_Real             lpobjval,           /**< objective value to get difference to as gain */
   SCIP_Real             baselpobjval        /**< base objective value to get difference to as gain */
   )
{
   SCIP_Real score;
   char scoringfunction;

   assert(scip != NULL);
   assert(config != NULL);
   assert(branchvar != NULL);
   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   if( config->inscoring )
      scoringfunction = config->scoringscoringfunction;
   else if( SCIPgetProbingDepth(scip) > 0 )
      scoringfunction = config->deeperscoringfunction;
   else
      scoringfunction = config->scoringfunction;

   switch( scoringfunction )
   {
   case 's':
      score = calculateScaledCutoffScore(downbranchingresult, upbranchingresult);
      break;
   case 'w':
      score = calculateWeightedCutoffScore(config, downbranchingresult, upbranchingresult);
      break;
   case 'f':
      score = calculateWeightedGain(scip, config, downbranchingresult, upbranchingresult, baselpobjval);
      break;
   case 'p':
      score = calculateScoreFromDeeperscore(scip, branchvar, downbranchingresult, upbranchingresult);
      break;
   case 'a':
      score = calculateScoreFromDeeperscoreAndCutoffs(scip, branchvar, downbranchingresult, upbranchingresult);
      break;
   case 'l':
      score = calculateScoreFromResult2(scip, branchvar, downbranchingresult, upbranchingresult, lpobjval);
      break;
   case 'c':
      score = calculateCutoffScore(scip, branchvar, downbranchingresult, upbranchingresult, lpobjval);
      break;
   case 'r':
      score = calculateRelCutoffScore(scip, branchvar, downbranchingresult, upbranchingresult, lpobjval);
      break;
   case 'x':
      score = calculateScoreFromResult(scip, branchvar, downbranchingresult, upbranchingresult, baselpobjval);
      break;
   default:
      assert(scoringfunction == 'd');
      score = calculateScoreFromResult(scip, branchvar, downbranchingresult, upbranchingresult, lpobjval);
   }

   return score;
}

/** calculates the score based on the pseudocosts of the given variable */
static
SCIP_Real calculateScoreFromPseudocosts(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATE*            lpcand              /**< candidate to get the score for */
   )
{
   SCIP_Real downpseudocost;
   SCIP_Real uppseudocost;
   SCIP_Real score;

   assert(scip != NULL);
   assert(lpcand != NULL);

   downpseudocost = SCIPgetVarPseudocostVal(scip, lpcand->branchvar, 0-lpcand->fracval);
   uppseudocost = SCIPgetVarPseudocostVal(scip, lpcand->branchvar, 1-lpcand->fracval);

   score = SCIPgetBranchScore(scip, lpcand->branchvar, downpseudocost, uppseudocost);

   return score;
}

#ifdef SCIP_DEBUG
/** prints the names of the candidates of the given candidate list with their corresponding scores */
static
void printCandidateList(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATELIST*        candidatelist,      /**< list to be printed */
   SCORECONTAINER*       scorecontainer      /**< container with all scores */
   )
{
   int i;

   assert(scip != NULL);
   assert(candidatelist != NULL);
   assert(scorecontainer != NULL);

   for( i = 0; i < candidatelist->ncandidates; i++ )
   {
      SCIP_VAR* var = candidatelist->candidates[i]->branchvar;
      SCIP_Real score = scorecontainer->scores[SCIPvarGetProbindex(var)];

      assert(var != NULL);

      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, " Index %2i: Var %s Score %.9g\n", i, SCIPvarGetName(var), score);
   }
}
#endif

/** sorts the best candidates (w.r.t. the score in the container) of the candidate list to the front of the list */
static
void sortFirstCandidatesByScore(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATELIST*        candidatelist,      /**< candidates to be sorted */
   SCORECONTAINER*       scorecontainer,     /**< container with the scores for each candidate */
   int                   nbestcandidates     /**< number of candidates that should be kept sorted at the start of the list*/
   )
{
   int i;

   assert(scip != NULL);
   assert(candidatelist != NULL);
   assert(scorecontainer != NULL);
   assert(candidatelist->ncandidates > 0);
   assert(nbestcandidates <= candidatelist->ncandidates);

   for( i = 1; i < candidatelist->ncandidates; i++ )
   {
      CANDIDATE* movecand = candidatelist->candidates[i];
      int moveprobindex;
      SCIP_Real movescore;
      int nsorted;
      int insertionindex;
      assert(movecand != NULL);

      moveprobindex = SCIPvarGetProbindex(movecand->branchvar);
      movescore = scorecontainer->scores[moveprobindex];

      /* the length of the sorted portion of the array, starting at 0 */
      nsorted = MIN(i, nbestcandidates);

      insertionindex = findInsertionPoint(scip, scorecontainer, movescore, candidatelist->candidates, nsorted);

      assert(insertionindex <= nsorted);

      /* if no change has to be made, skip the reordering;
       * if the insertionindex lies after the sorted block, skip the reordering
       */
      if( insertionindex != i && insertionindex < nsorted )
      {
         int j;
         CANDIDATE* reordercand = movecand;

         /* move everything inside the sorted block one place further */
         for( j = insertionindex; j < nsorted; j++ )
         {
            CANDIDATE* oldcand = candidatelist->candidates[j];
            assert(oldcand != NULL);

            candidatelist->candidates[j] = reordercand;
            reordercand = oldcand;
         }
         /* the dropped element gets placed in the position of the actually moved element */
         candidatelist->candidates[i] = reordercand;
      }
   }

#ifdef SCIP_DEBUG
   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "All %i candidates, with the first %i candidates sorted by their FSB score:"
      "\n", candidatelist->ncandidates, nbestcandidates);
   printCandidateList(scip, candidatelist, scorecontainer);
#endif
}

/** checks whether the given candidates is reliable, so that its pseudocosts may be used */
static
SCIP_Bool isCandidateReliable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar           /**< var to check for reliability */
   )
{
   SCIP_Real size;
   SCIP_Real downsize;
   SCIP_Real upsize;
   SCIP_Real reliable = 5;

   assert(scip != NULL);
   assert(branchvar != NULL);

   downsize = SCIPgetVarPseudocostCountCurrentRun(scip, branchvar, SCIP_BRANCHDIR_DOWNWARDS);
   upsize = SCIPgetVarPseudocostCountCurrentRun(scip, branchvar, SCIP_BRANCHDIR_UPWARDS);
   size = MIN(downsize, upsize);

   return size >= reliable;
}

/** checks whether the current problem is feasible or cutoff */
static
SCIP_Bool isCurrentNodeCutoff(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return (SCIPgetCutoffdepth(scip) <= SCIPgetDepth(scip));
}

/** Ensures that the scores are present in the scorecontainer for each of the candidates to consider */
static
SCIP_RETCODE ensureScoresPresent(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS*               status,             /**< current status */
   PERSISTENTDATA*       persistent,         /**< container to store data over multiple calls to the branching rule; or NULL */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   SCIP_SOL*             baselpsol,          /**< base lp solution */
   DOMAINREDUCTIONS*     domainreductions,   /**< container collecting all domain reductions found; or NULL */
   BINCONSDATA*          binconsdata,        /**< container collecting all binary constraints; or NULL */
   CANDIDATELIST*        allcandidates,      /**< list containing all candidates to consider */
   BRANCHINGDECISION*    decision,           /**< struct to store the final decision */
   SCORECONTAINER*       scorecontainer,     /**< container to retrieve already calculated scores; or NULL */
   LEVEL2DATA*           level2data,         /**< level 2 LP results data */
   SCIP_Real             lpobjval            /**< base LP objective value */
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics          /**< general statistical data */
   ,LOCALSTATISTICS*     localstats          /**< local statistics, may be disregarded */
#endif
   )
{
   int i;
   int nunscoredcandidates = 0;
   int* candidateunscored;

   assert(scip != NULL);
   assert(config != NULL);
   assert(status != NULL);
   assert(allcandidates != NULL);
   assert(scorecontainer != NULL);
   assert(allcandidates->candidates != NULL || allcandidates->ncandidates == 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &candidateunscored, allcandidates->ncandidates) );

   /* filter the candidates based on the presence of a score in the 'scorecontainer'. Only those without a score need a
    * new one.
    */
   for( i = 0; i < allcandidates->ncandidates; i++ )
   {
      CANDIDATE* lpcand = allcandidates->candidates[i];
      SCIP_VAR* branchvar = lpcand->branchvar;
      int probindex = SCIPvarGetProbindex(branchvar);
      SCIP_Real knownscore = scorecontainer->scores[probindex];

      assert(lpcand != NULL);
      assert(branchvar != NULL);

      if( SCIPisLT(scip, knownscore, 0.0) )
      {
         if( config->abbrevpseudo && isCandidateReliable(scip, branchvar) )
         {
            SCIP_Real score = calculateScoreFromPseudocosts(scip, lpcand);
            SCIP_CALL( scoreContainerSetScore(scip, scorecontainer, lpcand, score, 0.0, 0.0) );
         }
         else if( config->level2avgscore && SCIPgetProbingDepth(scip) > 0 )
         {
            assert(scorecontainer->nsetscores > 0);
            SCIP_CALL( scoreContainerSetScore(scip, scorecontainer, lpcand,
                  scorecontainer->scoresum / scorecontainer->nsetscores, 0.0, 0.0) );
         }
         else if( config->level2zeroscore && SCIPgetProbingDepth(scip) > 0 )
         {
            assert(scorecontainer->nsetscores > 0);
            SCIP_CALL( scoreContainerSetScore(scip, scorecontainer, lpcand,
                  -0.1, 0.0, 0.0) );
         }
         else
         {
            /* score is unknown and needs to be calculated */
            candidateunscored[nunscoredcandidates] = i;
            nunscoredcandidates++;
         }
      }
   }

   if( nunscoredcandidates > 0 )
   {
      CANDIDATELIST* unscoredcandidates;

      /* allocate the list of candidates without any score (gets updated further on) */
      SCIP_CALL( candidateListCreate(scip, &unscoredcandidates, nunscoredcandidates) );

      /* move the unscored candidates to the temp list */
      for( i = 0; i < nunscoredcandidates; i++ )
      {
         int candindex = candidateunscored[i];

         assert(allcandidates->candidates[candindex] != NULL);

         unscoredcandidates->candidates[i] = allcandidates->candidates[candindex];
      }

#ifdef SCIP_DEBUG
      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Of the given %i candidates, %i have no score: ",
         allcandidates->ncandidates, nunscoredcandidates);
      printCandidates(scip, SCIP_VERBLEVEL_HIGH, unscoredcandidates);
      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Calculating the FSB result to get a score for the remaining "
         "candidates.\n");
#endif

      /* Calculate all remaining FSB scores and collect the scores in the container */;
#ifdef SCIP_STATISTIC
      SCIP_CALL( getFSBResult(scip, status, persistent, config, baselpsol, domainreductions, binconsdata, unscoredcandidates,
            decision, scorecontainer, level2data, lpobjval, statistics, localstats) );
#else
      SCIP_CALL( getFSBResult(scip, status, persistent, config, baselpsol, domainreductions, binconsdata, unscoredcandidates,
            decision, scorecontainer, level2data, lpobjval) );
#endif

      /* move the now scored candidates back to the original list */
      for( i = 0; i < nunscoredcandidates; i++ )
      {
         assert(allcandidates->candidates[candidateunscored[i]] == unscoredcandidates->candidates[i]);

         assert(unscoredcandidates->candidates[i] != NULL);
         unscoredcandidates->candidates[i] = NULL;
      }

      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Calculated the scores for the remaining candidates\n");

      SCIP_CALL( candidateListFree(scip, &unscoredcandidates) );
   }

   /* reset the best sorted indices, as those are only valid on the FSB run already completed */
   scoreContainterResetBestSortedCands(scorecontainer);

   SCIPfreeBufferArray(scip, &candidateunscored);

   return SCIP_OKAY;
}

/** Get the candidates to temporarily branch on. In the LAB case this is the complete list of possible candidates. In the
 *  ALAB case only the 'best' candidates are returned. */
static
SCIP_RETCODE filterCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS*               status,             /**< current status */
   PERSISTENTDATA*       persistent,         /**< container to store data over multiple calls to the branching rule; or NULL */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   SCIP_SOL*             baselpsol,          /**< base lp solution */
   DOMAINREDUCTIONS*     domainreductions,   /**< container collecting all domain reductions found; or NULL */
   BINCONSDATA*          binconsdata,        /**< container collecting all binary constraints; or NULL */
   CANDIDATELIST*        candidatelist,      /**< list of candidates to branch on */
   BRANCHINGDECISION*    decision,           /**< struct to store the final decision */
   SCORECONTAINER*       scorecontainer,     /**< container to retrieve already calculated scores; or NULL */
   LEVEL2DATA*           level2data,         /**< level 2 LP results data */
   SCIP_Real             lpobjval            /**< base LP objective value */
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics          /**< general statistical data */
   ,LOCALSTATISTICS*     localstats          /**< local statistics, may be disregarded */
#endif
   )
{
   assert(scip != NULL);
   assert(config != NULL);
   assert(status != NULL);
   assert(candidatelist != NULL);
   assert(SCIPinProbing(scip));

   /* abbreviated LAB: only use the "best" candidates */
   if( config->abbreviated )
   {
      assert(scorecontainer != NULL);

      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Getting the best (at most) %i of the given %i candidates: ",
         config->maxncands, candidatelist->ncandidates);
#ifdef SCIP_DEBUG
      printCandidates(scip, SCIP_VERBLEVEL_HIGH, candidatelist);
#endif

      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "%s", "Ensuring that all candidates have a score.\n");
#ifdef SCIP_STATISTIC
      SCIP_CALL( ensureScoresPresent(scip, status, persistent, config, baselpsol, domainreductions, binconsdata, candidatelist,
            decision, scorecontainer, level2data, lpobjval, statistics, localstats) );
#else
      SCIP_CALL( ensureScoresPresent(scip, status, persistent, config, baselpsol, domainreductions, binconsdata, candidatelist,
            decision, scorecontainer, level2data, lpobjval) );
#endif

      /* if we didn't find any domreds or constraints during the FSB scoring, we branch on */
      if( isBranchFurther(status, SCIPgetProbingDepth(scip) == 0) )
      {
         int nusedcands;
         int i;

         if( SCIPgetProbingDepth(scip) == 0 || config->maxndeepercands == 0 )
            nusedcands = MIN(config->maxncands, candidatelist->ncandidates);
         else
            nusedcands = MIN(config->maxndeepercands, candidatelist->ncandidates);

         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "%s", "Filter the candidates by their score.\n");

         sortFirstCandidatesByScore(scip, candidatelist, scorecontainer, nusedcands);

         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Best candidate according to FSB scores: <%s>\n",
            SCIPvarGetName(candidatelist->candidates[0]->branchvar));

         if( config->worsefactor >= 0 )
         {
            for( i = 1; i < nusedcands; ++i )
            {
               if( scorecontainer->scores[SCIPvarGetProbindex(candidatelist->candidates[0]->branchvar)] >
                  config->worsefactor * scorecontainer->scores[SCIPvarGetProbindex(candidatelist->candidates[i]->branchvar)] )
                  break;
            }
            nusedcands = i;
         }

         if( config->filterbymaxgain && SCIPgetProbingDepth(scip) == 0 )
         {
            SCIP_Real maxgain;
            SCIP_Real bestmaxgain = MAX(scorecontainer->downgains[SCIPvarGetProbindex(candidatelist->candidates[0]->branchvar)],
               scorecontainer->upgains[SCIPvarGetProbindex(candidatelist->candidates[0]->branchvar)]); /*lint !e666*/

            if( bestmaxgain == 0.0 )
               nusedcands = 1;
            else
            {
               for( i = nusedcands - 1; i >= 1; --i )
               {
                  maxgain = MAX(scorecontainer->downgains[SCIPvarGetProbindex(candidatelist->candidates[i]->branchvar)],
                     scorecontainer->upgains[SCIPvarGetProbindex(candidatelist->candidates[i]->branchvar)]); /*lint !e666*/

                  if( SCIPisSumLE(scip, maxgain / bestmaxgain, 1.0) )
                  {
                     --nusedcands;

                     if( i < nusedcands )
                     {
                        CANDIDATE* tmp = candidatelist->candidates[i];
                        candidatelist->candidates[i] = candidatelist->candidates[nusedcands];
                        candidatelist->candidates[nusedcands] = tmp;
                     }
                  }
               }
            }
         }

         if( SCIPgetProbingDepth(scip) > 0 && scorecontainer->scores[SCIPvarGetProbindex(candidatelist->candidates[0]->branchvar)] > -0.05)
         {
            for( i = 1; i < nusedcands; ++i )
            {
               if( scorecontainer->scores[SCIPvarGetProbindex(candidatelist->candidates[i]->branchvar)] < -0.05 )
                  break;
            }
            nusedcands = i;
         }

         SCIP_CALL( candidateListKeep(scip, candidatelist, nusedcands) );
      }
#ifdef SCIP_DEBUG
      else
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Strong Branching would have stopped.\n");
      }
#endif

      if( isCurrentNodeCutoff(scip) )
         status->cutoff = TRUE;
   }
   else
   {
      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Getting the branching candidates by selecting all candidates.\n");
   }

   return SCIP_OKAY;
}

/** Executes the general branching on a variable in a given direction (up/down) and repeats the algorithm from the new node */
static
SCIP_RETCODE executeBranchingRecursive(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS*               status,             /**< current status */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   SCIP_SOL*             baselpsol,          /**< the base lp solution */
   CANDIDATE*            candidate,          /**< candidate to branch on */
   SCIP_Real             localbaselpsolval,  /**< the objective value of the current temporary problem */
   SCIP_Real             baselpobjval,       /**< LP objective value of focus node (not probing) */
   int                   recursiondepth,     /**< remaining recursion depth */
   DOMAINREDUCTIONS*     domainreductions,   /**< container collecting all domain reductions found; or NULL */
   BINCONSDATA*          binconsdata,        /**< container collecting all binary constraints; or NULL */
   LEVEL2DATA*           level2data,         /**< level 2 LP results data */
   BRANCHINGRESULTDATA*  branchingresult,    /**< container to store the result of the branching in */
   SCORECONTAINER*       scorecontainer,     /**< container to retrieve already calculated scores; or NULL */
   SCIP_Bool             downbranching       /**< should we branch up or down in here? */
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics          /**< general statistical data */
   ,LOCALSTATISTICS*     localstats          /**< local statistics, may be disregarded */
#endif
   )
{
   int probingdepth;
   SCIP_VAR* branchvar;
   SCIP_Real branchvalfrac;
   SCIP_Real branchval;
   SCIP_Bool varisbinary;
   SCIP_Bool solvedlp = TRUE;

   assert(scip != NULL);
   assert(status != NULL);
   assert(config != NULL);
   assert(candidate != NULL);
   assert(branchingresult != NULL);

   branchvar = candidate->branchvar;
   branchvalfrac = candidate->fracval;
   branchval = candidate->branchval;

   assert(branchvar != NULL);

   probingdepth = SCIPgetProbingDepth(scip);
   varisbinary = SCIPvarIsBinary(branchvar);

   if( binconsdata != NULL && varisbinary )
   {
      if( downbranching )
      {
         /* In case that the branch variable is binary, add the negated var to the list.
          * This list is used to generate a set packing constraint for cutoff branches which were reached by only using
          * binary variables.
          * DownBranching on a binary variable x means: x <= 0
          * When this cutoff occurs we have that: x >= 1 <=> 1-x <= 0
          */
         SCIP_VAR* negbranchvar;

         SCIP_CALL( SCIPgetNegatedVar(scip, branchvar, &negbranchvar) );

         assert(negbranchvar != NULL);

         binaryVarListAppend(scip, binconsdata->binaryvars, negbranchvar);
      }
      else
      {
         /* In case that the branch variable is binary, add the var to the list.
          * This list is used to generate a set packing constraint for cutoff branches which were reached by only using
          * binary variables.
          * UpBranching on a binary variable x means: x >= 1
          * When this cutoff occurs we have that: x <= 0
          */
         binaryVarListAppend(scip, binconsdata->binaryvars, branchvar);
      }
   }

   if( level2data != NULL )
   {
      SCIP_Real newbound = downbranching ? SCIPfeasFloor(scip, branchval) : SCIPfeasCeil(scip, branchval);

      if( SCIPgetProbingDepth(scip) == 0 )
      {
         assert(SCIPvarGetProbindex(branchvar) >= 0);
         level2data->branchvar1 = (unsigned int) SCIPvarGetProbindex(branchvar);
         level2data->branchdir1 = !downbranching;
         level2data->branchval1 = newbound;
      }
      else
      {
         LEVEL2RESULT* result;

         assert(SCIPgetProbingDepth(scip) == 1);
         assert(SCIPvarGetProbindex(branchvar) >= 0);

         level2data->branchvar2 = (unsigned int) SCIPvarGetProbindex(branchvar);
         level2data->branchdir2 = !downbranching;
         level2data->branchval2 = newbound;

         SCIP_CALL( level2dataGetResult(scip, level2data, &result) );

         /* we already processed a similar level 2 node */
         if( result != NULL )
         {
            solvedlp = FALSE;
#ifdef SCIP_STATISTIC
            statistics->nduplicatelps[probingdepth]++;
#endif
            branchingresult->objval = result->lpobjval;
            branchingresult->dualbound = result->lpobjval;
            branchingresult->dualboundvalid = result->valid;
            branchingresult->cutoff = result->cutoff;
            branchingresult->niterations = 0;

            if( !branchingresult->cutoff && branchingresult->dualboundvalid
               && SCIPisGE(scip, branchingresult->objval, SCIPgetCutoffbound(scip)) )
               branchingresult->cutoff = TRUE;

            LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH,
               "Use old %s branching result on var <%s> with 'val > %g' and bounds [<%g>..<%g>]: objval <%.9g>, cutoff <%d> "
               "(the parent objval was <%.9g>)\n",
               downbranching ? "down" : "up", SCIPvarGetName(branchvar), branchval, SCIPvarGetLbLocal(branchvar),
               SCIPvarGetUbLocal(branchvar), branchingresult->objval, branchingresult->cutoff, localbaselpsolval);
         }
      }
   }

   if( solvedlp )
   {
      SCIP_CALL( executeBranching(scip, config, downbranching, candidate, branchingresult, baselpsol, domainreductions,
            status) );

      assert(SCIPgetProbingDepth(scip) == 1 || SCIPgetProbingDepth(scip) == 2);

      if( level2data != NULL && SCIPgetProbingDepth(scip) == 2)
      {
         SCIP_Bool duplicate;

         SCIP_CALL( level2dataStoreResult(scip, level2data, branchingresult->objval, branchingresult->cutoff, branchingresult->dualboundvalid, &duplicate) );
         assert(!duplicate);
      }

#ifdef SCIP_STATISTIC
      statistics->nlpssolved[probingdepth]++;
      statistics->nlpiterations[probingdepth] += branchingresult->niterations;

      if( config->inscoring )
      {
         statistics->nlpssolvedfsb[probingdepth]++;
         statistics->nlpiterationsfsb[probingdepth] += branchingresult->niterations;
      }
#endif
      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Solving the LP took %" SCIP_LONGINT_FORMAT " iterations (status %d).\n",
         branchingresult->niterations, SCIPgetLPSolstat(scip));

#ifdef SCIP_DEBUG
      if( status->lperror )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The LP could not be solved.\n");
      }
      else if( branchingresult->cutoff )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The solved LP was infeasible and as such is cutoff\n");
      }
      else
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The solved LP was feasible and has an objval <%.9g> (the parent objval was "
            "<%.9g>)\n", branchingresult->objval, localbaselpsolval);
      }
#endif
   }

   if( !branchingresult->cutoff && !status->lperror && !status->limitreached )
   {
      SCIP_Real localgain;

      localgain = MAX(0, branchingresult->objval - localbaselpsolval);

      /* update pseudo costs */
      if( downbranching )
      {
         SCIP_CALL( SCIPupdateVarPseudocost(scip, branchvar, 0.0 - branchvalfrac, localgain, 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPupdateVarPseudocost(scip, branchvar, 1.0 - branchvalfrac, localgain, 1.0) );
      }
   }

   if( solvedlp && !branchingresult->cutoff && !status->lperror && !status->limitreached )
   {
      /* store the warm start information in the candidate, so that it can be reused in a later branching */
      if( config->reusebasis && config->inscoring )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Storing warm start information for %s branching on var <%s>\n",
            downbranching ? "down" : "up", SCIPvarGetName(branchvar));

         SCIP_CALL( candidateStoreWarmStartInfo(scip, candidate, downbranching) );
      }

      if( recursiondepth > 1 && !config->inscoring )
      {
         CANDIDATELIST* candidatelist;

         SCIP_CALL( candidateListGetAllFractionalCandidates(scip, &candidatelist) );
         assert(candidatelist != NULL);

         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "%sbranching has <%i> candidates.\n", downbranching ? "Down" : "Up",
            candidatelist->ncandidates);

         if( candidatelist->ncandidates > 0 )
         {
            BRANCHINGDECISION* deeperdecision;
            STATUS* deeperstatus;
            PERSISTENTDATA* deeperpersistent = NULL;
            SCIP_Real deeperlpobjval = branchingresult->objval;
#ifdef SCIP_STATISTIC
            LOCALSTATISTICS* deeperlocalstats;

            SCIP_CALL( localStatisticsAllocate(scip, &deeperlocalstats) );
#endif
            SCIP_CALL( statusCreate(scip, &deeperstatus) );

            SCIP_CALL( branchingDecisionCreate(scip, &deeperdecision) );

#ifdef SCIP_STATISTIC
            SCIP_CALL( filterCandidates(scip, deeperstatus, deeperpersistent, config, baselpsol, domainreductions, binconsdata, candidatelist,
               deeperdecision, scorecontainer, level2data, deeperlpobjval,
               statistics, localstats) );
#else
            SCIP_CALL( filterCandidates(scip, deeperstatus, deeperpersistent, config, baselpsol, domainreductions, binconsdata, candidatelist,
               deeperdecision, scorecontainer, level2data, deeperlpobjval) );
#endif
            if( deeperstatus->lperror )
            {
               LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "ignoring lperror in filtering call...\n");
               deeperstatus->lperror = FALSE;
            }
            if( deeperstatus->cutoff )
            {
               branchingresult->ndeepestnodes += 2;
               branchingresult->ndeepestcutoffs += 2;
            }

            /* the status may have changed because of FSB to get the best candidates */
            if( isBranchFurther(deeperstatus, FALSE) )
            {
               LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Now the objval is <%.9g>\n", branchingresult->objval);

#ifdef SCIP_STATISTIC
               deeperlocalstats->ncutoffproofnodes = 0;
               SCIP_CALL( selectVarRecursive(scip, deeperstatus, deeperpersistent, config, baselpsol, domainreductions,
                     binconsdata, candidatelist, deeperdecision, scorecontainer, level2data, recursiondepth - 1,
                     deeperlpobjval, baselpobjval, &branchingresult->niterations, &branchingresult->ndeepestcutoffs,
                     &branchingresult->bestgain, &branchingresult->totalgains, &branchingresult->ntotalgains,
                     &branchingresult->ndeepestnodes,
                     statistics, deeperlocalstats, NULL, NULL) );
#else
               SCIP_CALL( selectVarRecursive(scip, deeperstatus, deeperpersistent, config, baselpsol, domainreductions,
                     binconsdata, candidatelist, deeperdecision, scorecontainer, level2data, recursiondepth - 1,
                     deeperlpobjval, baselpobjval, &branchingresult->niterations, &branchingresult->ndeepestcutoffs,
                     &branchingresult->bestgain, &branchingresult->totalgains, &branchingresult->ntotalgains,
                     &branchingresult->ndeepestnodes) );
#endif

               assert(deeperstatus->cutoff || deeperstatus->domred || deeperstatus->lperror
                  || branchingresult->ndeepestnodes == 8
                  || branchingresult->ndeepestnodes == 2 * candidatelist->ncandidates
                  || SCIPisStopped(scip));

               /* the proved dual bound of the deeper branching cannot be less than the current dual bound, as every deeper
                * node has more/tighter constraints and as such cannot be better than the base LP. */
               assert(SCIPisGE(scip, deeperdecision->proveddb, branchingresult->dualbound));
               branchingresult->dualbound = deeperdecision->proveddb;
               branchingresult->deeperscore = deeperdecision->score;
               branchingresult->dualboundvalid = TRUE;
            }
#ifdef SCIP_STATISTIC
            else
            {
               assert(SCIPgetProbingDepth(scip) == probingdepth + 1);

               statistics->stopafterfsb[probingdepth+1]++;

               if( deeperstatus->cutoff )
               {
                  statistics->cutoffafterfsb[probingdepth+1]++;
               }
               else if( deeperstatus->domred )
               {
                  statistics->domredafterfsb[probingdepth+1]++;
               }
            }
#endif
            /* deeperstatus->cutoff is TRUE, if any up/down child pair of the up child were cutoff */
            if( deeperstatus->cutoff )
            {
               branchingresult->cutoff = TRUE;
#ifdef SCIP_STATISTIC
               localstats->ncutoffproofnodes += deeperlocalstats->ncutoffproofnodes;
#endif
               LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Both deeper children were cutoff, so the %s branch is "
                  "cutoff\n", downbranching ? "down" : "up");
            }

            branchingDecisionFree(scip, &deeperdecision);
            statusFree(scip, &deeperstatus);
#ifdef SCIP_STATISTIC
            localStatisticsFree(scip, &deeperlocalstats);
#endif
         }
         else
         {
            branchingresult->deeperscore = (branchingresult->dualbound - baselpobjval) * (branchingresult->dualbound - baselpobjval) * 10;
         }
         SCIP_CALL( candidateListFree(scip, &candidatelist) );
      }
   }

   if( recursiondepth == 1 && !config->inscoring )
   {
      branchingresult->ndeepestnodes++;
      /* this is a cutoff on the lowest tree level */
      if( branchingresult->cutoff )
      {
         branchingresult->ndeepestcutoffs++;
      }
   }

   if( binconsdata != NULL && varisbinary )
   {
      /* the current branching child is infeasible and we only branched on binary variables in lookahead branching */
      if( solvedlp && branchingresult->cutoff && !status->lperror && SCIPallColsInLP(scip)
         && binconsdata->binaryvars->nbinaryvars == (probingdepth + 1) )
      {
#ifdef SCIP_STATISTIC
         SCIP_CALL( addBinaryConstraint(scip, config, binconsdata, baselpsol, statistics) );
#else
         SCIP_CALL( addBinaryConstraint(scip, config, binconsdata, baselpsol) );
#endif
      }

      binaryVarListDrop(binconsdata->binaryvars);
   }

   /* reset the probing depth to undo the previous branching */
   SCIP_CALL( SCIPbacktrackProbing(scip, probingdepth) );

   return SCIP_OKAY;
}

/** branches recursively on all given candidates */
static
SCIP_RETCODE selectVarRecursive(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS*               status,             /**< current status */
   PERSISTENTDATA*       persistent,         /**< container to store data over multiple calls to the branching rule; or NULL */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   SCIP_SOL*             baselpsol,          /**< base lp solution */
   DOMAINREDUCTIONS*     domainreductions,   /**< container collecting all domain reductions found; or NULL */
   BINCONSDATA*          binconsdata,        /**< container collecting all binary constraints; or NULL */
   CANDIDATELIST*        candidatelist,      /**< list of candidates to branch on */
   BRANCHINGDECISION*    decision,           /**< struct to store the final decision */
   SCORECONTAINER*       scorecontainer,     /**< container to retrieve already calculated scores; or NULL */
   LEVEL2DATA*           level2data,         /**< level 2 LP results data */
   int                   recursiondepth,     /**< remaining recursion depth */
   SCIP_Real             lpobjval,           /**< LP objective value of current probing node*/
   SCIP_Real             baselpobjval,       /**< LP objective value of focus node (not probing) */
   SCIP_Longint*         niterations,        /**< pointer to store the total number of iterations for this variable; or NULL*/
   int*                  ndeepestcutoffs,    /**< pointer to store the total number of cutoffs on the deepest level; or NULL */
   SCIP_Real*            bestgain,           /**< pointer to store the best gain found with these candidates; or NULL */
   SCIP_Real*            totalgains,         /**< pointer to store the sum over all gains that are valid in both children;
                                              *   or NULL, if bestgain == NULL */
   int*                  ntotalgains,        /**< pointer to store the number of gains summed in totalgains;
                                              *   or NULL, if bestgain == NULL */
   int*                  ndeepestnodes       /**< pointer to store the number of nodes processed in the deepest level */
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics          /**< general statistical data */
   ,LOCALSTATISTICS*     localstats          /**< local statistics, may be disregarded */
   ,SCIP_Real*           firstscoreptr       /**< pointer to store score of first candidate, or NULL */
   ,SCIP_Real*           bestscoreptr        /**< pointer to store best score, or NULL */
#endif
   )
{
   BRANCHINGRESULTDATA* downbranchingresult = NULL;
   BRANCHINGRESULTDATA* upbranchingresult = NULL;
   BRANCHINGRESULTDATA* bestdownbranchingresult = NULL;
   BRANCHINGRESULTDATA* bestupbranchingresult = NULL;
   SCIP_LPI* lpi;
   SCIP_Real bestscore = -SCIPinfinity(scip);
   SCIP_Real bestscorelowerbound;
   SCIP_Real bestscoreupperbound;
   SCIP_Real bestscoringlpobjval = -SCIPinfinity(scip);
   int start = 0;
   int i;
   int c;
   int nlpcands;
   int probingdepth;
   SCIP_Bool stopafterinfeasible = FALSE;

   assert(scip != NULL);
   assert(status != NULL);
   assert(config != NULL);
   assert(!config->usedomainreduction || domainreductions != NULL);
   assert(candidatelist != NULL);
   assert(candidatelist->ncandidates > 0);
   assert(decision != NULL);
   assert(recursiondepth >= 1);
#ifdef SCIP_STATISTIC
   assert(statistics != NULL);

   if( firstscoreptr != NULL )
      *firstscoreptr = -1.0;
   if( bestscoreptr != NULL )
      *bestscoreptr = -1.0;
#endif

   nlpcands = candidatelist->ncandidates;
   probingdepth = SCIPgetProbingDepth(scip);
   assert(probingdepth >= 0 && probingdepth < config->recursiondepth);

   if( persistent != NULL && (!config->abbreviated || config->inscoring) && probingdepth == 0 )
      start = persistent->restartindex;

   /* init default decision */
   decision->branchvar = candidatelist->candidates[0]->branchvar;
   decision->branchval = candidatelist->candidates[0]->branchval;
   decision->downdb = lpobjval;
   decision->downdbvalid = TRUE;
   decision->updb = lpobjval;
   decision->updbvalid = TRUE;
   decision->proveddb = lpobjval;
   decision->score = 0.0;

   bestscorelowerbound = SCIPvarGetLbLocal(decision->branchvar);
   bestscoreupperbound = SCIPvarGetUbLocal(decision->branchvar);

   SCIP_CALL( branchingResultDataCreate(scip, &downbranchingresult) );
   SCIP_CALL( branchingResultDataCreate(scip, &upbranchingresult) );

   SCIP_CALL( branchingResultDataCreate(scip, &bestdownbranchingresult) );
   SCIP_CALL( branchingResultDataCreate(scip, &bestupbranchingresult) );

   assert(downbranchingresult != NULL);
   assert(upbranchingresult != NULL);

   if( config->inscoring )
   {
      SCIP_CALL( SCIPgetBoolParam(scip, "branching/forceallchildren", &stopafterinfeasible) );
      stopafterinfeasible = !stopafterinfeasible;
   }

   SCIP_CALL( SCIPgetLPI(scip, &lpi) );

#ifdef SCIP_DEBUG
   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Started selectVarRecursive with <%i> candidates: ", nlpcands);
   printCandidates(scip, SCIP_VERBLEVEL_HIGH, candidatelist);
#endif

   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Starting loop from index %d\n", start);

   /* iterate over all current branching candidates and evaluate their two potential child nodes by:
    * - potentially applying domain propagation at each node before
    * - solving the LP at the nodes to obtain a dual bound
    * - potentially evaluating branching candidates at the potential child node again by applying this method recursively
    *
    * some improvements of the general scheme:
    * - results obtained for a candidate in a previous lookahead branching call at this node may be re-used
    * - while i counts the number of candidates evaluated in this call, we do not always start at the front
    *   of the candidate array, but rather store at which index we stopped last time (e.g., because a domain reduction was
    *   found and applied) and start from that index next time. Even though the set of branching candidates is probably different
    *   it is often reasonably close and we avoid evaluating the same variables again and again.
    */
   for( i = 0, c = start;
        isBranchFurtherLoopDecrement(status, &c) && i < nlpcands && !SCIPisStopped(scip); i++, c++)
   {
      DOMAINREDUCTIONS* downdomainreductions = NULL;
      DOMAINREDUCTIONS* updomainreductions = NULL;
      SCIP_Bool useoldbranching = FALSE;
      SCIP_Real oldlpobjval = -SCIPinfinity(scip);
      CANDIDATE* candidate;
      SCIP_VAR* branchvar;
      SCIP_Real branchval;
      SCIP_Real branchlb;
      SCIP_Real branchub;

      c = c % nlpcands;

      candidate = candidatelist->candidates[c];

      assert(candidate != NULL);

      branchvar = candidate->branchvar;
      branchval = candidate->branchval;

      assert(branchvar != NULL);

      branchlb = SCIPvarGetLbLocal(branchvar);
      branchub = SCIPvarGetUbLocal(branchvar);

      if( SCIPisEQ(scip, branchlb, branchub) )
      {
         /* if both bounds are equal the variable is fixed and we cannot branch
          * this may happen if domain propagation on other candidates finds better bounds for the current candidate
          */
         status->domred = TRUE;
#ifdef SCIP_STATISTIC
         statistics->npropdomred[probingdepth]++;
#endif
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Domain Propagation changed the bounds of a branching candidate."
               "\n");
         continue;
      }

      /* @todo apply already found domainreductions for this candidate? */

#ifdef SCIP_STATISTIC
      /* Reset the cutoffproofnodes, as the number of proof nodes from previous branching vars (which where not
       * cutoff, as we didn't break the loop) is not relevant for the min total sum of proof nodes.
       */
      localstats->ncutoffproofnodes = 0;
#endif

      branchingResultDataInit(scip, downbranchingresult);
      branchingResultDataInit(scip, upbranchingresult);

      /* use old lookahead branching result, if last call on this variable is not too long ago */
      if( persistent != NULL && (config->inscoring || probingdepth == 0) && isUseOldBranching(scip, persistent, config, branchvar) )
      {
         SCIP_CALL( getOldBranching(scip, persistent, config, branchvar, downbranchingresult, upbranchingresult,
               &oldlpobjval) );
         useoldbranching = TRUE;
#ifdef SCIP_STATISTIC
         if( config->inscoring )
            statistics->noldbranchusedfsb[probingdepth]++;
         else
            statistics->noldbranchused[probingdepth]++;
#endif
      }
      else
      {
         SCIP_Bool down;
         int k;

         LABdebugMessage(scip, SCIP_VERBLEVEL_NORMAL, "Started branching on var <%s> with val <%g> and bounds "
            "[<%g>..<%g>]\n", SCIPvarGetName(branchvar), branchval, SCIPvarGetLbLocal(branchvar),
            SCIPvarGetUbLocal(branchvar));

         if( config->usedomainreduction )
         {
            SCIP_CALL( domainReductionsCreate(scip, &downdomainreductions) );
            SCIP_CALL( domainReductionsCreate(scip, &updomainreductions) );
         }

         down = SCIPisStrongbranchDownFirst(scip, branchvar);

         /* @todo break if result is infeasible (probably only in first layer)? */
         for( k = 0; k < 2; ++k )
         {
            DOMAINREDUCTIONS* localdomainreductions;
            BRANCHINGRESULTDATA* localbranchingresult;
            BRANCHINGRESULTDATA* otherbranchingresult;

            localdomainreductions = down ? downdomainreductions : updomainreductions;
            localbranchingresult = down ? downbranchingresult : upbranchingresult;
            otherbranchingresult = down ? upbranchingresult : downbranchingresult;

#ifdef SCIP_STATISTIC
            SCIP_CALL( executeBranchingRecursive(scip, status, config, baselpsol, candidate, lpobjval, baselpobjval,
                  recursiondepth, localdomainreductions, binconsdata, level2data, localbranchingresult, scorecontainer,
                  down, statistics, localstats) );
#else

            SCIP_CALL( executeBranchingRecursive(scip, status, config, baselpsol, candidate, lpobjval, baselpobjval,
                  recursiondepth, localdomainreductions, binconsdata, level2data, localbranchingresult, scorecontainer,
                  down) );
#endif

            /* check whether a new solutions rendered the previous child infeasible */
            if( SCIPallColsInLP(scip) && !otherbranchingresult->cutoff )
            {
               if( k == 1 && SCIPisGE(scip, otherbranchingresult->dualbound, SCIPgetCutoffbound(scip)) )
               {
                  otherbranchingresult->cutoff = TRUE;
                  LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH,
                     "The %s branching changed the cutoffbound and rendered the %s branching result infeasible.\n",
                     down ? "down" : "up", down ? "up" : "down");
               }
            }
            if( stopafterinfeasible &&  k == 0 && localbranchingresult->cutoff )
               break;

            /* the second iteration of the loop should branch in the other direction */
            down = !down;
         }

         LABdebugMessage(scip, SCIP_VERBLEVEL_NORMAL, "-> down=%.9g (gain=%.9g, valid=%u, inf=%u), up=%.9g "
               "(gain=%.9g, valid=%u, inf=%u)\n", downbranchingresult->dualbound,
               downbranchingresult->dualbound - lpobjval, downbranchingresult->dualboundvalid,
               downbranchingresult->cutoff, upbranchingresult->dualbound, upbranchingresult->dualbound - lpobjval,
               upbranchingresult->dualboundvalid, upbranchingresult->cutoff);

         if( niterations != NULL )
            *niterations += downbranchingresult->niterations + upbranchingresult->niterations;

         /* store results of branching call */
         if( persistent != NULL && !upbranchingresult->cutoff && !downbranchingresult->cutoff && (config->inscoring || probingdepth == 0) )
         {
            SCIP_CALL( updateOldBranching(scip, persistent, config, branchvar, branchval, downbranchingresult,
                  upbranchingresult, lpobjval) );
         }
      }

      if( ndeepestcutoffs != NULL )
         *ndeepestcutoffs += downbranchingresult->ndeepestcutoffs + upbranchingresult->ndeepestcutoffs;

      if( ndeepestnodes != NULL )
         *ndeepestnodes += downbranchingresult->ndeepestnodes + upbranchingresult->ndeepestnodes;

      if( !status->lperror && !status->limitreached )
      {
         SCIP_Real scoringlpobjval = useoldbranching ? oldlpobjval : lpobjval;
         SCIP_Real score = calculateScore(scip, config, branchvar, downbranchingresult, upbranchingresult,
            scoringlpobjval, baselpobjval);

#ifdef SCIP_STATISTIC
         if( i == 0 && firstscoreptr != NULL )
            *firstscoreptr = score;
#endif

         if( bestgain != NULL && !config->inscoring && SCIPgetProbingDepth(scip) == 1 && !useoldbranching )
         {
            assert(totalgains != NULL);
            assert(ntotalgains != NULL);

            *bestgain = MAX(*bestgain, score);

            if( !downbranchingresult->cutoff && !upbranchingresult->cutoff )
            {
               (*totalgains) += score;
               (*ntotalgains)++;
            }
         }

         /* both child nodes are infeasible -> the current node is infeasible */
         if( SCIPallColsInLP(scip) && upbranchingresult->cutoff && downbranchingresult->cutoff )
         {
            LABdebugMessage(scip, SCIP_VERBLEVEL_NORMAL, " -> variable <%s> is infeasible in both directions\n",
               SCIPvarGetName(branchvar));

            /* this cutoff may be transferred to a higher level as a domain reduction/valid bound */
            status->cutoff = TRUE;
#ifdef SCIP_STATISTIC
            statistics->nfullcutoffs[probingdepth]++;
            localstats->ncutoffproofnodes += 2;
#endif
         }
         /* up child is infeasible */
         else if( SCIPallColsInLP(scip) && upbranchingresult->cutoff )
         {
            LABdebugMessage(scip, SCIP_VERBLEVEL_NORMAL, " -> variable <%s> is infeasible in upward branch\n",
               SCIPvarGetName(branchvar));

            /* apply down branching bound change at current node if we proved that this node is really infeasible and
             * parameters are set accordingly
             */
            if( config->usedomainreduction && !useoldbranching )
            {
#ifdef SCIP_STATISTIC
               assert(localstats->ncutoffproofnodes == 0 || localstats->ncutoffproofnodes == 2);
               addUpperBound(scip, branchvar, branchval, baselpsol, TRUE, domainreductions,
               2 + localstats->ncutoffproofnodes, TRUE);
#else
               addUpperBound(scip, branchvar, branchval, baselpsol, TRUE, domainreductions);
#endif
            }

            /* the proved bound is given by the bound of the down child alone */
            if( downbranchingresult->dualboundvalid )
            {
               decision->proveddb = MAX(decision->proveddb, downbranchingresult->dualbound);
            }

#ifdef SCIP_STATISTIC
            statistics->nsinglecutoffs[probingdepth]++;
#endif
         }
         /* down child is infeasible */
         else if( SCIPallColsInLP(scip) && downbranchingresult->cutoff )
         {
            LABdebugMessage(scip, SCIP_VERBLEVEL_NORMAL, " -> variable <%s> is infeasible in downward branch\n",
               SCIPvarGetName(branchvar));

            /* apply up branching bound change at current node if we proved that this node is really infeasible and
             * parameters are set accordingly
             */
            if( config->usedomainreduction && !useoldbranching )
            {
#ifdef SCIP_STATISTIC
               assert(localstats->ncutoffproofnodes == 0 || localstats->ncutoffproofnodes == 2);
               addLowerBound(scip, branchvar, branchval, baselpsol, TRUE, domainreductions,
                  2 + localstats->ncutoffproofnodes, TRUE);
#else
               addLowerBound(scip, branchvar, branchval, baselpsol, TRUE, domainreductions);
#endif
            }

            /* the proved bound is given by the bound of the up child alone */
            if( upbranchingresult->dualboundvalid )
            {
               decision->proveddb = MAX(decision->proveddb, upbranchingresult->dualbound);
            }

#ifdef SCIP_STATISTIC
            statistics->nsinglecutoffs[probingdepth]++;
#endif
         }
         /* "normal" case: both child nodes are LP-feasible */
         else
         {
            LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Neither branch is cut off and no limit reached.\n");

            /* the proved dual bound is the minimum of the dual bounds of both child nodes */
            if( upbranchingresult->dualboundvalid && downbranchingresult->dualboundvalid )
            {
               decision->proveddb = MAX(decision->proveddb, MIN(upbranchingresult->dualbound,
                     downbranchingresult->dualbound));
            }
         }

         /* merge domain changes from the two child nodes */
         if( updomainreductions != NULL && config->usedomainreduction && SCIPallColsInLP(scip) )
         {
            int maxstoredomreds = INT_MAX;

            assert(downdomainreductions != NULL);

            if( config->enforcemaxdomreds && config->maxnviolateddomreds > 0)
               maxstoredomreds = config->maxnviolateddomreds;

            if( !upbranchingresult->cutoff && !downbranchingresult->cutoff && config->mergedomainreductions )
               applyDeeperDomainReductions(scip, baselpsol, maxstoredomreds, domainreductions, downdomainreductions,
                  updomainreductions);
            else if( upbranchingresult->cutoff && !downbranchingresult->cutoff )
               applySingleDeeperDomainReductions(scip, baselpsol, maxstoredomreds, domainreductions, downdomainreductions);
            else if( downbranchingresult->cutoff && !upbranchingresult->cutoff )
               applySingleDeeperDomainReductions(scip, baselpsol, maxstoredomreds, domainreductions, updomainreductions);
         }

         if( config->updatebranchingresults && bestscore > -1.0 &&
            (SCIPisGT(scip, decision->proveddb, bestdownbranchingresult->dualbound)
               || SCIPisGT(scip, decision->proveddb, bestupbranchingresult->dualbound)) )
         {
            SCIP_Real newscore;

            bestdownbranchingresult->dualbound = MAX(bestdownbranchingresult->dualbound, decision->proveddb);
            bestupbranchingresult->dualbound = MAX(bestupbranchingresult->dualbound, decision->proveddb);

            newscore = calculateScore(scip, config, decision->branchvar, bestdownbranchingresult, bestupbranchingresult,
               bestscoringlpobjval, baselpobjval);

            if( newscore > bestscore )
            {
               bestscore = newscore;

#ifdef SCIP_STATISTIC
               if( bestscoreptr != NULL )
                  *bestscoreptr = newscore;
#endif
               decision->score = newscore;
               decision->downdb = bestdownbranchingresult->dualbound;
               decision->updb = bestupbranchingresult->dualbound;
            }
         }

         /* the current candidate variable has a better score than the best candidate investigated so far */
         if( SCIPisRelGT(scip, score, bestscore) )
         {
            int nvars = SCIPgetNVars(scip);

            LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Old best var <%s> with bounds [<%g>..<%g>] and score %.9g\n",
               SCIPvarGetName(decision->branchvar), bestscorelowerbound, bestscoreupperbound, bestscore);

            bestscore = score;

#ifdef SCIP_STATISTIC
            if( bestscoreptr != NULL )
               *bestscoreptr = score;
#endif
            decision->branchvar = candidate->branchvar;
            decision->branchval = candidate->branchval;
            decision->downdb = downbranchingresult->dualbound;
            decision->downdbvalid = downbranchingresult->dualboundvalid;
            decision->updb = upbranchingresult->dualbound;
            decision->updbvalid = upbranchingresult->dualboundvalid;
            decision->score = score;

            branchingResultDataCopy(downbranchingresult, bestdownbranchingresult);
            branchingResultDataCopy(upbranchingresult, bestupbranchingresult);

            /* store domain reductions found at the child nodes */
            if( !config->inscoring && updomainreductions != NULL )
            {
               assert(downdomainreductions != NULL);

               SCIP_CALL( branchingDecisionEnsureBoundArraysSize(scip, decision, nvars) );

               BMScopyMemoryArray(decision->uplowerbounds, updomainreductions->lowerbounds, nvars);
               BMScopyMemoryArray(decision->upupperbounds, updomainreductions->upperbounds, nvars);
               BMScopyMemoryArray(decision->downlowerbounds, downdomainreductions->lowerbounds, nvars);
               BMScopyMemoryArray(decision->downupperbounds, downdomainreductions->upperbounds, nvars);
               decision->boundsvalid = TRUE;
            }
            else
            {
               decision->boundsvalid = FALSE;
            }

            bestscorelowerbound = branchlb;
            bestscoreupperbound = branchub;
            bestscoringlpobjval = scoringlpobjval;
            assert(!SCIPisEQ(scip, bestscorelowerbound, bestscoreupperbound));

            LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "New best var <%s> with bounds [<%g>..<%g>] and score %.9g\n",
               SCIPvarGetName(decision->branchvar), bestscorelowerbound, bestscoreupperbound, bestscore);
         }

#ifdef SCIP_DEBUG
         LABdebugMessage(scip, SCIP_VERBLEVEL_NORMAL, " -> cand %d/%d var <%s> (solval=%.9g, downgain=%.9g->%.9g, upgain=%.9g->%.9g,"
            " score=%.9g) -- best: <%s> (%.9g)\n", c, nlpcands, SCIPvarGetName(branchvar), branchval,
            MAX(downbranchingresult->objval - scoringlpobjval, 0), MAX(downbranchingresult->dualbound - scoringlpobjval, 0),
            MAX(upbranchingresult->objval - scoringlpobjval, 0), MAX(upbranchingresult->dualbound - scoringlpobjval, 0),
            score, SCIPvarGetName(decision->branchvar), bestscore);
#endif

         if( config->inscoring )
         {
            assert(scorecontainer != NULL);
            /* only for abbreviated lookahead branching: we are in the FSB filtering step and store the score for this
             * variable and the warm starting basis to reuse it in the subsequent lookahead evaluation of the best
             * candidates
             */
            SCIP_CALL( scoreContainerSetScore(scip, scorecontainer, candidate, score,
                  downbranchingresult->dualbound - scoringlpobjval, upbranchingresult->dualbound - scoringlpobjval) );
         }

         if( probingdepth == 0 && (binconsdata != NULL || domainreductions != NULL) && !useoldbranching
            && (config->maxnviolatedcons >= 0 || config->maxnviolatedbincons >= 0 || config->maxnviolateddomreds >= 0 ) )
         {
            int nbincons = 0;
            int ndomreds = 0;

            if( binconsdata != NULL )
            {
               assert(binconsdata != NULL); /* for lint */
               nbincons = binconsdata->conslist->nviolatedcons;
               LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Found %d binary constraints (%d violated by the LP solution)\n",
                  binconsdata->conslist->nelements, nbincons);

               if( (config->maxnviolatedbincons > 0) && (nbincons >= config->maxnviolatedbincons) )
               {
                  LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The max number of violated binary constraints <%i> is "
                     "exceeded.\n", config->maxnviolatedbincons);
                  status->maxnconsreached = TRUE;
               }
            }

            if( domainreductions != NULL )
            {
               assert(domainreductions != NULL); /* for lint */
               ndomreds = domainreductions->nviolatedvars;
               if( config->prefersimplebounds && ndomreds > domainreductions->nsimplebounds )
                  ndomreds = domainreductions->nsimplebounds;

               LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Found %d bound changes (%d violated by the LP solution)\n",
                  domainreductions->nchangedvars, ndomreds);

               if( (config->maxnviolateddomreds > 0) && (ndomreds >= config->maxnviolateddomreds) )
               {
                  LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The max number of violated bound changes <%i> is "
                     "exceeded.\n", config->maxnviolateddomreds);
                  status->maxnconsreached = TRUE;
               }
            }

            if( config->maxnviolatedcons > 0 && (nbincons + ndomreds >= config->maxnviolatedcons) )
            {
               LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The max number of violated binary constraints and bound "
                  "changes <%d> is exceeded.\n", config->maxnviolatedcons);
               status->maxnconsreached = TRUE;
            }
         }
      }

      if( !(status->domred && decision->branchvar == candidate->branchvar) && areBoundsChanged(scip, decision->branchvar, bestscorelowerbound, bestscoreupperbound) )
      {
         /* in case the bounds of the current highest scored solution have changed due to domain propagation during
          * the lookahead branching we can/should not branch on this variable but instead report the domain
          * reduction */
         status->domred = TRUE;
#ifdef SCIP_STATISTIC
         statistics->npropdomred[probingdepth]++;
#endif
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Domain Propagation changed the bounds of a branching candidate."
               "\n");
      }

      /* free domain reductions */
      if( updomainreductions != NULL )
      {
         assert(downdomainreductions != NULL);

         domainReductionsFree(scip, &updomainreductions);
         domainReductionsFree(scip, &downdomainreductions);
      }
   }

   branchingResultDataFree(scip, &bestupbranchingresult);
   branchingResultDataFree(scip, &bestdownbranchingresult);

   branchingResultDataFree(scip, &upbranchingresult);
   branchingResultDataFree(scip, &downbranchingresult);

   if( persistent != NULL && (!config->abbreviated || config->inscoring) && probingdepth == 0 )
   {
      persistent->restartindex = c;
   }

   return SCIP_OKAY;
}

/** checks whether the current decision should be stored. This is the case if we found domain reductions
 *  or constraints that will be applied, but none of them cuts off the current LP solution.
 *  Then our current decision still holds true for the next call and can be reused without further calculations
 */
static
SCIP_Bool isStoreDecision(
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   BINCONSDATA*          binconsdata,        /**< container collecting all binary constraints; or NULL */
   DOMAINREDUCTIONS*     domainreductions    /**< container collecting all domain reductions found; or NULL */
   )
{
   assert(config != NULL);

   if( !config->storeunviolatedsol )
      return FALSE;

   /* there are violating binary constraints */
   if( binconsdata != NULL && binconsdata->conslist->nviolatedcons > 0 )
      return FALSE;

   /* there are violating domain changes */
   if( domainreductions != NULL && domainreductions->nviolatedvars > 0 )
      return FALSE;

   /* return TRUE if there is at least one domain change or binary constraint */
   return (domainreductions != NULL && domainreductions->nchangedvars > 0)
      || (binconsdata != NULL && binconsdata->conslist->nelements > 0);
}

/** starting point to obtain a branching decision via LAB/ALAB. */
static
SCIP_RETCODE selectVarStart(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,             /**< the configuration of the branching rule */
   PERSISTENTDATA*       persistent,         /**< container to store data over multiple calls to the branching rule; or NULL */
   STATUS*               status,             /**< current status */
   BRANCHINGDECISION*    decision,           /**< struct to store the final decision */
   SCORECONTAINER*       scorecontainer,     /**< container to retrieve already calculated scores; or NULL */
   CANDIDATELIST*        candidatelist       /**< list of candidates to branch on */
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics          /**< general statistical data */
   ,LOCALSTATISTICS*     localstats          /**< local statistics, may be disregarded */
#endif
   )
{
   int recursiondepth;
   DOMAINREDUCTIONS* domainreductions = NULL;
   BINCONSDATA* binconsdata = NULL;
   LEVEL2DATA* level2data = NULL;
   SCIP_SOL* baselpsol = NULL;
   SCIP_Real lpobjval;
#ifdef SCIP_STATISTIC
   SCIP_Real firstscore = -1.0;
   SCIP_Real bestscore = -1.0;
   int chosencandnr = -1;
   SCIP_Bool performedlab = FALSE;
#endif

   assert(scip != NULL);
   assert(config != NULL);
   assert(status != NULL);
   assert(decision != NULL);
   assert(candidatelist != NULL);
#ifdef SCIP_STATISTIC
   assert(statistics != NULL);
#endif

   recursiondepth = config->recursiondepth;
   lpobjval = SCIPgetLPObjval(scip);

   assert(recursiondepth > 0);

   if( SCIP_MAXTREEDEPTH <= (SCIPgetDepth(scip) + recursiondepth) )
   {
      /* we need at least 'recursiondepth' space for the branching */
      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Cannot perform probing in selectVarRecursive, depth limit reached. "
         "Current:<%i>, Max:<%i>\n", SCIP_MAXTREEDEPTH, SCIPgetDepth(scip) + recursiondepth);
      status->depthtoosmall = TRUE;
#ifdef SCIP_STATISTIC
      statistics->ndepthreached++;
#endif
      return SCIP_OKAY;
   }

   assert(!config->inscoring);

   if( candidatelist->ncandidates == 1 )
   {
      decision->branchvar = candidatelist->candidates[0]->branchvar;
      decision->branchval = candidatelist->candidates[0]->branchval;
      decision->downdb = lpobjval;
      decision->downdbvalid = TRUE;
      decision->updb = lpobjval;
      decision->updbvalid = TRUE;
      decision->proveddb = lpobjval;

      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Only one candidate (<%s>) is given. This one is chosen without "
         "calculations.\n", SCIPvarGetName(decision->branchvar));

#ifdef SCIP_STATISTIC
      statistics->nsinglecandidate++;
#endif
      return SCIP_OKAY;
   }
   assert(!SCIPinProbing(scip));

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The objective value of the base lp is <%.9g>\n", lpobjval);

   if( config->usedomainreduction || config->usebincons )
   {
      /* we have to copy the current solution before getting the candidates, as we possibly solve some LPs during
       * the getter and as such would get a wrong LP copied */
      SCIP_CALL( copyCurrentSolution(scip, &baselpsol) );
   }

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "About to start probing.\n");
   SCIP_CALL( SCIPstartStrongbranch(scip, TRUE) );
   SCIPenableVarHistory(scip);

   /* create the binary constraint data */
   if( config->usebincons )
   {
      SCIP_CALL( binConsDataCreate(scip, &binconsdata, recursiondepth,
            (int)SCIPceil(scip, 0.5*candidatelist->ncandidates)) );
   }

   /* collect domain reductions in FSB scoring or LAB branching */
   if( config->usedomainreduction )
   {
      SCIP_CALL( domainReductionsCreate(scip, &domainreductions) );
   }

#ifdef SCIP_STATISTIC
   SCIP_CALL( filterCandidates(scip, status, persistent, config, baselpsol, domainreductions, NULL, candidatelist,
            decision, scorecontainer, level2data, lpobjval,
            statistics, localstats) );
#else
   SCIP_CALL( filterCandidates(scip, status, persistent, config, baselpsol, domainreductions, NULL, candidatelist,
            decision, scorecontainer, level2data, lpobjval) );
#endif

   if( candidatelist->ncandidates == 1 )
   {
      decision->branchvar = candidatelist->candidates[0]->branchvar;
      decision->branchval = candidatelist->candidates[0]->branchval;
      decision->downdb = lpobjval;
      decision->downdbvalid = TRUE;
      decision->updb = lpobjval;
      decision->updbvalid = TRUE;
      decision->proveddb = lpobjval;

      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Only one candidate (<%s>) is given. This one is chosen without "
         "calculations.\n", SCIPvarGetName(decision->branchvar));

#ifdef SCIP_STATISTIC
      statistics->nsingleafterfilter++;
#endif
      goto TERMINATE;
   }

   /* the status may have changed because of FSB to get the best candidates
    * if that is the case, we already changed the base node and should start again */
   if( isBranchFurther(status, TRUE) && candidatelist->ncandidates > 1 )
   {
      assert(candidatelist->ncandidates > 0);

      SCIPstatistic(performedlab = TRUE);

      /* we do not need the level 2 data for FSB scoring, so we do not need to create it before */
      if( recursiondepth == 2 && config->uselevel2data )
      {
         SCIP_CALL( level2dataCreate(scip, &level2data) );
      }

#ifdef SCIP_STATISTIC
      SCIP_CALL( selectVarRecursive(scip, status, persistent, config, baselpsol, domainreductions, binconsdata, candidatelist,
            decision, scorecontainer, level2data, recursiondepth, lpobjval, lpobjval, NULL, NULL, NULL, NULL, NULL, NULL,
            statistics, localstats, &firstscore, &bestscore) );
#else
      SCIP_CALL( selectVarRecursive(scip, status, persistent, config, baselpsol, domainreductions, binconsdata, candidatelist,
            decision, scorecontainer, level2data, recursiondepth, lpobjval, lpobjval, NULL, NULL, NULL, NULL, NULL, NULL) );
#endif

      if( level2data != NULL )
      {
         level2dataFree(scip, &level2data);
      }

      /* only unviolating constraints and domain changes: store branching decision */
      if( persistent != NULL && !status->lperror && isStoreDecision(config, binconsdata, domainreductions) )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "store decision: lpiters=%lld, cand <%s>[%g,%g - %g]\n",
            SCIPgetNNodeLPIterations(scip), SCIPvarGetName(decision->branchvar),
            SCIPvarGetLbLocal(decision->branchvar), SCIPvarGetUbLocal(decision->branchvar),
            SCIPgetSolVal(scip, NULL, decision->branchvar));

         persistent->oldntotalnodes = SCIPgetNTotalNodes(scip);
         persistent->oldnnodelpiterations = SCIPgetNNodeLPIterations(scip);
         persistent->oldnnodelps = SCIPgetNNodeLPs(scip) + SCIPgetNNodeZeroIterationLPs(scip);
         branchingDecisionCopy(decision, persistent->olddecision);
      }

#ifdef SCIP_STATISTIC
      if( config->abbreviated && !status->cutoff && !status->maxnconsreached
         && !status->addedbinconss && !status->domred)
      {
         if( candidatelist->ncandidates > 0 )
         {
            assert(candidatelist->ncandidates <= statistics->maxnbestcands);

            /* find the "FSB-index" of the decision */
            for( chosencandnr = 0; chosencandnr < candidatelist->ncandidates; ++chosencandnr )
            {
               if( decision->branchvar == candidatelist->candidates[chosencandnr]->branchvar )
               {
                  break;
               }
            }
            assert(chosencandnr < candidatelist->ncandidates);
         }
      }
   }
   else
   {
      int probingdepth = 0;
      if( SCIPinProbing(scip) )
         probingdepth = SCIPgetProbingDepth(scip);
      statistics->stopafterfsb[probingdepth]++;

      if( status->cutoff )
      {
         statistics->cutoffafterfsb[probingdepth]++;
      }
      else if( status->maxnconsreached )
      {
         statistics->domredafterfsb[probingdepth]++;
      }
#endif
   }

 TERMINATE:
   SCIP_CALL( SCIPendStrongbranch(scip) );
   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Ended probing.\n");

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Applying found data to the base node.\n");

   /* apply domain reductions */
   if( domainreductions != NULL )
   {
      assert(config->usedomainreduction);

      if( !status->cutoff )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Applying domain reductions to the base node.\n");
#ifdef SCIP_STATISTIC
         SCIP_CALL( applyDomainReductions(scip, config, baselpsol, domainreductions, &status->domredcutoff,
               &status->domred, statistics) );
#else
         SCIP_CALL( applyDomainReductions(scip, config, baselpsol, domainreductions, &status->domredcutoff,
               &status->domred) );
#endif
      }
      domainReductionsFree(scip, &domainreductions);
   }

   /* apply binary constraints */
   if( binconsdata != NULL )
   {
      assert(config->usebincons);
      assert(binconsdata->binaryvars->nbinaryvars == 0);

      if( !status->cutoff )
      {
         SCIP_NODE* basenode = SCIPgetCurrentNode(scip);

         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Applying %d binary constraints to the base node.\n", binconsdata->conslist->nelements);
#ifdef SCIP_STATISTIC
         SCIP_CALL( applyBinaryConstraints(scip, basenode, binconsdata->conslist, config,
               &status->addedbinconss, &status->cutoff, &status->domred, statistics) );
#else
         SCIP_CALL( applyBinaryConstraints(scip, basenode, binconsdata->conslist, config,
               &status->addedbinconss, &status->cutoff, &status->domred) );
#endif
      }
      else
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Discarding %d binary constraints because the base node is cut off.\n", binconsdata->conslist->nelements);
      }
      binConsDataFree(scip, &binconsdata);
   }
   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Applied found data to the base node.\n");

#if defined(SCIP_DEBUG) || defined(SCIP_STATISTIC)
   if( config->abbreviated )
   {
      if( status->domred )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Lookahead Branching has added domain reductions. LAB restarts.\n");

#ifdef SCIP_STATISTIC
         if( candidatelist->ncandidates == 1 )
            statistics->nsingleafterfilter--;
#endif
      }
      else if( status->addedbinconss )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Lookahead Branching has added binary constraints. LAB restarts.\n");

#ifdef SCIP_STATISTIC
         if( candidatelist->ncandidates == 1 )
            statistics->nsingleafterfilter--;
#endif
      }
      else if( status->cutoff )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Lookahead Branching cut this node off.\n");
      }
      else if( candidatelist->ncandidates > 0 )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Strong Branching would branch on variable <%s>\n",
            SCIPvarGetName(candidatelist->candidates[0]->branchvar));

         if( isBranchFurther(status, FALSE) && branchingDecisionIsValid(decision) )
         {
#ifdef SCIP_STATISTIC
            if( chosencandnr >= 0 )
            {
               ++statistics->chosenfsbcand[chosencandnr];

               LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "node %lld chose candidate %d score %16.9g vs %16.9g FSB: %16.9g vs %16.9g\n",
                  SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), chosencandnr,
                  scorecontainer->scores[SCIPvarGetProbindex(candidatelist->candidates[chosencandnr]->branchvar)],
                  scorecontainer->scores[SCIPvarGetProbindex(candidatelist->candidates[0]->branchvar)],
                  bestscore, firstscore);
            }
            else
               assert(!performedlab);
#endif

            LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Lookahead Branching branches on variable <%s>\n",
               SCIPvarGetName(decision->branchvar));
         }
      }
      else
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Something unexpected happened.");
         SCIPABORT();
      }
   }
#endif

   if( baselpsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &baselpsol) );
   }

   return SCIP_OKAY;
}

/**
 * We can use the previous result, stored in the branchruledata, if the branchingvariable (as an indicator) is set and
 * the current lp solution is equal to the previous lp solution.
 *
 * @return \ref TRUE, if we can branch on the previous decision, \ref FALSE, else.
 */
static
SCIP_Bool isUsePreviousResult(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   PERSISTENTDATA* persistent;

   assert(scip != NULL);
   assert(branchruledata != NULL);

   persistent = branchruledata->persistent;
   assert(persistent != NULL);

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "check if previous result should be used: valid=%d, "\
      "nodes=%lld (old=%lld), iterations=%lld (old=%lld), lps=%lld (old=%lld)\n",
      branchingDecisionIsValid(persistent->olddecision),
      SCIPgetNTotalNodes(scip), persistent->oldntotalnodes,
      SCIPgetNNodeLPIterations(scip), persistent->oldnnodelpiterations,
      SCIPgetNNodeLPs(scip) + SCIPgetNNodeZeroIterationLPs(scip), persistent->oldnnodelps);

   return branchingDecisionIsValid(persistent->olddecision)
      && (persistent->oldntotalnodes == SCIPgetNTotalNodes(scip))
      && (persistent->oldnnodelpiterations == SCIPgetNNodeLPIterations(scip))
      && (persistent->oldnnodelps == SCIPgetNNodeLPs(scip) + SCIPgetNNodeZeroIterationLPs(scip));
}

/**
 * Uses the results from the previous run saved in the branchruledata to branch.
 * This is the case, if in the previous run only non-violating constraints were added. In that case we can use the
 * branching decision we would have made then.
 * If everything worked, the result pointer contains SCIP_BRANCHED.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
static
SCIP_RETCODE usePreviousResult(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< branching rule data */
   SCIP_RESULT*          result              /**< the pointer to the branching result */
   )
{
   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(result != NULL);
   assert(branchruledata->config != NULL);
   assert(branchruledata->persistent != NULL);
   assert(branchruledata->persistent->olddecision != NULL);

   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Branching based on previous solution.\n");

   /* execute the actual branching */
   SCIP_CALL( branchOnVar(scip, branchruledata->config, branchruledata->persistent->olddecision) );
   *result = SCIP_BRANCHED;

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Result: Branched based on previous solution. Variable <%s>\n",
      SCIPvarGetName(branchruledata->persistent->olddecision->branchvar));

   /* reset the var pointer, as this is our indicator whether we should branch on prev data in the next call */
   branchruledata->persistent->olddecision->branchvar = NULL;

   return SCIP_OKAY;
}

/** free persistent data structure */
static
SCIP_RETCODE freePersistent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   PERSISTENTDATA* persistent;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(branchruledata != NULL);

   persistent = branchruledata->persistent;
   assert(persistent != NULL);

   nvars = persistent->nvars;

   for( i = nvars - 1; i >= 0; i--)
   {
      assert(persistent->lastbranchdownres[i] != NULL);
      assert(persistent->lastbranchupres[i] != NULL);

      SCIPfreeBlockMemory(scip, &persistent->lastbranchdownres[i]); /*lint !e866*/
      SCIPfreeBlockMemory(scip, &persistent->lastbranchupres[i]); /*lint !e866*/
   }

   SCIPfreeBlockMemory(scip, &branchruledata->persistent->olddecision);

   assert(persistent->lastbranchlpobjval != NULL);
   assert(persistent->lastbranchdownres != NULL);
   assert(persistent->lastbranchupres != NULL);
   assert(persistent->lastbranchnlps != NULL);
   assert(persistent->lastbranchid != NULL);

   SCIPfreeBlockMemoryArray(scip, &persistent->lastbranchlpobjval, nvars);
   SCIPfreeBlockMemoryArray(scip, &persistent->lastbranchdownres, nvars);
   SCIPfreeBlockMemoryArray(scip, &persistent->lastbranchupres, nvars);
   SCIPfreeBlockMemoryArray(scip, &persistent->lastbranchnlps, nvars);
   SCIPfreeBlockMemoryArray(scip, &persistent->lastbranchid, nvars);

   branchruledata->isinitialized = FALSE;

   return SCIP_OKAY;
}

/** initializes the branchruledata and the contained structs */
static
SCIP_RETCODE initBranchruleData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< the branch rule data to initialize */
   )
{
   int nvars;
   int i;

   assert(scip != NULL);
   assert(branchruledata != NULL);

   /* the branching rule data is already initialized and no new variables have been added in the meantime */
   if( branchruledata->isinitialized &&
      (SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) == branchruledata->persistent->nvars) )
      return SCIP_OKAY;

   if( branchruledata->isinitialized )
   {
      SCIP_CALL( freePersistent(scip, branchruledata) );
   }

   /* The variables given by the SCIPgetVars() array are sorted with the binaries at first and the integer variables
    * directly afterwards. With the SCIPvarGetProbindex() method we can access the index of a given variable in the
    * SCIPgetVars() array and as such we can use it to access our arrays which should only contain binary and integer
    * variables.
    */
   nvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->persistent->lastbranchid, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->persistent->lastbranchnlps, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->persistent->lastbranchupres, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->persistent->lastbranchdownres, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->persistent->lastbranchlpobjval, nvars) );
   branchruledata->persistent->nvars = nvars;
   branchruledata->persistent->oldntotalnodes = -1;
   branchruledata->persistent->oldnnodelpiterations = -1;
   branchruledata->persistent->oldnnodelps = -1;

   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata->persistent->olddecision) );
   branchingDecisionInit(scip, branchruledata->persistent->olddecision);

   for( i = 0; i < nvars; i++ )
   {
      branchruledata->persistent->lastbranchid[i] = -1;
      branchruledata->persistent->lastbranchnlps[i] = 0;

      SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata->persistent->lastbranchupres[i]) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata->persistent->lastbranchdownres[i]) ); /*lint !e866*/
   }

   branchruledata->isinitialized = TRUE;

   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Initialized the branchruledata\n");

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyLookahead)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->config != NULL);
   assert(branchruledata->persistent != NULL);

   SCIPfreeBlockMemory(scip, &branchruledata->persistent);
   SCIPfreeBlockMemory(scip, &branchruledata->config);
   SCIPfreeBlockMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL);
   assert(branchrule != NULL);

   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Entering branchInitLookahead\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->persistent != NULL);

   branchruledata->persistent->restartindex = 0;

#ifdef SCIP_STATISTIC
   {
      int recursiondepth;
      int maxncands;

      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Allocating space for the statistics struct.\n");

      recursiondepth = branchruledata->config->recursiondepth;
      maxncands = branchruledata->config->maxncands;
      maxncands = MIN(maxncands, SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip));

      SCIP_CALL( SCIPallocMemory(scip, &branchruledata->statistics) );
      /* RESULT enum is 1 based, so use MAXRESULT + 1 as array size with unused 0 element */
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nresults, MAXRESULT + 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nsinglecutoffs, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nfullcutoffs, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nlpssolved, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nlpssolvedfsb, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nlpiterations, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nlpiterationsfsb, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nduplicatelps, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->npropdomred, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->noldbranchused, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->noldbranchusedfsb, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->chosenfsbcand, maxncands) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->domredafterfsb, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->cutoffafterfsb, recursiondepth) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->stopafterfsb, recursiondepth) );

      branchruledata->statistics->recursiondepth = recursiondepth;
      branchruledata->statistics->maxnbestcands = maxncands;

      statisticsInit(branchruledata->statistics);
   }
#endif

   LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Leaving branchInitLookahead\n");

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitLookahead)
{  /*lint --e{715}*/
#ifdef SCIP_STATISTIC
   SCIP_BRANCHRULEDATA* branchruledata;
   STATISTICS* statistics;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   statistics = branchruledata->statistics;
   assert(statistics != NULL);

   statisticsPrint(scip, statistics);

   SCIPfreeMemoryArray(scip, &statistics->stopafterfsb);
   SCIPfreeMemoryArray(scip, &statistics->cutoffafterfsb);
   SCIPfreeMemoryArray(scip, &statistics->domredafterfsb);
   SCIPfreeMemoryArray(scip, &statistics->chosenfsbcand);
   SCIPfreeMemoryArray(scip, &statistics->noldbranchusedfsb);
   SCIPfreeMemoryArray(scip, &statistics->noldbranchused);
   SCIPfreeMemoryArray(scip, &statistics->npropdomred);
   SCIPfreeMemoryArray(scip, &statistics->nlpiterationsfsb);
   SCIPfreeMemoryArray(scip, &statistics->nlpiterations);
   SCIPfreeMemoryArray(scip, &statistics->nduplicatelps);
   SCIPfreeMemoryArray(scip, &statistics->nlpssolvedfsb);
   SCIPfreeMemoryArray(scip, &statistics->nlpssolved);
   SCIPfreeMemoryArray(scip, &statistics->nfullcutoffs);
   SCIPfreeMemoryArray(scip, &statistics->nsinglecutoffs);
   SCIPfreeMemoryArray(scip, &statistics->nresults);
   SCIPfreeMemory(scip, &statistics);
#endif

   return SCIP_OKAY;
}

/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
static
SCIP_DECL_BRANCHEXITSOL(branchExitSolLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( branchruledata->isinitialized )
   {
      SCIP_CALL( freePersistent(scip, branchruledata) );
   }

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   CONFIGURATION* config;
   SCIP_Bool userusebincons;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Entering branchExeclpLookahead at node %lld.\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   config = branchruledata->config;

   /* we are only allowed to add binary constraints, if the corresponding flag is given */
   userusebincons = config->usebincons;
   config->usebincons = config->usebincons && allowaddcons;

   SCIP_CALL( initBranchruleData(scip, branchruledata) );

   if( config->storeunviolatedsol
      && isUsePreviousResult(scip, branchruledata) )
   {
      /* in case we stopped the previous run without a branching decision, we have stored the decision and execute it
       * now */
      SCIP_CALL( usePreviousResult(scip, branchruledata, result) );

#ifdef SCIP_STATISTIC
      branchruledata->statistics->noldcandidate++;
#endif
   }
   else
   {
      BRANCHINGDECISION* decision;
      SCORECONTAINER* scorecontainer = NULL;
      CANDIDATELIST* candidatelist;
      STATUS* status;
#ifdef SCIP_STATISTIC
      LOCALSTATISTICS* localstats;
#endif

      /* create a struct to store the algorithm status */
      SCIP_CALL( statusCreate(scip, &status) );

      /* create a struct to store the branching decision (in case there is one) */
      SCIP_CALL( branchingDecisionCreate(scip, &decision) );
      if( config->abbreviated )
      {
         /* allocate and init the container used to store the FSB scores, later used to filter the candidates */
         SCIP_CALL( scoreContainerCreate(scip, &scorecontainer, config) );
      }

      SCIP_CALL( candidateListGetAllFractionalCandidates(scip, &candidatelist) );

      LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "The base lp has <%i> variables with fractional value.\n",
         candidatelist->ncandidates);

      /* execute the main logic */
#ifdef SCIP_STATISTIC
      /* create a struct to store the statistics needed for this single run */
      SCIP_CALL( localStatisticsAllocate(scip, &localstats) );
      SCIP_CALL( selectVarStart(scip, config, branchruledata->persistent, status, decision,
            scorecontainer, candidatelist, branchruledata->statistics, localstats) );
#else
      SCIP_CALL( selectVarStart(scip, config, branchruledata->persistent, status, decision,
            scorecontainer, candidatelist) );
#endif

      if( status->cutoff || status->domredcutoff )
      {
         *result = SCIP_CUTOFF;
#ifdef SCIP_STATISTIC
         branchruledata->statistics->ncutoffproofnodes += localstats->ncutoffproofnodes;
#endif
      }
      else if( status->addedbinconss )
      {
         *result = SCIP_CONSADDED;
      }
      else if( status->domred )
      {
         *result = SCIP_REDUCEDDOM;
      }
      else if( status->lperror )
      {
#ifdef SCIP_STATISTIC
         ++branchruledata->statistics->nlperrorcalls;
#endif
         if( !branchingDecisionIsValid(decision) )
         {
            LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "LP error with no valid candidate: select first candidate variable\n");

            assert(candidatelist->ncandidates > 0);
            decision->branchvar = candidatelist->candidates[0]->branchvar;
            decision->branchval = candidatelist->candidates[0]->branchval;
         }
      }
      else if( status->maxnconsreached )
      {
         /* this case may occure if the domain reductions that reached the limit were already applied via domain
          * propagation
          */
         *result = SCIP_REDUCEDDOM;
      }
#ifdef SCIP_STATISTIC
      else if( status->limitreached )
      {
         ++branchruledata->statistics->nlimitcalls;
      }
#endif

      LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Result before branching is %s\n", getStatusString(*result));

      if( *result != SCIP_CUTOFF /* a variable could not be branched in any direction or any of the calculated domain
                                  * reductions was infeasible */
         && *result != SCIP_REDUCEDDOM /* the domain of a variable was reduced by evaluating the calculated cutoffs */
         && *result != SCIP_CONSADDED /* implied binary constraints were already added */
         && !status->depthtoosmall /* branching depth wasn't high enough */
         && branchingDecisionIsValid(decision)
         /*&& (0 <= bestcand && bestcand < nlpcands)*/ /* no valid candidate index could be found */
         )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_NORMAL, " -> %d candidates, selected variable <%s> (solval=%g, down=%.9g, "
            "up=%.9g)\n", candidatelist->ncandidates, SCIPvarGetName(decision->branchvar), decision->branchval,
            decision->downdb, decision->updb);

         /* execute the branching as a result of the branching logic */
         SCIP_CALL( branchOnVar(scip, config, decision) );

         *result = SCIP_BRANCHED;
      }

#ifdef SCIP_DEBUG
      LABdebugMessage(scip, SCIP_VERBLEVEL_FULL, "Result after branching is %s\n", getStatusString(*result));

      if( *result == SCIP_BRANCHED )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Result: Finished LookaheadBranching by branching.\n");
      }
      else if( *result == SCIP_REDUCEDDOM )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Result: Finished LookaheadBranching by reducing domains.\n");
      }
      else if( *result == SCIP_CUTOFF )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Result: Finished LookaheadBranching by cutting off, as the current "
            "problem is infeasible.\n");
      }
      else if( *result == SCIP_CONSADDED )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Result: Finished LookaheadBranching by adding constraints.\n");
      }
      else if( status->depthtoosmall )
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Result: The remaining tree depth did not allow for multi level "
            "lookahead branching.\n");
      }
      else
      {
         LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Result: Could not find any variable to branch on.\n");
      }
#endif

#ifdef SCIP_STATISTIC
      localStatisticsFree(scip, &localstats);
#endif
      SCIP_CALL( candidateListFree(scip, &candidatelist) );

      /* scorecontainer != NULL iff branchruledata->config->abbreviated == TRUE */
      if( scorecontainer != NULL )
      {
         SCIP_CALL( scoreContainerFree(scip, &scorecontainer) );
      }
      branchingDecisionFree(scip, &decision);
      statusFree(scip, &status);
   }

#ifdef SCIP_STATISTIC
   assert(*result >= 1);
   assert(*result <= MAXRESULT);
   branchruledata->statistics->ntotalresults++;
   branchruledata->statistics->nresults[*result]++;

   if( config->abbreviated )
   {
      int sum;
      int i;

      sum = branchruledata->statistics->nsinglecandidate + branchruledata->statistics->nsingleafterfilter
         + branchruledata->statistics->noldcandidate + branchruledata->statistics->nlperrorcalls
         + branchruledata->statistics->nlimitcalls;

      for( i = 0; i < branchruledata->statistics->maxnbestcands; i++ )
      {
         sum += branchruledata->statistics->chosenfsbcand[i];
      }
      if( sum != branchruledata->statistics->nresults[SCIP_BRANCHED] )
      {
         printf("branched = %d != sum = %d (%d/%d/%d/%d/%d)\n",
            branchruledata->statistics->nresults[SCIP_BRANCHED], sum,
            branchruledata->statistics->nsinglecandidate, branchruledata->statistics->nsingleafterfilter,
            branchruledata->statistics->noldcandidate,
            branchruledata->statistics->nlperrorcalls, branchruledata->statistics->nlimitcalls);
         assert(SCIPisStopped(scip));
      }
      assert(sum == branchruledata->statistics->nresults[SCIP_BRANCHED]);
   }

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "#### ncutoffproofnodes: %d ndomredproofnodes: %d\n",
      branchruledata->statistics->ncutoffproofnodes, branchruledata->statistics->ndomredproofnodes);
#endif

   config->usebincons = userusebincons;

   LABdebugMessage(scip, SCIP_VERBLEVEL_HIGH, "Exiting branchExeclpLookahead.\n");
   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the lookahead branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleLookahead(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create lookahead branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata->config) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata->persistent) );
   branchruledata->persistent->restartindex = 0;
   branchruledata->isinitialized = FALSE;
   branchruledata->config->inscoring = FALSE;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyLookahead) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeLookahead) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitLookahead) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitLookahead) );
   SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitSolLookahead) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpLookahead) );

   /* add lookahead branching rule parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/useimpliedbincons",
         "should binary constraints be collected and applied?",
         &branchruledata->config->usebincons, TRUE, DEFAULT_USEBINARYCONSTRAINTS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/lookahead/addbinconsrow",
         "should binary constraints be added as rows to the base LP? (0: no, 1: separate, 2: as initial rows)",
         &branchruledata->config->addbinconsrow, TRUE, DEFAULT_ADDBINCONSROW, 0, 2, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead/maxnviolatedcons",
         "how many constraints that are violated by the base lp solution should be gathered until the rule is stopped and "\
         "they are added? [0 for unrestricted]",
         &branchruledata->config->maxnviolatedcons, TRUE, DEFAULT_MAXNVIOLATEDCONS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead/maxnviolatedbincons",
         "how many binary constraints that are violated by the base lp solution should be gathered until the rule is "\
         "stopped and they are added? [0 for unrestricted]",
         &branchruledata->config->maxnviolatedbincons, TRUE, DEFAULT_MAXNVIOLATEDBINCONS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead/maxnviolateddomreds",
         "how many domain reductions that are violated by the base lp solution should be gathered until the rule is "\
         "stopped and they are added? [0 for unrestricted]",
         &branchruledata->config->maxnviolateddomreds, TRUE, DEFAULT_MAXNVIOLATEDDOMREDS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "branching/lookahead/reevalage",
         "max number of LPs solved after which a previous prob branching results are recalculated",
         &branchruledata->config->reevalage, TRUE, DEFAULT_REEVALAGE, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "branching/lookahead/reevalagefsb",
         "max number of LPs solved after which a previous FSB scoring results are recalculated",
         &branchruledata->config->reevalagefsb, TRUE, DEFAULT_REEVALAGEFSB, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead/recursiondepth",
         "the max depth of LAB.",
         &branchruledata->config->recursiondepth, TRUE, DEFAULT_RECURSIONDEPTH, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/usedomainreduction",
         "should domain reductions be collected and applied?",
         &branchruledata->config->usedomainreduction, TRUE, DEFAULT_USEDOMAINREDUCTION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/mergedomainreductions",
         "should domain reductions of feasible siblings should be merged?",
         &branchruledata->config->mergedomainreductions, TRUE, DEFAULT_MERGEDOMAINREDUCTIONS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/prefersimplebounds",
         "should domain reductions only be applied if there are simple bound changes?",
         &branchruledata->config->prefersimplebounds, TRUE, DEFAULT_PREFERSIMPLEBOUNDS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/onlyvioldomreds",
         "should only domain reductions that violate the LP solution be applied?",
         &branchruledata->config->onlyvioldomreds, TRUE, DEFAULT_ONLYVIOLDOMREDS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/addnonviocons",
         "should binary constraints, that are not violated by the base LP, be collected and added?",
         &branchruledata->config->addnonviocons, TRUE, DEFAULT_ADDNONVIOCONS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/abbreviated",
         "toggles the abbreviated LAB.",
         &branchruledata->config->abbreviated, TRUE, DEFAULT_ABBREVIATED, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead/maxncands",
         "if abbreviated: The max number of candidates to consider at the node.",
         &branchruledata->config->maxncands, TRUE, DEFAULT_MAXNCANDS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead/maxndeepercands",
         "if abbreviated: The max number of candidates to consider per deeper node.",
         &branchruledata->config->maxndeepercands, TRUE, DEFAULT_MAXNDEEPERCANDS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/reusebasis",
         "if abbreviated: Should the information gathered to obtain the best candidates be reused?",
         &branchruledata->config->reusebasis, TRUE, DEFAULT_REUSEBASIS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/storeunviolatedsol",
         "if only non violating constraints are added, should the branching decision be stored till the next call?",
         &branchruledata->config->storeunviolatedsol, TRUE, DEFAULT_STOREUNVIOLATEDSOL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/abbrevpseudo",
         "if abbreviated: Use pseudo costs to estimate the score of a candidate.",
         &branchruledata->config->abbrevpseudo, TRUE, DEFAULT_ABBREVPSEUDO, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/level2avgscore",
         "should the average score be used for uninitialized scores in level 2?",
         &branchruledata->config->level2avgscore, TRUE, DEFAULT_LEVEL2AVGSCORE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/level2zeroscore",
         "should uninitialized scores in level 2 be set to 0?",
         &branchruledata->config->level2zeroscore, TRUE, DEFAULT_LEVEL2ZEROSCORE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/addclique",
         "add binary constraints with two variables found at the root node also as a clique",
         &branchruledata->config->addclique, TRUE, DEFAULT_ADDCLIQUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/propagate",
         "should domain propagation be executed before each temporary node is solved?",
         &branchruledata->config->propagate, TRUE, DEFAULT_PROPAGATE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/uselevel2data",
         "should branching data generated at depth level 2 be stored for re-using it?",
         &branchruledata->config->uselevel2data, TRUE, DEFAULT_USELEVEL2DATA, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/applychildbounds",
         "should bounds known for child nodes be applied?",
         &branchruledata->config->applychildbounds, TRUE, DEFAULT_APPLYCHILDBOUNDS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/enforcemaxdomreds",
         "should the maximum number of domain reductions maxnviolateddomreds be enforced?",
         &branchruledata->config->enforcemaxdomreds, TRUE, DEFAULT_ENFORCEMAXDOMREDS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/updatebranchingresults",
         "should branching results (and scores) be updated w.r.t. proven dual bounds?",
         &branchruledata->config->updatebranchingresults, TRUE, DEFAULT_UPDATEBRANCHINGRESULTS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/lookahead/maxproprounds",
         "maximum number of propagation rounds to perform at each temporary node (-1: unlimited, 0: SCIP default)",
         &branchruledata->config->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip,
         "branching/lookahead/scoringfunction",
         "scoring function to be used at the base level",
         &branchruledata->config->scoringfunction, TRUE, DEFAULT_SCORINGFUNCTION, "dfswplcra", NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip,
         "branching/lookahead/deeperscoringfunction",
         "scoring function to be used at deeper levels",
         &branchruledata->config->deeperscoringfunction, TRUE, DEFAULT_DEEPERSCORINGFUNCTION, "dfswlcrx", NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip,
         "branching/lookahead/scoringscoringfunction",
         "scoring function to be used during FSB scoring",
         &branchruledata->config->scoringscoringfunction, TRUE, DEFAULT_SCORINGSCORINGFUNCTION, "dfswlcr", NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/lookahead/minweight",
         "if scoringfunction is 's', this value is used to weight the min of the gains of two child problems in the convex combination",
         &branchruledata->config->minweight, TRUE, DEFAULT_MINWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/lookahead/worsefactor",
         "if the FSB score is of a candidate is worse than the best by this factor, skip this candidate (-1: disable)",
         &branchruledata->config->worsefactor, TRUE, DEFAULT_WORSEFACTOR, -1.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/filterbymaxgain",
         "should lookahead branching only be applied if the max gain in level 1 is not uniquely that of the best candidate?",
         &branchruledata->config->filterbymaxgain, TRUE, DEFAULT_FILTERBYMAXGAIN, NULL, NULL) );

   return SCIP_OKAY;
}
