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

/**@file   heur_padm.c
 * @brief  PADM primal heuristic based on ideas published in the paper
 *         "A Decomposition Heuristic for Mixed-Integer Supply Chain Problems"
 *         by Martin Schmidt, Lars Schewe, and Dieter Weninger
 * @author Dieter Weninger
 * @author Katrin Halbig
 *
 * The penalty alternating direction method (PADM) heuristic is a construction heuristic which additionally needs a
 * user decomposition with linking variables only.
 *
 * PADM splits the problem into several sub-SCIPs according to the decomposition, whereby the linking variables get
 * copied and the difference is penalized. Then the sub-SCIPs are solved on an alternating basis until they arrive at
 * the same values of the linking variables (ADM-loop). If they don't reconcile after a couple of iterations,
 * the penalty parameters are increased (penalty-loop) and the sub-SCIPs are solved again on an alternating basis.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "blockmemshell/memory.h"
#include "scip/cons_linear.h"
#include "scip/debug.h"
#include "scip/heur_padm.h"
#include "scip/heuristics.h"
#include "scip/pub_cons.h"
#include "scip/pub_tree.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_select.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/scipdefplugins.h"
#include "scip/scip_branch.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_dcmp.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nodesel.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_table.h"
#include "scip/scip_timing.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"

#define HEUR_NAME             "padm"
#define HEUR_DESC             "penalty alternating direction method primal heuristic"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_LNS
#define HEUR_PRIORITY         70000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE                  /**< does the heuristic use a secondary SCIP instance? */

#define COUPLINGSIZE          3
#define DEFAULT_MINNODES      50LL
#define DEFAULT_MAXNODES      5000LL
#define DEFAULT_NODEFAC       0.8
#define DEFAULT_ADMIT         4
#define DEFAULT_PENALTYIT     100
#define DEFAULT_GAP           2.0

/*
 * Data structures
 */

/** data related to one problem (see below) */
typedef struct Problem PROBLEM;

/** data related to one block */
typedef struct Block
{
   PROBLEM*              problem;            /**< the problem this block belongs to */
   SCIP*                 subscip;            /**< sub-SCIP representing this block */
   int                   number;             /**< component number */
   SCIP_VAR**            subvars;            /**< variables belonging to this block (without slack variables) */
   int                   nsubvars;           /**< number of variables belonging to this block (without slack variables) */
   SCIP_VAR**            slackspos;          /**< positive slack variables */
   SCIP_VAR**            slacksneg;          /**< negative slack variables */
   SCIP_CONS**           couplingcons;       /**< coupling contraints */
   int                   ncoupling;          /**< number of coupling contraints (equal to positive/negative slack variables) */
   SCIP_Real             size;               /**< share of total problem */
} BLOCK;

/** data related to one problem */
struct Problem
{
   SCIP*                 scip;               /**< the SCIP instance this problem belongs to */
   char*                 name;               /**< name of the problem */
   BLOCK*                blocks;             /**< blocks into which the problem will be divided */
   int                   nblocks;            /**< number of blocks */
};

/** set data structure */
typedef struct set
{
   int                   size;               /**< size of the set */
   int*                  indexes;            /**< set of indexes */
} SET;

/** data of one linking variable related to one block */
typedef struct blockinfo
{
   int                   block;              /**< index of this block */
   int                   otherblock;         /**< index of the other connected block */
   int                   linkvaridx;         /**< linking variable index */
   SCIP_Real             linkvarval;         /**< value of linking variable */
   SCIP_VAR*             linkvar;            /**< linking variable */
   SCIP_Real             slackposobjcoeff;   /**< penalty coefficient of positive slack variable */
   SCIP_VAR*             slackposvar;        /**< positive slack variable */
   SCIP_Real             slacknegobjcoeff;   /**< penalty coefficient of negative slack variable */
   SCIP_VAR*             slacknegvar;        /**< negative slack variable */
   SCIP_CONS*            couplingCons;       /**< coupling contraint (equation) */
} BLOCKINFO;

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(indexesEqual)
{  /*lint --e{715}*/
   BLOCKINFO* binfo1;
   BLOCKINFO* binfo2;

   binfo1 = (BLOCKINFO*) key1;
   binfo2 = (BLOCKINFO*) key2;

   if( binfo1->block != binfo2->block || binfo1->otherblock != binfo2->otherblock ||
         binfo1->linkvaridx != binfo2->linkvaridx )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(indexesHashval)
{  /*lint --e{715}*/
   BLOCKINFO* binfo;
   binfo = (BLOCKINFO*) key;

   return SCIPhashFour(SCIPrealHashCode((double)binfo->block), SCIPrealHashCode((double)binfo->otherblock),
                       SCIPrealHashCode((double)binfo->linkvaridx), SCIPrealHashCode((double)binfo->linkvaridx));
}

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in all subproblems */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in one subproblem */
   int                   admiterations;      /**< maximal number of ADM iterations in each penalty loop */
   int                   penaltyiterations;  /**< maximal number of penalty iterations */
   int                   timing;             /**< should the heuristic run before or after the processing of the node?
                                                  (0: before, 1: after, 2: both) */
   SCIP_Real             nodefac;            /**< factor to control nodelimits of subproblems */
   SCIP_Real             gap;                /**< mipgap at start */
   SCIP_Bool             reoptimize;         /**< should the problem get reoptimized with the original objective function? */
   SCIP_Bool             scaling;            /**< enable sigmoid rescaling of penalty parameters */
   SCIP_Bool             assignlinking;      /**< should linking constraints be assigned? */
   SCIP_Bool             original;           /**< should the original problem be used? */
};

/*
 * Local methods
 */

/** initializes one block */
static
SCIP_RETCODE initBlock(
   PROBLEM*              problem             /**< problem structure */
   )
{
   BLOCK* block;

   assert(problem != NULL);
   assert(problem->scip != NULL);

   block = &problem->blocks[problem->nblocks];

   block->problem = problem;
   block->subscip = NULL;
   block->subvars = NULL;
   block->nsubvars = 0;
   block->number = problem->nblocks;
   block->slackspos = NULL;
   block->slacksneg = NULL;
   block->couplingcons = NULL;
   block->ncoupling = 0;
   block->size = 0;

   ++problem->nblocks;

   return SCIP_OKAY;
}

/** frees component structure */
static
SCIP_RETCODE freeBlock(
   BLOCK*                block               /**< block structure */
   )
{
   assert(block != NULL);

   block->ncoupling = 0;

   if( block->subvars != NULL )
   {
      SCIPfreeBufferArray(block->problem->scip, &(block->subvars));
   }

   if( block->subscip != NULL )
   {
      SCIP_CALL( SCIPfree(&block->subscip) );
   }

   return SCIP_OKAY;
}

/** initializes subproblem structure */
static
SCIP_RETCODE initProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBLEM**             problem,            /**< pointer to problem structure */
   int                   nblocks             /**< number of blocks */
   )
{
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, problem) );
   assert(*problem != NULL);

   (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*problem)->name, name, strlen(name) + 1) );

   SCIPdebugMessage("initialized problem <%s>\n", (*problem)->name);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*problem)->blocks, nblocks) );

   (*problem)->scip = scip;
   (*problem)->nblocks = 0;

   return SCIP_OKAY;
}

/** frees subproblem structure */
static
SCIP_RETCODE freeProblem(
   PROBLEM**             problem,            /**< pointer to problem to free */
   int                   nblocks             /**< number of blocks in decomposition */
   )
{
   SCIP* scip;
   int c;

   assert(problem != NULL);
   assert(*problem != NULL);

   scip = (*problem)->scip;
   assert(scip != NULL);

   /* free all blocks */
   for( c = nblocks - 1; c >= 0; --c )
   {
      SCIP_CALL( freeBlock(&(*problem)->blocks[c]) );
   }
   if( (*problem)->blocks != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*problem)->blocks, nblocks);
   }

   /* free problem name */
   SCIPfreeMemoryArray(scip, &(*problem)->name);

   /* free PROBLEM struct and set the pointer to NULL */
   SCIPfreeBlockMemory(scip, problem);
   *problem = NULL;

   return SCIP_OKAY;
}

/** creates a sub-SCIP for the given variables and constraints */
static
SCIP_RETCODE createSubscip(
   SCIP*                 scip,               /**< main SCIP data structure */
   SCIP**                subscip             /**< pointer to store created sub-SCIP */
   )
{
   SCIP_Real infvalue;

   /* create a new SCIP instance */
   SCIP_CALL( SCIPcreate(subscip) );

   SCIP_CALL( SCIPincludeDefaultPlugins(*subscip) );

   /* copy value for infinity */
   SCIP_CALL( SCIPgetRealParam(scip, "numerics/infinity", &infvalue) );
   SCIP_CALL( SCIPsetRealParam(*subscip, "numerics/infinity", infvalue) );

   SCIP_CALL( SCIPcopyLimits(scip, *subscip) );

   /* avoid recursive calls */
   SCIP_CALL( SCIPsetSubscipsOff(*subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(*subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(*subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* disable expensive techniques */
   SCIP_CALL( SCIPsetIntParam(*subscip, "misc/usesymmetry", 0) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(*subscip, "misc/catchctrlc", FALSE) );

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/freq", 100000000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(*subscip, "timing/statistictiming", FALSE) );
#endif

   return SCIP_OKAY;
}

/** copies the given constraints and the corresponding variables to the given sub-SCIP */
static
SCIP_RETCODE copyToSubscip(
   SCIP*                 scip,               /**< source SCIP */
   SCIP*                 subscip,            /**< target SCIP */
   const char*           name,               /**< name for copied problem */
   SCIP_CONS**           conss,              /**< constraints to copy */
   SCIP_HASHMAP*         varmap,             /**< hashmap used for the copy process of variables */
   SCIP_HASHMAP*         consmap,            /**< hashmap used for the copy process of constraints */
   int                   nconss,             /**< number of constraints to copy */
   SCIP_Bool             useorigprob,        /**< do we use the original problem? */
   SCIP_Bool*            success             /**< pointer to store whether copying was successful */
   )
{
   SCIP_CONS* newcons;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(conss != NULL);
   assert(consmap != NULL);
   assert(success != NULL);

   *success = TRUE;
   newcons = NULL;

   /* create problem in sub-SCIP */
   SCIP_CALL( SCIPcopyProb(scip, subscip, varmap, consmap, FALSE, name) );

   /* copy constraints */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);

      /* do not check this if we use the original problem
       * Since constraints can be deleted etc. during presolving, these assertions would fail.
       */
      if( !useorigprob )
      {
         assert(!SCIPconsIsModifiable(conss[i]));
         assert(SCIPconsIsActive(conss[i]));
         assert(!SCIPconsIsDeleted(conss[i]));
      }

      /* copy the constraint */
      SCIP_CALL( SCIPgetConsCopy(scip, subscip, conss[i], &newcons, SCIPconsGetHdlr(conss[i]), varmap, consmap, NULL,
                                SCIPconsIsInitial(conss[i]), SCIPconsIsSeparated(conss[i]), SCIPconsIsEnforced(conss[i]),
                                SCIPconsIsChecked(conss[i]), SCIPconsIsPropagated(conss[i]), FALSE, FALSE,
                                SCIPconsIsDynamic(conss[i]), SCIPconsIsRemovable(conss[i]), FALSE, FALSE, success) );

      /* abort if constraint was not successfully copied */
      if( !(*success) || newcons == NULL)
         return SCIP_OKAY;

      SCIP_CALL( SCIPaddCons(subscip, newcons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &newcons) );
   }

   return SCIP_OKAY;
}

/** creates the subscip for a given block */
static
SCIP_RETCODE blockCreateSubscip(
   BLOCK*                block,              /**< block structure */
   SCIP_HASHMAP*         varmap,             /**< variable hashmap used to improve performance */
   SCIP_HASHMAP*         consmap,            /**< constraint hashmap used to improve performance */
   SCIP_CONS**           conss,              /**< constraints contained in this block */
   int                   nconss,             /**< number of constraints contained in this block */
   SCIP_Bool             useorigprob,        /**< do we use the original problem? */
   SCIP_Bool*            success             /**< pointer to store whether the copying process was successful */
   )
{
   char name[SCIP_MAXSTRLEN];
   PROBLEM* problem;
   SCIP* scip;
   SCIP_VAR** subscipvars;
   int nsubscipvars;
   int i;

   assert(block != NULL);
   assert(varmap != NULL);
   assert(consmap != NULL);
   assert(conss != NULL);
   assert(success != NULL);

   problem = block->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   (*success) = TRUE;

   SCIP_CALL( createSubscip(scip, &block->subscip) );

   if( block->subscip != NULL )
   {
      /* get name of the original problem and add "comp_nr" */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d", problem->name, block->number);

      SCIP_CALL( copyToSubscip(scip, block->subscip, name, conss, varmap, consmap, nconss, useorigprob, success) );

      if( !(*success) )
      {
         SCIP_CALL( SCIPfree(&block->subscip) );
         block->subscip = NULL;
      }
      else
      {
         /* save variables of subscip (without slack variables) */
         nsubscipvars = SCIPgetNOrigVars(block->subscip);
         subscipvars = SCIPgetOrigVars(block->subscip);
         SCIP_CALL( SCIPallocBufferArray(scip, &(block->subvars), nsubscipvars) );
         block->nsubvars = nsubscipvars;
         for( i = 0; i < nsubscipvars; i++ )
            block->subvars[i] = subscipvars[i];

         /* calculate size of sub-SCIP with focus on the number of integer variables
          * we use this value to determine the nodelimit
          */
         block->size = (SCIP_Real)(SCIPgetNOrigVars(block->subscip) + SCIPgetNOrigIntVars(block->subscip) + SCIPgetNOrigBinVars(block->subscip)) /
                       (SCIPgetNVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNBinVars(scip));
      }
   }
   else
      (*success) = FALSE;

   SCIPdebugMsg(scip, "created subscip of block %d\n", block->number);

   return SCIP_OKAY;
}

/** creates problem structure and split it into blocks */
static
SCIP_RETCODE createAndSplitProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           sortedconss,        /**< array of (checked) constraints sorted by blocks */
   int                   nconss,             /**< number of constraints */
   int*                  consssize,          /**< number of constraints per block (and border at index 0) */
   int                   nblocks,            /**< number of blocks */
   PROBLEM**             problem,            /**< pointer to store problem structure */
   SCIP_Bool             useorigprob,        /**< do we use the original problem? */
   SCIP_Bool*            success             /**< pointer to store whether the process was successful */
   )
{
   BLOCK* block;
   SCIP_HASHMAP* varmap;
   SCIP_HASHMAP* consmap;
   SCIP_CONS** blockconss;
   int nhandledconss;
   int nblockconss;
   int b;

   (*success) = TRUE;

   /* init subproblem data structure */
   SCIP_CALL( initProblem(scip, problem, nblocks) );
   assert((*problem)->blocks != NULL);

   /* hashmap mapping from original constraints to constraints in the sub-SCIPs (for performance reasons) */
   SCIP_CALL( SCIPhashmapCreate(&consmap, SCIPblkmem(scip), nconss) );

   for( b = 0; b < nblocks; b++ )
   {
      SCIP_CALL( initBlock(*problem) );
   }

   /* loop over all blocks and create subscips */
   nhandledconss = 0;
   for( b = 0; b < nblocks; b++ )
   {
      block = &(*problem)->blocks[b];

      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );

      /* get block constraints */
      blockconss = &(sortedconss[nhandledconss]);
      nblockconss = consssize[b + 1];

      /* build subscip for block */
      SCIP_CALL( blockCreateSubscip(block, varmap, consmap, blockconss, nblockconss, useorigprob, success) );

      SCIPhashmapFree(&varmap);
      nhandledconss += nblockconss;

      if( !(*success) )
         break;
   }

   SCIPhashmapFree(&consmap);

   if( !(*success) )
   {
      /* free subproblem data structure since not all blocks could be copied */
      SCIP_CALL( freeProblem(problem, nblocks) );
   }

   return SCIP_OKAY;
}

/** copies labels to newdecomp and assigns linking constraints if possible*/
static
SCIP_RETCODE assignLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          newdecomp,          /**< decomposition with (partially) assigned linking constraints */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_CONS**           sortedconss,        /**< sorted array of constraints */
   int*                  varlabels,          /**< array of variable labels */
   int*                  conslabels,         /**< sorted array of constraint labels */
   int                   nvars,              /**< number of variables */
   int                   nconss,             /**< number of constraints */
   int                   nlinkconss          /**< number of linking constraints */
   )
{
   assert(scip != NULL);
   assert(vars != NULL);
   assert(sortedconss != NULL);
   assert(varlabels != NULL);
   assert(conslabels != NULL);

   /* copy the labels */
   SCIP_CALL( SCIPdecompSetVarsLabels(newdecomp, vars, varlabels, nvars) );
   SCIP_CALL( SCIPdecompSetConsLabels(newdecomp, sortedconss, conslabels, nconss) );

   SCIPdebugMsg(scip, "try to assign %d linking constraints\n", nlinkconss);

   /* reassign linking constraints */
   SCIP_CALL( SCIPassignDecompLinkConss(scip, newdecomp, &sortedconss[0], nlinkconss, NULL) );

   SCIP_CALL( SCIPcomputeDecompVarsLabels(scip, newdecomp, sortedconss, nconss) );

   SCIP_CALL( SCIPcomputeDecompStats(scip, newdecomp, TRUE) );

   SCIPdecompGetConsLabels(newdecomp, sortedconss, conslabels, nconss);
   SCIPdecompGetVarsLabels(newdecomp, vars, varlabels, nvars);

   SCIPsortIntPtr(conslabels, (void**)sortedconss, nconss);

   return SCIP_OKAY;
}

/** computes feasible solution from last stored solution of the block*/
static
SCIP_RETCODE reuseSolution(
   SCIP*                 subscip,            /**< SCIP data structure */
   BLOCK*                block               /**< block structure*/
   )
{
   SCIP_SOL** sols;
   SCIP_SOL* sol; /* solution of block that will be repaired */
   SCIP_SOL* newsol;
   SCIP_VAR** blockvars;
   SCIP_VAR** consvars;
   SCIP_Real* blockvals;
   int nsols;
   int nvars;
   int c;
   SCIP_Bool success;

   assert(subscip != NULL);
   assert(block != NULL);

   nsols = SCIPgetNSols(subscip);

   /* no solution in solution candidate storage found */
   if( nsols == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(subscip, &consvars, COUPLINGSIZE) );

   sols = SCIPgetSols(subscip);
   sol = sols[nsols - 1];

   /* copy the solution */
   nvars = SCIPgetNVars(subscip);
   blockvars = SCIPgetVars(subscip);
   SCIP_CALL( SCIPallocBufferArray(subscip, &blockvals, nvars) );
   SCIP_CALL( SCIPgetSolVals(subscip, sol, nvars, blockvars, blockvals) );
   SCIP_CALL( SCIPcreateOrigSol(subscip, &newsol, NULL) );
   SCIP_CALL( SCIPsetSolVals(subscip, newsol, nvars, blockvars, blockvals) );

   /* correct each coupling constraint;
    * orig_var + slackpos - slackneg == side
    * adapt slack variables so that constraint is feasible */
   for( c = 0; c < block->ncoupling; c++ )
   {
      SCIP_Real solval; /* old solution values of variables; [0] original variable, [1] slackpos, [2] slackneg */
      SCIP_Real side; /* current right hand side */
      SCIP_Real diff;

      SCIP_CALL( SCIPgetConsVars(subscip, block->couplingcons[c], consvars, COUPLINGSIZE, &success) );
      solval = SCIPgetSolVal(subscip, sol, consvars[0]);

      side = SCIPgetRhsLinear(subscip, block->couplingcons[c]);
      assert(SCIPisEQ(subscip, SCIPgetRhsLinear(subscip, block->couplingcons[c]), SCIPgetLhsLinear(subscip, block->couplingcons[c])));

      diff = side - solval;

      /* slackpos is strict positiv */
      if( diff > 0 )
      {
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slackspos[c], diff) );
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slacksneg[c], 0.0) );
      }
      /* slackneg is strict positiv */
      else if( diff < 0 )
      {
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slacksneg[c], -diff) );
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slackspos[c], 0.0) );
      }
      /* no slack variable necessary */
      else
      {
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slackspos[c], 0.0) );
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slacksneg[c], 0.0) );
      }
   }

   SCIPdebugMsg(subscip, "Try adding solution with objective value %.2f\n", SCIPgetSolOrigObj(subscip, newsol));
   SCIP_CALL( SCIPaddSolFree(subscip, &newsol, &success) );

   if( !success )
      SCIPdebugMsg(subscip, "Correcting solution failed\n"); /* maybe not better than old solutions */
   else
   {
      SCIPdebugMsg(subscip, "Correcting solution successful\n");
   }

   SCIPfreeBufferArray(subscip, &blockvals);
   SCIPfreeBufferArray(subscip, &consvars);

   return SCIP_OKAY;
}

/** reoptimizes the heuristic solution with original objective function
 *
 * Since the main algorithm of padm ignores the objective function, this method can be called to obtain better solutions.
 * It copies the main scip, fixes the linking variables at the values of the already found solution
 * and solves the new problem with small limits.
 */
static
SCIP_RETCODE reoptimize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< pointer to heuristic*/
   SCIP_SOL*             sol,                /**< heuristic solution  */
   SCIP_VAR**            vars,               /**< pointer to variables */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            linkvars,           /**< pointer to linking variables */
   int                   nlinkvars,          /**< number of linking variables */
   SCIP_SOL**            newsol,             /**< pointer to store improved solution */
   SCIP_Bool*            success             /**< pointer to store whether reoptimization was successful */
   )
{
   SCIP* scipcopy;
   SCIP_SOL* startsol;
   SCIP_HASHMAP* varmap;
   SCIP_VAR** subvars;
   SCIP_STATUS status;
   SCIP_Real* linkvals;
   SCIP_Real time;
   int v;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(sol != NULL);
   assert(vars != NULL);
   assert(linkvars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &linkvals, nlinkvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   SCIP_CALL( SCIPgetSolVals(scip, sol, nlinkvars, linkvars, linkvals) );

   /* initializing the problem copy*/
   SCIP_CALL( SCIPcreate(&scipcopy) );

   /* - create the variable mapping hash map
    * - create a problem copy of main SCIP
    */
   if( SCIPheurGetData(heur)->original )
   {
      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scipcopy), SCIPgetNOrigVars(scip)) );
      SCIP_CALL( SCIPcopyOrigConsCompression(scip, scipcopy, varmap, NULL, "reopt_padm", linkvars, linkvals, nlinkvars,
               FALSE, FALSE, TRUE, success) );
   }
   else
   {
      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scipcopy), SCIPgetNVars(scip)) );
      SCIP_CALL( SCIPcopyConsCompression(scip, scipcopy, varmap, NULL, "reopt_padm", linkvars, linkvals, nlinkvars,
               TRUE, FALSE, FALSE, TRUE, success) );
   }
   for( v = 0; v < nvars; v++ )
   {
      subvars[v] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[v]);
   }

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(scipcopy, "misc/catchctrlc", FALSE) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(scipcopy, TRUE) );

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(scipcopy, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(scipcopy, "display/freq", 100000000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(scipcopy, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(scipcopy, "timing/statistictiming", FALSE) );
#endif

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(scipcopy, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving but enable components presolver */
   SCIP_CALL( SCIPsetPresolving(scipcopy, SCIP_PARAMSETTING_FAST, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scipcopy, "constraints/components/maxprerounds", 1) );
   SCIP_CALL( SCIPsetLongintParam(scipcopy, "constraints/components/nodelimit", 0LL) );

   /* disable expensive techniques */
   SCIP_CALL( SCIPsetIntParam(scipcopy, "misc/usesymmetry", 0) );

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(scipcopy, "lp/checkdualfeas", FALSE) );

   /* add heuristic solution as start solution */
   SCIP_CALL( SCIPtransformProb(scipcopy) );
   *success = FALSE;
   if( SCIPisLT(scip, SCIPgetSolOrigObj(scip, sol), SCIPinfinity(scip)) )
   {
      SCIP_CALL( SCIPcreateSol(scipcopy, &startsol, heur) );
      for( v = 0; v < nvars; v++ )
      {
         SCIP_VAR* subvar;
         subvar = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[v]);
         if( subvar != NULL )
         {
            SCIP_CALL( SCIPsetSolVal(scipcopy, startsol, subvar, SCIPgetSolVal(scip, sol, vars[v])) );
         }
      }

      SCIP_CALL( SCIPtrySolFree(scipcopy, &startsol, FALSE, FALSE, FALSE, FALSE, FALSE, success) );
      if( *success )
         SCIPdebugMsg(scip, "set start solution\n");
      else
         SCIPdebugMsg(scip, "start solution for reoptimizing is not feasible\n");
   }

   /* set limits; do not use more time than the heuristic has already used */
   SCIP_CALL( SCIPcopyLimits(scip, scipcopy) );
   SCIP_CALL( SCIPsetLongintParam(scipcopy, "limits/nodes", 1LL) );
   SCIP_CALL( SCIPgetRealParam(scipcopy, "limits/time", &time) );
   if( SCIPheurGetTime(heur) < time - 1.0 )
   {
      SCIP_CALL( SCIPsetRealParam(scipcopy, "limits/time", SCIPheurGetTime(heur) + 1.0) );
   }
   if( *success )
   {
      SCIP_CALL( SCIPsetIntParam(scipcopy, "limits/bestsol", 2) ); /* first solution is start solution */
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scipcopy, "limits/bestsol", 1) );
   }

   /* reoptimize problem */
   SCIP_CALL_ABORT( SCIPsolve(scipcopy) );
   status = SCIPgetStatus(scipcopy);

   /* copy solution to main scip */
   if( status == SCIP_STATUS_BESTSOLLIMIT || status == SCIP_STATUS_OPTIMAL )
   {
      SCIP_SOL* solcopy;

      solcopy = SCIPgetBestSol(scipcopy);
      SCIP_CALL( SCIPtranslateSubSol(scip, scipcopy, solcopy, heur, subvars, newsol) );

      SCIPdebugMsg(scip, "Objective value of reoptimized solution %.2f\n", SCIPgetSolOrigObj(scip, *newsol));
      *success = TRUE;
   }
   else
   {
      *success = FALSE;
   }

   /* free data */
   SCIPhashmapFree(&varmap);
   SCIP_CALL( SCIPfree(&scipcopy) );
   SCIPfreeBufferArray(scip, &subvars);
   SCIPfreeBufferArray(scip, &linkvals);

   return SCIP_OKAY;
}

/** rescales the penalty parameters
 *
 *  A sigmoid function is a function with an "S"-shaped graph, e.g. S(x) = x/(1+|x|).
 *  In order to avoid numerical instabilities due to large penalty parameters we rescale them
 *  using the sigmoid function
 *  S(x) = (x - shift)/(flatness + |x - shift|) * (range/2) + offset.
 *  The parameters are mapped into the more controllable interval [lowestslack, range + lowestslack].
 */
static
SCIP_RETCODE scalePenalties(
   PROBLEM*              problem,            /**< block structure */
   SET*                  linkvartoblocks,    /**< linking variable to blocks set */
   SET*                  blocktolinkvars,    /**< block to linking variable set */
   SCIP_HASHTABLE*       htable,             /**< hashtable containing blockinfo*/
   SCIP_Real             maxpenalty          /**< maximum penalty parameter */
   )
{
   SCIP_Real shift;
   SCIP_Real lowestslack;
   SCIP_Real range;
   SCIP_Real offset;
   SCIP_Real flatness;
   int b;
   int i;
   int k;

   shift = maxpenalty / 2.0;
   lowestslack = 0.1;
   range = 10.0;
   offset = range / 2.0 + lowestslack;
   flatness = maxpenalty / 10.0;

   for( b = 0; b < problem->nblocks; b++ )
   {
      for( i = 0; i < blocktolinkvars[b].size; i++ )
      {
         int linkvaridx;
         linkvaridx = blocktolinkvars[b].indexes[i];

         for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
         {
            int b2;
            b2 = linkvartoblocks[linkvaridx].indexes[k];

            if( b2 != b )
            {
               BLOCKINFO binfo;
               BLOCKINFO* binfoout;
               SCIP_Real oldcoeff;

               binfo.block = b;
               binfo.otherblock = b2;
               binfo.linkvaridx = linkvaridx;
               binfoout = (BLOCKINFO*) SCIPhashtableRetrieve(htable, (void*) &binfo);
               assert(binfoout != NULL);

               /* scale coefficient of positive slack variable */
               oldcoeff = binfoout->slackposobjcoeff;
               binfoout->slackposobjcoeff = ((oldcoeff - shift) / (flatness + REALABS(oldcoeff - shift))) * range / 2.0 + offset;

               /* scale coefficient of negative slack variable */
               oldcoeff = binfoout->slacknegobjcoeff;
               binfoout->slacknegobjcoeff = ((oldcoeff - shift) / (flatness + REALABS(oldcoeff - shift))) * range / 2.0 + offset;
            }
         }
      }
   }
   return SCIP_OKAY;
}

/** returns the available time limit that is left */
static
SCIP_RETCODE getTimeLeft(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            time                /**< pointer to store remaining time */
   )
{
   SCIP_Real timelim;
   SCIP_Real solvingtime;

   assert(scip != NULL);

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelim) );
   solvingtime = SCIPgetSolvingTime(scip);

   if( !SCIPisInfinity(scip, timelim) )
      *time = MAX(0.0, (timelim - solvingtime));
   else
      *time = SCIPinfinity(scip);

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyPADM)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurPADM(scip) );

   return SCIP_OKAY;
}


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreePADM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   SCIPheurSetData(heur, NULL);

   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static SCIP_DECL_HEUREXEC(heurExecPADM)
{  /*lint --e{715}*/
   char name[SCIP_MAXSTRLEN];
   char info[SCIP_MAXSTRLEN];
   SCIP_HEURDATA* heurdata;
   PROBLEM* problem;
   SCIP_DECOMP** alldecomps;
   SCIP_DECOMP* decomp;
   SCIP_DECOMP* assigneddecomp;
   SCIP_VAR** vars;
   SCIP_VAR** linkvars;
   SCIP_VAR** tmpcouplingvars;
   SCIP_CONS** conss;
   SCIP_CONS** sortedconss;
   SET* linkvartoblocks;
   SET* blocktolinkvars;
   BLOCKINFO* blockinfolist;
   SCIP_HASHTABLE* htable;
   int* varlabels;
   int* conslabels;
   int* consssize;
   int* alllinkvartoblocks; /* for efficient memory allocation */
   SCIP_Bool* varonlyobj;
   SCIP_Real* tmpcouplingcoef;
   SCIP_Real gap;
   SCIP_Real maxpenalty;
   SCIP_Real slackthreshold;
   SCIP_Real memory; /* in MB */
   SCIP_Real timeleft;
   SCIP_STATUS status;
   SCIP_Bool solutionsdiffer;
   SCIP_Bool solved;
   SCIP_Bool doscaling;
   SCIP_Bool istimeleft;
   SCIP_Bool success;
   SCIP_Bool avoidmemout;
   SCIP_Bool disablemeasures;
   int maxgraphedge;
   int ndecomps;
   int nconss;
   int nvars;
   int nblocks;
   int numlinkvars;
   int nentries;
   int aiter;
   int piter;
   int increasedslacks;
   int blockinfolistfill;
   SCIP_Longint nodesleft;
   int i;
   int b;
   int k;
   int j;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   problem = NULL;
   assigneddecomp = NULL;
   sortedconss = NULL;
   varlabels = NULL;
   conslabels = NULL;
   consssize = NULL;
   alllinkvartoblocks = NULL;
   linkvars = NULL;
   linkvartoblocks = NULL;
   blocktolinkvars = NULL;
   tmpcouplingvars = NULL;
   tmpcouplingcoef = NULL;
   varonlyobj = NULL;
   blockinfolist = NULL;
   htable = NULL;

   nodesleft = heurdata->maxnodes;
   gap = heurdata->gap;

   if( (heurtiming & SCIP_HEURTIMING_BEFORENODE) && heurdata->timing !=1 )
   {
      SCIPdebugMsg(scip, "Initialize padm heuristic before node\n");
   }
   else if( (heurtiming & SCIP_HEURTIMING_AFTERNODE) && heurdata->timing >=1 )
   {
      SCIPdebugMsg(scip, "Initialize padm heuristic after node\n");
   }
   else
   {
      return SCIP_OKAY;
   }

#ifdef PADM_WRITE_PROBLEMS
   SCIP_CALL( SCIPwriteOrigProblem(scip, "orig_problem", NULL, FALSE) );
   SCIP_CALL( SCIPwriteTransProblem(scip, "trans_problem", NULL, FALSE) );
#endif

   /* take original problem (This is only for testing and not recommended!) */
   if( heurdata->original )
   {
      /* multiaggregation of variables has to be switched off */
      if( !SCIPdoNotMultaggr(scip) )
      {
         SCIPwarningMessage(scip, "Heuristic %s does not support multiaggregation when the original problem is used.\nPlease turn multiaggregation off to use this feature.\n", HEUR_NAME);
         return SCIP_OKAY;
      }

      SCIPgetDecomps(scip, &alldecomps, &ndecomps, TRUE);
      if( ndecomps == 0)
         return SCIP_OKAY;

      /* it takes the first decomposition */
      decomp = alldecomps[0];
      SCIPdebugMsg(scip, "First original decomposition is selected\n");
      assert(decomp != NULL);

      nconss = SCIPgetNOrigConss(scip);
      conss = SCIPgetOrigConss(scip);
      nvars = SCIPgetNOrigVars(scip);
      vars = SCIPgetOrigVars(scip);
   }
   /* take transformed problem */
   else
   {
      SCIPgetDecomps(scip, &alldecomps, &ndecomps, FALSE);
      if( ndecomps == 0)
         return SCIP_OKAY;

      /* it takes the first decomposition */
      decomp = alldecomps[0];
      SCIPdebugMsg(scip, "First transformed decomposition is selected\n");
      assert(decomp != NULL);

      nconss = SCIPgetNConss(scip);
      conss = SCIPgetConss(scip);
      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);
   }

   nblocks = SCIPdecompGetNBlocks(decomp);

   /* if problem has no constraints, no variables or less than two blocks, return */
   if( nconss == 0 || nvars == 0 || nblocks <= 1 )
   {
      SCIPdebugMsg(scip, "problem has no constraints, no variables or less than two blocks\n");
      goto TERMINATE;
   }

   /* estimate required memory for all blocks and terminate if not enough memory is available */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memory) );
   SCIP_CALL( SCIPgetBoolParam(scip, "misc/avoidmemout", &avoidmemout) );
   if( avoidmemout && (((SCIPgetMemUsed(scip) + SCIPgetMemExternEstim(scip))/1048576.0) * nblocks >= memory) )
   {
      SCIPdebugMsg(scip, "The estimated memory usage for %d blocks is too large.\n", nblocks);
      goto TERMINATE;
   }

   /* we do not need the block decomposition graph and expensive measures of the decomposition statistics */
   SCIP_CALL( SCIPgetIntParam(scip, "decomposition/maxgraphedge", &maxgraphedge) );
   if( !SCIPisParamFixed(scip, "decomposition/maxgraphedge") )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "decomposition/maxgraphedge", 0) );
   }
   SCIP_CALL( SCIPgetBoolParam(scip, "decomposition/disablemeasures", &disablemeasures) );
   if( !SCIPisParamFixed(scip, "decomposition/disablemeasures") )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "decomposition/disablemeasures", TRUE) );
   }

   /* don't change problem by sorting constraints */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sortedconss, conss, nconss) );

   SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consssize, nblocks + 1) );

   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);
   SCIPdecompGetVarsLabels(decomp, vars, varlabels, nvars);

   /* sort constraints by blocks */
   SCIPsortIntPtr(conslabels, (void**)sortedconss, nconss);

   /* try to assign linking constraints */
   if( heurdata->assignlinking && conslabels[0] == SCIP_DECOMP_LINKCONS )
   {
      /* create new decomposition; don't change the decompositions in the decompstore */
      SCIP_CALL( SCIPcreateDecomp(scip, &assigneddecomp, nblocks, heurdata->original, SCIPdecompUseBendersLabels(decomp)) );

      SCIP_CALL( assignLinking(scip, assigneddecomp, vars, sortedconss, varlabels, conslabels, nvars, nconss, SCIPdecompGetNBorderConss(decomp)) );
      assert(SCIPdecompGetNBlocks(decomp) >= SCIPdecompGetNBlocks(assigneddecomp));
      decomp = assigneddecomp;

      /* number of blocks can get smaller (since assigning constraints can lead to empty blocks) */
      nblocks = SCIPdecompGetNBlocks(decomp);
   }
   else
   {
      /* The decomposition statistics were computed during transformation of the decomposition store.
       * Since propagators can have changed the number of constraints/variables,
       * the statistics are no longer up-to-date and have to be recomputed.
       */
      SCIP_CALL( SCIPcomputeDecompStats(scip, decomp, TRUE) );
      nblocks = SCIPdecompGetNBlocks(decomp);
   }

   /* reset parameters */
   SCIP_CALL( SCIPsetIntParam(scip, "decomposition/maxgraphedge", maxgraphedge) );
   SCIP_CALL( SCIPsetBoolParam(scip, "decomposition/disablemeasures", disablemeasures) );

   /* @note the terms 'linking' and 'border' (constraints/variables) are used interchangeably */

   if( SCIPdecompGetNBorderConss(decomp) != 0 )
   {
      SCIPdebugMsg(scip, "No support for linking contraints\n");
      goto TERMINATE;
   }

   /* get number of linking variables */
   numlinkvars = SCIPdecompGetNBorderVars(decomp);
   SCIPdebugMsg(scip, "%d linking variables\n", numlinkvars);

   if( numlinkvars == 0 )
   {
      SCIPdebugMsg(scip, "No linking variables\n");
      goto TERMINATE;
   }

   *result = SCIP_DIDNOTFIND;

   /* get for every block the number of constraints (first entry belongs to border) */
   SCIP_CALL( SCIPdecompGetConssSize(decomp, consssize, nblocks + 1) );

   /* create blockproblems */
   SCIP_CALL( createAndSplitProblem(scip, sortedconss, nconss, consssize, nblocks, &problem, heurdata->original, &success) );

   if( !success )
   {
      SCIPdebugMsg(scip, "Some subscips could not be created successfully.\n");
      goto TERMINATE;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &linkvartoblocks, numlinkvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocktolinkvars, problem->nblocks) );

   /* set pointer to NULL for safe memory release */
   for( i = 0; i < numlinkvars; i++ )
      linkvartoblocks[i].indexes = NULL;
   for( i = 0; i < problem->nblocks; i++ )
      blocktolinkvars[i].indexes = NULL;

   /* extract linking variables and init linking variable to blocks set */
   SCIP_CALL( SCIPallocBufferArray(scip, &alllinkvartoblocks, problem->nblocks * numlinkvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linkvars, numlinkvars) );

   b = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( varlabels[i] == SCIP_DECOMP_LINKVAR )
      {
         linkvars[b] = vars[i];
         linkvartoblocks[b].indexes = &alllinkvartoblocks[b * problem->nblocks];
         linkvartoblocks[b].size = 0;
         b++;
      }
   }

   /* fill linking variable to blocks set */
   for( i = 0; i < numlinkvars; i++ )
   {
      SCIP_VAR* var;
      const char* vname;

      vname = SCIPvarGetName(linkvars[i]);
      k = 0;
      for( b = 0; b < problem->nblocks; b++ )
      {
         var = SCIPfindVar((problem->blocks[b]).subscip, vname);
         if( var != NULL )
         {
            linkvartoblocks[i].indexes[k] = b;
            linkvartoblocks[i].size = k + 1;
            k++;
         }
      }
   }

   /* check whether there is enough time left */
   SCIP_CALL( getTimeLeft(scip, &timeleft) );
   if( timeleft <= 0 )
   {
      SCIPdebugMsg(scip, "no time left\n");
      goto TERMINATE;
   }

   /* init varonlyobj; true if variable is only part of the objective function */
   SCIP_CALL( SCIPallocBufferArray(scip, &varonlyobj, numlinkvars) );
   for( i = 0; i < numlinkvars; ++i)
      varonlyobj[i] = TRUE;

   /* init and fill block to linking variables set */
   for( b = 0; b < problem->nblocks; b++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(blocktolinkvars[b].indexes), numlinkvars) );
      blocktolinkvars[b].size = 0;

      k = 0;
      for( i = 0; i < numlinkvars; i++ )
      {
         SCIP_VAR* var;
         const char* vname;

         vname = SCIPvarGetName(linkvars[i]);
         var = SCIPfindVar((problem->blocks[b]).subscip, vname);
         if( var != NULL )
         {
            varonlyobj[i] = FALSE;
            blocktolinkvars[b].indexes[k] = i;
            blocktolinkvars[b].size = k + 1;
            k++;
         }
      }
   }

   /* init arrays for slack variables and coupling constraints */
   for( b = 0; b < problem->nblocks; b++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->blocks[b]).slackspos), blocktolinkvars[b].size * (nblocks - 1)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->blocks[b]).slacksneg), blocktolinkvars[b].size * (nblocks - 1)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->blocks[b]).couplingcons), blocktolinkvars[b].size * (nblocks - 1)) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcouplingvars, COUPLINGSIZE) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcouplingcoef, COUPLINGSIZE) );
   tmpcouplingcoef[0] = 1.0;
   tmpcouplingcoef[1] = 1.0;
   tmpcouplingcoef[2] = -1.0;

   /* count hashtable entries */
   nentries = 0;
   for( b = 0; b < problem->nblocks; b++ )
   {
      for( i = 0; i < blocktolinkvars[b].size; i++ )
      {
         int linkvaridx = blocktolinkvars[b].indexes[i];
         for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
         {
            if( linkvartoblocks[linkvaridx].indexes[k] != b )
               nentries++;
         }
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &blockinfolist, nentries) );
   SCIP_CALL( SCIPhashtableCreate(&htable, SCIPblkmem(scip), 1, SCIPhashGetKeyStandard, indexesEqual, indexesHashval, (void*) scip) );
   blockinfolistfill = 0;

   /* extend submips */
   SCIPdebugMsg(scip, "Extending %d block models\n", problem->nblocks);
   for( b = 0; b < problem->nblocks; b++ )
   {
      SCIP_VAR** blockvars;
      int nblockvars;

      blockvars = SCIPgetVars((problem->blocks[b]).subscip);
      nblockvars = SCIPgetNVars((problem->blocks[b]).subscip);

      /* set objective function of each block to zero */
      for( i = 0; i < nblockvars; i++ )
      {
         SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, blockvars[i], 0.0) );
      }

      /* add two slack variables for each linking variable in block */
      for( i = 0; i < blocktolinkvars[b].size; i++ )
      {
         int linkvaridx;
         linkvaridx = blocktolinkvars[b].indexes[i];

         for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
         {
            int b2;
            b2 = linkvartoblocks[linkvaridx].indexes[k];

            /* handle different blocks with common linking variable */
            if( b2 != b )
            {
               BLOCKINFO* binfo;
               binfo = &blockinfolist[blockinfolistfill];
               blockinfolistfill++;
               binfo->block = b;
               binfo->otherblock = b2;
               binfo->linkvaridx = linkvaridx;
               binfo->linkvar = SCIPfindVar((problem->blocks[b]).subscip, SCIPvarGetName(linkvars[linkvaridx]));
               j = (problem->blocks[b]).ncoupling;

               /* create positive slack variable */
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackpos_block_%d", SCIPvarGetName(linkvars[linkvaridx]), b2);
               (problem->blocks[b]).slackspos[j] = NULL;
               SCIP_CALL( SCIPcreateVarBasic((problem->blocks[b]).subscip,
                                            &((problem->blocks[b]).slackspos[j]), name,
                                            0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
               SCIP_CALL( SCIPaddVar((problem->blocks[b]).subscip, (problem->blocks[b]).slackspos[j]) );
               assert((problem->blocks[b]).slackspos[j] != NULL);
               binfo->slackposobjcoeff = 1.0;
               binfo->slackposvar = (problem->blocks[b]).slackspos[j];

               /* create negative slack variable */
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackneg_block_%d", SCIPvarGetName(linkvars[linkvaridx]), b2);
               (problem->blocks[b]).slacksneg[j] = NULL;
               SCIP_CALL( SCIPcreateVarBasic((problem->blocks[b]).subscip,
                                            &((problem->blocks[b]).slacksneg[j]), name,
                                            0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
               SCIP_CALL( SCIPaddVar((problem->blocks[b]).subscip, (problem->blocks[b]).slacksneg[j]) );
               assert((problem->blocks[b]).slacksneg[j] != NULL);
               binfo->slacknegobjcoeff = 1.0;
               binfo->slacknegvar = (problem->blocks[b]).slacksneg[j];

               /* fill variables for linking constraint */
               tmpcouplingvars[0] = binfo->linkvar;
               tmpcouplingvars[1] = binfo->slackposvar;
               tmpcouplingvars[2] = binfo->slacknegvar;

               /* create linking constraint */
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_coupling_block_%d",
                            SCIPvarGetName(linkvars[linkvaridx]), b2);
               (problem->blocks[b]).couplingcons[j] = NULL;

               /* create linking constraint with initial side equal to zero (or lower bound of linking variable) */
               if( heurtiming & SCIP_HEURTIMING_BEFORENODE )
               {
                  SCIP_Real initval;
                  SCIP_Real lb;

                  lb = SCIPvarGetLbOriginal(binfo->linkvar);
                  initval = MAX(lb, 0.0);

                  SCIP_CALL( SCIPcreateConsBasicLinear((problem->blocks[b]).subscip, &((problem->blocks[b]).couplingcons[j]),
                        name, COUPLINGSIZE, tmpcouplingvars, tmpcouplingcoef, initval, initval) );

                  /* set initial value of linking variable */
                  binfo->linkvarval = initval;
               }

               /* create linking constraint with initial side equal to LP solution (rounded if variable is integer) */
               if( heurtiming & SCIP_HEURTIMING_AFTERNODE )
               {
                  SCIP_Real initval;

                  initval = SCIPvarGetLPSol(linkvars[linkvaridx]);
                  if( SCIPvarGetType(binfo->linkvar) != SCIP_VARTYPE_CONTINUOUS )
                     initval = SCIPround(scip, initval);

                  SCIP_CALL( SCIPcreateConsBasicLinear((problem->blocks[b]).subscip, &((problem->blocks[b]).couplingcons[j]),
                        name, COUPLINGSIZE, tmpcouplingvars, tmpcouplingcoef, initval, initval) );

                  /* set initial value of linking variable */
                  binfo->linkvarval = initval;
               }

               SCIP_CALL( SCIPaddCons((problem->blocks[b]).subscip, (problem->blocks[b]).couplingcons[j]) );
               assert((problem->blocks[b]).couplingcons[j] != NULL);
               binfo->couplingCons = (problem->blocks[b]).couplingcons[j];

               (problem->blocks[b]).ncoupling++;

               /* feed hashtable */
               SCIP_CALL( SCIPhashtableSafeInsert(htable, (void*) binfo) );
            }
         }
      }
   }
   assert(nentries == SCIPhashtableGetNElements(htable));

#ifdef PADM_WRITE_PROBLEMS
   /* write extended submips */
   for( b = 0; b < problem->nblocks; b++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "extended_block_%d.lp", b);
      SCIP_CALL( SCIPwriteOrigProblem((problem->blocks[b]).subscip, name, NULL, FALSE) );
   }
#endif

   /* determine threshold for penalty coefficients via maximum norm */
   slackthreshold = SCIP_REAL_MIN;
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real obj;

      obj = REALABS(SCIPvarGetObj(vars[i]));
      if( obj > slackthreshold )
         slackthreshold = obj;
   }

   /* ------------------------------------------------------------------------------------------------- */

   /* check whether there is enough time left */
   SCIP_CALL( getTimeLeft(scip, &timeleft) );
   if( timeleft <= 0 )
   {
      SCIPdebugMsg(scip, "no time left\n");
      goto TERMINATE;
   }

   SCIPdebugMsg(scip, "Starting iterations\n");
   SCIPdebugMsg(scip, "PIt\tADMIt\tSlacks\tInfo\n");

   piter = 0;
   increasedslacks = 0;
   (void) SCIPsnprintf(info, SCIP_MAXSTRLEN, "-");
   solved = FALSE;
   istimeleft = TRUE;

   /* Penalty loop */
   while( !solved && piter < heurdata->penaltyiterations && istimeleft )
   {
      piter++;
      solutionsdiffer = TRUE;
      aiter = 0;

      /*  Alternating direction method loop */
      while( solutionsdiffer && aiter < heurdata->admiterations && istimeleft )
      {
         aiter++;
         solutionsdiffer = FALSE;
         SCIPdebugMsg(scip, "%d\t%d\t%d\t%s\n", piter, aiter, increasedslacks, info);

         /* Loop through the blocks and solve each sub-SCIP, potentially multiple times */
         for( b = 0; b < problem->nblocks; b++ )
         {
            for( i = 0; i < blocktolinkvars[b].size; i++ )
            {
               int linkvaridx;
               linkvaridx = blocktolinkvars[b].indexes[i];

               for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
               {
                  int b2;
                  b2 = linkvartoblocks[linkvaridx].indexes[k];

                  if( b2 != b )
                  {
                     BLOCKINFO binfo;
                     BLOCKINFO* binfoout;
                     BLOCKINFO binfo2;
                     BLOCKINFO* binfo2out;

                     SCIP_CONS* couplingcons;
                     SCIP_Real newrhs;

                     binfo.block = b;
                     binfo.otherblock = b2;
                     binfo.linkvaridx = linkvaridx;

                     binfoout = (BLOCKINFO*) SCIPhashtableRetrieve(htable, (void *)&binfo);
                     assert(binfoout != NULL);
                     couplingcons = binfoout->couplingCons;

                     /* interchange blocks b and b2 for getting new right hand side */
                     binfo2.block = b2;
                     binfo2.otherblock = b;
                     binfo2.linkvaridx = linkvaridx;
                     binfo2out = (BLOCKINFO*) SCIPhashtableRetrieve(htable, (void*) &binfo2);
                     assert(binfo2out != NULL);
                     newrhs = binfo2out->linkvarval;

                     /* change side of coupling constraint equation with linking variable value of the other block */
                     SCIP_CALL( SCIPchgLhsLinear((problem->blocks[b]).subscip, couplingcons, newrhs) );
                     SCIP_CALL( SCIPchgRhsLinear((problem->blocks[b]).subscip, couplingcons, newrhs) );

                     /* change penalty coefficients of slack variables */
                     SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slackposvar, binfoout->slackposobjcoeff) );
                     SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slacknegvar, binfoout->slacknegobjcoeff) );
                  }
               }
            }

            /* increase slack penalty coeffs until each subproblem can be solved to optimality */
            do
            {
               SCIP_Longint nnodes;
               int iteration;

#ifdef PADM_WRITE_PROBLEMS
               SCIPdebugMsg(scip, "write subscip of block %d in piter=%d and aiter=%d\n", b, piter, aiter);
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "blockproblem_%d_%d_%d.lp", b, piter, aiter);
               SCIP_CALL( SCIPwriteOrigProblem((problem->blocks[b]).subscip, name, NULL, FALSE) );
#endif

               SCIP_CALL( SCIPsetRealParam((problem->blocks[b]).subscip, "limits/gap", gap) );

               /* reuse old solution if available */
               SCIP_CALL( reuseSolution((problem->blocks[b]).subscip, &problem->blocks[b]) );

               /* update time and memory limit of subproblem */
               SCIP_CALL( SCIPcopyLimits(scip, (problem->blocks[b]).subscip) );

               /* stop if there are not enough nodes left */
               if( nodesleft < heurdata->minnodes )
               {
                  SCIPdebugMsg(scip, "Node limit reached.\n");
                  goto TERMINATE;
               }

               /* update node limit of subproblem
                * in the first iterations we have a smaller node limit
                */
               iteration = ((piter - 1) * heurdata->admiterations) + aiter;
               nnodes = (SCIP_Longint)SCIPceil(scip, (problem->blocks[b]).size * nodesleft * ( 1 - pow(heurdata->nodefac, (double)iteration) ));
               nnodes = MAX( heurdata->minnodes, nnodes );
               SCIP_CALL( SCIPsetLongintParam((problem->blocks[b]).subscip, "limits/nodes", nnodes) );

               /* solve block
                *
                * errors in solving the subproblem should not kill the overall solving process;
                * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
                */
               SCIP_CALL_ABORT( SCIPsolve((problem->blocks[b]).subscip) );
               status = SCIPgetStatus((problem->blocks[b]).subscip);

               /* subtract used nodes from the total nodelimit */
               nodesleft -= (SCIP_Longint)SCIPceil(scip, SCIPgetNNodes((problem->blocks[b]).subscip) * (problem->blocks[b]).size);

               /* check solution if one of the four cases occurs
                * - solution is optimal
                * - solution reached gaplimit
                * - node limit reached with at least one feasible solution
                * - time limit is reached but best solution needs no slack variables (no dual solution available)
                */
               if( status == SCIP_STATUS_OPTIMAL || status == SCIP_STATUS_GAPLIMIT ||
                     (status == SCIP_STATUS_NODELIMIT && SCIPgetNSols((problem->blocks[b]).subscip) > 0) ||
                     (status == SCIP_STATUS_TIMELIMIT && SCIPgetNSols((problem->blocks[b]).subscip) > 0 &&
                     SCIPisEQ(scip, SCIPgetSolOrigObj((problem->blocks[b]).subscip, SCIPgetBestSol((problem->blocks[b]).subscip)), 0.0) ) )
               {
                  SCIPdebugMsg(scip, "Block is optimal or reached gaplimit or nodelimit.\n");

                  if( status == SCIP_STATUS_TIMELIMIT )
                  {
                     SCIPdebugMsg(scip, "Block reached time limit with at least one feasible solution.\n");
                     istimeleft = FALSE;
                  }

                  for( i = 0; i < blocktolinkvars[b].size; i++ )
                  {
                     int linkvaridx;
                     linkvaridx = blocktolinkvars[b].indexes[i];

                     for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
                     {
                        int b2;
                        b2 = linkvartoblocks[linkvaridx].indexes[k];

                        if( b2 != b )
                        {
                           SCIP_SOL* sol;
                           BLOCKINFO binfo;
                           BLOCKINFO* binfoout;
                           SCIP_VAR* var;
                           SCIP_Real val;

                           binfo.block = b;
                           binfo.otherblock = b2;
                           binfo.linkvaridx = linkvaridx;
                           binfoout = (BLOCKINFO *)SCIPhashtableRetrieve(htable, (void *)&binfo);
                           assert(binfoout != NULL);

                           sol = SCIPgetBestSol((problem->blocks[b]).subscip);
                           assert(sol != NULL);
                           var = binfoout->linkvar;
                           val = SCIPgetSolVal((problem->blocks[b]).subscip, sol, var);

                           if( !EPSEQ(binfoout->linkvarval, val, SCIP_DEFAULT_EPSILON) )
                              solutionsdiffer = TRUE;

                           binfoout->linkvarval = val;
                        }
                     }
                  }
               }
               else if( status == SCIP_STATUS_UNBOUNDED )
               {
                  SCIPdebugMsg(scip, "Block is unbounded.\n");
                  for( i = 0; i < blocktolinkvars[b].size; i++ )
                  {
                     int linkvaridx;
                     linkvaridx = blocktolinkvars[b].indexes[i];

                     for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
                     {
                        int b2;
                        b2 = linkvartoblocks[linkvaridx].indexes[k];

                        if( b2 != b )
                        {
                           BLOCKINFO binfo;
                           BLOCKINFO* binfoout;

                           binfo.block = b;
                           binfo.otherblock = b2;
                           binfo.linkvaridx = linkvaridx;
                           binfoout = (BLOCKINFO*) SCIPhashtableRetrieve(htable, (void*) &binfo);
                           assert(binfoout != NULL);

                           /* increase penalty coefficients to obtain a bounded subproblem */
                           binfoout->slackposobjcoeff *= 10.0;
                           binfoout->slacknegobjcoeff *= 10.0;
                           SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slackposvar, binfoout->slackposobjcoeff) );
                           SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slacknegvar, binfoout->slacknegobjcoeff) );
                        }
                     }
                  }
               }
               else if( status == SCIP_STATUS_TIMELIMIT )
               {
                  SCIPdebugMsg(scip, "Block reached time limit. No optimal solution available.\n");
                  goto TERMINATE;
               }
               else
               {
                  SCIPdebugMsg(scip, "Block solving status %d not supported\n", status);
                  goto TERMINATE;
               }

               /* free solving data in order to change problem */
               SCIP_CALL( SCIPfreeTransform((problem->blocks[b]).subscip) );
            }
            while( status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_GAPLIMIT &&
                   !(status == SCIP_STATUS_NODELIMIT && SCIPgetNSols((problem->blocks[b]).subscip) > 0) &&
                   !(status == SCIP_STATUS_TIMELIMIT && SCIPgetNSols((problem->blocks[b]).subscip) > 0 &&
                   SCIPisEQ(scip, SCIPgetSolOrigObj((problem->blocks[b]).subscip, SCIPgetBestSol((problem->blocks[b]).subscip)), 0.0) ) );
         }
      }

      /* check wether problem has been solved and if not update penalty coeffs */
      doscaling = FALSE;
      solved = TRUE;
      increasedslacks = 0;
      maxpenalty = SCIP_REAL_MIN;
      for( b = 0; b < problem->nblocks; b++ )
      {
         for( i = 0; i < blocktolinkvars[b].size; i++ )
         {
            int linkvaridx;
            linkvaridx = blocktolinkvars[b].indexes[i];

            for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
            {
               int b2;
               b2 = linkvartoblocks[linkvaridx].indexes[k];

               if( b2 != b )
               {
                  SCIP_SOL* sol;
                  BLOCKINFO binfo;
                  BLOCKINFO* binfoout;
                  SCIP_Real slackposval;
                  SCIP_Real slacknegval;

                  binfo.block = b;
                  binfo.otherblock = b2;
                  binfo.linkvaridx = linkvaridx;
                  binfoout = (BLOCKINFO*) SCIPhashtableRetrieve(htable, (void*) &binfo);
                  assert(binfoout != NULL);

                  sol = SCIPgetBestSol((problem->blocks[b]).subscip);
                  slackposval = SCIPgetSolVal((problem->blocks[b]).subscip, sol, binfoout->slackposvar);
                  slacknegval = SCIPgetSolVal((problem->blocks[b]).subscip, sol, binfoout->slacknegvar);

                  /* increase penalty coefficient of positive slack variable */
                  if( SCIPisGT(scip, slackposval, 0.0) )
                  {
                     binfoout->slackposobjcoeff *= 10.0;

                     if( binfoout->slackposobjcoeff > slackthreshold )
                        doscaling = TRUE;

                     if( binfoout->slackposobjcoeff > maxpenalty )
                        maxpenalty = binfoout->slackposobjcoeff;

                     solved = FALSE;
                     increasedslacks++;
                  }

                  /* increase penalty coefficient of negative slack variable */
                  if( SCIPisGT(scip, slacknegval, 0.0) )
                  {
                     binfoout->slacknegobjcoeff *= 10.0;

                     if( binfoout->slacknegobjcoeff > slackthreshold )
                        doscaling = TRUE;

                     if( binfoout->slacknegobjcoeff > maxpenalty )
                        maxpenalty = binfoout->slacknegobjcoeff;

                     solved = FALSE;
                     increasedslacks++;
                  }
               }
            }
         }
      }

      /* should sigmoid scaling be applied to the penalty parameters? */
      if( doscaling && heurdata->scaling )
      {
         SCIPdebugMsg(scip, "rescale penalty parameters\n");

         /* reset counter */
         increasedslacks = 0;

         /* rescale penalty parameters */
         SCIP_CALL( scalePenalties(problem, linkvartoblocks, blocktolinkvars, htable, maxpenalty) );
      }

      /* adapt in some cases the gap parameter */
      if( (aiter == 1 && solutionsdiffer == FALSE) || (doscaling && heurdata->scaling) )
      {
         SCIP_Real mingap = 0.001; //todo
         SCIP_Real newgap = MAX(gap * 0.5, mingap);

         if( newgap >= mingap )
         {
            if( doscaling && heurdata->scaling )
               (void) SCIPsnprintf(info, SCIP_MAXSTRLEN, "scale, %f", newgap);
            else
               (void) SCIPsnprintf(info, SCIP_MAXSTRLEN, "%f", newgap);

            gap = newgap;
         }
      }

      /* free solution process data */
      if( !solved )
         for( b = 0; b < problem->nblocks; b++ )
         {
            SCIP_CALL( SCIPfreeTransform((problem->blocks[b]).subscip) );
         }
   }

   /* copy solution if present */
   if( solved )
   {
      SCIP_SOL* newsol;
      SCIP_Real* blocksolvals;

      assert(increasedslacks == 0);

      SCIP_CALL( SCIPallocBufferArray(scip, &blocksolvals, nvars) );
      SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

      for( b = 0; b < problem->nblocks; b++ )
      {
         SCIP_SOL* blocksol;
         SCIP_VAR** blockvars;
         int nblockvars;

         /* get solution of block variables (without slack variables) */
         blocksol = SCIPgetBestSol((problem->blocks[b]).subscip);
         assert(blocksol != NULL);
         blockvars = (problem->blocks[b]).subvars;
         nblockvars = (problem->blocks[b]).nsubvars;
         SCIP_CALL( SCIPgetSolVals((problem->blocks[b]).subscip, blocksol, nblockvars, blockvars, blocksolvals) );

         for( i = 0; i < nblockvars; i++ )
         {
            SCIP_VAR* origvar;
            SCIP_Real solval;

            origvar = SCIPfindVar(scip, SCIPvarGetName(blockvars[i]));
            solval = blocksolvals[i];
            SCIP_CALL_ABORT( SCIPsetSolVal(scip, newsol, origvar, solval) );
         }
      }

      /* treat variables with no constraints; set value of variable to bound */
      for( i = 0; i < numlinkvars; i++ )
      {
         if( varonlyobj[i] )
         {
            SCIP_Real fixedvalue;
            if( SCIPvarGetObj(linkvars[i]) < 0 )
            {
               fixedvalue = SCIPvarGetUbLocal(linkvars[i]);
               if( SCIPisInfinity(scip, fixedvalue) )
                  break; // todo: maybe we should return the status UNBOUNDED instead
            }
            else
            {
               fixedvalue = SCIPvarGetLbLocal(linkvars[i]);
               if( SCIPisInfinity(scip, fixedvalue) )
                  break; // todo: maybe we should return the status UNBOUNDED instead
            }
            SCIP_CALL_ABORT( SCIPsetSolVal(scip, newsol, linkvars[i], fixedvalue) );
         }
      }

      SCIPdebugMsg(scip, "Objective value %.2f\n", SCIPgetSolOrigObj(scip, newsol));

      /* fix linking variables and reoptimize with original objective function */
      if( heurdata->reoptimize )
      {
         SCIP_SOL* improvedsol = NULL;
         SCIP_CALL( reoptimize(scip, heur, newsol, vars, nvars, linkvars, numlinkvars, &improvedsol, &success) );
         assert(improvedsol != NULL || success == FALSE);

         if( success )
         {
            SCIP_CALL( SCIPtrySolFree(scip, &improvedsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
            if( !success )
            {
               SCIPdebugMsg(scip, "Reoptimizing solution failed\n");
            }
            else
            {
               SCIPdebugMsg(scip, "Reoptimizing solution successful\n");
               *result = SCIP_FOUNDSOL;
            }
         }
      }

      /* if reoptimization is turned off or reoptimization found no solution, try initial solution */
      if( *result != SCIP_FOUNDSOL )
      {
         SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
         if( !success )
         {
            SCIPdebugMsg(scip, "Solution copy failed\n");
         }
         else
         {
            SCIPdebugMsg(scip, "Solution copy successful\n");
            *result = SCIP_FOUNDSOL;
         }
      }
      else
      {
         SCIP_CALL( SCIPfreeSol(scip, &newsol) );
      }

      SCIPfreeBufferArray(scip, &blocksolvals);
   }
   else
   {
      SCIPdebugMsg(scip, "maximum number of penalty loops reached\n");
      *result = SCIP_DIDNOTFIND;
   }

TERMINATE:
   /* release variables, constraints and free memory */
   if( problem != NULL )
   {
      for( b = 0; b < problem->nblocks; b++ )
      {
         BLOCK curr_block = problem->blocks[b];
         for( i = 0; i < (problem->blocks[b]).ncoupling; i++ )
         {
            SCIP_CALL( SCIPreleaseCons(curr_block.subscip, &curr_block.couplingcons[i]) );
            SCIP_CALL( SCIPreleaseVar(curr_block.subscip, &curr_block.slackspos[i]) );
            SCIP_CALL( SCIPreleaseVar(curr_block.subscip, &curr_block.slacksneg[i]) );
         }
      }
   }

   if( htable != NULL )
      SCIPhashtableFree(&htable);

   if( blockinfolist != NULL )
      SCIPfreeBufferArray(scip, &blockinfolist);

   if( tmpcouplingcoef != NULL )
      SCIPfreeBufferArray(scip, &tmpcouplingcoef);

   if( tmpcouplingvars != NULL )
      SCIPfreeBufferArray(scip, &tmpcouplingvars);

   if( problem != NULL )
   {
      for( b = problem->nblocks - 1; b >= 0; b-- )
      {
         if( problem->blocks[b].couplingcons != NULL )
         {
            SCIPfreeBufferArray(scip, &problem->blocks[b].couplingcons);
            SCIPfreeBufferArray(scip, &problem->blocks[b].slacksneg);
            SCIPfreeBufferArray(scip, &problem->blocks[b].slackspos);
         }
      }
   }

   if( varonlyobj != NULL )
      SCIPfreeBufferArray(scip, &varonlyobj);

   if( problem != NULL && blocktolinkvars != NULL )
   {
      for( b = problem->nblocks -1; b >= 0; b-- )
      {
         if( blocktolinkvars[b].indexes != NULL )
            SCIPfreeBufferArray(scip, &(blocktolinkvars[b].indexes));
      }
   }

   if( linkvars != NULL )
      SCIPfreeBufferArray(scip, &linkvars);

   if( alllinkvartoblocks != NULL )
      SCIPfreeBufferArray(scip, &alllinkvartoblocks);

   if( blocktolinkvars != NULL )
      SCIPfreeBufferArray(scip, &blocktolinkvars);

   if( linkvartoblocks != NULL )
      SCIPfreeBufferArray(scip, &linkvartoblocks);

   if( assigneddecomp != NULL )
      SCIPfreeDecomp(scip, &assigneddecomp);

   if( consssize != NULL )
      SCIPfreeBufferArray(scip, &consssize);

   if( conslabels != NULL )
      SCIPfreeBufferArray(scip, &conslabels);

   if( varlabels != NULL )
      SCIPfreeBufferArray(scip, &varlabels);

   if( sortedconss != NULL )
      SCIPfreeBufferArray(scip, &sortedconss);

   if( problem != NULL )
   {
      SCIP_CALL( freeProblem(&problem, nblocks) );
   }

   SCIPdebugMsg(scip, "Leave padm heuristic\n");
   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the PADM primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurPADM(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur = NULL;

   /* create PADM primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */

   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
               HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
               HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecPADM, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyPADM) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreePADM) );

   /* add padm primal heuristic parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
      "maximum number of nodes to regard in all subproblems",
      &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
      "minimum number of nodes to regard in one subproblem",
      &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodefac",
      "factor to control nodelimits of subproblems", &heurdata->nodefac, TRUE, DEFAULT_NODEFAC, 0.0, 0.99, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/admiterations",
      "maximal number of ADM iterations in each penalty loop", &heurdata->admiterations, TRUE, DEFAULT_ADMIT, 1, 100, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/penaltyiterations",
      "maximal number of penalty iterations", &heurdata->penaltyiterations, TRUE, DEFAULT_PENALTYIT, 1, 100000, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/gap",
      "mipgap at start", &heurdata->gap, TRUE, DEFAULT_GAP, 0.0, 16.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/reoptimize",
      "should the problem get reoptimized with the original objective function?", &heurdata->reoptimize, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/scaling",
      "enable sigmoid rescaling of penalty parameters", &heurdata->scaling, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/assignlinking",
      "should linking constraints be assigned?", &heurdata->assignlinking, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/original",
      "should the original problem be used? This is only for testing and not recommended!", &heurdata->original, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/timing",
      "should the heuristic run before or after the processing of the node? (0: before, 1: after, 2: both)",
      &heurdata->timing, FALSE, 0, 0, 2, NULL, NULL) );

   return SCIP_OKAY;
}
