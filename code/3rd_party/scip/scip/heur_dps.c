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

/**@file   heur_dps.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  dynamic partition search
 * @author Katrin Halbig
 *
 * The dynamic partition search (DPS) is a construction heuristic which additionally needs a
 * user decomposition with linking constraints only.
 *
 * This heuristic splits the problem into several sub-SCIPs according to the given decomposition. Thereby the linking constraints
 * with their right-hand and left-hand sides are also split. DPS searches for a partition of the sides on the blocks
 * so that a feasible solution is obtained.
 * For each block the parts of the original linking constraints are extended by slack variables. Moreover, the objective function
 * is replaced by the sum of these additional variables weighted by penalty parameters lambda. If all blocks have an optimal solution
 * of zero, the algorithm terminates with a feasible solution for the main problem. Otherwise, the partition and the penalty parameters
 * are updated, and the sub-SCIPs are solved again.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/heur_dps.h"
#include "scip/pub_dcmp.h"
#include "scip/pub_heur.h"
#include "scip/pub_misc.h"
#include "scip/scipdefplugins.h"
#include "scip/scip_cons.h"
#include "scip/scip_dcmp.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"


#define HEUR_NAME             "dps"
#define HEUR_DESC             "primal heuristic for decomposable MIPs"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_LNS
#define HEUR_PRIORITY         75000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXIT         50    /**< maximum number of iterations */
#define DEFAULT_PENALTY       100.0 /**< multiplier for absolute increase of penalty parameters */

/* event handler properties */
#define EVENTHDLR_NAME        "Dps"
#define EVENTHDLR_DESC        "event handler for " HEUR_NAME " heuristic"

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_CONS**           linkingconss;       /**< linking constraints */
   int                   nlinking;           /**< number of linking constraints */
   int                   nblocks;            /**< number of blocks */
   int                   maxit;              /**< maximal number of iterations */
   SCIP_Real             maxlinkscore;       /**< maximal linking score of used decomposition (equivalent to percentage of linking constraints) */
   SCIP_Real             penalty;            /**< multiplier for absolute increase of penalty parameters */
   SCIP_Bool             reoptimize;         /**< should the problem get reoptimized with the original objective function? */
   SCIP_Bool             reuse;              /**< should solutions get reused in subproblems? */
};

/** data related to one block */
struct Blockproblem
{
   SCIP*                 blockscip;          /**< SCIP data structure */
   SCIP_VAR**            slackvars;          /**< slack variables */
   SCIP_CONS**           linkingconss;       /**< linking constraints */
   int*                  linkingindices;     /**< indices of linking constraints in original problem */
   int                   nlinking;           /**< number of linking constraints */
   int                   nblockvars;         /**< number of block variables */
   int                   nslackvars;         /**< number of slack variables */
   SCIP_Real*            origobj;            /**< original objective coefficients */
};
typedef struct Blockproblem BLOCKPROBLEM;

/** data related to one linking constraint */
struct Linking
{
   SCIP_CONS*            linkingcons;        /**< corresponding linking constraint of original problem */
   SCIP_CONS**           blockconss;         /**< linking constraints of the blocks */
   SCIP_VAR**            slacks;             /**< slackvars of the blocks */
   SCIP_Real*            minactivity;        /**< minimal activity of constraint for each block */
   SCIP_Real*            maxactivity;        /**< maximal activity of constraint for each block */
   SCIP_Real*            currentrhs;         /**< current partition of rhs */
   SCIP_Real*            currentlhs;         /**< current partition of lhs */
   int*                  blocknumbers;       /**< number of the blocks */
   int                   nblocks;            /**< number of blocks in which this linking constraint participates; dimension of arrays */
   int                   nslacks;            /**< number of slack variables */
   int                   nslacksperblock;    /**< 2, if ranged constraint; 1, if only rhs or lhs */
   int                   lastviolations;     /**< number of iterations in which the constraint was violated in succession */
   SCIP_Bool             hasrhs;             /**< has linking constraint finite right-hand side? */
   SCIP_Bool             haslhs;             /**< has linking constraint finite left-hand side? */
};
typedef struct Linking LINKING;

/*
 * Local methods
 */

/** assigns linking variables to last block
 *
 * The labels are copied to newdecomp and the linking variables are assigned to the last block (i.e., highest block label).
 * Constraint labels and statistics are recomputed.
 */
static
SCIP_RETCODE assignLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          newdecomp,          /**< decomposition with assigned linking variables */
   SCIP_VAR**            vars,               /**< sorted array of variables */
   SCIP_CONS**           conss,              /**< sorted array of constraints */
   int*                  varlabels,          /**< sorted array of variable labels */
   int*                  conslabels,         /**< sorted array of constraint labels */
   int                   nvars,              /**< number of variables */
   int                   nconss,             /**< number of constraints */
   int                   nlinkvars           /**< number of linking variables */
   )
{
   int newlabel;
   int v;

   assert(scip != NULL);
   assert(newdecomp != NULL);
   assert(vars != NULL);
   assert(conss != NULL);
   assert(varlabels != NULL);
   assert(conslabels != NULL);

   /* copy the labels */
   SCIP_CALL( SCIPdecompSetVarsLabels(newdecomp, vars, varlabels, nvars) );
   SCIP_CALL( SCIPdecompSetConsLabels(newdecomp, conss, conslabels, nconss) );

   /* assign linking variables */
   newlabel = varlabels[nvars - 1]; /* take always label of last block */
   assert(newlabel >= 0);
   for( v = 0; v < nlinkvars; v++ )
   {
      SCIP_CALL( SCIPdecompSetVarsLabels(newdecomp, &vars[v], &newlabel, 1) );
   }
   SCIPdebugMsg(scip, "assigned %d linking variables\n", nlinkvars);

   /* recompute constraint labels and statistics */
   SCIP_CALL( SCIPcomputeDecompConsLabels(scip, newdecomp, conss, nconss) );
   SCIP_CALL( SCIPcomputeDecompStats(scip, newdecomp, TRUE) );
   nlinkvars = SCIPdecompGetNBorderVars(newdecomp);

   /* get new labels and sort */
   SCIPdecompGetConsLabels(newdecomp, conss, conslabels, nconss);
   SCIPdecompGetVarsLabels(newdecomp, vars, varlabels, nvars);
   SCIPsortIntPtr(conslabels, (void**)conss, nconss);
   SCIPsortIntPtr(varlabels, (void**)vars, nvars);

   /* After assigning the linking variables, blocks can have zero constraints.
    * So the remaining variables are labeled as linking in SCIPcomputeDecompStats().
    * We assign this variables to the same label as above.
    */
   if( nlinkvars >= 1 )
   {
      assert(varlabels[0] == SCIP_DECOMP_LINKVAR);
      SCIPdebugMsg(scip, "assign again %d linking variables\n", nlinkvars);

      for( v = 0; v < nlinkvars; v++ )
      {
         SCIP_CALL( SCIPdecompSetVarsLabels(newdecomp, &vars[v], &newlabel, 1) );
      }
      SCIP_CALL( SCIPcomputeDecompConsLabels(scip, newdecomp, conss, nconss) );
      SCIP_CALL( SCIPcomputeDecompStats(scip, newdecomp, TRUE) );

      SCIPdecompGetConsLabels(newdecomp, conss, conslabels, nconss);
      SCIPdecompGetVarsLabels(newdecomp, vars, varlabels, nvars);
      SCIPsortIntPtr(conslabels, (void**)conss, nconss);
      SCIPsortIntPtr(varlabels, (void**)vars, nvars);
   }
   assert(varlabels[0] != SCIP_DECOMP_LINKVAR);

   return SCIP_OKAY;
}

/** creates a sub-SCIP and sets parameters */
static
SCIP_RETCODE createSubscip(
    SCIP*                scip,               /**< main SCIP data structure */
    SCIP**               subscip             /**< pointer to store created sub-SCIP */
   )
{
   assert(scip != NULL);
   assert(subscip != NULL);

   /* create a new SCIP instance */
   SCIP_CALL( SCIPcreate(subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(*subscip) );

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

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(*subscip, "lp/checkdualfeas", FALSE) );

   return SCIP_OKAY;
}

/** copies the given variables and constraints to the given sub-SCIP */
static
SCIP_RETCODE copyToSubscip(
   SCIP*                 scip,               /**< source SCIP */
   SCIP*                 subscip,            /**< target SCIP */
   const char*           name,               /**< name for copied problem */
   SCIP_VAR**            vars,               /**< array of variables to copy */
   SCIP_CONS**           conss,              /**< array of constraints to copy */
   SCIP_HASHMAP*         varsmap,            /**< hashmap for copied variables */
   SCIP_HASHMAP*         conssmap,           /**< hashmap for copied constraints */
   int                   nvars,              /**< number of variables to copy */
   int                   nconss,             /**< number of constraints to copy */
   SCIP_Bool*            success             /**< was copying successful? */
   )
{
   SCIP_CONS* newcons;
   SCIP_VAR* newvar;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(vars != NULL);
   assert(conss != NULL);
   assert(varsmap != NULL);
   assert(conssmap != NULL);
   assert(success != NULL);

   SCIPdebugMsg(scip, "copyToSubscip\n");

   /* create problem in sub-SCIP */
   SCIP_CALL( SCIPcreateProb(subscip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* copy variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPgetVarCopy(scip, subscip, vars[i], &newvar, varsmap, conssmap, FALSE, success) );

      /* abort if variable was not successfully copied */
      if( !(*success) )
      {
         SCIPwarningMessage(scip, "Abort heuristic dps since not all variables were successfully copied.\n");
         SCIPABORT();
         return SCIP_OKAY;
      }
   }
   assert(nvars == SCIPgetNOrigVars(subscip));

   /* copy constraints */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);
      assert(!SCIPconsIsModifiable(conss[i]));
      assert(SCIPconsIsActive(conss[i]));
      assert(!SCIPconsIsDeleted(conss[i]));

      /* copy the constraint */
      SCIP_CALL( SCIPgetConsCopy(scip, subscip, conss[i], &newcons, SCIPconsGetHdlr(conss[i]), varsmap, conssmap, NULL,
            SCIPconsIsInitial(conss[i]), SCIPconsIsSeparated(conss[i]), SCIPconsIsEnforced(conss[i]),
            SCIPconsIsChecked(conss[i]), SCIPconsIsPropagated(conss[i]), FALSE, FALSE,
            SCIPconsIsDynamic(conss[i]), SCIPconsIsRemovable(conss[i]), FALSE, FALSE, success) );

      /* abort if constraint was not successfully copied */
      if( !(*success) )
         return SCIP_OKAY;

      SCIP_CALL( SCIPaddCons(subscip, newcons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &newcons) );
   }

   /* block constraint contains variables which are not part of this block
    *
    * todo: maybe they are part of the block, but it is not recognized, because they are, for example, negated or aggregated.
    */
   if( nvars != SCIPgetNOrigVars(subscip) )
      *success = FALSE;

   return SCIP_OKAY;
}

/** creates the subscip for a given block */
static
SCIP_RETCODE createBlockproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   BLOCKPROBLEM*         blockproblem,       /**< blockproblem that should be created */
   LINKING**             linkings,           /**< linkings that will be (partially) initialized */
   SCIP_CONS**           conss,              /**< sorted array of constraints of this block */
   SCIP_VAR**            vars,               /**< sorted array of variables of this block */
   int                   nconss,             /**< number of constraints of this block */
   int                   nvars,              /**< number of variables of this block */
   SCIP_CONS**           linkingconss,       /**< linking constraints in the original problem */
   int                   nlinking,           /**< number of linking constraints in the original problem */
   int                   blocknumber,        /**< number of block that should be created */
   SCIP_Bool*            success             /**< pointer to store whether creation was successful */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_HASHMAP* varsmap;
   SCIP_HASHMAP* conssmap;
   SCIP_VAR** consvars; /* all vars in original linking cons */
   SCIP_Real* consvals;
   int nconsvars;
   SCIP_VAR** blockvars; /* vars of current linking cons of current block */
   SCIP_Real* blockvals;
   int nblockvars;
   SCIP_VAR** subvars; /* all vars of subscip */
   int maxnconsvars; /* current size of arrays */
   int c;
   int v;

   assert(scip != NULL);
   assert(blockproblem != NULL);
   assert(conss != NULL);
   assert(vars != NULL);
   assert(blockproblem->blockscip != NULL);

   maxnconsvars = 20; /* start size; increase size if necessary */

   SCIPdebugMsg(scip, "Create blockproblem %d\n", blocknumber);

   /* create the variable/constraint mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varsmap, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&conssmap, SCIPblkmem(scip), nconss) );

   /* get name of the original problem and add "comp_nr" */
   (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d", SCIPgetProbName(scip), blocknumber);

   SCIP_CALL( copyToSubscip(scip, blockproblem->blockscip, name, vars, conss, varsmap, conssmap,
                            nvars, nconss, success) );
   if( !(*success) )
   {
      SCIPdebugMsg(scip, "Copy to subscip failed\n");
      SCIPhashmapFree(&conssmap);
      SCIPhashmapFree(&varsmap);

      return SCIP_OKAY;
   }

   /* save number of variables that have a corresponding variable in original problem*/
   blockproblem->nblockvars = SCIPgetNVars(blockproblem->blockscip);

   /* save original objective and set objective to zero */
   subvars = SCIPgetVars(blockproblem->blockscip);
   for( v = 0; v < nvars; v++ )
   {
      blockproblem->origobj[v] = SCIPvarGetObj(subvars[v]);
      SCIP_CALL( SCIPchgVarObj(blockproblem->blockscip, subvars[v], 0.0) );
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(blockproblem->blockscip, &blockvars, nvars + 2) ); /* two entries for the slack variables */
   SCIP_CALL( SCIPallocBufferArray(blockproblem->blockscip, &blockvals, nvars + 2) );
   SCIP_CALL( SCIPallocBufferArray(blockproblem->blockscip, &consvars, maxnconsvars) );
   SCIP_CALL( SCIPallocBufferArray(blockproblem->blockscip, &consvals, maxnconsvars) );

   /* find and add parts of linking constraints */
   SCIPdebugMsg(scip, "add parts of linking constraints\n");
   for( c = 0; c < nlinking; c++ )
   {
      const char* conshdlrname;
      char consname[SCIP_MAXSTRLEN];
      SCIP_CONS* newcons;
      SCIP_Real rhs;
      SCIP_Real lhs;
      SCIP_Real minact;
      SCIP_Real maxact;
      SCIP_Bool mininfinite;
      SCIP_Bool maxinfinite;

      assert(linkingconss[c] != NULL);

      newcons = NULL;

#ifdef SCIP_MORE_DEBUG
      SCIPdebugMsg(scip, "consider constraint %s\n", SCIPconsGetName(linkingconss[c]));
      SCIPdebugPrintCons(scip, linkingconss[c], NULL);
#endif

      nblockvars = 0;

      /* every constraint with linear representation is allowed */
      conshdlrname = SCIPconshdlrGetName(SCIPconsGetHdlr(linkingconss[c]));
      if( !( (strcmp(conshdlrname, "linear") == 0) || (strcmp(conshdlrname, "setppc") == 0)
            || (strcmp(conshdlrname, "logicor") == 0) || (strcmp(conshdlrname, "knapsack") == 0)
            || (strcmp(conshdlrname, "varbound") == 0) ) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "Heuristic %s cannot handle linking constraints of type %s\n", HEUR_NAME, conshdlrname);
         /* TODO which other types can we handle/transform in a linear constraint? */

         *success = FALSE;
         break; /* releases memory and breaks heuristic */
      }

      SCIP_CALL( SCIPgetConsNVars(scip, linkingconss[c], &nconsvars, success) );

      /* reallocate memory if we have more variables than maxnconsvars */
      if( nconsvars > maxnconsvars )
      {
         int newsize;

         /* calculate new size */
         newsize = SCIPcalcMemGrowSize(scip, MAX(2 * maxnconsvars, nconsvars)); /* at least double size */

         SCIP_CALL( SCIPreallocBufferArray(blockproblem->blockscip, &consvars, newsize) );
         SCIP_CALL( SCIPreallocBufferArray(blockproblem->blockscip, &consvals, newsize) );
         maxnconsvars = newsize;
      }

      SCIP_CALL( SCIPgetConsVars(scip, linkingconss[c], consvars, nconsvars, success) );
      SCIP_CALL( SCIPgetConsVals(scip, linkingconss[c], consvals, nconsvars, success) );

      if( !(*success) )
      {
         SCIPdebugMsg(scip, "Create blockproblem failed\n");
         break; /* releases memory and breaks heuristic */
      }

      /* check if constraint contains variables of this block */
      for( v = 0; v < nconsvars; v++ )
      {
         if( SCIPhashmapExists(varsmap, (void*)consvars[v]) )
         {
            blockvars[nblockvars] = (SCIP_VAR*) SCIPhashmapGetImage(varsmap, (void*)consvars[v]);
            blockvals[nblockvars] = consvals[v];
            ++nblockvars;
         }
         /* handle negated variables*/
         else if( SCIPvarGetStatus(consvars[v]) == SCIP_VARSTATUS_NEGATED)
         {
            if( SCIPhashmapExists(varsmap, (void*)SCIPvarGetNegationVar(consvars[v])) ) /* negation exists in this block */
            {
               /* save negated variable */
               SCIP_VAR* origblockvar = (SCIP_VAR*) SCIPhashmapGetImage(varsmap, (void*)SCIPvarGetNegationVar(consvars[v]));
               SCIP_VAR* negblockvar = NULL;
               SCIP_CALL( SCIPgetNegatedVar(blockproblem->blockscip, origblockvar, &negblockvar) );
               blockvars[nblockvars] = negblockvar;
               blockvals[nblockvars] = consvals[v];
               ++nblockvars;
            }
         }
      }

      /* continue with next linking constraint if it has no part in current block */
      if( nblockvars == 0 )
         continue;

      /* get rhs and/or lhs */
      rhs = SCIPconsGetRhs(scip, linkingconss[c], success);
      if( !(*success) )
      {
         SCIPdebugMsg(scip, "Create blockproblem failed\n");
         return SCIP_OKAY;
      }
      lhs = SCIPconsGetLhs(scip, linkingconss[c], success);
      if( !(*success) )
      {
         SCIPdebugMsg(scip, "Create blockproblem failed\n");
         return SCIP_OKAY;
      }
      assert(!SCIPisInfinity(scip, rhs) || !SCIPisInfinity(scip, -lhs)); /* at least one side bounded */
      assert(SCIPisLE(scip, lhs, rhs));

      if( !SCIPisInfinity(scip, rhs) )
         linkings[c]->hasrhs = TRUE;
      if( !SCIPisInfinity(scip, -lhs) )
         linkings[c]->haslhs = TRUE;
      if( !SCIPisInfinity(scip, rhs) && !SCIPisInfinity(scip, -lhs))
         linkings[c]->nslacksperblock = 2;
      else
         linkings[c]->nslacksperblock = 1;

      /* add slack variable for rhs */
      if( linkings[c]->hasrhs )
      {
         /* slack variable z_r >= 0 */
         char varname[SCIP_MAXSTRLEN];
         (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z_r_%s", SCIPconsGetName(linkingconss[c]));
         SCIP_CALL( SCIPcreateVarBasic(blockproblem->blockscip, &blockvars[nblockvars], varname,
                                          0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
         blockvals[nblockvars] = -1.0;
         SCIP_CALL( SCIPaddVar(blockproblem->blockscip, blockvars[nblockvars]) );
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMsg(scip, "Add variable %s\n", SCIPvarGetName(blockvars[nblockvars]));
#endif
         linkings[c]->slacks[linkings[c]->nslacks] = blockvars[nblockvars];
         blockproblem->slackvars[blockproblem->nslackvars] = blockvars[nblockvars];
         ++blockproblem->nslackvars;
         ++linkings[c]->nslacks;
         ++nblockvars;
      }

      /* add slack variable for lhs */
      if( linkings[c]->haslhs )
      {
         /* slack variable z_l >= 0 */
         char varname[SCIP_MAXSTRLEN];
         (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z_l_%s", SCIPconsGetName(linkingconss[c]));
         SCIP_CALL( SCIPcreateVarBasic(blockproblem->blockscip, &blockvars[nblockvars], varname,
                                          0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
         blockvals[nblockvars] = 1.0;
         SCIP_CALL( SCIPaddVar(blockproblem->blockscip, blockvars[nblockvars]) );
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMsg(scip, "Add variable %s\n", SCIPvarGetName(blockvars[nblockvars]));
#endif
         linkings[c]->slacks[linkings[c]->nslacks] = blockvars[nblockvars];
         blockproblem->slackvars[blockproblem->nslackvars] = blockvars[nblockvars];
         ++blockproblem->nslackvars;
         ++linkings[c]->nslacks;
         ++nblockvars;
      }

      /* add linking constraint with slack variable */
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s", SCIPconsGetName(linkingconss[c]));
      SCIP_CALL( SCIPcreateConsBasicLinear(blockproblem->blockscip, &newcons, consname, nblockvars, blockvars, blockvals, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(blockproblem->blockscip, newcons) );
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMsg(blockproblem->blockscip, "add constraint %s\n", SCIPconsGetName(newcons));
      SCIPdebugPrintCons(blockproblem->blockscip, newcons, NULL);
#endif

      blockproblem->linkingconss[blockproblem->nlinking] = newcons;
      linkings[c]->blockconss[linkings[c]->nblocks] = newcons;
      linkings[c]->blocknumbers[linkings[c]->nblocks] = blocknumber;
      blockproblem->linkingindices[blockproblem->nlinking] = c;

      /* calculate minimal und maximal activity (exclude slack variables) */
      minact = 0;
      maxact = 0;
      mininfinite = FALSE;
      maxinfinite = FALSE;
      for( v = 0; v < nblockvars - linkings[c]->nslacksperblock && (!mininfinite || !maxinfinite); v++ )
      {
         SCIP_Real lb;
         SCIP_Real ub;
         lb = SCIPvarGetLbGlobal(blockvars[v]);
         ub = SCIPvarGetUbGlobal(blockvars[v]);

         if( blockvals[v] >= 0.0 )
         {
            mininfinite = (mininfinite || SCIPisInfinity(scip, -lb));
            maxinfinite = (maxinfinite || SCIPisInfinity(scip, ub));
            if( !mininfinite )
               minact += blockvals[v] * lb;
            if( !maxinfinite )
               maxact += blockvals[v] * ub;
         }
         else
         {
            mininfinite = (mininfinite || SCIPisInfinity(scip, ub));
            maxinfinite = (maxinfinite || SCIPisInfinity(scip, -lb));
            if( !mininfinite )
               minact += blockvals[v] * ub;
            if( !maxinfinite )
               maxact += blockvals[v] * lb;
         }
      }

      if( mininfinite )
         linkings[c]->minactivity[linkings[c]->nblocks] = -SCIPinfinity(scip);
      else
         linkings[c]->minactivity[linkings[c]->nblocks] = minact;
      if( maxinfinite )
         linkings[c]->maxactivity[linkings[c]->nblocks] = SCIPinfinity(scip);
      else
         linkings[c]->maxactivity[linkings[c]->nblocks] = maxact;
      assert(SCIPisLE(scip, linkings[c]->minactivity[linkings[c]->nblocks], linkings[c]->maxactivity[linkings[c]->nblocks]));

      linkings[c]->nblocks++;
      blockproblem->nlinking++;

      for( v = 1; v <= linkings[c]->nslacksperblock; v++ )
      {
         SCIP_CALL( SCIPreleaseVar(blockproblem->blockscip, &blockvars[nblockvars - v]) );
      }

      SCIP_CALL( SCIPreleaseCons(blockproblem->blockscip, &newcons) );
   }
   assert(blockproblem->nlinking <= nlinking);

   /* free memory */
   SCIPfreeBufferArray(blockproblem->blockscip, &consvals);
   SCIPfreeBufferArray(blockproblem->blockscip, &consvars);
   SCIPfreeBufferArray(blockproblem->blockscip, &blockvals);
   SCIPfreeBufferArray(blockproblem->blockscip, &blockvars);

   SCIPhashmapFree(&conssmap);
   SCIPhashmapFree(&varsmap);

   return SCIP_OKAY;
}

/** creates data structures and splits problem into blocks */
static
SCIP_RETCODE createAndSplitProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   BLOCKPROBLEM**        blockproblem,       /**< array of blockproblem data structures */
   LINKING**             linkings,           /**< array of linking data structures */
   SCIP_VAR**            vars,               /**< sorted array of variables */
   SCIP_CONS**           conss,              /**< sorted array of constraints */
   SCIP_Bool*            success             /**< pointer to store whether splitting was successful */
   )
{
   int* nconssblock;
   int* nvarsblock;
   int conssoffset;
   int varsoffset;
   int i;   /* blocknumber */

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(vars != NULL);
   assert(conss != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &nvarsblock, heurdata->nblocks + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nconssblock, heurdata->nblocks + 1) );
   SCIP_CALL( SCIPdecompGetVarsSize(decomp, nvarsblock, heurdata->nblocks + 1) );
   SCIP_CALL( SCIPdecompGetConssSize(decomp, nconssblock, heurdata->nblocks + 1) );
   assert(0 == nvarsblock[0]);

   varsoffset = 0;
   conssoffset = 0;

   for( i = 0; i < heurdata->nblocks; i++)
   {
      conssoffset += nconssblock[i];
      varsoffset += nvarsblock[i];

      SCIP_CALL( createBlockproblem(scip, blockproblem[i], linkings, &conss[conssoffset], &vars[varsoffset], nconssblock[i+1], nvarsblock[i+1],
                                    heurdata->linkingconss, heurdata->nlinking, i, success) );
      if( !(*success) )
         break;
   }

   SCIPfreeBufferArray(scip, &nconssblock);
   SCIPfreeBufferArray(scip, &nvarsblock);

   return SCIP_OKAY;
}

/** rounds partition for one linking constraint to integer value if variables and coefficients are integer
 *
 *  changes only currentrhs/currentlhs
 */
static
SCIP_RETCODE roundPartition(
   SCIP*                 scip,               /**< SCIP data structure */
   LINKING*              linking,            /**< one linking data structure */
   BLOCKPROBLEM**        blockproblem,       /**< array of blockproblem data structures */
   SCIP_Bool             roundbyrhs          /**< round by right hand side? */
   )
{
   SCIP_Real* fracPart;
   int* sorting;
   int* isinteger;
   SCIP_Real sumbefor; /* includes value at idx */
   SCIP_Real sumafter;
   SCIP_Real diff;
   int nnonintblocks; /* number of non integer blocks */
   int idx;
   int b;
   int i;
   int k;

   assert(scip != NULL);
   assert(linking != NULL);
   assert(blockproblem != NULL);

   nnonintblocks = 0;
   idx = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &fracPart, linking->nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sorting, linking->nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &isinteger, linking->nblocks) );

   /* get integer blocks and fractional parts */
   for( b = 0; b < linking->nblocks; b++ )
   {
      SCIP* subscip;
      SCIP_CONS* blockcons;
      SCIP_VAR** blockvars;
      SCIP_Real* blockvals;
      int nblockvars;
      int length; /* number of block variables without slack variables */
      SCIP_Bool success;

      subscip = blockproblem[linking->blocknumbers[b]]->blockscip;
      blockcons = linking->blockconss[b];
      sorting[b] = b; /* store current sorting to sort back */

      SCIP_CALL( SCIPgetConsNVars(subscip, blockcons, &nblockvars, &success) );
      assert(success);
      SCIP_CALL( SCIPallocBufferArray(scip, &blockvars, nblockvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &blockvals, nblockvars) );

      SCIP_CALL( SCIPgetConsVars(subscip, blockcons, blockvars, nblockvars, &success) );
      assert(success);
      SCIP_CALL( SCIPgetConsVals(subscip, blockcons, blockvals, nblockvars, &success) );
      assert(success);

      /* get number of block variables in this constraint without slack variables */
      length = nblockvars - linking->nslacksperblock;

      /* get blocks with integer value */
      isinteger[b] = 1;
      for( i = 0; i < length; i++ )
      {
         if( SCIPvarGetType(blockvars[i]) == SCIP_VARTYPE_CONTINUOUS ||  !SCIPisIntegral(scip, blockvals[i]) )
         {
            isinteger[b] = 0;
            nnonintblocks++;
            break;
         }
      }

      /* get fractional part of blockconstraints */
      if( roundbyrhs )
         fracPart[b] = linking->currentrhs[b] - floor(linking->currentrhs[b]); /* do not use SCIPfrac()! */
      else
         fracPart[b] = linking->currentlhs[b] - floor(linking->currentlhs[b]);

      SCIPfreeBufferArray(scip, &blockvals);
      SCIPfreeBufferArray(scip, &blockvars);
   }

   /* sort non-integer blocks to the front */
   SCIPsortIntIntReal(isinteger, sorting, fracPart, linking->nblocks);

   /* sort by fractional part */
   SCIPsortRealInt(fracPart, sorting, nnonintblocks);
   SCIPsortRealInt(&fracPart[nnonintblocks], &sorting[nnonintblocks], linking->nblocks - nnonintblocks);

   /* detect blocks for rounding down and rounding up:
    * integer blocks with small fractional parts are rounded down
    * integer blocks with big fractional parts are rounded up
    */

   sumbefor = 0.0;
   sumafter = 0.0;

   for( i = 0; i < linking->nblocks - nnonintblocks; i++ )
      sumafter += 1 - fracPart[nnonintblocks + i];

   for( i = 0; i < linking->nblocks - nnonintblocks; i++ )
   {
      sumbefor += fracPart[nnonintblocks + i];
      sumafter -= 1 - fracPart[nnonintblocks + i];

      if( sumbefor >= sumafter )
      {
         for( k = 0; k <= i; k++ )
            fracPart[nnonintblocks + k] = -fracPart[nnonintblocks + k];

         for( k = i + 1; k < linking->nblocks - nnonintblocks; k++ )
            fracPart[nnonintblocks + k] = 1 - fracPart[nnonintblocks + k];

         idx = i;
         break;
      }
   }
   diff = sumbefor - sumafter;
   assert(SCIPisGE(scip, diff, 0.0));

   /* add difference to last non integer block */
   for( i = nnonintblocks - 1; i >= 0; i-- )
   {
      if( SCIPisGT(scip, diff, 0.0) )
      {
         fracPart[i] = diff;
         diff = 0;
      }
      else
         fracPart[i] = 0;
   }

   /* add difference to last rounded down block if no non integer block exists */
   if( SCIPisGT(scip, diff, 0.0))
   {
      assert(nnonintblocks == 0);
      fracPart[idx] += diff;
   }

   /* sort back */
   SCIPsortIntReal(sorting, fracPart, linking->nblocks);

   /* round partition
    * if we have a ranged constraint, both sides get rounded in the same way
    */
   for( b = 0; b < linking->nblocks; b++ )
   {
      if( linking->hasrhs )
         linking->currentrhs[b] += fracPart[b];

      if( linking->haslhs )
         linking->currentlhs[b] += fracPart[b];
   }

   SCIPfreeBufferArray(scip, &isinteger);
   SCIPfreeBufferArray(scip, &sorting);
   SCIPfreeBufferArray(scip, &fracPart);

   return SCIP_OKAY;
}

/** calculates initial partition and sets rhs/lhs in blockproblems */
static
SCIP_RETCODE initCurrent(
   SCIP*                 scip,               /**< SCIP data structure of main scip */
   LINKING**             linkings,           /**< array of linking data structures */
   BLOCKPROBLEM**        blockproblem,       /**< array of blockproblem data structures */
   int                   nlinking,           /**< number of linking constraints */
   SCIP_Bool*            success             /**< pointer to store whether initialization was successful */
   )
{
   LINKING* linking;
   SCIP_Real rhs;
   SCIP_Real lhs;
   SCIP_Real residualrhs;
   SCIP_Real residuallhs;
   SCIP_Real goalvalue;
   int c;
   int b;

   assert(scip != NULL);
   assert(linkings != NULL);
   assert(blockproblem != NULL);
   assert(nlinking > 0);

   SCIPdebugMsg(scip, "initialize partition\n");

   /* for each linking constraint the rhs/lhs is split between the blocks */
   for( c = 0; c < nlinking; c++ )
   {
      linking = linkings[c];
      rhs = SCIPconsGetRhs(scip, linking->linkingcons, success);
      assert(*success);
      lhs = SCIPconsGetLhs(scip, linking->linkingcons, success);
      assert(*success);
      residualrhs = rhs;
      residuallhs = lhs;

      /* equal parts for each block with respect to minimal/maximal activity */
      if( linking->hasrhs || linking->haslhs )
      {
         if( linking->hasrhs )
         {
            for( b = 0; b < linking->nblocks; b++ )
            {
               goalvalue = residualrhs / (linking->nblocks - b);
               linking->currentrhs[b] = MIN(MAX(goalvalue, linking->minactivity[b]), linking->maxactivity[b]);
               residualrhs -= linking->currentrhs[b];
            }
            /* add residual partition to first block */
            linking->currentrhs[0] += residualrhs;
         }
         if( linking->haslhs )
         {
            for( b = 0; b < linking->nblocks; b++ )
            {
               goalvalue = residuallhs / (linking->nblocks - b);
               linking->currentlhs[b] = MIN(MAX(goalvalue, linking->minactivity[b]), linking->maxactivity[b]);
               residuallhs -= linking->currentlhs[b];
            }
            /* add residual partition to first block */
            linking->currentlhs[0] += residuallhs;
         }
      }
      else
      {
         assert(linking->nblocks == 0 && !SCIPconsIsChecked(linking->linkingcons));
      }

      SCIP_CALL( roundPartition(scip, linking, blockproblem, linking->hasrhs) );

      /* set sides in blockproblem at initial partition */
      for( b = 0; b < linking->nblocks; b++ )
      {
         if( linking->hasrhs )
         {
            SCIP_CALL( SCIPchgRhsLinear(blockproblem[linking->blocknumbers[b]]->blockscip,
                                       linking->blockconss[b], linking->currentrhs[b]) );
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMsg(scip, "change rhs of %s in block %d to %f\n",
                        SCIPconsGetName(linking->linkingcons), linking->blocknumbers[b], linking->currentrhs[b]);
#endif
         }
         if( linking->haslhs )
         {
            SCIP_CALL( SCIPchgLhsLinear(blockproblem[linking->blocknumbers[b]]->blockscip,
                                       linking->blockconss[b], linking->currentlhs[b]) );
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMsg(scip, "change lhs of %s in block %d to %f\n",
                        SCIPconsGetName(linking->linkingcons), linking->blocknumbers[b], linking->currentlhs[b]);
#endif
         }
      }
   }

   return SCIP_OKAY;
}

/** update partition */
static
SCIP_RETCODE updatePartition(
   SCIP*                 scip,               /**< SCIP data structure of main scip */
   LINKING*              linking,            /**< one linking data structure */
   BLOCKPROBLEM**        blockproblem,       /**< array of blockproblem data structures */
   int*                  nviolatedblocksrhs, /**< pointer to store number of blocks which violate rhs */
   int*                  nviolatedblockslhs, /**< pointer to store number of blocks which violate lhs */
   int                   iteration           /**< number of iteration */
   )
{
   SCIP_Real* shift; /* direction vector */
   int* nonviolatedblocksrhs;
   int* nonviolatedblockslhs;
   SCIP_Real sumviols = 0.0;
   int v;

   assert(scip != NULL);
   assert(linking != NULL);
   assert(blockproblem != NULL);
   assert(iteration >= 0);

   *nviolatedblocksrhs = 0;
   *nviolatedblockslhs = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &shift, linking->nblocks) );
   BMSclearMemoryArray(shift, linking->nblocks);
   if( linking->hasrhs )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &nonviolatedblocksrhs, linking->nblocks) );
   }
   if( linking->haslhs )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &nonviolatedblockslhs, linking->nblocks) );
   }

   /* get violated blocks and calculate shift */
   for( v = 0; v < linking->nblocks; v++ )
   {
      SCIP* subscip;
      SCIP_SOL* subsol;
      SCIP_Real slackval;

      subscip = blockproblem[linking->blocknumbers[v]]->blockscip;
      subsol = SCIPgetBestSol(subscip);

      /* if we have ranged constraints, the slack variables of the rhs are in front;
       * get slack variable of block; save violated blocks
       */
      if( linking->hasrhs )
      {
         slackval = SCIPgetSolVal(subscip, subsol, linking->slacks[v * linking->nslacksperblock]);

         /* block is violated */
         if( SCIPisPositive(scip, slackval) )
         {
            (*nviolatedblocksrhs)++;

            shift[v] += slackval;
            sumviols += slackval;
         }
         else
         {
            nonviolatedblocksrhs[v - *nviolatedblocksrhs] = v; /*lint !e644*/
         }
      }
      if( linking->haslhs )
      {
         slackval = SCIPgetSolVal(subscip, subsol, linking->slacks[(v * linking->nslacksperblock) + linking->nslacksperblock - 1]);

         /* block is violated */
         if( SCIPisPositive(scip, slackval) )
         {
            (*nviolatedblockslhs)++;

            shift[v] -= slackval;
            sumviols -= slackval;
         }
         else
         {
            nonviolatedblockslhs[v - *nviolatedblockslhs] = v; /*lint !e644*/
         }
      }
   }

   /* linking constraint is in no block violated or always violated -> do not update partition */
   if( *nviolatedblocksrhs + *nviolatedblockslhs == 0 ||
         linking->nblocks == *nviolatedblocksrhs || linking->nblocks == *nviolatedblockslhs )
   {
      /* free memory */
      if( linking->haslhs )
         SCIPfreeBufferArray(scip, &nonviolatedblockslhs);
      if( linking->hasrhs )
         SCIPfreeBufferArray(scip, &nonviolatedblocksrhs);
      SCIPfreeBufferArray(scip, &shift);
      return SCIP_OKAY;
   }

   /* set values of non violated blocks */
   if( SCIPisPositive(scip, sumviols) )
   {
      /* rhs of original linking constraint is violated */
      SCIP_Real residual = sumviols;
      SCIP_Real part;
      SCIP_Real shift_tmp;

      assert(linking->hasrhs);
      assert(*nviolatedblocksrhs != 0);

      /* substract from each non violated block the same amount with respect to minimal/maximal activity,
       * so that the shift is zero in sum
       */
      for( v = 0; v < linking->nblocks - *nviolatedblocksrhs; v++ )
      {
         part = linking->currentrhs[nonviolatedblocksrhs[v]] - residual/(linking->nblocks - *nviolatedblocksrhs - v);
         part = MIN(MAX(part, linking->minactivity[nonviolatedblocksrhs[v]]), linking->maxactivity[nonviolatedblocksrhs[v]]);
         shift_tmp = part - linking->currentrhs[nonviolatedblocksrhs[v]];
         residual += shift_tmp;
         shift[nonviolatedblocksrhs[v]] += shift_tmp;
      }
      if( !SCIPisZero(scip, residual) )
      {
         /* assign residual to first block */
         shift[nonviolatedblocksrhs[0]] -= residual;
      }
   }
   if( SCIPisNegative(scip, sumviols) )
   {
      /* lhs of original linking constraint is violated */
      SCIP_Real residual = sumviols;
      SCIP_Real part;
      SCIP_Real shift_tmp;

      assert(linking->haslhs);
      assert(*nviolatedblockslhs != 0);

      /* add to each non violated block the same amount with respect to minimal/maximal activity,
       * so that the shift is zero in sum
       */
      for( v = 0; v < linking->nblocks - *nviolatedblockslhs; v++ )
      {
         part = linking->currentlhs[nonviolatedblockslhs[v]] - residual/(linking->nblocks - *nviolatedblockslhs - v);
         part = MIN(MAX(part, linking->minactivity[nonviolatedblockslhs[v]]), linking->maxactivity[nonviolatedblockslhs[v]]);
         shift_tmp = part - linking->currentlhs[nonviolatedblockslhs[v]];
         residual += shift_tmp;
         shift[nonviolatedblockslhs[v]] += shift_tmp;
      }
      if( !SCIPisZero(scip, residual) )
      {
         /* assign residual to first block */
         shift[nonviolatedblockslhs[0]] -= residual;
      }
   }

#ifdef SCIP_DEBUG
   SCIP_Real sum = 0.0;
   /* check sum of shift; must be zero */
   for( v = 0; v < linking->nblocks; v++ )
      sum += shift[v];
   assert(SCIPisZero(scip, sum));
#endif

   /* add shift to both sides */
   for( v = 0; v < linking->nblocks; v++)
   {
      if( linking->hasrhs )
         linking->currentrhs[v] += shift[v];

      if( linking->haslhs )
         linking->currentlhs[v] += shift[v];
   }

   SCIP_CALL( roundPartition(scip, linking, blockproblem, ((linking->hasrhs && (*nviolatedblocksrhs != 0)) || !linking->haslhs)) );

   /* set sides in blockproblems to new partition */
   for( v = 0; v < linking->nblocks; v++ )
   {
      SCIP* subscip;
      subscip = blockproblem[linking->blocknumbers[v]]->blockscip;

      if( linking->hasrhs )
      {
         SCIP_CALL( SCIPchgRhsLinear(subscip, linking->blockconss[v], linking->currentrhs[v]) );
      }
      if( linking->haslhs )
      {
         SCIP_CALL( SCIPchgLhsLinear(subscip, linking->blockconss[v], linking->currentlhs[v]) );
      }
   }

   /* free memory */
   if( linking->haslhs )
      SCIPfreeBufferArray(scip, &nonviolatedblockslhs);
   if( linking->hasrhs )
      SCIPfreeBufferArray(scip, &nonviolatedblocksrhs);
   SCIPfreeBufferArray(scip, &shift);

   return SCIP_OKAY;
}

/** update penalty parameters lambda
 *
 * if a linking constraint is violated two times in succession, the corresponding penalty parameter is increased in each block
 */
static
SCIP_RETCODE updateLambda(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   LINKING**             linkings,           /**< array of linking data structures */
   BLOCKPROBLEM**        blockproblem,       /**< array of blockproblem data structures */
   int*                  nviolatedblocksrhs, /**< number of blocks which violate rhs */
   int*                  nviolatedblockslhs, /**< number of blocks which violate lhs */
   int                   nlinking            /**< number of linking constraints */
   )
{
   SCIP_VAR* slackvar;
   SCIP_Real old_obj;
   SCIP_Real new_obj;
   int c;
   int b;

   assert(scip != NULL);
   assert(linkings != NULL);
   assert(blockproblem != NULL);

   for( c = 0; c < nlinking; c++ )
   {
      assert(nviolatedblocksrhs[c] >= 0);
      assert(nviolatedblockslhs[c] >= 0);

      /* skip constraint if it is not violated */
      if( nviolatedblocksrhs[c] + nviolatedblockslhs[c] == 0 )
      {
         linkings[c]->lastviolations = 0; /* reset flag */
         continue;
      }

      /* add number of violated blocks multiplied with parameter "penalty" to lambda (initial value is 1) */
      for( b = 0; b < linkings[c]->nblocks; b++ )
      {
         SCIP* subscip;
         subscip = blockproblem[linkings[c]->blocknumbers[b]]->blockscip;
         assert(subscip != NULL);

         if( linkings[c]->hasrhs && (nviolatedblocksrhs[c] >= 1) && (linkings[c]->lastviolations >= 1) )
         {
            slackvar = linkings[c]->slacks[b * linkings[c]->nslacksperblock];
            old_obj = SCIPvarGetObj(slackvar);
            new_obj = old_obj + heurdata->penalty * nviolatedblocksrhs[c];

            SCIP_CALL( SCIPchgVarObj(subscip, slackvar, new_obj) );
         }
         if( linkings[c]->haslhs && (nviolatedblockslhs[c] >= 1) && (linkings[c]->lastviolations >= 1) )
         {
            slackvar = linkings[c]->slacks[b * linkings[c]->nslacksperblock + linkings[c]->nslacksperblock - 1];
            old_obj = SCIPvarGetObj(slackvar);
            new_obj = old_obj + heurdata->penalty * nviolatedblockslhs[c];

            SCIP_CALL( SCIPchgVarObj(subscip, slackvar, new_obj) );
         }
      }

      /* update number of violations in the last iterations */
      linkings[c]->lastviolations += 1;
   }

   return SCIP_OKAY;
}

/** computes feasible solution from last stored solution for each block and adds it to the solution storage */
static
SCIP_RETCODE reuseSolution(
   LINKING**             linkings,           /**< array of linking data structures */
   BLOCKPROBLEM**        blockproblem,       /**< array of blockproblem data structures */
   int                   nblocks             /**< number of blocks */
   )
{
   SCIP_SOL** sols;
   SCIP_SOL* sol; /* solution of block that will be repaired */
   SCIP_SOL* newsol;
   SCIP_VAR** blockvars;
   SCIP_Real* blockvals;
   int nsols;
   int nvars;
   int b;
   int c;
   int i;
   SCIP_Bool success;

   assert(linkings != NULL);
   assert(blockproblem != NULL);

   for( b = 0; b < nblocks; b++ )
   {
      SCIP* subscip;

      subscip = blockproblem[b]->blockscip;
      nsols = SCIPgetNSols(subscip);

      /* no solution in solution candidate storage found */
      if( nsols == 0 )
         return SCIP_OKAY;

      /* take last solution */
      sols = SCIPgetSols(subscip);
      sol = sols[nsols - 1];

      /* copy the solution */
      nvars = SCIPgetNVars(subscip);
      blockvars = SCIPgetVars(subscip);
      SCIP_CALL( SCIPallocBufferArray(subscip, &blockvals, nvars) );
      SCIP_CALL( SCIPgetSolVals(subscip, sol, nvars, blockvars, blockvals) );
      SCIP_CALL( SCIPcreateOrigSol(subscip, &newsol, NULL) );
      SCIP_CALL( SCIPsetSolVals(subscip, newsol, nvars, blockvars, blockvals) );

      /* correct each coupling constraint:
      * lhs <= orig_var - z_r + z_l <= rhs
      * adapt slack variables so that constraint is feasible
      */
      for( c = 0; c < blockproblem[b]->nlinking; c++ )
      {
         LINKING* linking; /* linking data structure of this constraint */
         SCIP_VAR* rvar; /* slack variable z_r */
         SCIP_VAR* lvar; /* slack variable z_l */
         SCIP_Real rval; /* value of slack variable z_r */
         SCIP_Real lval; /* value of slack variable z_l */
         SCIP_Real activitycons; /* activity of constraint*/
         SCIP_Real activity; /* activity of constraint without slack variables */
         SCIP_Real rhs; /* current right hand side */
         SCIP_Real lhs; /* current left hand side */

         linking = linkings[blockproblem[b]->linkingindices[c]];
         rhs = SCIPgetRhsLinear(subscip, blockproblem[b]->linkingconss[c]);
         lhs = SCIPgetLhsLinear(subscip, blockproblem[b]->linkingconss[c]);
         assert(SCIPisGE(subscip, rhs, lhs));

         activitycons = SCIPgetActivityLinear(subscip, blockproblem[b]->linkingconss[c], sol);

         /* get slack variables and subtract their value from the activity;
          * calculate and set values of slack variables
          */
         for( i = 0; i < linking->nblocks; i++ )
         {
            if( linking->blocknumbers[i] == b )
            {
               if( linking->hasrhs && linking->haslhs )
               {
                  rvar = linking->slacks[2 * i];
                  lvar = linking->slacks[2 * i + 1];
                  rval = SCIPgetSolVal(subscip, sol, rvar);
                  lval = SCIPgetSolVal(subscip, sol, lvar);
                  activity = activitycons + rval - lval;
                  SCIP_CALL( SCIPsetSolVal(subscip, newsol, rvar, MAX(0.0, activity - rhs)) );
                  SCIP_CALL( SCIPsetSolVal(subscip, newsol, lvar, MAX(0.0, lhs - activity)) );
               }
               else if( linking->hasrhs )
               {
                  rvar = linking->slacks[i];
                  rval = SCIPgetSolVal(subscip, sol, rvar);
                  activity = activitycons + rval;
                  SCIP_CALL( SCIPsetSolVal(subscip, newsol, rvar, MAX(0.0, activity - rhs)) );
               }
               else /* linking->haslhs */
               {
                  assert(linking->haslhs);
                  lvar = linking->slacks[i];
                  lval = SCIPgetSolVal(subscip, sol, lvar);
                  activity = activitycons - lval;
                  SCIP_CALL( SCIPsetSolVal(subscip, newsol, lvar, MAX(0.0, lhs - activity)) );
               }
               break;
            }
         }
      }

      SCIPdebugMsg(subscip, "Try adding solution with objective value %.2f\n", SCIPgetSolOrigObj(subscip, newsol));
      SCIP_CALL( SCIPaddSolFree(subscip, &newsol, &success) );

      if( !success )
         SCIPdebugMsg(subscip, "Correcting solution failed\n"); /* maybe not better than old solutions */
      else
         SCIPdebugMsg(subscip, "Correcting solution successful\n");

      SCIPfreeBufferArray(subscip, &blockvals);
   }

   return SCIP_OKAY;
}

/** reoptimizes the heuristic solution with original objective function */
static
SCIP_RETCODE reoptimize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< pointer to heuristic */
   SCIP_SOL*             sol,                /**< heuristic solution */
   BLOCKPROBLEM**        blockproblem,       /**< array of blockproblem data structures */
   int                   nblocks,            /**< number of blockproblems */
   SCIP_SOL**            newsol,             /**< pointer to store improved solution */
   SCIP_Bool*            success             /**< pointer to store whether reoptimization was successful */
   )
{
   SCIP_Real time;
   SCIP_Real timesubscip;
   SCIP_Bool check;
   int b;
   int v;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(sol != NULL);
   assert(blockproblem != NULL);

   *success = FALSE;
   check = FALSE;

   /* for each blockproblem:
    * - change back to original objective function
    * - fix slack variables to zero
    * - set limits and solve problem
    */
   for( b = 0; b < nblocks; b++ )
   {
      SCIP* subscip;
      SCIP_VAR** blockvars;
      int nvars;

      subscip = blockproblem[b]->blockscip;
      timesubscip = SCIPgetTotalTime(subscip);
      blockvars = SCIPgetOrigVars(subscip);
      nvars = SCIPgetNOrigVars(subscip);

      /* in order to change objective function */
      SCIP_CALL( SCIPfreeTransform(subscip) );

      /* change back to original objective function */
      for( v = 0; v < blockproblem[b]->nblockvars; v++ )
      {
         SCIP_CALL( SCIPchgVarObj(subscip, blockvars[v], blockproblem[b]->origobj[v]) );
      }

      /* fix slack variables to zero */
      for( v = blockproblem[b]->nblockvars; v < nvars; v++ )
      {
         SCIP_CALL( SCIPchgVarUb(subscip, blockvars[v], 0.0) );
         SCIP_CALL( SCIPchgVarLb(subscip, blockvars[v], 0.0) );
      }

      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

      /* forbid recursive call of heuristics and separators solving sub-SCIPs */
      SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

#ifdef SCIP_DEBUG
      /* for debugging, enable full output */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
      SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#else
      /* disable statistic timing inside sub SCIP and output to console */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
      SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

      /* disable cutting plane separation */
      SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

      /* disable expensive presolving */
      SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

      /* disable expensive techniques */
      SCIP_CALL( SCIPsetIntParam(subscip, "misc/usesymmetry", 0) );

      /* speed up sub-SCIP by not checking dual LP feasibility */
      SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

      /* set limits; do not use more time than the heuristic has already used for first solution */
      SCIP_CALL( SCIPcopyLimits(scip, subscip) );
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", 1LL) );
      SCIP_CALL( SCIPgetRealParam(subscip, "limits/time", &time) );
      if( timesubscip <  time - 1.0 )
      {
         SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timesubscip + 1.0) );
      }
      SCIP_CALL( SCIPtransformProb(subscip) );
      SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", SCIPgetNSols(subscip) + 1) );

      /* reoptimize problem */
      SCIP_CALL_ABORT( SCIPsolve(subscip) );

      if( SCIPgetNSols(subscip) == 0 )
      {
         /* we found no solution */
         return SCIP_OKAY;
      }
      else if( SCIPgetStatus(subscip) == SCIP_STATUS_BESTSOLLIMIT || SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
      {
         check = TRUE;
      }
   }

   if( !check )
   {
      /* we have no better solution */
      return SCIP_OKAY;
   }

   /* create sol of main scip */
   SCIP_CALL( SCIPcreateSol(scip, newsol, heur) );

   /* copy solution to main scip */
   for( b = 0; b < nblocks; b++ )
   {
      SCIP_SOL* blocksol;
      SCIP_VAR** blockvars;
      SCIP_Real* blocksolvals;
      int nblockvars;

      /* get solution of block variables (without slack variables) */
      blocksol = SCIPgetBestSol(blockproblem[b]->blockscip);
      blockvars = SCIPgetOrigVars(blockproblem[b]->blockscip);
      nblockvars = blockproblem[b]->nblockvars;

      SCIP_CALL( SCIPallocBufferArray(scip, &blocksolvals, nblockvars) );
      SCIP_CALL( SCIPgetSolVals(blockproblem[b]->blockscip, blocksol, nblockvars, blockvars, blocksolvals) );

      for( v = 0; v < nblockvars; v++ )
      {
         SCIP_VAR* origvar;

         origvar = SCIPfindVar(scip, SCIPvarGetName(blockvars[v]));
         SCIP_CALL( SCIPsetSolVal(scip, *newsol, origvar, blocksolvals[v]) );
      }

      SCIPfreeBufferArray(scip, &blocksolvals);
   }

   *success = TRUE;

   return SCIP_OKAY;
}


/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * interrupt solution process of sub-SCIP if dual bound is greater than zero and a solution is available
 */
static
SCIP_DECL_EVENTEXEC(eventExecDps)
{
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_LPSOLVED);

   SCIPdebugMsg(scip, "dual bound: %.2f\n", SCIPgetDualbound(scip));

   if( SCIPisFeasGT(scip, SCIPgetDualbound(scip), 0.0) && SCIPgetNSols(scip) >= 1 )
   {
      SCIPdebugMsg(scip, "DPS: interrupt subscip\n");
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyDps)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurDps(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeDps)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecDps)
{  /*lint --e{715}*/
   SCIP_DECOMP** alldecomps;
   SCIP_DECOMP* decomp;
   SCIP_DECOMP* assigneddecomp;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_VAR** sortedvars;
   SCIP_CONS** sortedconss;
   SCIP_HEURDATA* heurdata;
   BLOCKPROBLEM** blockproblem;
   LINKING** linkings;
   int* sortedvarlabels;
   int* sortedconslabels;
   SCIP_EVENTHDLR* eventhdlr; /* event handler */
   SCIP_Real memory; /* in MB */
   SCIP_Real timelimit;
   SCIP_Real allslacksval;
   SCIP_Real blocksolval;
   SCIP_STATUS status;
   SCIP_Bool avoidmemout;
   SCIP_Bool disablemeasures;
   int maxgraphedge;
   int ndecomps;
   int nvars;
   int nconss;
   int nblocks;
   SCIP_Bool success;
   int b;
   int c;
   int k;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   assigneddecomp = NULL;
   blockproblem = NULL;
   linkings = NULL;
   eventhdlr = NULL;

   *result = SCIP_DIDNOTRUN;

   /* -------------------------------------------------------------------- */
   SCIPdebugMsg(scip, "initialize dps heuristic\n");

   /* take the first transformed decomposition */
   SCIPgetDecomps(scip, &alldecomps, &ndecomps, FALSE);
   if( ndecomps == 0)
      return SCIP_OKAY;

   decomp = alldecomps[0];
   assert(decomp != NULL);
   SCIPdebugMsg(scip, "First transformed decomposition is selected\n");

   nblocks = SCIPdecompGetNBlocks(decomp);
   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);

   /* if problem has no constraints, no variables or less than two blocks, return */
   if( nconss == 0 || nvars == 0 || nblocks <= 1 )
   {
      SCIPdebugMsg(scip, "problem has no constraints, no variables or less than two blocks\n");
      return SCIP_OKAY;
   }

   /* estimate required memory for all blocks and terminate if not enough memory is available */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memory) );
   SCIP_CALL( SCIPgetBoolParam(scip, "misc/avoidmemout", &avoidmemout) );
   if( avoidmemout && (((SCIPgetMemUsed(scip) + SCIPgetMemExternEstim(scip))/1048576.0) * (nblocks/4.0 + 2) >= memory) )
   {
      SCIPdebugMsg(scip, "The estimated memory usage for %d blocks is too large.\n", nblocks);
      return SCIP_OKAY;
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

   vars = SCIPgetVars(scip);
   conss = SCIPgetConss(scip);
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sortedvars, vars, nvars) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sortedconss, conss, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedvarlabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedconslabels, nconss) );

   /* get labels and sort in increasing order */
   SCIPdecompGetVarsLabels(decomp, sortedvars, sortedvarlabels, nvars);
   SCIPdecompGetConsLabels(decomp, sortedconss, sortedconslabels, nconss);
   SCIPsortIntPtr(sortedvarlabels, (void**)sortedvars, nvars);
   SCIPsortIntPtr(sortedconslabels, (void**)sortedconss, nconss);

   if( sortedvarlabels[0] == SCIP_DECOMP_LINKVAR )
   {
      /* create new decomposition; don't change the decompositions in the decompstore */
      SCIP_CALL( SCIPdecompCreate(&assigneddecomp, SCIPblkmem(scip), nblocks, SCIPdecompIsOriginal(decomp), SCIPdecompUseBendersLabels(decomp)) );

      SCIP_CALL( assignLinking(scip, assigneddecomp, sortedvars, sortedconss, sortedvarlabels, sortedconslabels, nvars, nconss, SCIPdecompGetNBorderVars(decomp)) );
      assert(SCIPdecompGetNBlocks(decomp) >= SCIPdecompGetNBlocks(assigneddecomp));
      decomp = assigneddecomp;

      /* number of blocks can get smaller */
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

#ifdef SCIP_DEBUG
      char buffer[SCIP_MAXSTRLEN];
      SCIPdebugMsg(scip, "DPS used decomposition:\n%s\n", SCIPdecompPrintStats(decomp, buffer));
#endif

   /* check properties of used decomposition */
   if( heurdata->maxlinkscore != 1.0 )
   {
      SCIP_Real linkscore;
      int nlinkconss;

      nlinkconss = SCIPdecompGetNBorderConss(decomp);

      linkscore = (SCIP_Real)nlinkconss / (SCIP_Real)nconss;

      if( linkscore > heurdata->maxlinkscore )
      {
         SCIPdebugMsg(scip, "decomposition has not the required properties\n");
         goto TERMINATE;
      }
   }

   if( sortedvarlabels[0] == SCIP_DECOMP_LINKVAR ||
         sortedconslabels[0] != SCIP_DECOMP_LINKCONS ||
         nblocks <= 1 )
   {
      SCIPdebugMsg(scip, "Problem has linking variables or no linking constraints or less than two blocks\n");
      goto TERMINATE;
   }

   /* initialize heurdata */
   heurdata->linkingconss = sortedconss;
   heurdata->nlinking = SCIPdecompGetNBorderConss(decomp);
   heurdata->nblocks = nblocks;

   /* allocate memory for blockproblems and initialize partially */
   SCIP_CALL( SCIPallocBufferArray(scip, &blockproblem, nblocks) );
   for( b = 0; b < nblocks; b++ )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &(blockproblem[b])) ); /*lint !e866*/
      SCIP_CALL( createSubscip(scip, &blockproblem[b]->blockscip) );

      SCIP_CALL( SCIPallocBufferArray(scip, &blockproblem[b]->linkingconss, heurdata->nlinking) );
      SCIP_CALL( SCIPallocBufferArray(scip, &blockproblem[b]->linkingindices, heurdata->nlinking) );
      SCIP_CALL( SCIPallocBufferArray(scip, &blockproblem[b]->slackvars, heurdata->nlinking * 2) ); /* maximum two slacks per linking constraint */
      SCIP_CALL( SCIPallocBufferArray(scip, &blockproblem[b]->origobj, nvars) );
      blockproblem[b]->nblockvars = 0;
      blockproblem[b]->nlinking = 0;
      blockproblem[b]->nslackvars = 0;
   }

   /* allocate memory for linkings and initialize partially */
   SCIP_CALL( SCIPallocBufferArray(scip, &linkings, heurdata->nlinking) );
   for( c = 0; c < heurdata->nlinking; c++ )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &(linkings[c])) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &(linkings[c])->blockconss, heurdata->nblocks) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(linkings[c])->slacks, heurdata->nblocks*2) ); /* maximum two slacks per block */
      SCIP_CALL( SCIPallocBufferArray(scip, &(linkings[c])->blocknumbers, heurdata->nblocks) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(linkings[c])->minactivity, heurdata->nblocks) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(linkings[c])->maxactivity, heurdata->nblocks) );

      linkings[c]->linkingcons = heurdata->linkingconss[c];
      linkings[c]->currentrhs = NULL;
      linkings[c]->currentlhs = NULL;
      linkings[c]->nblocks = 0;
      linkings[c]->nslacks = 0;
      linkings[c]->nslacksperblock = 0;
      linkings[c]->lastviolations = 0;
      linkings[c]->hasrhs = FALSE;
      linkings[c]->haslhs = FALSE;
   }

   SCIP_CALL( createAndSplitProblem(scip, heurdata, decomp, blockproblem, linkings, sortedvars, sortedconss, &success) );
   if( !success )
   {
      SCIPdebugMsg(scip, "Create and split problem failed\n");
      goto TERMINATE;
   }

   /* allocate memory for current partition*/
   for( c = 0; c < heurdata->nlinking; c++ )
   {
      if( linkings[c]->hasrhs )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &(linkings[c])->currentrhs, linkings[c]->nblocks ) );
      }

      if( linkings[c]->haslhs )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &(linkings[c])->currentlhs, linkings[c]->nblocks ) );
      }
   }

   /* initialize partition */
   SCIP_CALL( initCurrent(scip, linkings, blockproblem, heurdata->nlinking, &success) );

   /** ------------------------------------------------------------------------ */
   SCIPdebugMsg(scip, "Start heuristik DPS\n");
   *result = SCIP_DIDNOTFIND;

   for( k = 0; k < heurdata->maxit; k++ )
   {
      /* do not exceed the timelimit */
      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      if( (timelimit - SCIPgetSolvingTime(scip)) <= 0 )
         goto TERMINATE;

      /* solve the subproblems */
      allslacksval = 0.0;
      for( b = 0; b < heurdata->nblocks; b++ )
      {
         SCIP* subscip;
         subscip = blockproblem[b]->blockscip;

         /* update time and memory limit of subproblem */
         SCIP_CALL( SCIPcopyLimits(scip, subscip) );

         /* create event handler for LP events */
         if( k == 0 )
         {
            SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecDps, NULL) );
            if( eventhdlr == NULL )
            {
               SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
               return SCIP_PLUGINNOTFOUND;
            }
         }

         /* catch LP events of sub-SCIP */
         SCIP_CALL( SCIPtransformProb(subscip) );
         SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

         SCIPdebugMsg(scip, "Solve blockproblem %d\n", b);
         SCIP_CALL_ABORT( SCIPsolve(subscip) );

         /* drop LP events of sub-SCIP */
         SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

         /* get status and objective value if available */
         status = SCIPgetStatus(subscip);
         if( status == SCIP_STATUS_INFEASIBLE )
         {
            SCIPdebugMsg(scip, "Subproblem is infeasible\n");
            goto TERMINATE;
         }
         else if( status == SCIP_STATUS_UNBOUNDED )
         {
            SCIPdebugMsg(scip, "Subproblem is unbounded\n");
            goto TERMINATE;
         }
         else if( SCIPgetNSols(subscip) >= 1 )
         {
            blocksolval = SCIPgetPrimalbound(subscip);

            if( status == SCIP_STATUS_TIMELIMIT && !SCIPisZero(scip, blocksolval) )
            {
               SCIPdebugMsg(scip, "Subproblem reached timelimit without optimal solution\n");
               goto TERMINATE;
            }
            SCIPdebugMsg(scip, "Solution value: %f\n", blocksolval);
            allslacksval += blocksolval;
         }
         else
         {
            SCIPdebugMsg(scip, "No subproblem solution available\n");
            goto TERMINATE;
         }
      }

      /* all slack variables are zero -> we found a feasible solution */
      if( SCIPisZero(scip, allslacksval) )
      {
         SCIP_SOL* newsol;

         SCIPdebugMsg(scip, "Feasible solution found after %i iterations\n", k);

         /* create new solution */
         SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
         for( b = 0; b < heurdata->nblocks; b++ )
         {
            SCIP_SOL* blocksol;
            SCIP_VAR** blockvars;
            SCIP_Real* blocksolvals;
            int nblockvars;

            /* get solution of block variables (without slack variables) */
            blocksol = SCIPgetBestSol(blockproblem[b]->blockscip);
            blockvars = SCIPgetOrigVars(blockproblem[b]->blockscip);
            nblockvars = blockproblem[b]->nblockvars;

            SCIP_CALL( SCIPallocBufferArray(scip, &blocksolvals, nblockvars) );
            SCIP_CALL( SCIPgetSolVals(blockproblem[b]->blockscip, blocksol, nblockvars, blockvars, blocksolvals) );

            for( c = 0; c < nblockvars; c++ )
            {
               SCIP_VAR* origvar;

               origvar = SCIPfindVar(scip, SCIPvarGetName(blockvars[c]));
               SCIP_CALL( SCIPsetSolVal(scip, newsol, origvar, blocksolvals[c]) );
            }

            SCIPfreeBufferArray(scip, &blocksolvals);
         }

         /* if reoptimization is activated, fix partition and reoptimize with original objective function */
         if( heurdata->reoptimize )
         {
            SCIP_SOL* improvedsol = NULL;
            SCIP_CALL( reoptimize(scip, heur, newsol, blockproblem, heurdata->nblocks, &improvedsol, &success) );
            assert(improvedsol != NULL || success == FALSE);

            if( success )
            {
               SCIP_CALL( SCIPtrySolFree(scip, &improvedsol, TRUE, FALSE, TRUE, TRUE, TRUE, &success) );
               if( success )
               {
                  SCIPdebugMsg(scip, "Reoptimizing solution successful\n");
                  *result = SCIP_FOUNDSOL;
               }
            }
         }

         /* if reoptimization is turned off or reoptimization found no solution, try initial solution */
         if( *result != SCIP_FOUNDSOL )
         {
            SCIPdebugMsg(scip, "Solution has value: %.2f\n", SCIPgetSolOrigObj(scip, newsol));
            SCIP_CALL( SCIPtrySolFree(scip, &newsol, TRUE, FALSE, TRUE, TRUE, TRUE, &success) );
            if( success )
            {
               SCIPdebugMsg(scip, "Solution copy successful\n");
               *result = SCIP_FOUNDSOL;
            }
         }
         else
         {
            SCIP_CALL( SCIPfreeSol(scip, &newsol) );
         }

         goto TERMINATE;
      }
      /* current partition is not feasible:
       * - update partition
       * - update penalty parameters lambda
       * - reuse last solution
       */
      else
      {
         int* nviolatedblocksrhs; /* number of blocks which violate rhs for each linking constraint */
         int* nviolatedblockslhs; /* number of blocks which violate lhs for each linking constraint */

         SCIP_CALL( SCIPallocBufferArray(scip, &nviolatedblocksrhs, heurdata->nlinking) );
         SCIP_CALL( SCIPallocBufferArray(scip, &nviolatedblockslhs, heurdata->nlinking) );
         BMSclearMemoryArray(nviolatedblocksrhs, heurdata->nlinking);
         BMSclearMemoryArray(nviolatedblockslhs, heurdata->nlinking);

         SCIPdebugMsg(scip, "update partition\n");
         for( c = 0; c < heurdata->nlinking; c++ )
         {
            SCIP_CALL( updatePartition(scip, linkings[c], blockproblem, &nviolatedblocksrhs[c], &nviolatedblockslhs[c], k) );
         }

         /* in order to change objective function */
         for( b = 0; b < heurdata->nblocks; b++ )
         {
            SCIP_CALL( SCIPfreeTransform(blockproblem[b]->blockscip) );
         }

         SCIPdebugMsg(scip, "update lambda\n");
         SCIP_CALL( updateLambda(scip, heurdata, linkings, blockproblem, nviolatedblocksrhs, nviolatedblockslhs, heurdata->nlinking) );

         if( heurdata->reuse )
         {
            /* reuse old solution in each block if available */
            SCIP_CALL( reuseSolution(linkings, blockproblem, nblocks) );
         }

         SCIPfreeBufferArray(scip, &nviolatedblockslhs);
         SCIPfreeBufferArray(scip, &nviolatedblocksrhs);
      }
   }
   SCIPdebugMsg(scip, "maximum number of iterations reached\n");

   /** ------------------------------------------------------------------------ */
   /** free memory */
TERMINATE:
   if( linkings != NULL )
   {
      for( c = heurdata->nlinking - 1; c >= 0; c-- )
      {
         if( linkings[c]->currentlhs != NULL )
            SCIPfreeBufferArray(scip, &(linkings[c])->currentlhs);

         if( linkings[c]->currentrhs != NULL )
            SCIPfreeBufferArray(scip, &(linkings[c])->currentrhs);
      }

      for( c = heurdata->nlinking - 1; c >= 0; c-- )
      {
         linkings[c]->linkingcons = NULL;
         SCIPfreeBufferArray(scip, &(linkings[c])->maxactivity);
         SCIPfreeBufferArray(scip, &(linkings[c])->minactivity);
         SCIPfreeBufferArray(scip, &(linkings[c])->blocknumbers);
         SCIPfreeBufferArray(scip, &(linkings[c])->slacks);
         SCIPfreeBufferArray(scip, &(linkings[c])->blockconss);
         SCIPfreeBlockMemory(scip, &(linkings[c])); /*lint !e866*/
      }
      SCIPfreeBufferArray(scip, &linkings);
   }

   if( blockproblem != NULL )
   {
      for( b = nblocks - 1; b >= 0; b-- )
      {
         SCIPfreeBufferArray(scip, &(blockproblem[b])->origobj);
         SCIPfreeBufferArray(scip, &(blockproblem[b])->slackvars);
         SCIPfreeBufferArray(scip, &(blockproblem[b])->linkingindices);
         SCIPfreeBufferArray(scip, &(blockproblem[b])->linkingconss);
         SCIP_CALL( SCIPfree(&blockproblem[b]->blockscip) );
         SCIPfreeBlockMemory(scip, &(blockproblem[b])); /*lint !e866*/
      }
      SCIPfreeBufferArray(scip, &blockproblem);
   }

   if( assigneddecomp != NULL )
      SCIPdecompFree(&assigneddecomp, SCIPblkmem(scip));

   if( sortedconslabels != NULL )
      SCIPfreeBufferArray(scip, &sortedconslabels);

   if( sortedvarlabels != NULL )
      SCIPfreeBufferArray(scip, &sortedvarlabels);

   if( sortedconss != NULL )
      SCIPfreeBufferArray(scip, &sortedconss);

   if( sortedvars != NULL )
      SCIPfreeBufferArray(scip, &sortedvars);

   SCIPdebugMsg(scip, "Leave DPS heuristic\n");

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the dps primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurDps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create dps primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   heur = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecDps, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyDps) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeDps) );

   /* add dps primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxiterations",
   "maximal number of iterations", &heurdata->maxit, FALSE, DEFAULT_MAXIT, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxlinkscore",
   "maximal linking score of used decomposition (equivalent to percentage of linking constraints)",
   &heurdata->maxlinkscore, FALSE, 1.0, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/penalty",
   "multiplier for absolute increase of penalty parameters (0: no increase)",
   &heurdata->penalty, FALSE, DEFAULT_PENALTY, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/reoptimize",
   "should the problem get reoptimized with the original objective function?", &heurdata->reoptimize, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/reuse",
   "should solutions get reused in subproblems?", &heurdata->reuse, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
