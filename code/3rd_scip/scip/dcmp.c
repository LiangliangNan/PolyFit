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

/**@file   dcmp.c
 * @ingroup OTHER_CFILES
 * @brief  internal methods for decompositions and the decomposition store
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include "scip/dcmp.h"
#include "scip/mem.h"
#include "scip/pub_cons.h"
#include "scip/pub_dcmp.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_dcmp.h"
#include "scip/scip_mem.h"
#include "scip/scip_prob.h"
#include "scip/scip_var.h"
#include "scip/scip_general.h"
#include "scip/scip_var.h"
#include "scip/scip_datastructures.h"
#include "scip/scip_message.h"
#include "scip/struct_dcmp.h"
#include "scip/struct_scip.h"

/* create and free a decomposition */
#define INIT_MAP_SIZE 2000

/** creates a decomposition */
SCIP_RETCODE SCIPdecompCreate(
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   nblocks,            /**< the number of blocks (without the linking block) */
   SCIP_Bool             original,           /**< is this a decomposition in the original (TRUE) or transformed space? */
   SCIP_Bool             benderslabels       /**< should the variables be labeled for the application of Benders' decomposition */
   )
{
   int memsize;

   assert(decomp != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, decomp) );
   SCIP_CALL( SCIPhashmapCreate(&(*decomp)->var2block, blkmem, INIT_MAP_SIZE) );
   SCIP_CALL( SCIPhashmapCreate(&(*decomp)->cons2block, blkmem, INIT_MAP_SIZE) );

   /* we allocate one extra slot for the linking block */
   memsize = nblocks + 1;
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*decomp)->varssize, memsize) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*decomp)->consssize, memsize) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*decomp)->labels, memsize) );

   (*decomp)->memsize = memsize;
   (*decomp)->nblocks = nblocks;
   (*decomp)->modularity = -1.0;
   (*decomp)->idxsmallestblock = -1;
   (*decomp)->idxlargestblock = -1;
   (*decomp)->original = original;
   (*decomp)->benderslabels = benderslabels;
   (*decomp)->areascore = -1.0;
   (*decomp)->nedges = 0;
   (*decomp)->mindegree = 0;
   (*decomp)->maxdegree = 0;
   (*decomp)->ncomponents = 0;
   (*decomp)->narticulations = 0;
   (*decomp)->statscomplete = FALSE;

   return SCIP_OKAY;
}

/** free a decomposition */
void SCIPdecompFree(
   SCIP_DECOMP**         decomp,             /**< pointer to free the decomposition data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(decomp != NULL);
   assert(blkmem != NULL);

   if( *decomp == NULL )
      return;

   assert((*decomp)->var2block != NULL);
   SCIPhashmapFree(&(*decomp)->var2block);
   SCIPhashmapFree(&(*decomp)->cons2block);

   BMSfreeBlockMemoryArray(blkmem, &(*decomp)->varssize, (*decomp)->memsize);
   BMSfreeBlockMemoryArray(blkmem, &(*decomp)->consssize, (*decomp)->memsize);
   BMSfreeBlockMemoryArray(blkmem, &(*decomp)->labels, (*decomp)->memsize);

   BMSfreeBlockMemory(blkmem, decomp);
}

/* getter and setter for variable labels */

/** sets labels for an array of variables */
SCIP_RETCODE SCIPdecompSetVarsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int*                  labels,             /**< array of labels, one per variable */
   int                   nvars               /**< length of variables array */
   )
{
   int i;

   assert(decomp != NULL);
   assert(vars != NULL);
   assert(labels != NULL);

   /* store each label in hash map */
   for( i = 0; i < nvars; ++i )
   {
      assert(labels[i] == SCIP_DECOMP_LINKVAR || labels[i] >= 0);

      SCIP_CALL( SCIPhashmapSetImageInt(decomp->var2block, (void *)vars[i], labels[i]) );
   }

   return SCIP_OKAY;
}

/** queries labels for an array of variables */
void SCIPdecompGetVarsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int*                  labels,             /**< buffer to store labels, one per variable */
   int                   nvars               /**< length of variables array */
   )
{
   int i;

   assert(decomp != NULL);
   assert(vars != NULL);
   assert(labels != NULL);

   /* store variable labels in buffer array */
   for( i = 0; i < nvars; ++i )
   {
      if( SCIPhashmapExists(decomp->var2block, (void *)vars[i]) )
         labels[i] = SCIPhashmapGetImageInt(decomp->var2block, (void *)vars[i]);
      else
         labels[i] = SCIP_DECOMP_LINKVAR;
   }
}

/** sets labels for an array of constraints */
SCIP_RETCODE SCIPdecompSetConsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int*                  labels,             /**< array of labels, one per constraint */
   int                   nconss              /**< length of constraints array */
   )
{
   int i;

   assert(decomp != NULL);
   assert(conss != NULL);
   assert(labels != NULL);

   /* store each label in hash map */
   for( i = 0; i < nconss; ++i )
   {
      assert(labels[i] == SCIP_DECOMP_LINKCONS || labels[i] >= 0);

      SCIP_CALL( SCIPhashmapSetImageInt(decomp->cons2block, (void *)conss[i], labels[i]) );
   }

   return SCIP_OKAY;
}

/** queries labels for an array of constraints */
void SCIPdecompGetConsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int*                  labels,             /**< array of labels, one per constraint */
   int                   nconss              /**< length of constraints array */
   )
{
   int i;

   assert(decomp != NULL);
   assert(conss != NULL);
   assert(labels != NULL);

   /* store variable labels in buffer array */
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPhashmapExists(decomp->cons2block, (void *)conss[i]) )
      {
         labels[i] = SCIPhashmapGetImageInt(decomp->cons2block, (void *)conss[i]);
      }
      else
         labels[i] = SCIP_DECOMP_LINKCONS;
   }
}

/** clears the corresponding labeling (constraints, variables, or both) of this decomposition */
SCIP_RETCODE SCIPdecompClear(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Bool             clearvarlabels,     /**< should the variable labels be cleared? */
   SCIP_Bool             clearconslabels     /**< should the constraint labels be cleared? */
   )
{
   assert(decomp != NULL);

   if( clearvarlabels )
   {
      SCIP_CALL( SCIPhashmapRemoveAll(decomp->var2block) );
   }

   if( clearconslabels )
   {
      SCIP_CALL( SCIPhashmapRemoveAll(decomp->cons2block) );
   }

   return SCIP_OKAY;
}

/** returns TRUE if decomposition is in the original space */
SCIP_Bool SCIPdecompIsOriginal(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->original;
}

/** sets the parameter that indicates whether the variables must be labeled for the application of Benders'
 * decomposition
 */
void SCIPdecompSetUseBendersLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Bool             benderslabels       /**< whether Benders' variable labels should be used */
   )
{
   assert(decomp != NULL);

   decomp->benderslabels = benderslabels;
}

/** returns TRUE if the variables must be labeled for the application of Benders' decomposition */
SCIP_Bool SCIPdecompUseBendersLabels(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->benderslabels;
}

/** gets number of blocks of this decomposition */
int SCIPdecompGetNBlocks(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->nblocks;
}

/** gets area score of this decomposition */
SCIP_Real SCIPdecompGetAreaScore(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->areascore;
}

/** gets modularity of this decomposition */
SCIP_Real SCIPdecompGetModularity(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->modularity;
}

/** gets variable size for each block, sorted by increasing block label
 *
 * To get all variable sizes, set nlabels to SCIPdecompGetNBlocks() + 1.
 * The first entry corresponds to the number of border variables.
 *
 * @note Ensure that SCIPcomputeDecompStats() has been called before.
 *       If the decomposition was read from a file, this was done automatically.
 */
SCIP_RETCODE SCIPdecompGetVarsSize(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   int*                  varssize,           /**< array to store variable sizes of blocks*/
   int                   nlabels             /**< length of variable sizes array */
   )
{
   int i;
   int nsizes;

   assert(decomp != NULL);
   assert(decomp->labels[0] == SCIP_DECOMP_LINKVAR);
   assert(varssize != NULL);
   assert(nlabels >= 0);

   nsizes = MIN(nlabels, decomp->nblocks + 1);

   /* store variable sizes */
   for( i = 0; i < nsizes; ++i )
   {
      varssize[i] = decomp->varssize[i];
   }

   return SCIP_OKAY;
}

/** gets constraint size for each block, sorted by increasing block label
 *
 * To get all constraint sizes, set nlabels to SCIPdecompGetNBlocks() + 1.
 * The first entry corresponds to the number of border constraints.
 *
 * @note Ensure that SCIPcomputeDecompStats() has been called before.
 *       If the decomposition was read from a file, this was done automatically.
 */
SCIP_RETCODE SCIPdecompGetConssSize(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   int*                  consssize,          /**< array to store constraint sizes of blocks*/
   int                   nlabels             /**< length of constraint sizes array */
   )
{
   int i;
   int nsizes;

   assert(decomp != NULL);
   assert(decomp->labels[0] == SCIP_DECOMP_LINKVAR);
   assert(consssize != NULL);
   assert(nlabels >= 0);

   nsizes = MIN(nlabels, decomp->nblocks + 1);

   /* store constraint sizes */
   for( i = 0; i < nsizes; ++i )
   {
      consssize[i] = decomp->consssize[i];
   }

   return SCIP_OKAY;
}

/** gets number of border variables of this decomposition
 *
 * @note Ensure that SCIPcomputeDecompStats() has been called before.
 *       If the decomposition was read from a file, this was done automatically.
 */
int SCIPdecompGetNBorderVars(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);
   assert(decomp->labels[0] == SCIP_DECOMP_LINKVAR);

   return decomp->varssize[0];
}

/** gets number of border constraints of this decomposition
 *
 * @note Ensure that SCIPcomputeDecompStats() has been called before.
 *       If the decomposition was read from a file, this was done automatically.
 */
int SCIPdecompGetNBorderConss(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);
   assert(decomp->labels[0] == SCIP_DECOMP_LINKVAR);

   return decomp->consssize[0];
}

/** gets number of edges in the block-decomposition graph of this decomposition */
int SCIPdecompGetNBlockGraphEdges(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->nedges;
}

/** gets number of connected components in the block-decomposition graph of this decomposition */
int SCIPdecompGetNBlockGraphComponents(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->ncomponents;
}

/** gets number of articulation points in the block-decomposition graph of this decomposition */
int SCIPdecompGetNBlockGraphArticulations(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->narticulations;
}

/** gets the maximum degree of the block-decomposition graph of this decomposition */
int SCIPdecompGetBlockGraphMaxDegree(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->maxdegree;
}

/** gets the minimum degree of the block-decomposition graph of this decomposition */
int SCIPdecompGetBlockGraphMinDegree(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->mindegree;
}

/** prints decomposition statistics into string buffer */
char* SCIPdecompPrintStats(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   char*                 strbuf              /**< string buffer storage */
   )
{
   char* ptr;

   assert(decomp != NULL);
   assert(strbuf != NULL);

   ptr = strbuf;

   ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN,
            "Decomposition with %d blocks.\n",
            decomp->nblocks);
   ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN,
            "Largest block: Block %d with %d constraints and %d variables\n",
            decomp->nblocks == 0 ? -1 : decomp->labels[decomp->idxlargestblock],
            decomp->nblocks == 0 ? 0 : decomp->consssize[decomp->idxlargestblock],
            decomp->nblocks == 0 ? 0 : decomp->varssize[decomp->idxlargestblock]);
   ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN,
            "Smallest block: Block %d with %d constraints and %d variables\n",
            decomp->nblocks == 0 ? -1 : decomp->labels[decomp->idxsmallestblock],
            decomp->nblocks == 0 ? 0 : decomp->consssize[decomp->idxsmallestblock],
            decomp->nblocks == 0 ? 0 : decomp->varssize[decomp->idxsmallestblock]);
   ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN,
            "Border has %d constraints and %d variables\n",
            decomp->labels[0] == SCIP_DECOMP_LINKVAR ? decomp->consssize[0] : 0,
            decomp->labels[0] == SCIP_DECOMP_LINKVAR ? decomp->varssize[0] : 0
            );

   ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN,
            "Modularity: %.3f, Area Score: %.3f\n",
            decomp->modularity, decomp->areascore);

   (void) SCIPsnprintf(ptr, SCIP_MAXSTRLEN,
      "Constraint Block Graph: %d edges, %d articulation points, %d connected components, %d min., %d max. degree%s\n",
      decomp->nedges, decomp->narticulations, decomp->ncomponents, decomp->mindegree, decomp->maxdegree,
      decomp->statscomplete ? "" :
               "(approximately: graph construction hit size limit)");

   return strbuf;
}

/** creates a decomposition storage */
SCIP_RETCODE SCIPdecompstoreCreate(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to store decomposition storage */
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   int                   nslots              /**< maximum number of decomposition slots in storage */
   )
{
   assert(decompstore != NULL);
   assert(blkmem != NULL);
   assert(nslots > 0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, decompstore) );

   (*decompstore)->ndecomps = 0;
   (*decompstore)->norigdecomps = 0;
   (*decompstore)->decompssize = nslots;

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*decompstore)->decomps, nslots) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*decompstore)->origdecomps, nslots) );

   return SCIP_OKAY;
}

/** frees array of decompositions */
static
void freeDecompositions(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_DECOMP**         decomps,            /**< decomposition array */
   int*                  ndecomps            /**< pointer for initial number of decompositions, will be set to 0 */
   )
{
   int d;

   assert(decomps != NULL);
   assert(ndecomps != NULL);

   /* delete all remaining decompositions from this store */
   for( d = 0; d < *ndecomps; ++d )
      SCIPdecompFree(&decomps[d], blkmem);

   *ndecomps = 0;
}

/** frees all decompositions in transformed space */
void SCIPexitSolveDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DECOMPSTORE* decompstore = scip->decompstore;

   assert(decompstore != NULL);

   freeDecompositions(SCIPblkmem(scip), decompstore->decomps, &decompstore->ndecomps);
}

/** frees a decomposition storage */
void SCIPdecompstoreFree(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to store decomposition storage */
   BMS_BLKMEM*           blkmem              /**< block memory data structure */
   )
{
   assert(decompstore != NULL);

   if( *decompstore == NULL )
      return;

   freeDecompositions(blkmem, (*decompstore)->origdecomps, &(*decompstore)->norigdecomps);
   freeDecompositions(blkmem, (*decompstore)->decomps, &(*decompstore)->ndecomps);

   BMSfreeBlockMemoryArray(blkmem, &(*decompstore)->decomps, (*decompstore)->decompssize);
   BMSfreeBlockMemoryArray(blkmem, &(*decompstore)->origdecomps, (*decompstore)->decompssize);

   BMSfreeBlockMemory(blkmem, decompstore);
}

/** adds decomposition to storage */
SCIP_RETCODE SCIPdecompstoreAdd(
   SCIP_DECOMPSTORE*     decompstore,        /**< decomposition storage */
   SCIP_DECOMP*          decomp              /**< decomposition to add */
   )
{
   SCIP_DECOMP** decomps;
   int* ndecompsptr;

   assert(decompstore != NULL);
   assert(decomp != NULL);

   /* distinguish between storage for original or transformed decompositions */
   if( SCIPdecompIsOriginal(decomp) )
   {
      decomps = decompstore->origdecomps;
      ndecompsptr = &decompstore->norigdecomps;
   }
   else
   {
      decomps = decompstore->decomps;
      ndecompsptr = &decompstore->ndecomps;
   }

   /* ensure that storage capacity is not exceeded */
   if( *ndecompsptr == decompstore->decompssize )
   {
      SCIPerrorMessage("Error: Decomposition storage size exceeded, maximum is %d decompositions\n", decompstore->decompssize);
      return SCIP_ERROR;
   }

   decomps[(*ndecompsptr)++] = decomp;

   return SCIP_OKAY;
}

/** gets decompositions from storage */
SCIP_DECOMP** SCIPdecompstoreGetDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   )
{
   assert(decompstore != NULL);

   return decompstore->decomps;
}

/** gets number of decompositions in storage */
int SCIPdecompstoreGetNDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   )
{
   assert(decompstore != NULL);
   return decompstore->ndecomps;
}

/** gets decompositions from storage */
SCIP_DECOMP** SCIPdecompstoreGetOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   )
{
   assert(decompstore != NULL);

   return decompstore->origdecomps;
}

/** gets number of decompositions in storage */
int SCIPdecompstoreGetNOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   )
{
   assert(decompstore != NULL);
   return decompstore->norigdecomps;
}

/** transforms all available original decompositions into transformed space */
SCIP_RETCODE SCIPtransformDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   int d;
   int v;
   SCIP_DECOMPSTORE* decompstore;
   SCIP_VAR** vars;
   SCIP_VAR** origvars;
   SCIP_VAR** varssorted;
   SCIP_CONS** conss;
   int nconss;
   int nvars;
   int nvarsoriginal;
   int nvarsintroduced;
   int* varslabels;
   SCIP_Bool original = FALSE;

   assert(scip != NULL);
   assert(scip->decompstore != NULL);

   decompstore = scip->decompstore;
   assert(decompstore->ndecomps == 0);

   assert(SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVED);

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &varssorted, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &origvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, nvars) );

   /* determine if variable has an original counterpart or not, and put it into varssorted array at the front or back */
   nvarsoriginal = nvarsintroduced = 0;
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Real scalar;
      SCIP_Real constant;
      SCIP_VAR* origvar;

      origvar = vars[v];
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

      /* the variable has no original counterpart and is therefore put at the end of the array */
      if( origvar == NULL )
      {
         varssorted[nvars - 1 - nvarsintroduced] = vars[v];
         ++nvarsintroduced;
      }
      else
      {
         varssorted[nvarsoriginal] = vars[v];
         origvars[nvarsoriginal] = origvar;
         ++nvarsoriginal;
      }

      assert(nvarsoriginal + nvarsintroduced <= nvars);
   }

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   /* loop over available, original decompositions, transform and add them to the storage */
   for( d = 0; d < decompstore->norigdecomps; ++d )
   {
      SCIP_DECOMP* origdecomp = decompstore->origdecomps[d];
      SCIP_DECOMP* decomp;
      char strbuf[SCIP_MAXSTRLEN];

      /* 1. query the decomposition labels of the original variables and set them for the transformed variables
       * that have original counterparts
       */
      SCIP_CALL( SCIPcreateDecomp(scip, &decomp, SCIPdecompGetNBlocks(origdecomp), original, SCIPdecompUseBendersLabels(origdecomp)) );

      SCIPdecompGetVarsLabels(origdecomp, origvars, varslabels, nvarsoriginal);

      SCIP_CALL( SCIPdecompSetVarsLabels(decomp, varssorted, varslabels, nvarsoriginal) );

      /* 2. compute the constraint labels based on the preliminary variable labels */
      SCIP_CALL( SCIPcomputeDecompConsLabels(scip, decomp, conss, nconss) );

      /* 3. remove the variable labels now that we have constraint labels */
      SCIP_CALL( SCIPdecompClear(decomp, TRUE, FALSE) );

      /* 4. use the constraint labels for the final variable labeling */
      SCIP_CALL( SCIPcomputeDecompVarsLabels(scip, decomp, conss, nconss) );

      SCIP_CALL( SCIPcomputeDecompStats(scip, decomp, TRUE) );

      SCIP_CALL( SCIPdecompstoreAdd(decompstore, decomp) );

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Transformed Decomposition statistics %d\n%s", d, SCIPdecompPrintStats(decomp, strbuf));
   }

   SCIPfreeBufferArray(scip, &varslabels);
   SCIPfreeBufferArray(scip, &origvars);
   SCIPfreeBufferArray(scip, &varssorted);

   return SCIP_OKAY;
}
