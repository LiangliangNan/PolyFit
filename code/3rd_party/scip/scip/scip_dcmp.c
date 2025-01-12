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

/**@file   scip_dcmp.c
 * @ingroup OTHER_CFILES
 * @brief  methods for working with decompositions
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include "scip/struct_dcmp.h"
#include "scip/debug.h"
#include "scip/dcmp.h"
#include "scip/mem.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_prob.h"
#include "scip/scip_var.h"
#include "scip/scip_mem.h"
#include "scip/struct_scip.h"
#include "scip/pub_cons.h"
#include "scip/pub_dcmp.h"
#include "scip/scip_dcmp.h"
#include "scip/scip_general.h"
#include "scip/scip_param.h"
#include "scip/scip_var.h"
#include "scip/scip_datastructures.h"
#include "scip/scip_message.h"


#define LABEL_UNASSIGNED INT_MIN /* label constraints or variables as unassigned. Only for internal use */

/** count occurrences of label in array, starting from pos */
static
int countLabelFromPos(
   int*                  labels,             /**< array of labels */
   int                   pos,                /**< position to start counting from */
   int                   nlabels             /**< the number of labels */
   )
{
   int endpos = pos;
   int currlabel;

   assert(labels != NULL);
   assert(pos < nlabels);

   currlabel = labels[pos];

   do
   {
      endpos++;
   }
   while( endpos < nlabels && labels[endpos] == currlabel );

   return endpos - pos;
}

/** raises an error if the condition is not TRUE */
static
SCIP_RETCODE ensureCondition(
   SCIP_Bool             condition           /**< some condition that must hold */
   )
{
   return condition ? SCIP_OKAY : SCIP_ERROR;
}

/** get variable buffer storage size for the buffer to be large enough to hold all variables */
static
int getVarbufSize(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   int norigvars;
   int ntransvars;

   norigvars = SCIPgetNOrigVars(scip);
   ntransvars = SCIPgetNVars(scip);

   return 2 * MAX(norigvars, ntransvars);
}

/** get variables and constraints of the original or transformed problem, to which the decomposition corresponds */
static
void getDecompVarsConssData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_VAR***           vars,               /**< pointer to store original/transformed variables array, or NULL */
   SCIP_CONS***          conss,              /**< pointer to store original/transformed constraints array, or NULL */
   int*                  nvars,              /**< pointer to store original/transformed variables array's length, or NULL */
   int*                  nconss              /**< pointer to store original/transformed constraints array's length, or NULL */
   )
{
   SCIP_Bool original;
   assert(scip != NULL);
   assert(decomp != NULL);

   original = SCIPdecompIsOriginal(decomp);

   if( vars )
      *vars = original ? SCIPgetOrigVars(scip) : SCIPgetVars(scip);

   if( nvars )
         *nvars = original ? SCIPgetNOrigVars(scip) : SCIPgetNVars(scip);

   if( conss )
      *conss = original ? SCIPgetOrigConss(scip) : SCIPgetConss(scip);

   if( nconss )
      *nconss = original ? SCIPgetNOrigConss(scip) : SCIPgetNConss(scip);
}

/** query the constraints active variables and their labels */
static
SCIP_RETCODE decompGetConsVarsAndLabels(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS*            cons,               /**< the constraint */
   SCIP_VAR**            varbuf,             /**< variable buffer array */
   int*                  labelbuf,           /**< buffer to store labels, or NULL if not needed */
   int                   bufsize,            /**< size of buffer arrays */
   int*                  nvars,              /**< pointer to store number of variables */
   int*                  requiredsize,       /**< pointer to store required size */
   SCIP_Bool*            success             /**< pointer to store whether variables and labels were successfully inserted */
   )
{
   SCIP_Bool success2;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(cons != NULL);
   assert(varbuf != NULL);
   assert(nvars != NULL);
   assert(requiredsize != NULL);
   assert(success != NULL);

   *success = FALSE;
   *requiredsize = 0;
   *nvars = 0;
   SCIP_CALL( SCIPgetConsNVars(scip, cons, nvars, &success2) );

   /* the constraint does not have the corresponding callback */
   if( ! success2 )
   {
      return SCIP_OKAY;
   }

   if( bufsize < *nvars )
   {
      *requiredsize = *nvars;

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPgetConsVars(scip, cons, varbuf, bufsize, &success2) );
   /* the constraint does not have the corresponding callback */
   if( ! success2 )
   {
      return SCIP_OKAY;
   }

   if( ! SCIPdecompIsOriginal(decomp) )
   {
      SCIP_CALL( SCIPgetActiveVars(scip, varbuf, nvars, bufsize, requiredsize) );

      if( *requiredsize > bufsize )
         return SCIP_OKAY;
   }
   else
   {
      int v;
      for( v = 0; v < *nvars; ++v )
      {
         assert(SCIPvarIsActive(varbuf[v]) || SCIPvarIsNegated(varbuf[v]));

         /* some constraint handlers such as indicator may already return inactive variables */
         if( SCIPvarIsNegated(varbuf[v]) )
            varbuf[v] = SCIPvarGetNegatedVar(varbuf[v]);
      }
   }

   /* get variables labels, if requested */
   if( labelbuf != NULL )
   {
      SCIPdecompGetVarsLabels(decomp, varbuf, labelbuf, *nvars);
   }

   *success = TRUE;

   return SCIP_OKAY;
}

/** creates a decomposition */
SCIP_RETCODE SCIPcreateDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   int                   nblocks,            /**< the number of blocks (without the linking block) */
   SCIP_Bool             original,           /**< is this a decomposition in the original (TRUE) or transformed space? */
   SCIP_Bool             benderslabels       /**< should the variables be labeled for the application of Benders' decomposition */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPdecompCreate(decomp, SCIPblkmem(scip), nblocks, original, benderslabels) );

   return SCIP_OKAY;
}

/** frees a decomposition */
void SCIPfreeDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP**         decomp              /**< pointer to free the decomposition data structure */
   )
{
   assert(scip != NULL);

   SCIPdecompFree(decomp, SCIPblkmem(scip));
}

/** adds decomposition to SCIP */
SCIP_RETCODE SCIPaddDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp              /**< decomposition to add */
   )
{
   assert(scip != NULL);
   assert(decomp != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddDecomp", FALSE, SCIPdecompIsOriginal(decomp), SCIPdecompIsOriginal(decomp),
      SCIPdecompIsOriginal(decomp), SCIPdecompIsOriginal(decomp), TRUE, TRUE, TRUE, !SCIPdecompIsOriginal(decomp),
      !SCIPdecompIsOriginal(decomp), !SCIPdecompIsOriginal(decomp), FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPdecompstoreAdd(scip->decompstore, decomp) );

   return SCIP_OKAY;
}

/** gets available user decompositions for either the original or transformed problem */
void SCIPgetDecomps(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP***        decomps,            /**< pointer to store decompositions array */
   int*                  ndecomps,           /**< pointer to store number of decompositions */
   SCIP_Bool             original            /**< should the decompositions for the original problem be returned? */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetDecomps", FALSE, original, original, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   if( decomps != NULL )
      *decomps = original ? SCIPdecompstoreGetOrigDecomps(scip->decompstore) : SCIPdecompstoreGetDecomps(scip->decompstore);

   if( ndecomps != NULL )
      *ndecomps = original ? SCIPdecompstoreGetNOrigDecomps(scip->decompstore) : SCIPdecompstoreGetNDecomps(scip->decompstore);
}

/** returns TRUE if the constraint \p cons contains only linking variables in decomposition \p decomp */
SCIP_RETCODE SCIPhasConsOnlyLinkVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS*            cons,               /**< the constraint */
   SCIP_Bool*            hasonlylinkvars     /**< will be set to TRUE if this constraint has only linking variables */
   )
{
   SCIP_VAR** consvars;
   int nvars;
   int i;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(decomp != NULL);
   assert(hasonlylinkvars != NULL);

   SCIP_CALL( SCIPgetConsNVars(scip, cons, &nvars, &success) );
   SCIP_CALL( ensureCondition(success) );

   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );

   SCIP_CALL( SCIPgetConsVars(scip, cons, consvars, nvars, &success) );
   SCIP_CALL( ensureCondition(success) );

   if( ! SCIPdecompIsOriginal(decomp) )
   {
      int requiredsize;
      SCIP_CALL( SCIPgetActiveVars(scip, consvars, &nvars, nvars, &requiredsize) );
      assert(requiredsize <= nvars);
   }

   *hasonlylinkvars = TRUE;
   /* check if variables are all linking variables */
   for( i = 0; i < nvars && *hasonlylinkvars; ++i )
   {
      int label;

      SCIPdecompGetVarsLabels(decomp, &consvars[i], &label, 1);

      *hasonlylinkvars = (label == SCIP_DECOMP_LINKVAR);
   }

   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}


/** computes constraint labels from variable labels
 *
 *  Existing labels for the constraints are simply overridden
 *
 *  The computed labels depend on the flag SCIPdecompUseBendersLabels() of the decomposition. If the flag is set
 *  to FALSE, the labeling assigns
 *
 *  - label i, if only variables labeled i are present in the constraint (and optionally linking variables)
 *  - SCIP_DECOMP_LINKCONS, if there are either only variables labeled with SCIP_DECOMP_LINKVAR present, or
 *    if there are variables with more than one block label.
 *
 *  If the flag is set to TRUE, the assignment is the same, unless variables from 2 named blocks occur in the same
 *  constraint, which is an invalid labeling for the Benders case.
 */
SCIP_RETCODE SCIPcomputeDecompConsLabels(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_VAR** varbuffer;
   int c;
   int varbufsize;
   int* varlabels;
   int* conslabels;
   SCIP_Bool benderserror;
   SCIP_Bool benderslabels;

   assert(decomp != NULL);

   if( nconss == 0 )
      return SCIP_OKAY;

   varbufsize = getVarbufSize(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &varbuffer, varbufsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, varbufsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );

   benderslabels = SCIPdecompUseBendersLabels(decomp);
   benderserror = FALSE;

   /* assign label to each individual constraint */
   for( c = 0; c < nconss && ! benderserror; ++c )
   {
      int nconsvars;
      int v;
      int nlinkingvars = 0;
      int conslabel;
      int requiredsize;

      SCIP_Bool success;

      SCIP_CALL( decompGetConsVarsAndLabels(scip, decomp, conss[c], varbuffer, varlabels,
            varbufsize, &nconsvars, &requiredsize, &success) );
      SCIP_CALL( ensureCondition(success) );

      /* loop over variable labels to compute the constraint label */
      conslabel = LABEL_UNASSIGNED;
      for( v = 0; v < nconsvars; ++v )
      {
         int varlabel = varlabels[v];

         /* count the number of linking variables, and keep track if there are two variables with different labels */
         if( varlabel == SCIP_DECOMP_LINKVAR )
            ++nlinkingvars;
         else if( conslabel == LABEL_UNASSIGNED )
            conslabel = varlabel;
         else if( conslabel != varlabel )
         {
            /* there must not be two variables from different named blocks in a single constraint, since the presence
             * of named block variables forbids this constraint from the master (linking) block
             */
            if( benderslabels )
               benderserror = TRUE;

            conslabel = SCIP_DECOMP_LINKCONS;
            break;
         }
      }

      assert(nlinkingvars == nconsvars || conslabel != LABEL_UNASSIGNED);
      assert(v == nconsvars || conslabel == SCIP_DECOMP_LINKCONS);

      /* if there are only linking variables, the constraint is unassigned */
      if( conslabel == LABEL_UNASSIGNED )
         conslabel = SCIP_DECOMP_LINKCONS;
      conslabels[c] = conslabel;
   }

   SCIP_CALL( SCIPdecompSetConsLabels(decomp, conss, conslabels, nconss) );

   SCIPfreeBufferArray(scip, &conslabels);
   SCIPfreeBufferArray(scip, &varlabels);
   SCIPfreeBufferArray(scip, &varbuffer);

   /* throw an error and inform the user if the variable block decomposition does not allow a benders constraint labeling */
   if( benderserror )
   {
      SCIPerrorMessage("Error in constraint label computation; variables from multiple named blocks in a single constraint\n");

      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** creates a decomposition of the variables from a labeling of the constraints
 *
 *  NOTE: by default, the variable labeling is based on a Dantzig-Wolfe decomposition. This means that constraints in named
 *  blocks have have precedence over linking constraints. If a variable exists in constraints from
 *  two or more named blocks, then this variable is marked as a linking variable.
 *  If a variable occurs in exactly one named block i>=0, it is assigned label i.
 *  Variables which are only in linking constraints are unlabeled. However, SCIPdecompGetVarsLabels() will
 *  label them as linking variables.
 *
 *  If the variables should be labeled for the application of Benders' decomposition, the decomposition must be
 *  flagged explicitly via SCIPdecompSetUseBendersLabels().
 *  With this setting, the presence in linking constraints takes precedence over the presence in named blocks.
 *  Now, a variable is considered linking if it is present in at least one linking constraint and an arbitrary
 *  number of constraints from named blocks.
 */
SCIP_RETCODE SCIPcomputeDecompVarsLabels(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   )
{
   int c;
   int* conslabels;
   SCIP_VAR** varbuffer;
   int varbufsize;
   SCIP_Bool benderslabels;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(conss != NULL);

   /* make the buffer array large enough such that we do not have to reallocate buffer */
   varbufsize = getVarbufSize(scip);

   /* allocate buffer to store constraint variables and labels */
   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varbuffer, varbufsize) );

   /* query constraint labels */
   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);

   benderslabels = SCIPdecompUseBendersLabels(decomp);
   /* iterate over constraints and query the corresponding constraint labels */
   for( c = 0; c < nconss; ++c )
   {
      int conslabel;
      int v;
      int nconsvars;
      SCIP_Bool success;
      int requiredsize;
      int newvarlabel;

      conslabel = conslabels[c];

      if( conslabel == SCIP_DECOMP_LINKCONS )
      {
         /* skip linking constraints unless Benders labeling is used */
         if( ! benderslabels )
            continue;
         else
            newvarlabel = SCIP_DECOMP_LINKVAR;
      }
      else
         newvarlabel = conslabel;

      SCIP_CALL( decompGetConsVarsAndLabels(scip, decomp, conss[c], varbuffer, NULL,
            varbufsize, &nconsvars, &requiredsize, &success) );
      SCIP_CALL( ensureCondition(success) );

      /* each variable in this constraint gets the constraint label unless it already has a different label -> make it a linking variable */
      for( v = 0; v < nconsvars; ++v )
      {
         SCIP_VAR* var = varbuffer[v];

         assert(SCIPvarIsActive(var) || (SCIPdecompIsOriginal(decomp) && SCIPvarIsNegated(var)));

         /* some constraint handlers such as indicator may already return inactive variables */
         if( SCIPvarIsNegated(var) )
            var = SCIPvarGetNegatedVar(var);

         if( SCIPhashmapExists(decomp->var2block, (void *)var) )
         {
            int varlabel = SCIPhashmapGetImageInt(decomp->var2block, (void *)var);

            /* store the label linking variable explicitly to distinguish it from the default */
            if( varlabel != SCIP_DECOMP_LINKVAR && varlabel != newvarlabel )
               SCIP_CALL( SCIPhashmapSetImageInt(decomp->var2block, (void *)var, SCIP_DECOMP_LINKVAR) );
         }
         else
         {
            SCIP_CALL( SCIPhashmapInsertInt(decomp->var2block, (void *)var, newvarlabel) );
         }
      }
   }

   SCIPfreeBufferArray(scip, &varbuffer);
   SCIPfreeBufferArray(scip, &conslabels);

   return SCIP_OKAY;
}

/** assigns linking constraints to blocks
 *
 * Each linking constraint is assigned to the most frequent block among its variables.
 * Variables of other blocks are relabeled as linking variables.
 * Constraints that have only linking variables are skipped.
 *
 * @note: In contrast to SCIPcomputeDecompConsLabels(), this method potentially relabels variables.
 */
SCIP_RETCODE SCIPassignDecompLinkConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of linking constraints that should be reassigned */
   int                   nconss,             /**< number of constraints */
   int*                  nskipconss          /**< pointer to store the number of constraints that were skipped, or NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** allvars;
   int* varslabels;
   int requiredsize;
   int nvars;
   int nconsvars;
   int varbufsize;
   int c;
   int nskipconsslocal;
   int defaultlabel;

   assert(scip != NULL);
   assert(decomp != NULL);

   nvars = SCIPgetNVars(scip);
   varbufsize = getVarbufSize(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, varbufsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, varbufsize) );

   /* get one label as default label */
   allvars = SCIPdecompIsOriginal(decomp) ? SCIPgetOrigVars(scip) : SCIPgetVars(scip);
   SCIPdecompGetVarsLabels(decomp, allvars, varslabels, nvars);
   for( c = 0; c < nvars; c++ )
   {
      if( varslabels[c] != SCIP_DECOMP_LINKVAR )
      {
         defaultlabel = varslabels[c];
         break;
      }
   }

   nskipconsslocal = 0;
   for( c = 0; c < nconss; c++ )
   {
      SCIP_Bool success;

      SCIP_CALL( decompGetConsVarsAndLabels(scip, decomp, conss[c], vars, varslabels, varbufsize,
               &nconsvars, &requiredsize, &success) );
      SCIP_CALL( ensureCondition(success) );

      SCIPsortIntPtr(varslabels, (void **)vars, nconsvars);
      /* constraint contains no (active) variables */
      if( nconsvars == 0 )
      {
         SCIP_CALL( SCIPdecompSetConsLabels(decomp, &conss[c], &defaultlabel, 1) );
      }
      /* constraint contains only linking variables */
      else if( varslabels[nconsvars - 1] == SCIP_DECOMP_LINKVAR )
      {
         nskipconsslocal++;

         continue;
      }
      else
      {
         int startposs[2];
         int endposs[2];
         int nlinkvars;
         int block;
         int maxnblockvars;
         int curr;
         int v;
         int p;

         /* count linking variables */
         if( varslabels[0] == SCIP_DECOMP_LINKVAR )
         {
            nlinkvars = countLabelFromPos(varslabels, 0, nconsvars);
         }
         else
         {
            nlinkvars = 0;
         }

         assert(nlinkvars < nconsvars);

         curr = nlinkvars;
         /* find the most frequent block label among the nonlinking variables */
         maxnblockvars = 0;
         block = -1;
         do
         {
            int nblockvars = countLabelFromPos(varslabels, curr, nconsvars);
            if( nblockvars > maxnblockvars )
            {
               maxnblockvars = nblockvars;
               block = curr;
            }
            curr += nblockvars;
         }
         while( curr < nconsvars );

         /* reassign all variables from other blocks as linking variables */
         startposs[0] = nlinkvars;
         endposs[0] = block;
         startposs[1] = block + maxnblockvars;
         endposs[1] = nconsvars;

         /* loop over all variables before (p==0) and after (p==1) the most frequent block label */
         for( p = 0; p < 2; ++p )
         {
            /* relabel */
            for( v = startposs[p]; v < endposs[p]; ++v )
               varslabels[v] = SCIP_DECOMP_LINKVAR;

            /* set labels in the decomposition */
            SCIP_CALL( SCIPdecompSetVarsLabels(decomp, &vars[startposs[p]], &varslabels[startposs[p]], endposs[p] - startposs[p]) );
         }

         SCIP_CALL( SCIPdecompSetConsLabels(decomp, &conss[c], &varslabels[block], 1) );
      }
   }

   if( nskipconss != NULL )
      *nskipconss = nskipconsslocal;

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &varslabels);

   return SCIP_OKAY;
}

/** return position of a label in decomp array */
static
int findLabelIdx(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   int                   label               /**< the label */
   )
{
   int pos;

   (void)SCIPsortedvecFindInt(decomp->labels, label, decomp->nblocks + 1, &pos);

   return pos;
}

/** compute decomposition modularity (comparison of within block edges against a random decomposition) */
static
SCIP_RETCODE computeModularity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Real*            modularity          /**< pointer to store modularity value */
   )
{
   SCIP_CONS** conss;
   SCIP_VAR** varbuf;
   int* varslabels;
   int* conslabels;
   int* totaldegrees; /* the total degree for every block */
   int* withinedges; /* the number of edges within each block */
   int nnonzeroes = 0;
   int varbufsize;
   int nconss;
   int c;
   int b;

   /* allocate buffer arrays to hold constraint and variable labels, and store within-edges and total community degrees */
   getDecompVarsConssData(scip, decomp, NULL, &conss, NULL, &nconss);
   varbufsize = getVarbufSize(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, varbufsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varbuf, varbufsize) );

   SCIP_CALL( SCIPallocClearBufferArray(scip, &totaldegrees, decomp->nblocks + 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &withinedges, decomp->nblocks + 1) );

   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);

   /*
    * loop over all nonzeros, consider the labels of the incident nodes (cons and variable)
    * and increase the corresponding counters
    */
   for( c = 0; c < nconss; ++c )
   {
      int nconsvars;
      int conslabel = conslabels[c];
      int blockpos;
      int varblockstart;
      int requiredsize;
      SCIP_Bool success;

      /* linking constraints do not contribute to the modularity */
      if( conslabel == SCIP_DECOMP_LINKCONS )
         continue;

      /* find the position of the constraint label. Constraints of the border always belong to the first block at index 0 */
      blockpos = findLabelIdx(decomp, conslabel);

      SCIP_CALL( decompGetConsVarsAndLabels(scip, decomp, conss[c], varbuf, varslabels,
               varbufsize, &nconsvars, &requiredsize, &success) );
      SCIP_CALL( ensureCondition(success) );

      SCIPsortInt(varslabels, nconsvars);

      /* count occurrences of labels (blocks) in the sorted labels array */
      varblockstart = 0;
      while( varblockstart < nconsvars )
      {
         int varblockpos;
         int nblockvars = countLabelFromPos(varslabels, varblockstart, nconsvars);

         varblockpos = findLabelIdx(decomp, varslabels[varblockstart]);

         /* don't consider linking variables for modularity statistics */
         if( varslabels[varblockstart] != SCIP_DECOMP_LINKVAR )
         {
            /* increase the number of within edges for variable and constraints from the same block */
            if( varblockpos == blockpos )
               withinedges[varblockpos] += nblockvars;

            /* increase the total degrees and nonzero (edge) counts; it is intended that the total degrees sum up
             * to twice the number of edges
             */
            totaldegrees[blockpos] += nblockvars;
            totaldegrees[varblockpos] += nblockvars;
            nnonzeroes += nblockvars;
         }

         varblockstart += nblockvars;
      }
   }

/* ensure that total degrees sum up to twice the number of edges */
#ifndef NDEBUG
   {
      int totaldegreesum = 0;
      for( b = 1; b < decomp->nblocks + 1; ++b )
         totaldegreesum += totaldegrees[b];

      assert(totaldegreesum == 2 * nnonzeroes);
   }
#endif

   /* compute modularity */
   *modularity = 0.0;
   nnonzeroes = MAX(nnonzeroes, 1);
   for( b = 1; b < decomp->nblocks + 1; ++b )
   {
      SCIP_Real expectedval;
      expectedval = totaldegrees[b] / (2.0 * nnonzeroes);
      expectedval = SQR(expectedval);
      *modularity += (withinedges[b] / (SCIP_Real)nnonzeroes) - expectedval;
   }

   SCIPfreeBufferArray(scip, &withinedges);
   SCIPfreeBufferArray(scip, &totaldegrees);
   SCIPfreeBufferArray(scip, &varbuf);
   SCIPfreeBufferArray(scip, &varslabels);
   SCIPfreeBufferArray(scip, &conslabels);

   return SCIP_OKAY;
}

/** compute the area score */
static
void computeAreaScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   SCIP_Real areascore = 1.0;
   int nvars;
   int nconss;

   getDecompVarsConssData(scip, decomp, NULL, NULL, &nvars, &nconss);

   if( nvars > 0 && nconss > 0 )
   {
      int nlinkconss = decomp->consssize[0];
      int nlinkvars = decomp->varssize[0];
      SCIP_Real factor = 1.0 / ((SCIP_Real)nvars * nconss);

      int i;

      /* compute diagonal contributions to the area score */
      for( i = 1; i < decomp->nblocks + 1; ++i )
      {
         areascore -= (factor * decomp->consssize[i]) * decomp->varssize[i];
      }

      areascore -= ((SCIP_Real)nlinkconss * nvars + (SCIP_Real)nconss * nlinkvars - (SCIP_Real)nlinkconss * nlinkvars) * factor;
   }

   decomp->areascore = areascore;
}

/** build the block decomposition graph */
static
SCIP_RETCODE buildBlockGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   int                   maxgraphedge        /**< maximum number of edges in block graph computation (-1: no limit, 0: disable block graph computation) */
   )
{
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_VAR** consvars;
   SCIP_DIGRAPH* blocklinkingvargraph;
   SCIP_DIGRAPH* blockgraph = NULL;
   int* varlabels;
   int* conslabels;
   SCIP_CONS** consscopy; /* working copy of the constraints */
   int* linkvaridx;
   int* succnodesblk;
   int* succnodesvar;
   SCIP_Bool success;
   int varbufsize;
   int nvars;
   int nconss;
   int nblocks;
   int nlinkingvars = 0;
   int nconsvars;
   int conspos;
   int tempmin;
   int tempmax;
   int nblockgraphedges;
   int blocknodeidx;
   int i;
   int j;
   int v;
   int n;

   if( maxgraphedge == -1 )
      maxgraphedge = INT_MAX;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(decomp->statscomplete == FALSE);

   /* capture the trivial case that no linking variables are present */
   if( decomp->varssize[0] == 0 || decomp->nblocks == 0 )
   {
      decomp->mindegree = 0;
      decomp->maxdegree = 0;
      decomp->nedges = 0;
      decomp->ncomponents = SCIPdecompGetNBlocks(decomp);
      decomp->narticulations = 0;
      decomp->statscomplete = TRUE;

      return SCIP_OKAY;
   }

   getDecompVarsConssData(scip, decomp, &vars, &conss, &nvars, &nconss);
   varbufsize = getVarbufSize(scip);
   nblocks = SCIPdecompGetNBlocks(decomp);

   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, varbufsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linkvaridx, varbufsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, varbufsize) );

   /* store variable and constraint labels in buffer arrays */
   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);
   SCIPdecompGetVarsLabels(decomp, vars, varlabels, nvars);

   /* create a mapping of all linking variables to 0,..,nlinkingvars -1 and store it in array linkvaridx */
   for( v = 0; v < nvars; ++v )
   {
      if( varlabels[v] == SCIP_DECOMP_LINKVAR )
      {
         linkvaridx[v] = nlinkingvars;
         assert(SCIPvarGetProbindex(vars[v]) == v);
         ++nlinkingvars;
      }
      else
         linkvaridx[v] = -1;
   }

   /* create a bipartite graph composed of block and linking var nodes */
   SCIP_CALL( SCIPcreateDigraph(scip, &blocklinkingvargraph, nblocks + nlinkingvars) );/* nblocks does not include the linking constraints block */

   /* initialize position to start after the linking constraints, which we skip anyway */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &consscopy, conss, nconss) );
   SCIPsortIntPtr(conslabels, (void**)consscopy, nconss);
   if( conslabels[0] == SCIP_DECOMP_LINKCONS )
      conspos = countLabelFromPos(conslabels, 0, nconss);
   else
      /* no linking constraints present */
      conspos = 0;

   blocknodeidx = -1;
   /* loop over each block */
   while( conspos < nconss )
   {
      SCIP_Bool* adjacent;
      int* adjacentidxs;
      int nblockconss = countLabelFromPos(conslabels, conspos, nconss);
      int nblocklinkingvars = 0;
      int c;

      ++blocknodeidx;
      /* allocate buffer storage to store all linking variable indices adjacent to this block */
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &adjacent, nlinkingvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &adjacentidxs, nlinkingvars) );

      /* loop over the constraints of this block; stop if the block vertex has maximum degree */
      for( c = conspos; c < conspos + nblockconss && nblocklinkingvars < nlinkingvars; ++c )
      {
         int requiredsize;
         SCIP_CONS* cons = consscopy[c];
         assert(conslabels[c] != SCIP_DECOMP_LINKCONS);

         SCIP_CALL( decompGetConsVarsAndLabels(scip, decomp, cons, consvars, varlabels,
               varbufsize, &nconsvars, &requiredsize, &success) );
         SCIP_CALL( ensureCondition(success) );

         /* search for linking variables that are not connected so far; stop as soon as block vertex has max degree */
         for( j = 0; j < nconsvars && nblocklinkingvars < nlinkingvars; ++j )
         {
            int linkingvarnodeidx;
            /* consider only linking variables */
            if( varlabels[j] != SCIP_DECOMP_LINKVAR )
               continue;

            linkingvarnodeidx = linkvaridx[SCIPvarGetProbindex(consvars[j])];
            assert(linkingvarnodeidx >= 0);

            if( !adjacent[linkingvarnodeidx] )
            {
               adjacent[linkingvarnodeidx] = TRUE;
               adjacentidxs[nblocklinkingvars++] = linkingvarnodeidx;
            }
         }
      }

      /* connect block and linking variables in the digraph */
      assert(blocknodeidx == findLabelIdx(decomp, conslabels[conspos]) - 1);
      for( i = 0; i < nblocklinkingvars; ++i )
      {
         SCIP_CALL( SCIPdigraphAddArc(blocklinkingvargraph, blocknodeidx, nblocks + adjacentidxs[i], NULL) );
         SCIP_CALL( SCIPdigraphAddArc(blocklinkingvargraph, nblocks + adjacentidxs[i], blocknodeidx, NULL) );
      }

      /* clean up the adjacent array before freeing */
      for( i = 0; i < nblocklinkingvars; ++i )
         adjacent[adjacentidxs[i]] = FALSE;

      /* check that adjacent has been entirely cleaned up */
#ifndef NDEBUG
      for( i = 0; i < nlinkingvars; ++i )
         assert(adjacent[i] == FALSE);
#endif

      SCIPfreeBufferArray(scip, &adjacentidxs);
      SCIPfreeCleanBufferArray(scip, &adjacent);

      conspos += nblockconss;
   }

   SCIPfreeBufferArray(scip, &consscopy);

   assert(SCIPdigraphGetNNodes(blocklinkingvargraph) > 0);
   /* check first if any of the linking variables is connected with all blocks -> block graph is complete and connected */
   for( n = nblocks; n < SCIPdigraphGetNNodes(blocklinkingvargraph); ++n )
   {
      int nsuccvar;
      nsuccvar = (int) SCIPdigraphGetNSuccessors(blocklinkingvargraph, n);

      if( nsuccvar == nblocks )
      {
         decomp->nedges = nblocks * (nblocks - 1) / 2;
         decomp->mindegree = decomp->maxdegree = nblocks - 1;
         decomp->narticulations = 0;
         decomp->ncomponents = 1;
         decomp->statscomplete = TRUE;

         goto TERMINATE;
      }
   }

   /* from the information of the above bipartite graph, build the block-decomposition graph: nodes -> blocks and double-direction arcs -> linking variables */
   SCIP_CALL( SCIPcreateDigraph(scip, &blockgraph, nblocks) );

   /* we count the number of block graph edges manually, because SCIPdigraphGetNArcs() iterates over all nodes */
   nblockgraphedges = 0;
   for( n = 0; n < nblocks - 1 && nblockgraphedges < maxgraphedge; ++n )
   {
      SCIP_Bool* adjacent; /* an array to mark the adjacent blocks to the current block */
      int* adjacentidxs;
      int nsuccblk;
      int nadjacentblks = 0;
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &adjacent, nblocks) );
      SCIP_CALL( SCIPallocBufferArray(scip, &adjacentidxs, nblocks) );

      assert(n < SCIPdigraphGetNNodes(blocklinkingvargraph));

      /* loop over the connected linking variables to the current block and their connected blocks to update the adjacency array */
      nsuccblk = SCIPdigraphGetNSuccessors(blocklinkingvargraph, n);
      succnodesblk = SCIPdigraphGetSuccessors(blocklinkingvargraph, n);
      for( i = 0; i < nsuccblk && nadjacentblks < nblocks - (n + 1); ++i )
      {
         int startpos;
         int nsuccvar;

         assert(succnodesblk[i] < SCIPdigraphGetNNodes(blocklinkingvargraph));

         nsuccvar = SCIPdigraphGetNSuccessors(blocklinkingvargraph, succnodesblk[i]);
         succnodesvar = SCIPdigraphGetSuccessors(blocklinkingvargraph, succnodesblk[i]);

         /* previously visited blocks can be skipped in this step */
         if( ! SCIPsortedvecFindInt(succnodesvar, n, nsuccvar, &startpos) )
            SCIPABORT();
         for( j = startpos + 1; j < nsuccvar; ++j )
         {
            assert( succnodesvar[j] > n );
            if( !adjacent[succnodesvar[j]] )
            {
               adjacent[succnodesvar[j]] = TRUE;
               adjacentidxs[nadjacentblks] = succnodesvar[j];
               ++nadjacentblks;
            }
         }
      }

      /* double-direction arcs are added in this step between the current block and its adjacent block nodes */
      for( i = 0; i < nadjacentblks && nblockgraphedges < maxgraphedge; ++i )
      {
          SCIP_CALL( SCIPdigraphAddArc(blockgraph, n, adjacentidxs[i], NULL) );
          SCIP_CALL( SCIPdigraphAddArc(blockgraph, adjacentidxs[i], n, NULL) );

          ++nblockgraphedges;
      }

      /* clean up the adjacent array and free it */
      for( i = 0; i < nadjacentblks; ++i )
         adjacent[adjacentidxs[i]] = FALSE;

      SCIPfreeBufferArray(scip, &adjacentidxs);
      SCIPfreeCleanBufferArray(scip, &adjacent);
   }

   assert(SCIPdigraphGetNNodes(blockgraph) > 0);

   /* Get the number of edges in the block-decomposition graph.*/
   decomp->nedges = nblockgraphedges;
   decomp->statscomplete = nblockgraphedges < maxgraphedge;

   /* Get the minimum and maximum degree of the block-decomposition graph */
   tempmin = (int) SCIPdigraphGetNSuccessors(blockgraph, 0);
   tempmax = (int) SCIPdigraphGetNSuccessors(blockgraph, 0);
   for( n = 1; n < SCIPdigraphGetNNodes(blockgraph); ++n )
   {
      int nsuccblk = SCIPdigraphGetNSuccessors(blockgraph, n);

      if( nsuccblk < tempmin )
         tempmin = nsuccblk;
      else if( nsuccblk > tempmax )
         tempmax = nsuccblk;
   }

   decomp->mindegree = tempmin;
   decomp->maxdegree = tempmax;

   /* Get the number of connected components in the block-decomposition graph.*/
   SCIP_CALL( SCIPdigraphComputeUndirectedComponents(blockgraph, -1, NULL, NULL) );
   decomp->ncomponents = SCIPdigraphGetNComponents(blockgraph);

   /* Get the number of articulation points in the block-decomposition graph using DFS.*/
   SCIP_CALL( SCIPdigraphGetArticulationPoints(blockgraph, NULL, &decomp->narticulations) );

TERMINATE:
   SCIPfreeBufferArray(scip, &consvars);
   SCIPfreeBufferArray(scip, &linkvaridx);
   SCIPfreeBufferArray(scip, &varlabels);
   SCIPfreeBufferArray(scip, &conslabels);

   /* blockgraph has probably not been allocated */
   if( blockgraph != NULL)
      SCIPdigraphFree(&blockgraph);

   SCIPdigraphFree(&blocklinkingvargraph);

   return SCIP_OKAY;
}

/** computes decomposition statistics and store them in the decomposition object */
SCIP_RETCODE SCIPcomputeDecompStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Bool             uselimits           /**< respect user limits on potentially expensive graph statistics? */
   )
{
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_CONS** conssarray;
   SCIP_VAR** varsarray;
   int* varslabels;
   int* conslabels;
   int nvars;
   int nconss;
   int varblockstart;
   int consblockstart;
   int currlabelidx;
   int varidx;
   int considx;
   int i;
   int maxgraphedge;
   SCIP_Bool disablemeasures;

   assert(scip != NULL);
   assert(decomp != NULL);

   getDecompVarsConssData(scip, decomp, &vars, &conss, &nvars, &nconss);

   /* return if problem is empty
    *
    * TODO ensure that statistics reflect this correctly
    */
   if( nvars == 0 || nconss == 0 )
   {
      decomp->nblocks = 0;
      decomp->varssize[0] = nvars;
      decomp->consssize[0] = nconss;
      decomp->labels[0] = SCIP_DECOMP_LINKVAR;
      return SCIP_OKAY;
   }

   decomp->statscomplete = FALSE;

   /* store variable and constraint labels in buffer arrays */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &conssarray, conss, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &varsarray, vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, nvars) );

   SCIPdecompGetVarsLabels(decomp, varsarray, varslabels, nvars);
   SCIPdecompGetConsLabels(decomp, conssarray, conslabels, nconss);

   /* sort both buffer arrays for quick counting */
   SCIPsortIntPtr(varslabels, (void**)varsarray, nvars);
   SCIPsortIntPtr(conslabels, (void**)conssarray, nconss);

   /* the first label is always LINKVAR, even if Benders' variable labels are used. We can ignore the variables
    * labelled as LINKCONS since this label is only required when computing the variable labels for Benders'
    * decomposition.
    */
   decomp->labels[0] = SCIP_DECOMP_LINKVAR;

   /* treating the linking variables first */
   if( varslabels[0] == SCIP_DECOMP_LINKVAR )
      decomp->varssize[0] = countLabelFromPos(varslabels, 0, nvars);
   else
      decomp->varssize[0] = 0;

   /* count border constraints and store their number */
   if( conslabels[0] == SCIP_DECOMP_LINKCONS )
      decomp->consssize[0] = countLabelFromPos(conslabels, 0, nconss);
   else
      decomp->consssize[0] = 0;

   /* merge labels (except for border at position 0) since neither variable nor constraint labels by themselves need to be complete */
   currlabelidx = 1;
   varidx = decomp->varssize[0];
   considx = decomp->consssize[0];

   while( varidx < nvars || considx < nconss )
   {
      int varlabel;
      int conslabel;

      varlabel = varidx < nvars ? varslabels[varidx] : INT_MAX;
      conslabel = considx < nconss ? conslabels[considx] : INT_MAX;

      assert(currlabelidx < decomp->memsize);
      /* store the smaller of the two current labels */
      decomp->labels[currlabelidx] = MIN(varlabel, conslabel);

      /* a strictly larger variable label means that there are no variables for the current label */
      if( varlabel <= conslabel )
         decomp->varssize[currlabelidx] = countLabelFromPos(varslabels, varidx, nvars);
      else
         decomp->varssize[currlabelidx] = 0;

      /* the same for constraint labels */
      if( conslabel <= varlabel )
         decomp->consssize[currlabelidx] = countLabelFromPos(conslabels, considx, nconss);
      else
         decomp->consssize[currlabelidx] = 0;

      /* increase indices appropriately */
      varidx += decomp->varssize[currlabelidx];
      considx += decomp->consssize[currlabelidx];

      currlabelidx++;
   }

   SCIPdebugMsg(scip, "Counted %d different labels (should be %d)\n", currlabelidx, decomp->nblocks + 1);

   /* strip the remaining, unused blocks */
   if( currlabelidx < decomp->nblocks + 1 )
      decomp->nblocks = currlabelidx - 1;

   /* delete empty blocks from statistics, relabel the corresponding constraints/variables as linking */
   varblockstart = decomp->varssize[0];
   consblockstart = decomp->consssize[0];

   for( i = 1; i < decomp->nblocks + 1; ++i )
   {
      assert(MAX(decomp->varssize[i], decomp->consssize[i]) > 0);
      /* relabel constraint blocks as linking, if there are no corresponding variables */
      if( decomp->varssize[i] == 0 )
      {
         int nblockconss = decomp->consssize[i];
         int c;
         /* relabel these constraints as linking */
         for( c = consblockstart; c < consblockstart + nblockconss; ++c )
            conslabels[c] = SCIP_DECOMP_LINKCONS;

         SCIP_CALL( SCIPdecompSetConsLabels(decomp, &conssarray[consblockstart], &conslabels[consblockstart], nblockconss) );

         /* increase number of linking constraints */
         decomp->consssize[0] += nblockconss;
      }

      /* same for constraints */
      if( decomp->consssize[i] == 0 )
      {
         int nblockvars = decomp->varssize[i];
         int v;

         /* relabel the variables as linking variables */
         for( v = varblockstart; v < varblockstart + nblockvars; ++v )
            varslabels[v] = SCIP_DECOMP_LINKVAR;

         SCIP_CALL( SCIPdecompSetVarsLabels(decomp, &varsarray[varblockstart], &varslabels[varblockstart], nblockvars) );

         /* increase number of linking variables */
         decomp->varssize[0] += nblockvars;
      }

      varblockstart += decomp->varssize[i];
      consblockstart += decomp->consssize[i];
   }

   currlabelidx = 1;

   /* delete empty blocks; they are no longer present */
   for( i = 1; i < decomp->nblocks + 1; ++i )
   {
      /* keep only nonempty blocks */
      if( decomp->varssize[i] > 0 && decomp->consssize[i] > 0 )
      {
         decomp->labels[currlabelidx] = decomp->labels[i];
         decomp->varssize[currlabelidx] = decomp->varssize[i];
         decomp->consssize[currlabelidx] = decomp->consssize[i];

         currlabelidx++;
      }
   }

   decomp->nblocks = currlabelidx - 1;

   decomp->idxsmallestblock = decomp->idxlargestblock = -1;
   /* now that indices are fixed, store indices with largest and smallest number of constraints */
   for( i = 1; i < decomp->nblocks + 1; ++i )
   {
      if( decomp->idxsmallestblock == -1 )
         decomp->idxsmallestblock = decomp->idxlargestblock = i;
      else if( decomp->consssize[decomp->idxsmallestblock] > decomp->consssize[i] )
         decomp->idxsmallestblock = i;
      else if( decomp->consssize[decomp->idxlargestblock] < decomp->consssize[i] )
         decomp->idxlargestblock = i;
   }

   /* compute more involved statistics such as the area score, the modularity, and the block graph statistics */
   SCIP_CALL( SCIPgetBoolParam(scip, "decomposition/disablemeasures", &disablemeasures) );
   if( !disablemeasures )
   {
      SCIP_CALL( computeModularity(scip, decomp, &decomp->modularity) );
      computeAreaScore(scip, decomp);
   }

   if( uselimits )
   {
      SCIP_CALL( SCIPgetIntParam(scip, "decomposition/maxgraphedge", &maxgraphedge) );
   }
   else
      maxgraphedge = -1;

   /* do not start computation of the block graph if maxgraphedge is set to 0 */
   if( maxgraphedge != 0 )
   {
      SCIP_CALL( buildBlockGraph(scip, decomp, maxgraphedge) );
   }

   SCIPfreeBufferArray(scip, &varslabels);
   SCIPfreeBufferArray(scip, &varsarray);
   SCIPfreeBufferArray(scip, &conslabels);
   SCIPfreeBufferArray(scip, &conssarray);

   return SCIP_OKAY;
}
