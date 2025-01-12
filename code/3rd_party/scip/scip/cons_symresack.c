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

/**@file   cons_symresack.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for symresack constraints
 * @author Christopher Hojny
 * @author Jasper van Doornmalen
 *
 * The type of constraints of this constraint handler is described in cons_symresack.h.
 *
 * The details of the method implemented here are described in the following papers:
 *
 * Fundamental Domains for Integer Programs with Symmetries@n
 * Eric J. Friedman,@n
 * Combinatorial Optimization, volume 4616 of LNCS, 146-153 (2007)
 *
 * This paper describes an inequality to handle symmetries of a single permutation. This
 * so-called FD-inequality is the basic for the propagation routine of our implementation.
 *
 * Polytopes Associated with Symmetry Handling@n
 * Christopher Hojny and Marc E. Pfetsch,@n
 * Mathematical Programming 175, No. 1, 197-240, 2019
 *
 * This paper describes an almost linear time separation routine for so-called cover
 * inequalities of symresacks. A slight modification of this algorithm allows for a linear
 * running time, which is used in this implementation.
 *
 * Packing, Partitioning, and Covering Symresacks@n
 * Christopher Hojny,@n
 * (2020), available at https://doi.org/10.1016/j.dam.2020.03.002
 * Discrete Applied Mathematics, volume 283, 689-717 (2020)
 *
 * This paper introduces linearly many inequalities with ternary coefficients that suffice to
 * characterize the binary points contained in a packing and partitioning symresack completely.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_orbisack.h"
#include "scip/cons_setppc.h"
#include "scip/cons_symresack.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip.h"
#include "scip/scip_branch.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cons.h"
#include "scip/scip_cut.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_var.h"
#include <ctype.h>
#include <string.h>

/* constraint handler properties */
#define CONSHDLR_NAME          "symresack"
#define CONSHDLR_DESC          "symmetry breaking constraint handler relying on symresacks"
#define CONSHDLR_SEPAPRIORITY    +40100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -1005200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -1005200 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             5 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             5 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP
#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_EXHAUSTIVE

#define DEFAULT_PPSYMRESACK        TRUE /**< whether we allow upgrading to packing/partitioning symresacks */
#define DEFAULT_CHECKMONOTONICITY  TRUE /**< check whether permutation is monotone when upgrading to packing/partitioning symresacks */
#define DEFAULT_FORCECONSCOPY     FALSE /**< whether symresack constraints should be forced to be copied to sub SCIPs */

/* Constants to store fixings */
#define FIXED0    1                     /* When a variable is fixed to 0. */
#define FIXED1    2                     /* When a variable is fixed to 1. */
#define UNFIXED   3                     /* When a variable is neither fixed to 0 or to 1. */
#define NOINIT    0                     /* A dummy entry for non-initialized variables.
                                         * Must have value 0 because of SCIPallocCleanBufferArray. */
/* A macro for checking if a variable was fixed during a bound change */
#define ISFIXED(x, bdchgidx)   (SCIPvarGetUbAtIndex(x, bdchgidx, FALSE) - SCIPvarGetLbAtIndex(x, bdchgidx, FALSE) < 0.5)



/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             checkppsymresack;   /**< whether we allow upgrading to packing/partitioning symresacks */
   SCIP_Bool             checkmonotonicity;  /**< check whether permutation is monotone when upgrading to packing/partitioning symresacks */
   int                   maxnvars;           /**< maximal number of variables in a symresack constraint */
   SCIP_Bool             forceconscopy;      /**< whether symresack constraints should be forced to be copied to sub SCIPs */
};


/** constraint data for symresack constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables */
   int                   nvars;              /**< number of variables */
   int*                  perm;               /**< permutation associated to the symresack */
   int*                  invperm;            /**< inverse permutation */
   SCIP_Bool             ppupgrade;          /**< whether constraint is upgraded to packing/partitioning symresack */
   SCIP_Bool             ismodelcons;        /**< whether the symresack is a model constraint */
#ifdef SCIP_DEBUG
   int                   debugcnt;           /**< counter to store number of added cover inequalities */
#endif

   /* data for upgraded symresack constraints */
   int                   ncycles;            /**< number of cycles in permutation */
   int**                 cycledecomposition; /**< cycle decomposition */
   int                   ndescentpoints;     /**< number of descent points in perm (only used if perm is not monotone) */
   int*                  descentpoints;      /**< descent points in perm (only used if perm is not monotone) */
};


/*
 * Local methods
 */

/** frees a symresack constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to symresack constraint data */
   )
{
   int nvars;
   int i;

   assert( consdata != NULL );
   assert( *consdata != NULL );

   nvars = (*consdata)->nvars;

   if ( nvars == 0 )
   {
      assert( (*consdata)->vars == NULL );
      assert( (*consdata)->perm == NULL );
      assert( (*consdata)->invperm == NULL );
      assert( (*consdata)->ncycles == 0 );
      assert( (*consdata)->cycledecomposition == NULL );

      SCIPfreeBlockMemory(scip, consdata);

      return SCIP_OKAY;
   }

   if ( (*consdata)->ndescentpoints > 0 )
   {
      assert( (*consdata)->descentpoints != NULL );

      SCIPfreeBlockMemoryArray(scip, &((*consdata)->descentpoints), (*consdata)->ndescentpoints);
   }

   if ( (*consdata)->ppupgrade )
   {
      for (i = 0; i < (*consdata)->ncycles; ++i)
      {
         SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cycledecomposition[i]), nvars + 1);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cycledecomposition), (*consdata)->ncycles);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->invperm), nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->perm), nvars);

   for (i = 0; i < nvars; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars), nvars);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** check whether constraint can be upgraded to packing/partitioning symresack */
static
SCIP_RETCODE packingUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables affected by permutation */
   int                   nvars,              /**< length of permutation */
   SCIP_Bool             checkmonotonicity,  /**< check whether permutation is monotone */
   SCIP_Bool*            upgrade             /**< pointer to store whether upgrade was successful */
   )
{
   SCIP_Bool* covered;
   SCIP_Bool descent;
   SCIP_CONSHDLR* setppcconshdlr;
   SCIP_CONS** setppcconss;
   SCIP_VAR* var;
   SCIP_Bool terminated = FALSE;
   int** cycledecomposition;
   int* indicesincycle;
   int nsetppcconss;
   int curcycle;
   int maxcyclelength;
   int ncycles = 0;
   int c;
   int i;
   int j;
   int ndescentpoints = 0;
   int* descentpoints;

   assert( scip != NULL );
   assert( perm != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );
   assert( upgrade != NULL );

   *upgrade = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &covered, nvars) );

   for (i = 0; i < nvars; ++i)
      covered[i] = FALSE;

   /* get number of cycles in permutation */
   for (i = 0; i < nvars; ++i)
   {
      /* skip checked indices */
      if ( covered[i] )
         continue;

      ++ncycles;
      j = i;
      descent = FALSE;

      do
      {
         covered[j] = TRUE;

         if ( perm[j] < j )
         {
            ++ndescentpoints;

            if ( ! descent )
               descent = TRUE;
            else if ( checkmonotonicity )
               break;
         }

         j = perm[j];
      }
      while ( j != i );

      /* if cycle is not monotone and we require the cycle to be monotone */
      if ( j != i )
      {
         assert( checkmonotonicity );
         SCIPfreeBufferArray(scip, &covered);

         return SCIP_OKAY;
      }
   }
   assert( ncycles <= nvars / 2 );

   /* check for packing/partitioning type */
   for (i = 0; i < nvars; ++i)
      covered[i] = FALSE;

   /* compute cycle decomposition: row i stores in entry 0 the length of the cycle,
    * the remaining entries are the coordinates in the cycle;
    * store descent points as well if permutation is not monotone */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cycledecomposition, ncycles) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &descentpoints, ndescentpoints) );
   for (i = 0; i < ncycles; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cycledecomposition[i], nvars + 1) );
   }

   curcycle = 0;
   maxcyclelength = 0;
   c = 0;
   for (i = 0; i < nvars; ++i)
   {
      int cyclelength = 0;

      /* skip checked indices */
      if ( covered[i] )
         continue;

      j = i;
      do
      {
         if ( perm[j] < j )
            descentpoints[c++] = j;

         covered[j] = TRUE;
         cycledecomposition[curcycle][++cyclelength] = j;
         j = perm[j];
      }
      while ( j != i );

      cycledecomposition[curcycle][0] = cyclelength;
      ++curcycle;

      if ( maxcyclelength < cyclelength )
         maxcyclelength = cyclelength;
   }
   assert( c == ndescentpoints );

   /* permutation can be upgraded -> check whether the symresack is of packing/partitioning type */
   setppcconshdlr = SCIPfindConshdlr(scip, "setppc");
   if ( setppcconshdlr == NULL )
   {
      SCIPerrorMessage("Setppc constraint handler not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   setppcconss = SCIPconshdlrGetConss(setppcconshdlr);
   nsetppcconss = SCIPconshdlrGetNConss(setppcconshdlr);

   /* Check whether each cycle of the symresack is contained in a set packing/partitioning constraint.
    * To this end, we have to guarantee that all affected variables are not negated since permutations
    * are given w.r.t. original variables. */
   *upgrade = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &indicesincycle, maxcyclelength) );

   for (i = 0; i < ncycles && *upgrade && ! terminated; ++i)
   {
      int cyclelength;

      /* get indices of variables in current cycle */
      for (j = 0; j < cycledecomposition[i][0]; ++ j)
      {
         var = vars[cycledecomposition[i][j + 1]];

         if ( SCIPvarIsNegated(var) )
         {
            terminated = TRUE;
            break;
         }

         indicesincycle[j] = SCIPvarGetProbindex(var);
      }

      cyclelength = cycledecomposition[i][0];

      /* iterate over constraints
       *
       * @todo Improve the check by sorting the constraints in the setppcconss array
       * by type and number of contained variables. */
      for (c = 0; c < nsetppcconss; ++c)
      {
         int nsetppcvars;
         SCIP_VAR** setppcvars;
         int varidx;
         int nfound = 0;
         int k;

         /* check type */
         if ( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_COVERING )
            continue;
         assert( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_PARTITIONING || SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_PACKING );

         /* get set packing/partitioning variables */
         nsetppcvars = SCIPgetNVarsSetppc(scip, setppcconss[c]);

         /* skip empty constraints (might not have been removed by presolving yet) */
         if ( nsetppcvars == 0 )
            continue;
         assert( nsetppcvars > 0 );

         setppcvars = SCIPgetVarsSetppc(scip, setppcconss[c]);
         assert( setppcvars != NULL );

         /* check whether all variables of the cycle are contained in setppc constraint */
         for (j = 0; j < nsetppcvars && nfound < cyclelength; ++j)
         {
            var = setppcvars[j];

            if ( SCIPvarIsNegated(var) )
               continue;

            varidx = SCIPvarGetProbindex(var);

            for (k = 0; k < cyclelength; ++k)
            {
               if ( varidx == indicesincycle[k] )
               {
                  ++nfound;
                  break;
               }
            }
         }
         assert( nfound <= cyclelength );

         if ( nfound == cyclelength )
            break;
      }

      /* row is not contained in a set packing/partitioning constraint */
      if ( c >= nsetppcconss )
         *upgrade = FALSE;
   }

   if ( *upgrade )
   {
      (*consdata)->ncycles = ncycles;
      (*consdata)->cycledecomposition = cycledecomposition;
      (*consdata)->ndescentpoints = ndescentpoints;
      (*consdata)->descentpoints = descentpoints;
      SCIPdebugMsg(scip, "added monotone PP symresack.\n");

      SCIPfreeBufferArray(scip, &indicesincycle);
      SCIPfreeBufferArray(scip, &covered);
   }
   else
   {
      SCIPfreeBlockMemoryArray(scip, &descentpoints, ndescentpoints);
      SCIPfreeBufferArray(scip, &indicesincycle);
      SCIPfreeBufferArray(scip, &covered);
      for (i = 0; i < ncycles; ++i)
      {
         SCIPfreeBlockMemoryArray(scip, &cycledecomposition[i], nvars + 1);
      }
      SCIPfreeBlockMemoryArray(scip, &cycledecomposition, ncycles);
   }

   return SCIP_OKAY;
}


/** creates symresack constraint data
 *
 *  If the input data contains non-binary variables or fixed
 *  points, we delete these variables in a preprocessing step.
 */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< symresack constraint handler */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   SCIP_VAR*const*       inputvars,          /**< input variables of the constraint handler */
   int                   inputnvars,         /**< input number of variables of the constraint handler*/
   int*                  inputperm,          /**< input permutation of the constraint handler */
   SCIP_Bool             ismodelcons         /**< whether the symresack is a model constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR** vars;
   SCIP_Bool upgrade;
   int* indexcorrection;
   int* invperm;
   int* perm;
   int naffectedvariables;
   int i;
   int j = 0;

   assert( consdata != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

#ifdef SCIP_DEBUG
   (*consdata)->debugcnt = 0;
#endif

   (*consdata)->ndescentpoints = 0;
   (*consdata)->descentpoints = NULL;
   (*consdata)->ismodelcons = ismodelcons;

   /* count the number of binary variables which are affected by the permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &indexcorrection, inputnvars) );
   indexcorrection[0] = -1;
   for (i = 0; i < inputnvars; ++i)
   {
      if ( inputperm[i] != i && SCIPvarIsBinary(inputvars[i]) )
      {
         if ( i == 0 )
            indexcorrection[i] = 0;
         else
            indexcorrection[i] = indexcorrection[i - 1] + 1;
      }
      else
      {
         if ( i > 0 )
            indexcorrection[i] = indexcorrection[i - 1];
      }
   }
   naffectedvariables = indexcorrection[inputnvars - 1] + 1;

   (*consdata)->nvars = naffectedvariables;

   /* Stop if we detect that the permutation fixes each binary point. */
   if ( naffectedvariables == 0 )
   {
      SCIPfreeBufferArrayNull(scip, &indexcorrection);

      (*consdata)->vars = NULL;
      (*consdata)->perm = NULL;
      (*consdata)->invperm = NULL;
      (*consdata)->ppupgrade = FALSE;
      (*consdata)->ncycles = 0;
      (*consdata)->cycledecomposition = NULL;
      return SCIP_OKAY;
   }

   /* remove fixed points from permutation representation */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, naffectedvariables) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &perm, naffectedvariables) );
   for (i = 0; i < inputnvars; ++i)
   {
      if ( i == 0 )
      {
         if ( indexcorrection[i] > -1 )
         {
            vars[j] = inputvars[i];
            perm[j++] = indexcorrection[inputperm[i]];
         }
      }
      else
      {
         if ( indexcorrection[i] > indexcorrection[i - 1] )
         {
            vars[j] = inputvars[i];
            perm[j++] = indexcorrection[inputperm[i]];
         }
      }
   }
   SCIPfreeBufferArrayNull(scip, &indexcorrection);

   (*consdata)->vars = vars;
   (*consdata)->perm = perm;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &invperm, naffectedvariables) );
   for (i = 0; i < naffectedvariables; ++i)
   {
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[i]) );
      invperm[perm[i]] = i;
   }
   (*consdata)->invperm = invperm;

   /* check for upgrade to packing/partitioning symresacks*/
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   upgrade = FALSE;
   if ( conshdlrdata->checkppsymresack )
   {
      SCIP_CALL( packingUpgrade(scip, consdata, perm, vars, naffectedvariables, conshdlrdata->checkmonotonicity, &upgrade) );
   }

   (*consdata)->ppupgrade = upgrade;

   /* get transformed variables, if we are in the transformed problem */
   if ( SCIPisTransformed(scip) )
   {
      /* Make sure that all variables cannot be multiaggregated (cannot be handled by cons_symresack, since one cannot
       * easily eliminate single variables from a symresack constraint. */
      for (i = 0; i < naffectedvariables; ++i)
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->vars[i], &(*consdata)->vars[i]) );
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*consdata)->vars[i]) );
      }
   }

   return SCIP_OKAY;
}


/** generate initial LP cut
 *
 *  We generate the ordering inequality for the pair \f$(1, \gamma^{-1}(1))\f$, i.e.,
 *  the inequality \f$-x_{1} + x_{\gamma^{-1}(1)} \leq 0\f$. This inequality is valid,
 *  because we guaranteed in a preprocessing step that all variables are binary.
 *
 *  Furthermore, we add facet inequalities of packing/partitioning symresacks if
 *  we deal with packing/partitioning symresacks.
 */
static
SCIP_RETCODE initLP(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             checkmonotonicity,  /**< has it been checked whether permutation is monotone for packing/partitioning symresacks? */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_ROW* row;
   int nvars;
#ifdef SCIP_DEBUG
   char name[SCIP_MAXSTRLEN];
#endif
   int i;
   int j;
   int k;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   nvars = consdata->nvars;

   /* avoid stupid problems */
   if ( nvars <= 1 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   vars = consdata->vars;

   /* there are no fixed points */
   assert( consdata->invperm[0] != 0 );

   /* add ordering inequality */
#ifdef SCIP_DEBUG
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_init_%s", SCIPconsGetName(cons));
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#else
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#endif
   SCIP_CALL( SCIPaddVarToRow(scip, row, vars[0], -1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, vars[consdata->invperm[0]], 1.0) );

   SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   /* check whether we have a packing/partioning symresack */
   if ( consdata->ppupgrade && ! *infeasible )
   {
      if ( checkmonotonicity )
      {
         SCIP_VAR** varsincons;
         SCIP_Real* coeffs;
         int** cycledecomposition;
         int ncycles;
         int nvarsincons;
         int nvarsincycle;
         int firstelemincycle;

         ncycles = consdata->ncycles;
         cycledecomposition = consdata->cycledecomposition;

         SCIP_CALL( SCIPallocBufferArray(scip, &varsincons, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &coeffs, nvars) );

         coeffs[0] = 1.0;

         /* add packing/partitioning symresack constraints */
         for (i = 0; i < ncycles; ++i)
         {
            assert( cycledecomposition[i][0] > 0 );

            nvarsincycle = cycledecomposition[i][0];
            varsincons[0] = vars[cycledecomposition[i][nvarsincycle]];
            firstelemincycle = cycledecomposition[i][1];

            assert( firstelemincycle == consdata->perm[cycledecomposition[i][nvarsincycle]] );

            nvarsincons = 1;

            /* add variables of other cycles to the constraint */
            for (j = 0; j < i; ++j)
            {
               nvarsincycle = cycledecomposition[j][0];
               for (k = 1; k <= nvarsincycle; ++k)
               {
                  if ( cycledecomposition[j][k] < firstelemincycle )
                  {
                     varsincons[nvarsincons] = vars[cycledecomposition[j][k]];
                     coeffs[nvarsincons++] = -1.0;
                  }
                  else
                     continue;
               }
            }

#ifdef SCIP_DEBUG
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ppSymresack_%d_%s", i, SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#else
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#endif
            SCIP_CALL( SCIPaddVarsToRow(scip, row, nvarsincons, varsincons, coeffs) );

            SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );

            if ( *infeasible )
               break;
         }

         SCIPfreeBufferArray(scip, &coeffs);
         SCIPfreeBufferArray(scip, &varsincons);
      }
      else
      {
         SCIP_Real* coeffs;
         SCIP_VAR** varsincons;
         int* imgdescentpoints;
         int* descentpoints;
         int* perm;
         int ndescentpoints;
         int lastascent = 0;
         int newlastascent = 0;
         int nvarsincons = 1;

         descentpoints = consdata->descentpoints;
         ndescentpoints = consdata->ndescentpoints;
         perm = consdata->perm;

         assert( descentpoints != NULL );
         assert( ndescentpoints > 0 );
         assert( perm != NULL );
         assert( vars != NULL );
         assert( nvars > 0 );

         SCIP_CALL( SCIPallocBufferArray(scip, &imgdescentpoints, ndescentpoints) );

         /* get images of descentpoints */
         for (j = 0; j < ndescentpoints; ++j)
            imgdescentpoints[j] = perm[descentpoints[j]];

         /* sort descent points increasingly w.r.t. the corresponding image */
         SCIPsortIntInt(imgdescentpoints, descentpoints, ndescentpoints);

         /* iteratively generate coefficient vector: the first entry is the descent point j and the remaining entries
          * are the corresponding ascent points less than perm[j]
          */
         SCIP_CALL( SCIPallocClearBufferArray(scip, &coeffs, nvars) );
         SCIP_CALL( SCIPallocClearBufferArray(scip, &varsincons, nvars) );
         coeffs[0] = 1.0;
         for (j = 0; j < ndescentpoints; ++j)
         {
            varsincons[0] = vars[descentpoints[j]];
            for (i = lastascent; i < imgdescentpoints[j]; ++i)
            {
               if ( perm[i] > i )
               {
                  coeffs[nvarsincons] = -1.0;
                  varsincons[nvarsincons++] = vars[i];
                  newlastascent = i;
               }
            }
            lastascent = newlastascent;

#ifdef SCIP_DEBUG
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ppSymresack_%d_%s", j, SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#else
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#endif
            SCIP_CALL( SCIPaddVarsToRow(scip, row, nvarsincons, varsincons, coeffs) );

            SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );

            if ( *infeasible )
               break;
         }

         SCIPfreeBufferArray(scip, &varsincons);
         SCIPfreeBufferArray(scip, &coeffs);
         SCIPfreeBufferArray(scip, &imgdescentpoints);
      }
   }

   return SCIP_OKAY;
}


/** Determines if a vector with additional fixings could exist that is lexicographically larger than its image.
 *
 * Given a vector of variables, a permutation, and a set of additional (virtual) fixings.
 * If a vector adhering to the local variable bounds (local fixings) and to the virtual fixings exists,
 * then infeasible is FALSE, otherwise TRUE.
 */
static
SCIP_RETCODE checkFeasible(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR**            vars,               /**< array of variables affected by permutation */
   int*                  invperm,            /**< inverse of permutation */
   int                   nvars,              /**< number of variables */
   int                   start,              /**< at which position to start (assuming previous positions are equal) */
   int*                  tempfixings,        /**< array with at entry i the virtual fixing of variable vars[i] */
   int*                  tempfixentries,     /**< the entries i that are virtually fixed until numfixentriesinit */
   int                   numfixentriesinit,  /**< the number of virtually fixed entries */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is detected in these fixings */
   int*                  infeasibleentry     /**< pointer to store at which entry a (0, 1) pattern is found */
)
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int var1fix;
   int var2fix;

   int i;
   int numfixentries;

   /* avoid trivial problems */
   if ( nvars < 2 )
      return SCIP_OKAY;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( invperm != NULL );
   assert( tempfixings != NULL );
   assert( tempfixentries != NULL );
   assert( infeasible != NULL );

   /* A counter for how many virtual fixings we have. */
   numfixentries = numfixentriesinit;

   *infeasible = FALSE;

   for (i = start; i < nvars; ++i)
   {
      /* there are no fixed points */
      assert( invperm[i] != i );

      /* get variables of first and second column */
      var1 = vars[i];
      var2 = vars[invperm[i]];

      assert( var1 != NULL );
      assert( var2 != NULL );

      /* Get virtual fixing of variable in left column */
      var1fix = tempfixings[i];
      if ( var1fix == NOINIT )
      {
         if ( SCIPvarGetUbLocal(var1) < 0.5 )
         {
            var1fix = FIXED0;
            assert( SCIPvarGetLbLocal(var1) <= 0.5 );
         }
         else if ( SCIPvarGetLbLocal(var1) > 0.5 )
            var1fix = FIXED1;
         else
            var1fix = UNFIXED;
      }
      assert( var1fix != NOINIT );

      /* Get virtual fixing of variable in right column */
      var2fix = tempfixings[invperm[i]];
      if ( var2fix == NOINIT )
      {
         if ( SCIPvarGetUbLocal(var2) < 0.5 )
         {
            var2fix = FIXED0;
            assert( SCIPvarGetLbLocal(var2) <= 0.5 );
         }
         else if ( SCIPvarGetLbLocal(var2) > 0.5 )
            var2fix = FIXED1;
         else
            var2fix = UNFIXED;
      }
      assert( var2fix != NOINIT );

      /* Encounter one of (_, _), (_, 0), (1, _), (1, 0). In all cases (1, 0) can be constructed. Thus feasible. */
      if ( var1fix != FIXED0 && var2fix != FIXED1 )
         break;
      /* Encounter (0, 1). Infeasible. */
      else if ( var1fix == FIXED0 && var2fix == FIXED1 )
      {
         *infeasible = TRUE;
         *infeasibleentry = i;
         break;
      }
      /* Encounter (0, _). Virtually fix var2 to 0. */
      else if ( var1fix == FIXED0 && var2fix == UNFIXED )
      {
         tempfixings[invperm[i]] = FIXED0;
         /* Mark that we have fixed invperm[i]. */
         tempfixentries[numfixentries++] = invperm[i];
      }
      /* Encounter (_, 1). Virtually fix var1 to 1. */
      else if(var1fix == UNFIXED && var2fix == FIXED1 )
      {
         tempfixings[i] = FIXED0;
         /* Mark that we have fixed invperm[i]. */
         tempfixentries[numfixentries++] = i;
      }
      /* Remaining cases are (0, 0) and (1, 1). In both cases: continue. */
   }

   /* Undo virtual fixings made in this function */
   for (i = numfixentriesinit; i < numfixentries; ++i)
   {
      tempfixings[tempfixentries[i]] = NOINIT;
      tempfixentries[i] = 0;
   }

   return SCIP_OKAY;
}


/** perform propagation of symresack constraint */
static
SCIP_RETCODE propVariables(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to be propagated */
   SCIP_Bool*            infeasible,         /**< pointer to store whether it was detected that the node is infeasible */
   int*                  ngen                /**< pointer to store number of generated bound strengthenings */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* invperm;
   int nvars;
   int i;
   int r;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int var1fix;
   int var2fix;
   SCIP_Bool tightened;
   SCIP_Bool peekinfeasible;
   int peekinfeasibleentry;
   int* tempfixings;
   int* tempfixentries;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   SCIPdebugMsg(scip, "Propagating variables of constraint <%s>.\n", SCIPconsGetName(cons));

   *ngen = 0;
   *infeasible = FALSE;

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   nvars = consdata->nvars;

   /* avoid trivial problems */
   if ( nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->invperm != NULL );
   vars = consdata->vars;
   invperm = consdata->invperm;

   /* loop through all variables */
   for (i = 0; i < nvars; ++i)
   {
      /* there are no fixed points */
      assert( invperm[i] != i );

      /* get variables of first and second column */
      var1 = vars[i];
      var2 = vars[invperm[i]];
      assert( var1 != NULL );
      assert( var2 != NULL );

      /* Get the fixing status of the left column variable var1 */
      if ( SCIPvarGetUbLocal(var1) < 0.5 )
      {
         var1fix = FIXED0;
         assert( SCIPvarGetLbLocal(var1) <= 0.5 );
      }
      else if ( SCIPvarGetLbLocal(var1) > 0.5 )
         var1fix = FIXED1;
      else
         var1fix = UNFIXED;

      /* Get the fixing status of the right column variable var2 */
      if ( SCIPvarGetUbLocal(var2) < 0.5 )
      {
         var2fix = FIXED0;
         assert( SCIPvarGetLbLocal(var2) <= 0.5 );
      }
      else if ( SCIPvarGetLbLocal(var2) > 0.5 )
         var2fix = FIXED1;
      else
         var2fix = UNFIXED;

      /* Encounter one of (_, _), (_, 0), (1, _), (1, 0). Check if (1, 1) or (0, 0) are possible, otherwise fix. */
      if ( var1fix != FIXED0 && var2fix != FIXED1 )
      {
         assert( SCIPvarGetUbLocal(var1) > 0.5 );
         assert( SCIPvarGetLbLocal(var2) < 0.5 );

         SCIPdebugMsg(scip, "Check variable pair (%d,%d).\n", i, invperm[i]);
         SCIPdebugMsg(scip, " -> node is feasible (could set pair to (1,0) and every earlier pair is constant).\n");

         if ( var1fix == UNFIXED || var2fix == UNFIXED )
         {
            /* Create arrays tempfixings and tempfixentries to store virtual fixings. */
            SCIP_CALL( SCIPallocCleanBufferArray(scip, &tempfixings, nvars) );
            SCIP_CALL( SCIPallocCleanBufferArray(scip, &tempfixentries, nvars) );

            if ( var1fix == UNFIXED )
            {
               assert( SCIPvarGetLbLocal(var1) < 0.5 );

               /* Peek whether a lexicographical larger-or-equal vector can be created with var1 fixed to 0 */
               SCIPdebugMsg(scip, " -> First entry is not fixed. Check if 0 is feasible.\n");
               tempfixings[i] = FIXED0;
               tempfixentries[0] = i;
               SCIP_CALL( checkFeasible(scip, vars, invperm, nvars, i, tempfixings, tempfixentries, 1,
                     &peekinfeasible, &peekinfeasibleentry) );

               if ( peekinfeasible )
               {
                  /* No feasible vector exists with var1 set to 0, so it must be a 1-fixing. */
                  SCIPdebugMsg(scip, " -> First entry is not fixed. 0 is not feasible. Fixing to 1.\n");
                  SCIP_CALL( SCIPinferVarLbCons(scip, var1, 1.0, cons, i + nvars * peekinfeasibleentry,
                     FALSE, infeasible, &tightened) ); /*lint !e713*/
                  assert( ! *infeasible );

                  if ( tightened )
                     ++(*ngen);
               }

               tempfixings[i] = NOINIT;
               tempfixentries[0] = 0;
            }

            if ( var2fix == UNFIXED )
            {
               assert( SCIPvarGetUbLocal(var2) > 0.5 );

               /* Peek whether a lexicographical larger-or-equal vector can be created with var2 fixed to 1 */
               SCIPdebugMsg(scip, " -> Second entry is not fixed. Check if 1 is feasible.\n");
               tempfixings[invperm[i]] = FIXED1;
               tempfixentries[0] = invperm[i];
               SCIP_CALL( checkFeasible(scip, vars, invperm, nvars, i, tempfixings, tempfixentries, 1,
                     &peekinfeasible, &peekinfeasibleentry) );

               if ( peekinfeasible )
               {
                  /* No feasible vector exists with var2 set to 1, so it must be a 1-fixing. */
                  SCIPdebugMsg(scip, " -> Second entry is not fixed. 1 is not feasible. Fixing to 0.\n");
                  SCIP_CALL( SCIPinferVarUbCons(scip, var2, 0.0, cons, i + nvars * peekinfeasibleentry,
                     FALSE, infeasible, &tightened) ); /*lint !e713*/
                  assert( ! *infeasible );

                  if ( tightened )
                     ++(*ngen);
               }

               tempfixings[invperm[i]] = NOINIT;
               tempfixentries[0] = 0;
            }

            SCIPfreeCleanBufferArray(scip, &tempfixentries);
            SCIPfreeCleanBufferArray(scip, &tempfixings);
         }

         /* Can stop here, because this row can become (1, 0). Therefore all next rows can take arbitrary values. */
         break;
      }
      /* Encounter (0, 1): If first part of variable pair fixed to 0 and second part is fixed to 1 */
      else if ( var1fix == FIXED0 && var2fix == FIXED1 )
      {
         SCIPdebugMsg(scip, "Check variable pair (%d,%d).\n", i, invperm[i]);

         SCIPdebugMsg(scip, " -> node infeasible (pair was fixed to (0,1) but there was no pair of type (1,0) before) ---> lexicographical order violated, infeasible.\n");

         /* perform conflict analysis */
         if ( SCIPisConflictAnalysisApplicable(scip) )
         {
            SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

            for (r = 0; r <= i; ++r)
            {
               /* there are no fixed points */
               assert( invperm[r] != r );

               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[r]) );
               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[invperm[r]]) );
            }

            SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
         }

         *infeasible = TRUE;
         break;
      }
      /* Encounter (0, _): Fix second part to 0 */
      else if ( var1fix == FIXED0 && var2fix != FIXED0 )
      {
         assert( SCIPvarGetUbLocal(var1) < 0.5 );
         assert( SCIPvarGetLbLocal(var2) < 0.5 );
         assert( SCIPvarGetUbLocal(var2) > 0.5 );

         SCIPdebugMsg(scip, "Check variable pair (%d,%d).\n", i, invperm[i]);

         assert( SCIPvarGetLbLocal(var2) < 0.5 );
         SCIP_CALL( SCIPinferVarUbCons(scip, var2, 0.0, cons, i, FALSE, infeasible, &tightened) ); /*lint !e713*/
         assert( ! *infeasible );

         if ( tightened )
            ++(*ngen);
      }
      /* Encounter (_, 1): fix first part to 1 */
      else if ( var1fix != FIXED1 && var2fix == FIXED1 )
      {
         assert( SCIPvarGetLbLocal(var1) < 0.5 );
         assert( SCIPvarGetUbLocal(var1) > 0.5 );
         assert( SCIPvarGetLbLocal(var2) > 0.5 );

         SCIPdebugMsg(scip, "Check variable pair (%d,%d).\n", i, invperm[i]);

         assert( SCIPvarGetUbLocal(var1) > 0.5 );
         SCIP_CALL( SCIPinferVarLbCons(scip, var1, 1.0, cons, i, FALSE, infeasible, &tightened) ); /*lint !e713*/
         assert( ! *infeasible );

         if ( tightened )
            ++(*ngen);
      }
      /* Remaining cases are (0, 0) and (1, 1). In these cases we can continue! */
   }

   return SCIP_OKAY;
}


/** add symresack cover inequality */
static
SCIP_RETCODE addSymresackInequality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars,               /**< variables */
   int*                  coeffs,             /**< coefficient vector of inequality to be added */
   SCIP_Real             rhs,                /**< right-hand side of inequality to be added */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_ROW* row;
   int i;
#ifdef SCIP_DEBUG
   SCIP_CONSDATA* consdata;
   char name[SCIP_MAXSTRLEN];
#endif

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nvars > 0 );
   assert( vars != NULL );
   assert( coeffs != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

#ifdef SCIP_DEBUG
   consdata = SCIPconsGetData(cons);
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_cover_%s_%d", SCIPconsGetName(cons), consdata->debugcnt);
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
   ++consdata->debugcnt;
#else
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
#endif
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for (i = 0; i < nvars; ++i)
   {
      if ( coeffs[i] == 1 || coeffs[i] == -1 )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i], (SCIP_Real) coeffs[i]) );
      }
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );
   SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}


/** Maximize a linear function on a "strict" symresack,
 *  that is a symresack where we do not allow the solution x = gamma(x).
 */
static
SCIP_RETCODE maximizeObjectiveSymresackStrict(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   nvars,              /**< number of variables in symresack */
   SCIP_Real*            objective,          /**< the objective vector */
   int*                  perm,               /**< the permutation (without fixed points) as an array */
   int*                  invperm,            /**< the inverse permutation as an array */
   int*                  maxcrit,            /**< pointer to the critical entry where optimality is found at */
   SCIP_Real*            maxsoluval          /**< pointer to store the optimal objective value */
   )
{
   /* The maximal objective in every iteration. */
   SCIP_Real tmpobj;
   /* The new value of componentobj when combining two components. */
   SCIP_Real tmpnewcompobj;
   /* helperobj is the sum of all positive objective-sums for all components. */
   SCIP_Real helperobj = 0.0;

   int crit;
   int critinv;
   int i;

   /* For every vertex of degree < 2 we maintain componentends and componentobj. */
   int* componentends;
   SCIP_Real* componentobj;

   assert( scip != NULL );
   assert( nvars > 0 );
   assert( objective != NULL );
   assert( perm != NULL );
   assert( invperm != NULL );
   assert( maxcrit != NULL );
   assert( maxsoluval != NULL );

   /* The current best known critical entry and objective */
   *maxcrit = -1;
   *maxsoluval = -SCIP_DEFAULT_INFINITY;

   SCIP_CALL( SCIPallocBufferArray(scip, &componentends, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &componentobj, nvars) );

   /* Initialization: Every entry is a component in the graph,
    * having the corresponding objective
    */
   for (i = 0; i < nvars; ++i)
   {
      componentends[i] = i;
      componentobj[i] = objective[i];
      if ( SCIPisGT(scip, objective[i], 0.0) )
         helperobj += objective[i];
   }

   /* Iterate over all critical rows, and of the graph maintain the components on the vertices of degree < 2. */
   for (crit = 0; crit < nvars; ++crit)
   {
      critinv = invperm[crit];

      /* Do not allow fixed points. */
      assert( crit != critinv );

      /* If the other end of the component of crit is critinv, then crit cannot be a critical entry. */
      if ( componentends[crit] == critinv )
         continue;

      /* Compute objective for crit as critical entry. Update if it is better than the best found objective */
      tmpobj = helperobj;
      if ( SCIPisLT(scip, componentobj[crit], 0.0) )
         tmpobj += componentobj[crit];
      if ( SCIPisGT(scip, componentobj[critinv], 0.0) )
         tmpobj -= componentobj[critinv];
      if ( SCIPisGT(scip, tmpobj, *maxsoluval) )
      {
         *maxsoluval = tmpobj;
         *maxcrit = crit;
      }

      /* Update helperobj */
      tmpnewcompobj = componentobj[crit] + componentobj[critinv];
      if ( SCIPisGT(scip, componentobj[crit], 0.0) )
         helperobj -= componentobj[crit];
      if ( SCIPisGT(scip, componentobj[critinv], 0.0) )
         helperobj -= componentobj[critinv];
      if ( SCIPisGT(scip, tmpnewcompobj, 0.0) )
         helperobj += tmpnewcompobj;

      /* Update the objective of a component */
      componentobj[componentends[crit]] = tmpnewcompobj;
      componentobj[componentends[critinv]] = tmpnewcompobj;

      /* Connect the endpoints of the newly created path */
      if ( componentends[crit] == crit )
      {
         componentends[crit] = componentends[critinv];
         componentends[componentends[critinv]] = crit;
      }
      else
      {
         componentends[componentends[crit]] = componentends[critinv];
         componentends[componentends[critinv]] = componentends[crit];
      }

      /* Early termination criterion. helperobj is upper bound to tmpobj for every next iteration,
         * so if helperobj <= maxsoluval then we can terminate earlier.
         */
      if ( SCIPisGE(scip, *maxsoluval, helperobj) )
         break;
   }

   /* It is always possible to make the first entry critical. */
   assert( *maxcrit >= 0 );

   SCIPfreeBufferArray(scip, &componentobj);
   SCIPfreeBufferArray(scip, &componentends);

   return SCIP_OKAY;
}


/** For a symresack, determine a maximizer for optimizing linear function
 *  over a symresack, where the critical entry is fixed.
 */
static
SCIP_RETCODE maximizeObjectiveSymresackCriticalEntry(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   nvars,              /**< number of variables in symresack */
   SCIP_Real*            objective,          /**< the objective vector */
   int*                  perm,               /**< the permutation (without fixed points) as an array */
   int*                  invperm,            /**< the inverse permutation as an array */
   int                   crit,               /**< critical entry where optimality is found at */
   int*                  maxsolu             /**< pointer to the optimal objective array */
   )
{
   /* Compute to which components all entries belong. */
   int* entrycomponent;
   SCIP_Real* componentobjective;

   int i;
   int c;

   assert( scip != NULL );
   assert( nvars > 0 );
   assert( objective != NULL );
   assert( perm != NULL );
   assert( invperm != NULL );
   assert( maxsolu != NULL );
   assert( crit >= 0 );
   assert( crit <= nvars );

   SCIP_CALL( SCIPallocBufferArray(scip, &entrycomponent, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &componentobjective, nvars) );

   /* Initially: Everything forms its own component */
   for (i = 0; i < nvars; ++i)
   {
      entrycomponent[i] = i;
      componentobjective[i] = objective[i];
   }
   for (i = 0; i < crit; ++i)
   {
      /* The graph with arcs {i, invperm[i]} if i < c is a collection of paths, cycles and singletons.
       * Label the vertices to the lowest entry in the component,  and store the value of that in this component.
       * Every inner while-loop labels one new vertex per iteration, and a vertex is relabeled exactly once.
       */
      if ( entrycomponent[i] < i )
      {
         /* This entry is already included in a component. */
         continue;
      }

      /* Follow the path forward: Take edges {c, invperm[c]} until c >= crit, or a cycle is found. */
      c = i;
      while( c < crit )
      {
         /* c < crit, so edge {c, invperm[c]} exists. Label invperm[c] as part of component of i */
         c = invperm[c];

         /* Stop if we find a cycle. */
         if ( entrycomponent[c] != c )
            break;

         entrycomponent[c] = i;
         componentobjective[i] += objective[c];
      }

      /* Follow the path backward: Take edges {c, perm[c]} until perm[c] >= crit, or a cycle is found. */
      c = perm[i];
      while( c < crit )
      {
         /* c < crit, so edge {c, invperm[c]} exists. Label c as part of component of i */

         /* Stop if we find a cycle. */
         if ( entrycomponent[c] != c )
            break;

         entrycomponent[c] = i;
         componentobjective[i] += objective[c];
         /* For next iteration: We do another step back */
         c = perm[c];
      }
   }

   /* Now fill the objective vector.
    * For the component containing crit, set the value to 1.
    * For the component contraining invperm[crit], set the value to 0.
    * For the other components, set the value to 1 if the objective sum is positive.
    * Otherwise to 0.
    */
   for (i = 0; i < nvars; ++i)
   {
      if ( entrycomponent[i] == entrycomponent[crit] )
         maxsolu[i] = 1;
      else if ( entrycomponent[i] == entrycomponent[invperm[crit]] )
         maxsolu[i] = 0;
      else if ( SCIPisGT(scip, componentobjective[entrycomponent[i]], 0.0) )
         maxsolu[i] = 1;
      else
         maxsolu[i] = 0;
   }

   SCIPfreeBufferArray(scip, &componentobjective);
   SCIPfreeBufferArray(scip, &entrycomponent);

   return SCIP_OKAY;
}


/** separate symresack cover inequalities
 *
 *  We currently do NOT enter cuts into the pool.
 */
static
SCIP_RETCODE separateSymresackCovers(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   const SCIP_CONSDATA*  consdata,           /**< constraint data */
   SCIP_Real*            vals,               /**< solution values of variables */
   int*                  ngen,               /**< pointer to store the number of separated covers */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_Real constobjective;
   SCIP_Real* sepaobjective;
   SCIP_Real maxsoluobj = 0.0;
   int* maxsolu;
   int* invperm;
   int* perm;
   int nvars;
   int maxcrit;
   int i;

   *infeasible = FALSE;
   *ngen = 0;

   assert( scip != NULL );
   assert( consdata != NULL );

   /* we do not have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->perm != NULL );
   assert( consdata->invperm != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   nvars = consdata->nvars;
   perm = consdata->perm;
   invperm = consdata->invperm;

   /* initialize objective */
   SCIP_CALL( SCIPallocBufferArray(scip, &sepaobjective, nvars) );

   constobjective = 1.0; /* constant part of separation objective */
   for (i = 0; i < nvars; ++i)
   {
      if ( i < perm[i] )
         sepaobjective[i] = - vals[i];
      else
      {
         sepaobjective[i] = 1.0 - vals[i];
         constobjective += vals[i] - 1.0;
      }
   }

   /* allocate memory for temporary and global solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &maxsolu, nvars) );

   /* Find critical row of a maximally violated cover */
   SCIP_CALL( maximizeObjectiveSymresackStrict(scip, nvars, sepaobjective, perm, invperm, &maxcrit, &maxsoluobj) );
   assert( maxcrit >= 0 );
   SCIPdebugMsg(scip, "Critical row %d found; Computing maximally violated cover.\n", maxcrit);
   SCIP_CALL( maximizeObjectiveSymresackCriticalEntry(scip, nvars, sepaobjective, perm, invperm, maxcrit, maxsolu) );

   /* Add constant to maxsoluobj to get the real objective */
   maxsoluobj += constobjective;

   /* Check whether the separation objective is positive, i.e., a violated cover was found. */
   if ( SCIPisEfficacious(scip, maxsoluobj) )
   {
      /* Now add the cut. Reuse array maxsolu as coefficient vector for the constraint. */
      SCIP_Real rhs = -1.0;
      for (i = 0; i < nvars; ++i)
      {
         if ( i < perm[i] )
            maxsolu[i] = -maxsolu[i];
         else
         {
            if ( maxsolu[i] == 0 )
               rhs += 1.0;
            maxsolu[i] = 1 - maxsolu[i];
         }
      }

      /* add cover inequality */
      SCIP_CALL( addSymresackInequality(scip, cons, nvars, consdata->vars, maxsolu, rhs, infeasible) );

      if ( ! *infeasible )
         ++(*ngen);
   }

   SCIPfreeBufferArrayNull(scip, &maxsolu);
   SCIPfreeBufferArrayNull(scip, &sepaobjective);

   return SCIP_OKAY;
}


/** check whether solution is feasible for symresacks */
static
SCIP_RETCODE checkSymresackSolution(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constrained for which we check the solution */
   SCIP_SOL*             sol,                /**< solution to be checked */
   SCIP_RESULT*          result,             /**< pointer to store whether we detected infeasibility */
   SCIP_Bool             printreason         /**< whether reason for infeasibility should be printed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* invperm;
   int nvars;
   int i;

   assert( cons != NULL );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);

   /* we do not have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->invperm != NULL );

   SCIPdebugMsg(scip, "Check method for symresack constraint <%s> (%d rows) ...\n", SCIPconsGetName(cons), consdata->nvars);

   nvars = consdata->nvars;
   vars = consdata->vars;
   invperm = consdata->invperm;

   /* detect first non-constant pair of variables */
   for (i = 0; i < nvars; ++i)
   {
      SCIP_Real solval;
      int val1;
      int val2;

      /* there are no fixed points */
      assert( invperm[i] != i );

      /* get value of variable i and its inverse */
      solval = SCIPgetSolVal(scip, sol, vars[i]);
      assert( SCIPisFeasIntegral(scip, solval) );
      if ( solval > 0.5 )
         val1 = 1;
      else
         val1 = 0;

      solval = SCIPgetSolVal(scip, sol, vars[invperm[i]]);
      assert( SCIPisFeasIntegral(scip, solval) );
      if ( solval > 0.5 )
         val2 = 1;
      else
         val2 = 0;

      /* if we detected a constant pair */
      if ( val1 == val2 )
         continue;
      /* pair is (1,0) --> lexicographically maximal */
      else if ( val1 > val2 )
         break;

      /* pair is (0,1) --> solution is infeasible */
      assert( val2 > val1 );
      SCIPdebugMsg(scip, "Solution is infeasible.\n");
      *result = SCIP_INFEASIBLE;

      if ( printreason )
         SCIPinfoMessage(scip, NULL, "First non-constant pair (%d, %d) of variables has pattern (0,1).\n", i, invperm[i]);

      break;
   }

   return SCIP_OKAY;
}


/** Upgrade symresack constraints to orbisacks */
static
SCIP_RETCODE orbisackUpgrade(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            inputvars,          /**< permuted variables array */
   int                   nvars,              /**< size of perm array */
   SCIP_Bool*            upgrade,            /**< whether constraint was upgraded */
   SCIP_Bool             ismodelcons,        /**< whether the symresack is a model constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   int nrows = 0;
   int i;

   assert( scip != NULL );
   assert( perm != NULL );
   assert( nvars > 0 );
   assert( inputvars != NULL );
   assert( upgrade != NULL );

   *upgrade = TRUE;

   /* check whether orbisack conshdlr is available */
   conshdlr = SCIPfindConshdlr(scip, "orbisack");
   if ( conshdlr == NULL )
   {
      *upgrade = FALSE;
      SCIPdebugMsg(scip, "Cannot check whether symresack constraint can be upgraded to orbisack constraint. ");
      SCIPdebugMsg(scip, "---> Orbisack constraint handler not found.\n");

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &vars1, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars2, nvars) );

   /* check whether permutation is a composition of 2-cycles */
   for (i = 0; i < nvars; ++i)
   {
      /* ignore non-binary variables */
      if ( ! SCIPvarIsBinary(inputvars[i]) )
         continue;

      if ( perm[perm[i]] != i )
      {
         *upgrade = FALSE;
         break;
      }

      if ( perm[i] > i )
      {
         vars1[nrows] = inputvars[i];
         vars2[nrows++] = inputvars[perm[i]];

         assert( nrows <= nvars );
      }
   }

   /* if permutation can be upgraded to an orbisack */
   if ( nrows == 0 )
      *upgrade = FALSE;
   else if ( *upgrade )
   {
      SCIP_CALL( SCIPcreateConsOrbisack(scip, cons, name, vars1, vars2, nrows, FALSE, FALSE, ismodelcons,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   SCIPfreeBufferArray(scip, &vars2);
   SCIPfreeBufferArray(scip, &vars1);

   return SCIP_OKAY;
}


/** creates a symmetry breaking constraint
 *
 * Depending on the given permutation, either an orbisack or symresack constraint
 * is created.
 */
SCIP_RETCODE SCIPcreateSymbreakCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables in vars array */
   SCIP_Bool             ismodelcons,        /**< whether the added constraint is a model constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_Bool upgrade = FALSE;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( perm != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );

   SCIP_CALL( orbisackUpgrade(scip, cons, name, perm, vars, nvars, &upgrade, ismodelcons,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   if ( ! upgrade )
   {
      SCIP_CALL( SCIPcreateConsSymresack(scip, cons, name, perm, vars, nvars, ismodelcons,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   return SCIP_OKAY;
}


/*--------------------------------------------------------------------------------------------
 *--------------------------------- SCIP functions -------------------------------------------
 *--------------------------------------------------------------------------------------------*/

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySymresack)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSymresack(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSymresack)
{   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( consdata != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** frees constraint handler */
static
SCIP_DECL_CONSFREE(consFreeSymresack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSymresack)
{
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* consdata = NULL;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   SCIPdebugMsg(scip, "Transforming constraint.\n");

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL);

   /* constraint might be empty and not deleted if no presolving took place */
   assert( sourcedata->nvars == 0 || sourcedata->vars != NULL );
   assert( sourcedata->nvars == 0 || sourcedata->perm != NULL );
   assert( sourcedata->nvars == 0 || sourcedata->invperm != NULL );
#ifndef NDEBUG
   if ( sourcedata->ppupgrade )
   {
      assert( sourcedata->nvars > 0 );
      assert( sourcedata->ncycles != 0 );
      assert( sourcedata->cycledecomposition != NULL );
      for (i = 0; i < sourcedata->ncycles; ++i)
      {
         assert( sourcedata->cycledecomposition[i] != NULL );
         assert( sourcedata->cycledecomposition[i][0] != 0 );
      }
   }
#endif

   /* create transformed constraint data
    *
    * do NOT call consdataCreate() again to avoid doing the packing-upgrade check twice
    */
   nvars = sourcedata->nvars;

   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->vars = NULL;
   consdata->nvars = nvars;
   consdata->perm = NULL;
   consdata->invperm = NULL;
   consdata->ppupgrade = sourcedata->ppupgrade;
   consdata->ismodelcons = sourcedata->ismodelcons;
#ifdef SCIP_DEBUG
   consdata->debugcnt = 0;
#endif
   consdata->ncycles = 0;
   consdata->cycledecomposition = NULL;
   consdata->ndescentpoints = 0;
   consdata->descentpoints = NULL;

   if ( nvars > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, nvars) );
      SCIP_CALL( SCIPgetTransformedVars(scip, nvars, sourcedata->vars, consdata->vars) );
      for (i = 0; i < nvars; ++i)
      {
         SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
      }

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->perm, sourcedata->perm, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->invperm, sourcedata->invperm, nvars) );

      if ( sourcedata->ppupgrade )
      {
         consdata->ncycles = sourcedata->ncycles;
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->cycledecomposition, sourcedata->cycledecomposition, sourcedata->ncycles) );
         for (i = 0; i < sourcedata->ncycles; ++i)
         {
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->cycledecomposition[i], sourcedata->cycledecomposition[i], nvars + 1) ); /*lint !e866*/
         }

         consdata->ndescentpoints = sourcedata->ndescentpoints;
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->descentpoints, sourcedata->descentpoints, sourcedata->ndescentpoints) );
      }

      /* Make sure that all variables cannot be multiaggregated (this cannot be handled by cons_symresack, since one cannot
       * easily eliminate single variables from a symresack constraint).
       *
       * We need to call this again to ensure that multiaggregation is forbidden also if the constraint was part
       * of the original problem.
       */
      for (i = 0; i < sourcedata->nvars; ++i)
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, consdata->vars[i], &consdata->vars[i]) );
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, consdata->vars[i]) );
      }
   }

   /* create transformed constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpSymresack)
{
   int c;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( infeasible != NULL );
   *infeasible = FALSE;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      assert( conss[c] != NULL );

      SCIPdebugMsg(scip, "Generating initial symresack cut for constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( initLP(scip, conss[c], conshdlrdata->checkmonotonicity, infeasible) );
      if ( *infeasible )
         break;
   }
   SCIPdebugMsg(scip, "Generated initial symresack cuts.\n");

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolSymresack)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* determine maximum number of vars in a symresack constraint */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   conshdlrdata->maxnvars = 0;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      assert( conss[c] != NULL );

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* update conshdlrdata if necessary */
      if ( consdata->nvars > conshdlrdata->maxnvars )
         conshdlrdata->maxnvars = consdata->nvars;
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpSymresack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real* vals;
   int maxnvars;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation method for symresack constraints\n");

   *result = SCIP_DIDNOTRUN;

   /* if solution is not integer */
   if ( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   if ( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   maxnvars = conshdlrdata->maxnvars;
   assert( maxnvars > 0 );

   SCIP_CALL( SCIPallocBufferArray(scip, &vals, maxnvars) );

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      int ngen = 0;

      SCIPdebugMsg(scip, "Separating symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);

      if ( consdata->nvars == 0 )
         continue;

      /* get solution */
      assert( consdata->nvars <= maxnvars );
      SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nvars, consdata->vars, vals) );
      SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, vals, &ngen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         SCIPfreeBufferArray(scip, &vals);

         return SCIP_OKAY;
      }

      if ( ngen > 0 )
         *result = SCIP_SEPARATED;

      if ( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;
   }
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solution */
static
SCIP_DECL_CONSSEPASOL(consSepasolSymresack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real* vals;
   int maxnvars;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation method for symresack constraints\n");

   *result = SCIP_DIDNOTRUN;

   if ( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   maxnvars = conshdlrdata->maxnvars;
   assert( maxnvars > 0 );

   SCIP_CALL( SCIPallocBufferArray(scip, &vals, maxnvars) );

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      int ngen = 0;

      SCIPdebugMsg(scip, "Separating symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);

      if ( consdata->nvars == 0 )
         continue;

      /* get solution */
      assert( consdata->nvars <= maxnvars );
      SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nvars, consdata->vars, vals) );
      SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, vals, &ngen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         SCIPfreeBufferArray(scip, &vals);

         return SCIP_OKAY;
      }

      if ( ngen > 0 )
         *result = SCIP_SEPARATED;

      if ( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;
   }
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions.
 *
 *  To check feasibility, we separate cover inequalities.
 *
 *  @pre It is assumed that the solution is integral (this can be ensured by appropriate priorities).
 */
static
SCIP_DECL_CONSENFOLP(consEnfolpSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Enforcing method for symresack constraints (lp solutions) ...\n");

   /* we have a negative priority, so we should come after the integrality conshdlr. */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   *result = SCIP_FEASIBLE;

   if ( nconss > 0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_Real* vals;
      int maxnvars;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );

      maxnvars = conshdlrdata->maxnvars;
      assert( maxnvars > 0 );

      SCIP_CALL( SCIPallocBufferArray(scip, &vals, maxnvars) );

      /* loop through constraints */
      for (c = 0; c < nconss; ++c)
      {
         SCIP_Bool infeasible = FALSE;
         int ngen = 0;

         SCIPdebugMsg(scip, "Enforcing symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

         /* get data of constraint */
         assert( conss[c] != NULL );
         consdata = SCIPconsGetData(conss[c]);
         assert( consdata != NULL );

         /* do not enforce non-model constraints */
         if ( !consdata->ismodelcons )
            continue;

         if ( consdata->nvars == 0 )
            continue;

         /* get solution */
         assert( consdata->nvars <= maxnvars );
         SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nvars, consdata->vars, vals) );
         SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, vals, &ngen, &infeasible) );

         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            SCIPfreeBufferArray(scip, &vals);

            return SCIP_OKAY;
         }

         /* SCIPdebugMsg(scip, "Generated symresack inequalities for <%s>: %d\n", SCIPconsGetName(conss[c]), ngen); */

         if ( ngen > 0 )
            *result = SCIP_SEPARATED;
      }
      SCIPfreeBufferArray(scip, &vals);
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Enforcing method for symresack constraints (pseudo solutions) ...\n");

   *result = SCIP_FEASIBLE;

   if ( objinfeasible || solinfeasible )
      return SCIP_OKAY;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* do not enforce non-model constraints */
      if ( !consdata->ismodelcons )
         continue;

      SCIP_CALL( checkSymresackSolution(scip, conss[c], NULL, result, FALSE) );

      if ( *result == SCIP_INFEASIBLE )
         break;
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions
 *
 *  To check feasibility, we separate cover inequalities.
 *
 */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSymresack)
{   /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Enforcing method for symresack constraints (relaxation solutions) ...\n");

   /* we have a negative priority, so we should come after the integrality conshdlr. */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   *result = SCIP_FEASIBLE;

   if ( nconss > 0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_Real* vals;
      int maxnvars;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );

      maxnvars = conshdlrdata->maxnvars;
      assert( maxnvars > 0 );

      SCIP_CALL( SCIPallocBufferArray(scip, &vals, maxnvars) );

      /* loop through constraints */
      for (c = 0; c < nconss; ++c)
      {
         SCIP_Bool infeasible = FALSE;
         int ngen = 0;

         SCIPdebugMsg(scip, "Enforcing symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

         /* get data of constraint */
         assert( conss[c] != NULL );
         consdata = SCIPconsGetData(conss[c]);
         assert( consdata != NULL );

         /* do not enforce non-model constraints */
         if ( !consdata->ismodelcons )
            continue;

         if ( consdata->nvars == 0 )
            continue;

          /* get solution */
         assert( consdata->nvars <= maxnvars );
         SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nvars, consdata->vars, vals) );
         SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, vals, &ngen, &infeasible) );

         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            SCIPfreeBufferArray(scip, &vals);

            return SCIP_OKAY;
         }

         if ( ngen > 0 )
            *result = SCIP_SEPARATED;
      }
      SCIPfreeBufferArray(scip, &vals);
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckSymresack)
{   /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* do not check non-model constraints */
      if ( !consdata->ismodelcons )
         continue;

      SCIP_CALL( checkSymresackSolution(scip, conss[c], sol, result, printreason) );

      if ( *result == SCIP_INFEASIBLE )
         break;
   }

   if ( *result == SCIP_FEASIBLE )
      SCIPdebugMsg(scip, "Solution is feasible.\n");

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSymresack)
{  /*lint --e{715}*/
   int c;
   SCIP_Bool success = FALSE;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "Propagation method of symresack constraint handler.\n");

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      int ngen = 0;

      assert( conss[c] != NULL );

      SCIP_CALL( propVariables(scip, conss[c], &infeasible, &ngen) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      success = success || ( ngen > 0 );

      *result = SCIP_DIDNOTFIND;
   }

   if ( success )
   {
      *result = SCIP_REDUCEDDOM;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSymresack)
{  /*lint --e{715}*/
   int c;
   SCIP_Bool success = FALSE;
   int oldndelconss;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   oldndelconss = *ndelconss;

   SCIPdebugMsg(scip, "Presolving method of symresack constraint handler. Propagating symresack inequalities.\n");
   *result = SCIP_DIDNOTRUN;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_CONSDATA* consdata;
      int ngen = 0;

      assert( conss[c] != NULL );

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* avoid trivial problems */
      if ( consdata->nvars == 0 )
      {
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         (*ndelconss)++;
      }
      else
      {
         SCIP_CALL( propVariables(scip, conss[c], &infeasible, &ngen) );
      }

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if ( ngen > 0 )
      {
         *nfixedvars += ngen;
         success = TRUE;
      }

      *result = SCIP_DIDNOTFIND;
   }

   if ( *ndelconss > oldndelconss ||  success )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** Propagation resolution for conflict analysis */
static
SCIP_DECL_CONSRESPROP(consRespropSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* perm;
   int* invperm;
   int nvars;
   int i;
   int varrow;
   int infrow;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Propagation resolution method of symresack constraint handler.\n");

   *result = SCIP_DIDNOTFIND;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* we do not have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->invperm != NULL );

   vars = consdata->vars;
   nvars = consdata->nvars;
   perm = consdata->perm;
   invperm = consdata->invperm;

   /* inferinfo == varrow + infrow * nvars.
    * infrow is 0 if the fixing is not caused by a lookahead.
    */
   varrow = inferinfo % nvars;
   infrow = inferinfo / nvars;

   assert( varrow >= 0 );
   assert( varrow < nvars );
   assert( infrow >= 0 );
   assert( infrow < nvars );
   assert( vars[varrow] == infervar || vars[invperm[varrow]] == infervar );

   /* Up to entry varrow the vectors x and perm[x] are equal. */
   for (i = 0; i < varrow; ++i)
   {
      /* Conflict caused by bounds of x[i] and perm(x)[i] = x[invperm[i]]. */

      /* No fixed points in the permutation. */
      assert( i != invperm[i] );

      /* Up to entry varrow the vectors x and perm[x] are fixed to the same value. */
      assert( ISFIXED(vars[i], bdchgidx) );
      assert( ISFIXED(vars[invperm[i]], bdchgidx) );
      assert( REALABS(SCIPvarGetUbAtIndex(vars[i], bdchgidx, FALSE) -
         SCIPvarGetUbAtIndex(vars[invperm[i]], bdchgidx, FALSE)) < 0.5 );
      assert( REALABS(SCIPvarGetLbAtIndex(vars[i], bdchgidx, FALSE) -
         SCIPvarGetLbAtIndex(vars[invperm[i]], bdchgidx, FALSE)) < 0.5 );

      /* At iteration i the vars x[i] and x[invperm[i]] are fixed.
       * So only new information is received if i < perm[i] (i.e. there is no j < i with j = invperm[i])
       * Or if invperm[i] > i.
       */
      if ( i < perm[i] )
      {
         assert( vars[i] != infervar );
         SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
      }
      if ( invperm[i] > i )
      {
         assert( vars[invperm[i]] != infervar );
         SCIP_CALL( SCIPaddConflictUb(scip, vars[invperm[i]], bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, vars[invperm[i]], bdchgidx) );
      }
   }

   /* Case distinction: Fixing due to propagation or due to lookahead */
   if ( infrow > 0 )
   {
      /* The fixing of infervar is caused by a lookahead (checkFeasible)
       * Up to row "varrow" the entries x[i] and perm(x)[i] are forced to be equal
       * If x[varrow] = perm(x)[varrow] is assumed, then until infrow we find x[i] = perm(x)[i] ( = x[invperm[i]] )
       * and (x[infrow], perm(x)[infrow]) = (0, 1).
       */

      /* Everything after varrow to infrow is forced to a constant, and row infrow is (0, 1) */
      for (i = varrow + 1; i <= infrow; ++i)
      {
         /* Conflict caused by bounds of x[i] and perm(x)[i] = x[invperm[i]]. */

         /* No fixed points in the permutation. */
         assert( i != invperm[i] );

         /* The fixing are applied 'virtually', i.e. if varrow is considered constant, then fixings will follow.
          * Thus, between entries varrow and infrow of vectorx x and gamma(x) the entries do not have to be fixed.
          * For conflict analysis, only the fixed entries matter.
          */
         if ( ( i < perm[i] || i == invperm[varrow] ) && ISFIXED(vars[i], bdchgidx) )
         {
            assert( vars[i] != infervar );
            SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
         }
         if ( ( invperm[i] > i || invperm[i] == varrow ) && ISFIXED(vars[invperm[i]], bdchgidx) )
         {
            assert( vars[invperm[i]] != infervar );
            SCIP_CALL( SCIPaddConflictUb(scip, vars[invperm[i]], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[invperm[i]], bdchgidx) );
         }
      }
   }
   else
   {
      /* This is not a fixing caused by lookahead (checkFeasible),
       * so row "varrow" was (0, _) or (_, 1) and for i < varrow x[i] = perm(x)[i].
       */
      if ( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         /* Changed the lower bound of infervar to 1. That means that this fixing is due to (_, 1) */
         assert( infervar == vars[varrow] );
         assert( ISFIXED(vars[invperm[varrow]], bdchgidx) );

         if ( invperm[varrow] > varrow )
         {
            SCIP_CALL( SCIPaddConflictUb(scip, vars[invperm[varrow]], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[invperm[varrow]], bdchgidx) );
         }
      }
      else
      {
         /* Changed the lower bound of infervar to 0. That means that this fixing is due to (0, _) */
         assert( infervar == vars[invperm[varrow]] );
         assert( ISFIXED(vars[varrow], bdchgidx) );

         if ( varrow < perm[varrow] )
         {
            SCIP_CALL( SCIPaddConflictUb(scip, vars[varrow], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[varrow], bdchgidx) );
         }
      }
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** lock variables
 *
 *  We assume we have only one global (void) constraint and lock all binary variables
 *  which do not correspond to fixed points of the permutation.
 *
 * - Symresack constraints may get violated if the variables with a negative coefficient
 *   in the FD inequality are rounded down, we therefor call
 *   SCIPaddVarLocksType(..., nlockspos, nlocksneg).
 * - Symresack constraints may get violated if the variables with a positive coefficient
 *   in the FD inequality are rounded up, we therefor call
 *   SCIPaddVarLocksType(..., nlocksneg, nlockspo ).
 */
static
SCIP_DECL_CONSLOCK(consLockSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* perm;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMsg(scip, "Locking method for symresack constraint handler.\n");

   /* get data of original constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* we do not have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->perm != NULL );

   nvars = consdata->nvars;
   vars = consdata->vars;
   perm = consdata->perm;

   for (i = 0; i < nvars; ++i)
   {
      /* due to clean-up in consdataCreate, there are no fixed points */
      assert( perm[i] != i );

      if ( perm[i] > i )
      {
         SCIP_CALL( SCIPaddVarLocksType(scip, vars[i], locktype, nlockspos, nlocksneg) );
      }
      else
      {
         SCIP_CALL( SCIPaddVarLocksType(scip, vars[i], locktype, nlocksneg, nlockspos) );
      }
   }

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySymresack)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_VAR** sourcevars;
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( sourcescip != NULL );
   assert( sourceconshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(sourceconshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( varmap != NULL );
   assert( valid != NULL );

   *valid = TRUE;

   SCIPdebugMsg(scip, "Copying method for symresack constraint handler.\n");

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->vars != NULL );
   assert( sourcedata->perm != NULL );
   assert( sourcedata->nvars > 0 );

   conshdlrdata = SCIPconshdlrGetData(sourceconshdlr);
   assert( conshdlrdata != NULL );

   /* do not copy non-model constraints */
   if ( !sourcedata->ismodelcons && !conshdlrdata->forceconscopy )
   {
      *valid = FALSE;

      return SCIP_OKAY;
   }

   sourcevars = sourcedata->vars;
   nvars = sourcedata->nvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   for (i = 0; i < nvars && *valid; ++i)
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[i], &(vars[i]), varmap, consmap, global, valid) );
      assert( !(*valid) || vars[i] != NULL );
   }

   /* only create the target constraint, if all variables could be copied */
   if ( *valid )
   {
      /* create copied constraint */
      if ( name == NULL )
         name = SCIPconsGetName(sourcecons);

      SCIP_CALL( SCIPcreateConsSymresack(scip, cons, name, sourcedata->perm, vars, nvars, sourcedata->ismodelcons,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseSymresack)
{  /*lint --e{715}*/
   const char* s;
   char* endptr;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   int* perm;
   int val;
   int nvars = 0;
   int cnt = 0;
   int nfoundpermidx = 0;
   int maxnvars = 128;

   assert( success != NULL );

   *success = TRUE;
   s = str;

   /* skip white space */
   while ( *s != '\0' && isspace((unsigned char)*s) )
      ++s;

   if ( strncmp(s, "symresack(", 10) != 0 )
   {
      SCIPerrorMessage("Syntax error - expected \"symresack(\", but got '%s'", s);
      *success = FALSE;
      return SCIP_OKAY;
   }
   s += 10;

   /* loop through string */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, maxnvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, maxnvars) );

   do
   {
      if ( cnt > 1 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected two arrays, but got more\n");
         *success = FALSE;

         SCIPfreeBufferArray(scip, &perm);
         SCIPfreeBufferArray(scip, &vars);
      }

      /* skip whitespace and ',' */
      while ( *s != '\0' && ( isspace((unsigned char)*s) ||  *s == ',' ) )
         ++s;

      /* if we could not find starting indicator of array */
      if ( *s != '[' )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected '[' to start new array\n");
         *success = FALSE;

         SCIPfreeBufferArray(scip, &perm);
         SCIPfreeBufferArray(scip, &vars);
      }
      ++s;

      /* read array, cnt = 0: variables; cnt = 1: permutation*/
      if ( cnt == 0 )
      {
         do
         {
            /* skip whitespace and ',' */
            while ( *s != '\0' && ( isspace((unsigned char)*s) ||  *s == ',' ) )
               ++s;

            /* parse variable name */
            SCIP_CALL( SCIPparseVarName(scip, s, &var, &endptr) );
            if ( var == NULL )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable name at '%s'\n", str);
               *success = FALSE;
               SCIPfreeBufferArray(scip, &vars);
               return SCIP_OKAY;
            }
            s = endptr;
            assert( s != NULL );

            vars[nvars++] = var;

            if ( nvars >= maxnvars )
            {
               int newsize;

               newsize = SCIPcalcMemGrowSize(scip, nvars + 1);
               SCIP_CALL( SCIPreallocBufferArray(scip, &vars, newsize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &perm, newsize) );
               maxnvars = newsize;
            }
         }
         while ( *s != ']' );
      }
      else
      {
         do
         {
            /* skip whitespace and ',' */
            while ( *s != '\0' && ( isspace((unsigned char)*s) ||  *s == ',' ) )
               ++s;

            /* parse integer value */
            if ( ! SCIPstrToIntValue(s, &val, &endptr) )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "could not extract int from string '%s'\n", str);
               *success = FALSE;
               SCIPfreeBufferArray(scip, &perm);
               SCIPfreeBufferArray(scip, &vars);
               return SCIP_OKAY;
            }
            s = endptr;
            assert( s != NULL );

            perm[nfoundpermidx++] = val;

            if ( nfoundpermidx > nvars )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "permutation is longer than vars array\n");
               *success = FALSE;
               SCIPfreeBufferArray(scip, &perm);
               SCIPfreeBufferArray(scip, &vars);
               return SCIP_OKAY;
            }
         }
         while ( *s != ']' );
      }
      ++s;
      ++cnt;
   }
   while ( *s != ')' );

   if ( nfoundpermidx == nvars )
   {
      SCIP_CALL( SCIPcreateConsBasicSymresack(scip, cons, name, perm, vars, nvars, TRUE) );
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL,
         "Length of permutation is not equal to number of given variables.\n");
      *success = FALSE;
   }

   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint display method of constraint handler
 *
 *  The constraint handler should output a representation of the constraint into the given text file.
 */
static
SCIP_DECL_CONSPRINT(consPrintSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* perm;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIPdebugMsg(scip, "Printing method for symresack constraint handler\n");

   /* we do not have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->perm != NULL );

   vars = consdata->vars;
   nvars = consdata->nvars;
   perm = consdata->perm;

   SCIPinfoMessage(scip, file, "symresack([");
   SCIP_CALL( SCIPwriteVarName(scip, file, vars[0], TRUE) );

   for (i = 1; i < nvars; ++i)
   {
      SCIPinfoMessage(scip, file, ",");
      SCIP_CALL( SCIPwriteVarName(scip, file, vars[i], TRUE) );
   }
   SCIPinfoMessage(scip, file, "],[%d", perm[0]);
   for (i = 1; i < nvars; ++i)
      SCIPinfoMessage(scip, file, ",%d", perm[i]);
   SCIPinfoMessage(scip, file, "])");

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( success != NULL );
   assert( vars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      int cnt = 0;
      int i;

      for (i = 0; i < consdata->nvars; ++i)
         vars[cnt++] = consdata->vars[i];
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( success != NULL );
   assert( nvars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}


/** creates the handler for symresack constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSymresack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSymresack, consEnfopsSymresack, consCheckSymresack, consLockSymresack,
         conshdlrdata) );
   assert( conshdlr != NULL );

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySymresack, consCopySymresack) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSymresack) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSymresack) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSymresack) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSymresack) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSymresack) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSymresack) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSymresack, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSymresack) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSymresack, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSymresack) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSymresack, consSepasolSymresack, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSymresack) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpSymresack) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolSymresack) );

   /* whether we allow upgrading to packing/partioning symresack constraints*/
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/ppsymresack",
         "Upgrade symresack constraints to packing/partioning symresacks?",
         &conshdlrdata->checkppsymresack, TRUE, DEFAULT_PPSYMRESACK, NULL, NULL) );

   /* whether we check for monotonicity of perm when upgrading to packing/partioning symresacks */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/checkmonotonicity",
         "Check whether permutation is monotone when upgrading to packing/partioning symresacks?",
         &conshdlrdata->checkmonotonicity, TRUE, DEFAULT_CHECKMONOTONICITY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/forceconscopy",
         "Whether symresack constraints should be forced to be copied to sub SCIPs.",
         &conshdlrdata->forceconscopy, TRUE, DEFAULT_FORCECONSCOPY, NULL, NULL) );

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates and captures a symresack constraint
 *
 *  In a presolving step, we check whether the permutation acts only on binary points. Otherwise, we eliminate
 *  the non-binary variables from the permutation.
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
SCIP_RETCODE SCIPcreateConsSymresack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables in vars array */
   SCIP_Bool             ismodelcons,        /**< whether the symresack is a model constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( nvars > 0 );

   /* find the symresack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("Symresack constraint handler not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, conshdlr, &consdata, vars, nvars, perm, ismodelcons) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate && (! consdata->ppupgrade), enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}


/** creates and captures a symresack constraint
 *  in its most basic variant, i.e., with all constraint flags set to their default values
 *
 *  In a presolving step, we remove all fixed points and cycles that act on non-binary variables of the permutation
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
SCIP_RETCODE SCIPcreateConsBasicSymresack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables in vars array */
   SCIP_Bool             ismodelcons         /**< whether the symresack is a model constraint */
   )
{
   SCIP_CALL( SCIPcreateConsSymresack(scip, cons, name, perm, vars, nvars, ismodelcons,
         TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
