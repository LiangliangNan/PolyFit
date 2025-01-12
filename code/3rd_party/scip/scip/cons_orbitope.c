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

/**@file   cons_orbitope.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for (partitioning/packing/full) orbitope constraints w.r.t. the full symmetric group
 * @author Timo Berthold
 * @author Marc Pfetsch
 * @author Christopher Hojny
 *
 * The type of constraints of this constraint handler is described in cons_orbitope.h.
 *
 * The details of the method implemented here are described in the following papers.
 *
 * Packing and Partitioning Orbitopes@n
 * Volker Kaibel and Marc E. Pfetsch,@n
 * Math. Program. 114, No. 1, 1-36 (2008)
 *
 * Among other things, this paper describes so-called shifted column inequalities of the following
 * form \f$x(S) \leq x(B)\f$, where \f$S\f$ is a so-called shifted column and \f$B\f$ is a so-called
 * bar. These inequalities can be used to handle symmetry and they are separated in this constraint
 * handler. We use the linear time separation algorithm of the paper.@par
 *
 * Orbitopal Fixing@n
 * Volker Kaibel, Matthias Peinhardt, and Marc E. Pfetsch,@n
 * Discrete Optimization 8, No. 4, 595-610 (2011)
 * (A preliminary version appears in Proc. IPCO 2007.)
 *
 * In this paper a linear time propagation algorithm is described, a variant of which is implemented
 * here. The implemented variant does not run in linear time, but is very fast in practice.
 *
 * <table>
 *   <caption>translation table</caption>
 *   <tr><td>here</td><td>paper</td></tr>
 *   <tr><td></td><td></td></tr>
 *   <tr><td>nspcons      </td><td>p       </td></tr>
 *   <tr><td>nblocks      </td><td>q       </td></tr>
 *   <tr><td>vars         </td><td>x       </td></tr>
 *   <tr><td>vals         </td><td>A^\\star</td></tr>
 *   <tr><td>weights      </td><td>\\omega </td></tr>
 *   <tr><td>cases        </td><td>\\tau   </td></tr>
 *   <tr><td>fixtriangle  </td><td>--      </td></tr>
 *   <tr><td>resolveprop  </td><td>--      </td></tr>
 *   <tr><td>firstnonzeros</td><td>\\mu    </td></tr>
 *   <tr><td>lastones     </td><td>\\alpha </td></tr>
 *   <tr><td>frontiersteps</td><td>\\Gamma </td></tr>
 * </table>
 *
 * Orbitopal fixing for the full (sub-)orbitope and application to the Unit Commitment Problem@n
 * Pascale Bendotti, Pierre Fouilhoux, and Cecile Rottner,@n
 * Optimization Online: http://www.optimization-online.org/DB_HTML/2017/10/6301.html
 *
 * Two linear time propagation algorithms for full orbitopes are described in this paper, a static
 * version and a dynamic one. While the static version uses a fixed variable order, the dynamic
 * version determines the variable order during the solving process via branching descisions.
 * We implemented the static version as well as a modified version of the dynamic one. The reason
 * for the latter is to simplify the compatibility with full orbitope cutting planes.
 *
 * Note, however, that the dynamic version may lead to conflicts if orbitopes are copied to subSCIPs.
 * Since the dynamic version is based on branching decisions, which may be different in main SCIP
 * and subSCIPs, orbitopes using the dynamic algorithm are not allowed to be copied. However, as a
 * user might use orbitopes to enforce a certain variable ordering in a solution, we distinguish
 * whether an orbitope is a model constraint or not. If it is a model constraint, we assume that
 * a variable order has already been fixed and disable the dynamic algorithm. In this case, orbitope
 * constraints are copied to subSCIPs. If it is not a model constraint, the orbitope was added to
 * handle symmetries but not to enforce a solution to have a certain structure. In this case, the
 * dynamic algorithm can be used and we do not copy orbitope constraints to subSCIPs.
 *
 * Polytopes associated with symmetry handling@n
 * Christopher Hojny and Marc E. Pfetsch,@n
 * Math. Program. (2018)
 *
 * In this paper, a linear time separation algorithm for orbisacks (full orbitopes with two columnes)
 * is described. We use this algorithm for every pair of adjacent columns within the orbitope as well
 * as a version that is adapted to the reordering based on the dynamic full orbitope propagation
 * algorithm to ensure validity of binary points via cutting planes.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_orbisack.h"
#include "scip/cons_orbitope.h"
#include "scip/cons_setppc.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/scip.h"
#include "scip/scip_branch.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_cut.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_var.h"
#include "scip/symmetry.h"
#include <ctype.h>
#include <string.h>
#include <symmetry/type_symmetry.h>

/* constraint handler properties */
#define CONSHDLR_NAME          "orbitope"
#define CONSHDLR_DESC          "symmetry breaking constraint handler relying on (partitioning/packing) orbitopes"
#define CONSHDLR_SEPAPRIORITY    +40100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -1005200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -1005200 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler */
#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */

#define DEFAULT_PPORBITOPE         TRUE /**< whether we check if full orbitopes can be strengthened to packing/partitioning orbitopes */
#define DEFAULT_SEPAFULLORBITOPE  FALSE /**< whether we separate inequalities for full orbitopes */
#define DEFAULT_FORCECONSCOPY     FALSE /**< whether orbitope constraints should be forced to be copied to sub SCIPs */

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             checkpporbitope;    /**< whether we allow upgrading to packing/partitioning orbitopes */
   SCIP_Bool             sepafullorbitope;   /**< whether we separate inequalities for full orbitopes orbitopes */
   SCIP_Bool             forceconscopy;      /**< whether orbitope constraints should be forced to be copied to sub SCIPs */
};

/** constraint data for orbitope constraints */
struct SCIP_ConsData
{
   SCIP_VAR***           vars;               /**< matrix of variables on which the symmetry acts            */
   SCIP_VAR**            tmpvars;            /**< temporary storage for variables                           */
   SCIP_HASHMAP*         rowindexmap;        /**< map of variables to row index in orbitope matrix */
   SCIP_Real**           vals;               /**< LP-solution for those variables                           */
   SCIP_Real*            tmpvals;            /**< temporary storage for values                              */
   SCIP_Real**           weights;            /**< SC weight table                                           */
   int**                 cases;              /**< indicator of the SC cases                                 */
   int                   nspcons;            /**< number of set partitioning/packing constraints  <=> p     */
   int                   nblocks;            /**< number of symmetric variable blocks             <=> q     */
   SCIP_ORBITOPETYPE     orbitopetype;       /**< type of orbitope constraint                               */
   SCIP_Bool             resolveprop;        /**< should propagation be resolved?                           */
   SCIP_Bool             istrianglefixed;    /**< has the upper right triangle already globally been fixed to zero?  */
   int*                  roworder;           /**< order of orbitope rows if dynamic propagation for full orbitopes
                                              *   is used. */
   SCIP_Bool*            rowused;            /**< whether a row has been considered in roworder */
   int                   nrowsused;          /**< number of rows that have already been considered in roworder */
   SCIP_Bool             ismodelcons;        /**< whether the orbitope is a model constraint */
   SCIP_Bool             mayinteract;        /**< whether symmetries corresponding to orbitope might interact
                                              *   with symmetries handled by other routines */
   SCIP_Bool             usedynamicprop;     /**< whether we use a dynamic version of the propagation routine */
};


/*
 * Local methods
 */

/** frees an orbitope constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to orbitope constraint data */
   )
{
   int p;
   int q;

   assert( consdata != NULL );
   assert( *consdata != NULL );

   if ( (*consdata)->usedynamicprop && (*consdata)->rowindexmap != NULL )
   {
      SCIPhashmapFree(&((*consdata)->rowindexmap));
   }

   p = (*consdata)->nspcons;
   q = (*consdata)->nblocks;
   for (int i = 0; i < p; ++i)
   {
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cases[i]), q);    /*lint !e866*/
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars[i]), q);     /*lint !e866*/
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->weights[i]), q);  /*lint !e866*/
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vals[i]), q);     /*lint !e866*/
   }

   if ( (*consdata)->usedynamicprop )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->rowused), p);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->roworder), p);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cases), p);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars), p);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->weights), p);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vals), p);

   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->tmpvals), p + q);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->tmpvars), p + q);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** creates orbitope constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure                                     */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data                        */
   SCIP_VAR***           vars,               /**< variables array, must have size nspcons x nblocks       */
   int                   nspcons,            /**< number of set partitioning (packing) constraints  <=> p */
   int                   nblocks,            /**< number of symmetric variable blocks               <=> q */
   SCIP_ORBITOPETYPE     orbitopetype,       /**< type of orbitope constraint                             */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved?                         */
   SCIP_Bool             usedynamicprop,     /**< whether we use a dynamic version of the propagation routine */
   SCIP_Bool             ismodelcons,        /**< whether the orbitope is a model constraint */
   SCIP_Bool             mayinteract         /**< whether symmetries corresponding to orbitope might interact
                                              *   with symmetries handled by other routines */
   )
{
   int i;
   int j;

   assert(consdata != NULL);
#ifndef NDEBUG
   if ( usedynamicprop )
   {
      assert( ! mayinteract );
   }
#endif

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vals, nspcons) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->weights, nspcons) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vars, nspcons) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->cases, nspcons) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->roworder, nspcons) );

   /* if orbitope might interact with other symmetries, we cannot adapt the row order of orbitopes dynamically */
   if ( usedynamicprop )
   {
      SCIP_CALL( SCIPhashmapCreate(&(*consdata)->rowindexmap, SCIPblkmem(scip), nspcons) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->rowused, nspcons) );
   }

   for (i = 0; i < nspcons; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vals[i], nblocks) );                 /*lint !e866*/
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->weights[i], nblocks) );              /*lint !e866*/
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars[i], vars[i], nblocks) );    /*lint !e866*/
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->cases[i], nblocks) );                /*lint !e866*/
      (*consdata)->roworder[i] = i;

      if ( usedynamicprop )
      {
         (*consdata)->rowused[i] = FALSE;
      }
   }
   (*consdata)->nrowsused = 0;

   (*consdata)->tmpvals = NULL;
   (*consdata)->tmpvars = NULL;
   (*consdata)->nspcons = nspcons;
   (*consdata)->nblocks = nblocks;
   (*consdata)->orbitopetype = orbitopetype;
   (*consdata)->resolveprop = resolveprop;
   (*consdata)->istrianglefixed = FALSE;
   (*consdata)->ismodelcons = ismodelcons;
   (*consdata)->mayinteract = mayinteract;
   (*consdata)->usedynamicprop = usedynamicprop;

   /* get transformed variables, if we are in the transformed problem */
   if ( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->tmpvals, nspcons + nblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->tmpvars, nspcons + nblocks) );

      for (i = 0; i < nspcons; ++i)
      {
         /* make sure that no variable gets multiaggregated (cannot be handled by cons_orbitope, since one cannot easily
          * eliminate single variables from an orbitope constraint).
          */
         for (j = 0; j < nblocks; ++j)
         {
            SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->vars[i][j], &(*consdata)->vars[i][j]) );
            SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*consdata)->vars[i][j]) );
            if ( usedynamicprop )
            {
               SCIP_CALL( SCIPhashmapInsert((*consdata)->rowindexmap, (*consdata)->vars[i][j], (void*) (size_t) i) );
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** strengthen full orbitopes to packing/partitioning orbitopes if possible */
static
SCIP_RETCODE strengthenOrbitopeConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< variable matrix of orbitope constraint */
   int*                  nrows,              /**< pointer to number of rows of variable matrix */
   int                   ncols,              /**< number of columns of variable matrix */
   SCIP_ORBITOPETYPE*    type,               /**< pointer to store type of orbitope constraint after strengthening */
   SCIP_Bool             mayinteract         /**< whether symmetries corresponding to orbitope might interact
                                              *   with symmetries handled by other routines */
   )
{
   SCIP_Bool* pprows = NULL;
   int npprows;
   int nrowsorig;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( vars != NULL );
   assert( *nrows > 0 );
   assert( ncols > 0 );
   assert( type != NULL );

   nrowsorig = *nrows;
   SCIP_CALL( SCIPisPackingPartitioningOrbitope(scip, vars, *nrows, ncols, &pprows, &npprows, type) );

   /* If only some rows are contained in set packing/partitioning constraints, it may still be worth it
    * to exploit the packing/partitioning structure on these rows, because packing/partitioning orbitopes
    * or more restrictive than full orbitopes. If at least three rows have this property, we discard
    * all rows not contained in set packing/partitioning constraints and add the smaller sub packing orbitope.
    * This is only possible if the orbitope's symmetries do not interact with other symmetry handling
    * methods (otherwise, dropping rows might change the variable order).
    */
   if ( npprows >= 3 && ! mayinteract )
   {
      int r = *nrows - 1;
      int i;

      assert( pprows != NULL );

      while ( r >= 0 )
      {
         if ( ! pprows[r] )
         {
            for (i = r; i < *nrows - 1; ++i)
            {
               SCIP_VAR** row;
               row = vars[i];
               vars[i] = vars[i+1];
               vars[i+1] = row;
            }
            *nrows -= 1;
         }
         --r;
      }
      *type = SCIP_ORBITOPETYPE_PACKING;
   }

   /* pprows might not have been initialized if there are no setppc conss */
   if ( pprows != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &pprows, nrowsorig);
   }

   return SCIP_OKAY;
}

#ifdef PRINT_MATRIX
/** debug method, prints variable matrix */
static
void printMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< the constraint data */
   )
{
   int i;
   int j;

   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vars != NULL );

   for (j = 0; j < consdata->nblocks; ++j)
      SCIPinfoMessage(scip, NULL, "-");

   SCIPinfoMessage(scip, NULL, "\n");
   for (i = 0; i < consdata->nspcons; ++i)
   {
      for (j = 0; j < consdata->nblocks; ++j)
      {
         if ( SCIPvarGetUbLocal(consdata->vars[i][j]) - SCIPvarGetLbLocal(consdata->vars[i][j]) < 0.5 )
            SCIPinfoMessage(scip, NULL, "%1.0f", REALABS(SCIPvarGetUbLocal(consdata->vars[i][j])));
         else
            SCIPinfoMessage(scip, NULL, " ");
      }
      SCIPinfoMessage(scip, NULL, "|\n");
   }
   for (j = 0; j < consdata->nblocks; ++j)
      SCIPinfoMessage(scip, NULL, "-");
   SCIPinfoMessage(scip, NULL, "\n");
}
#endif


#ifdef SHOW_SCI
/** Print SCI in nice form for debugging */
static
SCIP_RETCODE printSCI(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   p,                  /**< number of rows */
   int                   q,                  /**< number of columns */
   int**                 cases,              /**< SCI dynamic programming table */
   int                   i,                  /**< row position of bar */
   int                   j                   /**< column position of bar */
   )
{
   int k;
   int l;
   int** M;
   int p1;
   int p2;

   SCIP_CALL( SCIPallocBufferArray(scip, &M, p) );
   for (k = 0; k < p; ++k)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &M[k], q) ); /*lint !e866*/
      for (l = 0; l < q; ++l)
         M[k][l] = 0;
   }

   /* first add bar */
   for (l = j; l < q; ++l)
   {
      assert( M[i][l] == 0 );
      M[i][l] = 1;
   }

   /* then add shifted column */
   p1 = i-1;
   p2 = j-1;
   do
   {
      assert( cases[p1][p2] != -1 );
      assert( p1 >= 0 && p1 < i );
      assert( p2 >= 0 && p2 < j );

      /* if case 1 */
      if ( cases[p1][p2] == 1 )
         --p2;   /* decrease column */
      else
      {
         /* case 2 or 3: */
         assert( cases[p1][p2] == 2 || cases[p1][p2] == 3 );
         assert( M[p1][p2] == 0 );
         M[p1][p2] = -1;
         if ( cases[p1][p2] == 3 )
            break;
      }
      --p1;  /* decrease row */
   }
   while ( p1 >= 0 );   /* should always be true, i.e. the break should end the loop */
   assert( cases[p1][p2] == 3 );

   /* now output matrix M */
   for (l = 0; l < q; ++l)
      SCIPinfoMessage(scip, NULL, "-");
   SCIPinfoMessage(scip, NULL, "\n");

   for (k = 0; k < p; ++k)
   {
      for (l = 0; l < q; ++l)
      {
         if ( l > k )
            SCIPinfoMessage(scip, NULL, "*");
         else
         {
            switch (M[k][l])
            {
            case 1:
               SCIPinfoMessage(scip, NULL, "+");
               break;
            case -1:
               SCIPinfoMessage(scip, NULL, "-");
               break;
            case 0:
               SCIPinfoMessage(scip, NULL, "#");
               break;
            default:
               SCIPerrorMessage("unexpected matrix entry <%d>: should be -1, 0 or +1\n", M[k][l]);
               SCIPABORT();
            }
         }
      }
      SCIPinfoMessage(scip, NULL, "\n");
   }

   for (l = 0; l < q; ++l)
      SCIPinfoMessage(scip, NULL, "-");
   SCIPinfoMessage(scip, NULL, "\n");

   for (k = 0; k < p; ++k)
      SCIPfreeBufferArray(scip, &M[k]);
   SCIPfreeBufferArray(scip, &M);

   return SCIP_OKAY;
}
#endif


/** copies the variables values from the solution to the constraint data structure */
static
void copyValues(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< the constraint data     */
   SCIP_SOL*             sol                 /**< a primal solution or NULL for the current LP optimum */
   )
{
   int i;
   int j;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->vals != NULL );

   for (i = 0; i < consdata->nspcons; ++i)
   {
      for (j = 0; j < consdata->nblocks; ++j)
         consdata->vals[i][j] = SCIPgetSolVal(scip, sol, consdata->vars[i][j]);
   }
}


/** compute the dynamic programming table for SC
 *
 *  Build up dynamic programming table in order to find SCs with minimum weight.
 *
 *  The values of the minimal SCIs are stored in @a weights.
 *  The array @a cases[i][j] stores which of the cases were applied to get @a weights[i][j].
 *  Here, 3 means that we have reached the upper limit.
 *
 *  We assume that the upper right triangle is fixed to 0. Hence we can perform the computation a
 *  bit more efficient.
 */
static
void computeSCTable(
   SCIP*                 scip,               /**< SCIP pointer                                              */
   int                   nspcons,            /**< number of set partitioning (packing) constraints  <=> p   */
   int                   nblocks,            /**< number of symmetric variable blocks               <=> q   */
   SCIP_Real**           weights,            /**< SC weight table                                           */
   int**                 cases,              /**< indicator of the SC cases                                 */
   SCIP_Real**           vals                /**< current solution                                          */
   )
{
   SCIP_Real minvalue;
   int diagsize;
   int i;
   int j;

   assert( weights != NULL );
   assert( cases != NULL );
   assert( vals != NULL );

#ifndef NDEBUG
   /* for debugging */
   for (i = 0; i < nspcons; ++i)
   {
      for (j = 0; j < nblocks; ++j)
      {
         if ( i >= j )
         {
            weights[i][j] = -1.0;
            cases[i][j] = -1;
         }
      }
   }
#endif

   /* initialize diagonal */
   minvalue = vals[0][0];
   weights[0][0] = minvalue;
   cases[0][0] = 3;

   /* get last row of triangle */
   diagsize = nblocks;
   if ( nspcons < nblocks )
      diagsize = nspcons;

   for (j = 1; j < diagsize; ++j)
   {
      /* use LT to move entry as far to the left as possible */
      if ( SCIPisLT(scip, vals[j][j], minvalue) )
      {
         minvalue = vals[j][j];
         cases[j][j] = 3;
      }
      else
         cases[j][j] = 1;
      weights[j][j] = minvalue;
   }

   /* initialize first column */
   for (i = 1; i < nspcons; ++i)
   {
      weights[i][0] = weights[i-1][0] + vals[i][0];
      cases[i][0] = 2;  /* second case */
   }

   /* build the table */
   for (i = 2; i < nspcons; ++i)
   {
      for (j = 1; j < nblocks && j < i; ++j)
      {
         SCIP_Real weightleft;
         SCIP_Real weightright;

         assert( cases[i-1][j] != -1 );
         assert( cases[i-1][j-1] != -1 );

         weightleft = weights[i-1][j-1];
         weightright = vals[i][j] + weights[i-1][j];

         /* For first column: cannot take left possibility */
         if ( SCIPisLT(scip, weightleft, weightright) )
         {
            weights[i][j] = weightleft;
            cases[i][j] = 1;
         }
         else
         {
            weights[i][j] = weightright;
            cases[i][j] = 2;
         }
      }
   }
}


/** fix upper right triangle if necessary */
static
SCIP_RETCODE fixTriangle(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be processed */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   SCIP_Bool fixedglobal;
   SCIP_Bool fixed;
   int diagsize;
   int nspcons;
   int nblocks;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vars != NULL );

   *infeasible = FALSE;
   *nfixedvars = 0;

   if ( consdata->istrianglefixed )
      return SCIP_OKAY;

   nspcons = consdata->nspcons;
   nblocks = consdata->nblocks;
   vars = consdata->vars;
   fixedglobal = TRUE;

   /* get last row of triangle */
   diagsize = nblocks;
   if ( nspcons < nblocks )
      diagsize = nspcons;

   /* fix variables to 0 */
   for (i = 0; i < diagsize; ++i)
   {
      for (j = i+1; j < nblocks; ++j)
      {
         /* fix variable, if not in the root the fixation is local */
         SCIP_CALL( SCIPfixVar(scip, vars[i][j], 0.0, infeasible, &fixed) );

         if ( *infeasible )
         {
            SCIPdebugMsg(scip, "The problem is infeasible: some variable in the upper right triangle is fixed to 1.\n");
            return SCIP_OKAY;
         }

         if ( fixed )
            ++(*nfixedvars);

         if ( SCIPvarGetUbGlobal(vars[i][j]) > 0.5 )
            fixedglobal = FALSE;
      }
   }
   if ( *nfixedvars > 0 )
   {
      SCIPdebugMsg(scip, "<%s>: %s fixed upper right triangle to 0 (fixed vars: %d).\n", SCIPconsGetName(cons), fixedglobal ? "globally" : "locally", *nfixedvars);
   }
   else
   {
      SCIPdebugMsg(scip, "<%s>: Upper right triangle already fixed to 0.\n", SCIPconsGetName(cons));
   }

   if ( fixedglobal )
      consdata->istrianglefixed = TRUE;

   return SCIP_OKAY;
}


/** separates shifted column inequalities according to the solution stored in consdata->vals */
static
SCIP_RETCODE separateSCIs(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< the constraint data */
   SCIP_Bool*            infeasible,         /**< whether we detected infeasibility */
   int*                  nfixedvars,         /**< pointer to store the number of variables fixed */
   int*                  ncuts               /**< pointer to store number of separated SCIs */
   )
{
   SCIP_Real** vals;
   SCIP_Real** weights;
   SCIP_Real* tmpvals;
   SCIP_VAR*** vars;
   SCIP_VAR** tmpvars;
   int** cases;
   int nspcons;
   int nblocks;
   int i;
   int j;
   int l;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL);
   assert( nfixedvars != NULL );
   assert( ncuts != NULL );

   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->vals != NULL );
   assert( consdata->tmpvars != NULL );
   assert( consdata->tmpvals != NULL );
   assert( consdata->weights != NULL );
   assert( consdata->cases != NULL );

   *infeasible = FALSE;
   *nfixedvars = 0;
   *ncuts = 0;

   nspcons = consdata->nspcons;
   nblocks = consdata->nblocks;
   vars = consdata->vars;
   vals = consdata->vals;
   tmpvars = consdata->tmpvars;
   tmpvals = consdata->tmpvals;
   weights = consdata->weights;
   cases = consdata->cases;

   /* check for upper right triangle */
   if ( ! consdata->istrianglefixed )
   {
      SCIP_CALL( fixTriangle(scip, cons, infeasible, nfixedvars) );
      if ( *infeasible )
         return SCIP_OKAY;
      if ( *nfixedvars > 0 )
         return SCIP_OKAY;
   }

   /* compute table if necessary (i.e., not computed before) */
   computeSCTable(scip, nspcons, nblocks, weights, cases, vals);

   /* loop through rows */
   for (i = 1; i < nspcons && ! (*infeasible); ++i)
   {
      SCIP_Real bar;       /* value of bar: */
      int lastcolumn;      /* last column considered as part of the bar */

      bar = 0.0;
      lastcolumn = nblocks - 1;
      if ( lastcolumn > i )
         lastcolumn = i;

      /* traverse row from right to left: */
      /* j >= 1, since for j = 0, i.e., the bar is a complete row, there does not exist an SCI */
      for (j = lastcolumn; j > 0; --j)
      {
         bar += vals[i][j];

         /* check whether weights[i-1][j-1] < bar  (<=> bar - weights[i-1][j-1] > 0), i.e. cut is violated) */
         if ( SCIPisEfficacious(scip, bar - weights[i-1][j-1]) )
         {
            SCIP_Real weight;
            SCIP_ROW* row;
#ifdef SCIP_DEBUG
            char name[SCIP_MAXSTRLEN];
#endif
            int nvars;
            int p1;
            int p2;

            nvars = 0;
            p1 = i-1;
            p2 = j-1;
            weight = 0.0;

            /* first add bar */
            for (l = j; l <= lastcolumn; ++l)
            {
               tmpvars[nvars] = vars[i][l];
               tmpvals[nvars] = 1.0;
               nvars++;
            }

            /* then add shifted column */
            do
            {
               assert( cases[p1][p2] != -1 );
               assert( p1 >= 0 && p1 < i );
               assert( p2 >= 0 && p2 < j );

               /* if case 1 */
               if (cases[p1][p2] == 1)
                  p2--;   /* decrease column */
               else
               {
                  /* case 2 or 3: */
                  assert( cases[p1][p2] == 2 || cases[p1][p2] == 3 );
                  tmpvars[nvars] = vars[p1][p2];
                  tmpvals[nvars] = -1.0;
                  nvars++;
                  weight += vals[p1][p2];
                  if ( cases[p1][p2] == 3 )
                     break;
               }
               p1--;  /* decrease row */
            }
            while ( p1 >= 0 );   /* should always be true, i.e. the break should end the loop */
            assert( cases[p1][p2] == 3 );

            /* generate cut */
#ifdef SCIP_DEBUG
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sci_%d_%d", i, j);
            SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, name, -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#else
            SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#endif
            SCIP_CALL( SCIPaddVarsToRow(scip, row, nvars, tmpvars, tmpvals) );
            /*SCIP_CALL( SCIPprintRow(scip, row, NULL) ); */
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
            ++(*ncuts);

#ifdef SHOW_SCI
            SCIP_CALL( printSCI(scip, nspcons, nblocks, cases, i, j) );
#endif

            assert( SCIPisSumEQ(scip, weights[i-1][j-1], weight) );
         }
      }
   }
   return SCIP_OKAY;
}


/** propagation method for a single packing or partitioning orbitope constraint */
static
SCIP_RETCODE propagatePackingPartitioningCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be processed */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   SCIP_ORBITOPETYPE orbitopetype;
   int* firstnonzeros;
   int* lastones;
   int* frontiersteps;
   int lastoneprevrow;
   int nspcons;
   int nblocks;
   int nsteps;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   *nfixedvars = 0;

   if( !SCIPallowStrongDualReds(scip) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vars != NULL );

   nspcons = consdata->nspcons;
   nblocks = consdata->nblocks;
   vars = consdata->vars;
   orbitopetype = consdata->orbitopetype;

   assert( orbitopetype == SCIP_ORBITOPETYPE_PACKING || orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING );

   /* fix upper right triangle if still necessary */
   if ( ! consdata->istrianglefixed )
   {
      int nfixed = 0;
      SCIP_CALL( fixTriangle(scip, cons, infeasible, &nfixed) );
      *nfixedvars += nfixed;
   }

   /* prepare further propagation */
   SCIP_CALL( SCIPallocBufferArray(scip, &firstnonzeros, nspcons) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lastones, nspcons) );
   SCIP_CALL( SCIPallocBufferArray(scip, &frontiersteps, nblocks) );

#ifdef PRINT_MATRIX
   SCIPdebugMsg(scip, "Matrix:\n");
   printMatrix(scip, consdata);
#endif

   /* propagate */
   lastoneprevrow = 0;
   lastones[0] = 0;

   if ( orbitopetype == SCIP_ORBITOPETYPE_PACKING )
   {
      /* packing case: if entry (0,0) is fixed to 0 */
      if ( SCIPvarGetUbLocal(vars[0][0]) < 0.5 )
      {
         lastoneprevrow = -1;
         lastones[0] = -1;
      }
   }
   nsteps = 0;

   for (i = 1; i < nspcons; ++i)
   {
      int lastcolumn;
      int firstnonzeroinrow;
      int lastoneinrow;
      SCIP_Bool infrontier;

      /* last column considered as part of the bar: */
      lastcolumn = nblocks - 1;
      if ( lastcolumn > i )
         lastcolumn = i;

      /* find first position not fixed to 0 (partitioning) or fixed to 1 (packing) */
      firstnonzeroinrow = -1;
      for (j = 0; j <= lastcolumn; ++j)
      {
         if ( orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
         {
            /* partitioning case: if variable is not fixed to 0 */
            if ( SCIPvarGetUbLocal(vars[i][j]) > 0.5 )
            {
               firstnonzeroinrow = j;
               break;
            }
         }
         else
         {
            /* packing case: if variable is fixed to 1 */
            if ( SCIPvarGetLbLocal(vars[i][j]) > 0.5 )
            {
               firstnonzeroinrow = j;
               break;
            }
         }
      }
      /* if all variables are fixed to 0 in the partitioning case - should not happen */
      if ( firstnonzeroinrow == -1 && orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
      {
         SCIPdebugMsg(scip, " -> Infeasible node: all variables in row %d are fixed to 0.\n", i);
         *infeasible = TRUE;
         /* conflict should be analyzed by setppc constraint handler */
         goto TERMINATE;
      }
      firstnonzeros[i] = firstnonzeroinrow;
      assert( orbitopetype == SCIP_ORBITOPETYPE_PACKING || firstnonzeroinrow >= 0 );
      assert( -1 <= firstnonzeroinrow && firstnonzeroinrow <= lastcolumn );

      /* compute rightmost possible position for a 1 */
      assert( orbitopetype == SCIP_ORBITOPETYPE_PACKING || 0 <= lastoneprevrow );
      assert( lastoneprevrow <= lastcolumn );

      /* if we are at right border or if entry in column lastoneprevrow+1 is fixed to 0 */
      infrontier = FALSE;
      assert( lastoneprevrow + 1 >= 0 );
      if ( lastoneprevrow == nblocks-1 || SCIPvarGetUbLocal(vars[i][lastoneprevrow+1]) < 0.5 ) /*lint !e679*/
         lastoneinrow = lastoneprevrow;
      else
      {
         lastoneinrow = lastoneprevrow + 1;
         frontiersteps[nsteps++] = i;
         infrontier = TRUE;
      }

      /* store lastoneinrow */
      assert( orbitopetype == SCIP_ORBITOPETYPE_PACKING || 0 <= lastoneinrow );
      assert( lastoneinrow <= lastcolumn );
      lastones[i] = lastoneinrow;

      /* check whether we are infeasible */
      if ( firstnonzeroinrow > lastoneinrow )
      {
         int k;

#ifdef SCIP_DEBUG
         if ( orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
         {
            SCIPdebugMsg(scip, " -> Infeasible node: row %d, leftmost nonzero at %d, rightmost 1 at %d\n",
               i, firstnonzeroinrow, lastoneinrow);
         }
         else
         {
            SCIPdebugMsg(scip, " -> Infeasible node: row %d, 1 at %d, rightmost position for 1 at %d\n",
               i, firstnonzeroinrow, lastoneinrow);
         }
#endif
         /* check if conflict analysis is applicable */
         if ( SCIPisConflictAnalysisApplicable(scip) )
         {
            /* conflict analysis only applicable in SOLVING stage */
            assert( SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPinProbing(scip) );

            /* perform conflict analysis */
            SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

            if ( orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
            {
               /* add bounds (variables fixed to 0) that result in the first nonzero entry */
               for (j = 0; j <= lastcolumn; ++j)
               {
                  /* add varaibles in row up to the first variable fixed to 0 */
                  if ( SCIPvarGetUbLocal(vars[i][j]) > 0.5 )
                     break;

                  assert( SCIPvarGetUbLocal(vars[i][j]) < 0.5 );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) );
               }
            }
            else
            {
               /* add bounds that result in the last one - check top left entry for packing case */
               if ( lastones[0] == -1 )
               {
                  assert( SCIPvarGetUbLocal(vars[0][0]) < 0.5 );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[0][0]) );
               }

               /* mark variable fixed to 1 */
               assert( SCIPvarGetLbLocal(vars[i][firstnonzeroinrow]) > 0.5 );
               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][firstnonzeroinrow]) );
            }

            /* add bounds that result in the last one - pass through rows */
            for (k = 1; k < i; ++k)
            {
               int l;
               l = lastones[k] + 1;

               /* if the frontier has not moved and we are not beyond the matrix boundaries */
               if ( l <= nblocks-1 && l <= k && lastones[k-1] == lastones[k] )
               {
                  assert( SCIPvarGetUbLocal(vars[k][l]) < 0.5 );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[k][l]) );
               }
            }
            SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
         }

         *infeasible = TRUE;
         goto TERMINATE;
      }

      /* fix entries beyond the last possible position for a 1 in the row to 0 (see Lemma 1 in the paper) */
      for (j = lastoneinrow+1; j <= lastcolumn; ++j)
      {
         /* if the entry is not yet fixed to 0 */
         if ( SCIPvarGetUbLocal(vars[i][j]) > 0.5 )
         {
            SCIP_Bool tightened;
            int inferInfo;

            SCIPdebugMsg(scip, " -> Fixing entry (%d,%d) to 0.\n", i, j);

            tightened = FALSE;

            /* fix variable to 0 and store position of (i,lastoneinrow+1) for conflict resolution */
            inferInfo = i * nblocks + lastoneinrow + 1;
            /* correction according to Lemma 1 in the paper (second part): store (i,lastoneinrow+2) */
            if ( !infrontier )
               ++inferInfo;
            SCIP_CALL( SCIPinferBinvarCons(scip, vars[i][j], FALSE, cons, inferInfo, infeasible, &tightened) );

            /* if entry is fixed to one -> infeasible node */
            if ( *infeasible )
            {
               SCIPdebugMsg(scip, " -> Infeasible node: row %d, 1 in column %d beyond rightmost position %d\n", i, j, lastoneinrow);
               /* check if conflict analysis is applicable */
               if( SCIPisConflictAnalysisApplicable(scip) )
               {
                  int k;

                  /* conflict analysis only applicable in SOLVING stage */
                  assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPinProbing(scip));

                  /* perform conflict analysis */
                  SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

                  /* add current bound */
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) );

                  /* add bounds that result in the last one - check top left entry for packing case */
                  if ( orbitopetype == SCIP_ORBITOPETYPE_PACKING && lastones[0] == -1 )
                  {
                     assert( SCIPvarGetUbLocal(vars[0][0]) < 0.5 );
                     SCIP_CALL( SCIPaddConflictBinvar(scip, vars[0][0]) );
                  }

                  /* add bounds that result in the last one - pass through rows */
                  for (k = 1; k < i; ++k)
                  {
                     int l;
                     l = lastones[k] + 1;

                     /* if the frontier has not moved and we are not beyond the matrix boundaries */
                     if ( l <= nblocks-1 && l <= k && lastones[k-1] == lastones[k] )
                     {
                        assert( SCIPvarGetUbLocal(vars[k][l]) < 0.5 );
                        SCIP_CALL( SCIPaddConflictBinvar(scip, vars[k][l]) );
                     }
                  }
                  SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
               }

               goto TERMINATE;
            }
            if ( tightened )
               ++(*nfixedvars);
         }
      }

      lastoneprevrow = lastoneinrow;
   }

   /* check whether fixing any entry to 0 results in a contradiction -> loop through rows in frontiersteps (a.k.a. gamma) */
   for (j = 0; j < nsteps; ++j)
   {
      int s;
      int lastoneinrow;

      s = frontiersteps[j];
      lastoneinrow = lastones[s];
      /* note for packing case: if we are in a frontier step then lastoneinrow >= 0 */
      assert( 0 <= lastoneinrow && lastoneinrow < nblocks );

      /* if entry is not fixed */
      if ( SCIPvarGetLbLocal(vars[s][lastoneinrow]) < 0.5 && SCIPvarGetUbLocal(vars[s][lastoneinrow]) > 0.5 )
      {
         int betaprev;
         betaprev = lastoneinrow - 1;

         /* loop through rows below s */
         for (i = s+1; i < nspcons; ++i)
         {
            int beta;

            assert( betaprev + 1 >= 0 );
            if ( betaprev == nblocks-1 || SCIPvarGetUbLocal(vars[i][betaprev+1]) < 0.5 ) /*lint !e679*/
               beta = betaprev;
            else
               beta = betaprev + 1;
            assert( -1 <= beta && beta < nblocks );

            if ( firstnonzeros[i] > beta )
            {
               SCIP_Bool tightened;
               int inferInfo;

               /* can fix (s,lastoneinrow) (a.k.a (s,alpha)) to 1
                * (do not need to fix other entries to 0, since they will be
                * automatically fixed by SCIPtightenVarLb.)
                */
               assert( SCIPvarGetLbLocal(vars[s][lastoneinrow]) < 0.5 );
               SCIPdebugMsg(scip, " -> Fixing entry (%d,%d) to 1.\n", s, lastoneinrow);

               tightened = FALSE;

               /* store position (i,firstnonzeros[i]) */
               inferInfo = nblocks * nspcons + i * nblocks + firstnonzeros[i];
               SCIP_CALL( SCIPinferBinvarCons(scip, vars[s][lastoneinrow], TRUE, cons, inferInfo, infeasible, &tightened) );

               assert( !(*infeasible) );
               if ( tightened )
                  ++(*nfixedvars);
               break;
            }
            betaprev = beta;
         }
      }
   }

 TERMINATE:
   SCIPfreeBufferArray(scip, &frontiersteps);
   SCIPfreeBufferArray(scip, &lastones);
   SCIPfreeBufferArray(scip, &firstnonzeros);

   return SCIP_OKAY;
}


/* Compute dynamic order of rows based on the branching decisions, i.e., the row of the first branching variable
 * determines the first row in the new ordering, the row of the second branching variable determines the second
 * row in the new ordering if it differs from the row of the first branching variable, and so on.
 *
 * The roworder array stores this reordering, where acutally only the first maxrowlabel entries encode the
 * reordering.
 */
static
SCIP_RETCODE computeDynamicRowOrder(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_HASHMAP*         rowindexmap,        /**< map of variables to indices in orbitope vars matrix */
   SCIP_Bool*            rowused,            /**< bitset marking whether a row has been considered in the new order */
   int*                  roworder,           /**< reordering of the rows w.r.t. branching decisions */
   int                   nrows,              /**< number of rows in orbitope */
   int                   ncols,              /**< number of columns in orbitope */
   int*                  maxrowlabel         /**< maximum row label in ordering */
   )
{
   int i;
   SCIP_NODE* node;
   int* branchdecisions;
   int nbranchdecision;

   assert( scip != NULL );
   assert( rowindexmap != NULL );
   assert( roworder != NULL );
   assert( nrows > 0 );
   assert( maxrowlabel != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &branchdecisions, nrows * ncols) );
   nbranchdecision = 0;

   /* get current node */
   node = SCIPgetCurrentNode(scip);

   /* follow path to the root (in the root no domains were changed due to branching) */
   while ( SCIPnodeGetDepth(node) != 0 )
   {
      SCIP_BOUNDCHG* boundchg;
      SCIP_DOMCHG* domchg;
      SCIP_VAR* branchvar;
      int nboundchgs;

      /* get domain changes of current node */
      domchg = SCIPnodeGetDomchg(node);
      assert( domchg != NULL );

      /* loop through all bound changes */
      nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
      for (i = 0; i < nboundchgs; ++i)
      {
         /* get bound change info */
         boundchg = SCIPdomchgGetBoundchg(domchg, i);
         assert( boundchg != NULL );

         /* branching decisions have to be in the beginning of the bound change array */
         if ( SCIPboundchgGetBoundchgtype(boundchg) != SCIP_BOUNDCHGTYPE_BRANCHING )
            break;

         /* get corresponding branching variable */
         branchvar = SCIPboundchgGetVar(boundchg);

         /* we only consider binary variables */
         if ( SCIPvarGetType(branchvar) == SCIP_VARTYPE_BINARY )
         {
            int rowidx;

            /* make sure that branching variable is present in the orbitope */
            if ( ! SCIPhashmapExists(rowindexmap, (void*) branchvar) )
               continue;

            rowidx = (int) (size_t) SCIPhashmapGetImage(rowindexmap, (void*) branchvar);
            branchdecisions[nbranchdecision++] = rowidx;
         }
      }

      node = SCIPnodeGetParent(node);
   }

   /* Insert branching decisions of current path into global row order.
    * Iterate in reverse order over branching decisions to get the path
    * from the root to the current node.
    */
   for (i = nbranchdecision - 1; i >= 0; --i)
   {
      if ( ! rowused[branchdecisions[i]] )
      {
         roworder[*maxrowlabel] = branchdecisions[i];
         rowused[branchdecisions[i]] = TRUE;
         *maxrowlabel += 1;
      }
   }

   SCIPfreeBufferArray(scip, &branchdecisions);

   return SCIP_OKAY;
}


/* Compute lexicographically minimal face of the hypercube w.r.t. some coordinate fixing */
static
SCIP_RETCODE findLexMinFace(
   SCIP_VAR***           vars,               /**< variable matrix */
   int**                 lexminfixes,        /**< fixings characterzing lex-min face */
   int*                  minfixedrowlexmin,  /**< index of minimum fixed row for each column or
                                              *   NULL (if in prop) */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility has been
                                              *   detected or NULL (if in resprop) */
   int                   m,                  /**< number of rows in vars */
   int                   n,                  /**< number of columns in vars */
   int                   nrowsused,          /**< number of rows considered in propagation */
   SCIP_Bool             resprop             /**< whether we are in resprop (TRUE) or prop (FALSE) */
   )
{
   int i;
   int j;
   *infeasible = FALSE;

   assert( vars != NULL );
   assert( lexminfixes != NULL );
   assert( !resprop || minfixedrowlexmin != NULL );
   assert( m > 0 );
   assert( n > 0 );
   assert( nrowsused > 0 );
   assert( nrowsused <= m );
   assert( infeasible != NULL );

   /* iterate over columns in reverse order and find the lexicographically minimal face
    * of the hypercube containing lexminfixes
    */
   for (j = n - 2; j >= 0; --j)
   {
      int maxdiscriminating = m;
      int minfixed = -1;

      /* fix free entries in column j to the corresponding value in column j + 1 and collect some information */
      for (i = 0; i < nrowsused; ++i)
      {
         /* is row i j-discriminating? */
         if ( minfixed == -1 && lexminfixes[i][j] != 0 && lexminfixes[i][j + 1] != 1 )
         {
            assert( lexminfixes[i][j + 1] == 0 );

            maxdiscriminating = i;
         }

         /* is row i j-fixed? */
         if ( minfixed == -1 && lexminfixes[i][j] != lexminfixes[i][j + 1] && lexminfixes[i][j] != 2 )
         {
            assert( lexminfixes[i][j + 1] != 2 );

            minfixed = i;

            /* detect infeasibility */
            if ( maxdiscriminating > minfixed )
            {
               *infeasible = TRUE;

               return SCIP_OKAY;
            }
         }
      }

      /* ensure that column j is lexicographically not smaller than column j + 1 */
      for (i = 0; i < nrowsused; ++i)
      {
         if ( lexminfixes[i][j] == 2 )
         {
            if ( i < maxdiscriminating || minfixed == -1 )
               lexminfixes[i][j] = lexminfixes[i][j + 1];
            else if ( i == maxdiscriminating )
               lexminfixes[i][j] = 1;
            else
               lexminfixes[i][j] = 0;
         }
      }

      if ( resprop )
      {
         assert( minfixedrowlexmin != NULL );

         /* store minimum fixed row */
         if ( minfixed == -1 )
            minfixedrowlexmin[j] = nrowsused - 1;
         else
            minfixedrowlexmin[j] = minfixed;

         /* columns 1, ..., n-2 are contained in two columns (take the minimum) and
          * the minimum fixed row of column n-1 is determined by column n-2 */
         if ( minfixedrowlexmin[j + 1] < minfixedrowlexmin[j] )
            minfixedrowlexmin[j + 1] = minfixedrowlexmin[j];
      }
   }

   return SCIP_OKAY;
}


/* Compute lexicographically maximal face of the hypercube w.r.t. some coordinate fixing */
static
SCIP_RETCODE findLexMaxFace(
   SCIP_VAR***           vars,               /**< variable matrix */
   int**                 lexmaxfixes,        /**< fixings characterzing lex-max face */
   int*                  minfixedrowlexmax,  /**< index of minimum fixed row for each column or
                                              *   NULL (if in prop) */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility has been
                                              *   detected or NULL (if in resprop) */
   int                   m,                  /**< number of rows in vars */
   int                   n,                  /**< number of columns in vars */
   int                   nrowsused,          /**< number of rows considered in propagation */
   SCIP_Bool             resprop             /**< whether we are in resprop (TRUE) or prop (FALSE) */
   )
{
   int i;
   int j;
   *infeasible = FALSE;

   assert( vars != NULL );
   assert( lexmaxfixes != NULL );
   assert( !resprop || minfixedrowlexmax != NULL );
   assert( m > 0 );
   assert( n > 0 );
   assert( nrowsused > 0 );
   assert( nrowsused <= m );
   assert( infeasible != NULL );

   for (j = 1; j < n; ++j)
   {
      int maxdiscriminating = m;
      int minfixed = -1;

      /* fix free entries in column j to the corresponding value in column j - 1 and collect some information */
      for (i = 0; i < nrowsused; ++i)
      {
         /* is row i j-discriminating? */
         if ( minfixed == -1 && lexmaxfixes[i][j - 1] != 0 && lexmaxfixes[i][j] != 1 )
         {
            assert( lexmaxfixes[i][j - 1] == 1 );

            maxdiscriminating = i;
         }

         /* is row i j-fixed? */
         if ( minfixed == -1 && lexmaxfixes[i][j - 1] != lexmaxfixes[i][j] && lexmaxfixes[i][j] != 2 )
         {
            assert( lexmaxfixes[i][j - 1] != 2 );

            minfixed = i;

            /* detect infeasibility */
            if ( maxdiscriminating > minfixed )
            {
               *infeasible = TRUE;

               return SCIP_OKAY;
            }
         }
      }

      /* ensure that column j is lexicographically not greater than column j - 1 */
      for (i = 0; i < nrowsused; ++i)
      {
         if ( lexmaxfixes[i][j] == 2 )
         {
            if ( i < maxdiscriminating || minfixed == -1 )
               lexmaxfixes[i][j] = lexmaxfixes[i][j - 1];
            else if ( i == maxdiscriminating )
               lexmaxfixes[i][j] = 0;
            else
               lexmaxfixes[i][j] = 1;
         }
      }

      if ( resprop )
      {
         assert( minfixedrowlexmax != NULL );

         /* store minimum fixed row */
         if ( minfixed == -1 )
            minfixedrowlexmax[j] = nrowsused - 1;
         else
            minfixedrowlexmax[j] = minfixed;

         /* columns 1, ..., n-2 are contained in two columns (take the minimum) and
          * the minimum fixed row of column 0 is determined by column 1 */
         if ( minfixedrowlexmax[j - 1] < minfixedrowlexmax[j] )
            minfixedrowlexmax[j - 1] = minfixedrowlexmax[j];
      }
   }

   return SCIP_OKAY;
}


/** propagation method for a single packing or partitioning orbitope constraint */
static
SCIP_RETCODE propagateFullOrbitopeCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be processed */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars,         /**< pointer to add up the number of found domain reductions */
   SCIP_Bool             dynamic             /**< whether we use a dynamic propagation routine */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int** lexminfixes;
   int** lexmaxfixes;
   int* roworder;
   int nrowsused;
   int i;
   int j;
   int m;
   int n;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   *nfixedvars = 0;
   *infeasible = FALSE;

   /* @todo Can the following be removed? */
   if ( ! SCIPallowStrongDualReds(scip) )
      return SCIP_OKAY;

   /* do nothing if we use dynamic propagation and if we are in a probing node */
   if ( dynamic && SCIPinProbing(scip) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->orbitopetype == SCIP_ORBITOPETYPE_FULL );

   m = consdata->nspcons;
   n = consdata->nblocks;
   vars = consdata->vars;

   /* determine order of orbitope rows dynamically by branching decisions */
   if ( dynamic )
   {
      SCIP_CALL( computeDynamicRowOrder(scip, consdata->rowindexmap, consdata->rowused,
            consdata->roworder, m, n, &(consdata->nrowsused)) );

      /* if no branching variable is contained in the full orbitope */
      if ( consdata->nrowsused == 0 )
         return SCIP_OKAY;

      nrowsused = consdata->nrowsused;
   }
   else
      nrowsused = m;
   roworder = consdata->roworder;

   /* Initialize lexicographically minimal matrix by fixed entries at the current node.
    * Free entries in the last column are set to 0.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexminfixes, nrowsused) );
   for (i = 0; i < nrowsused; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lexminfixes[i], n) ); /*lint !e866*/
   }

   for (i = 0; i < nrowsused; ++i)
   {
      int origrow;

      origrow = roworder[i];

      for (j = 0; j < n; ++j)
      {
         if ( SCIPvarGetLbLocal(vars[origrow][j]) > 0.5 )
            lexminfixes[i][j] = 1;
         else if ( SCIPvarGetUbLocal(vars[origrow][j]) < 0.5 || j == n - 1 )
            lexminfixes[i][j] = 0;
         else
            lexminfixes[i][j] = 2;
      }
   }

   /* find lexicographically minimal face of hypercube containing lexmin fixes */
   SCIP_CALL( findLexMinFace(vars, lexminfixes, NULL, infeasible, m, n, nrowsused, FALSE) );

   if ( *infeasible == TRUE )
      goto FREELEXMIN;

   /* Initialize lexicographically maximal matrix by fixed entries at the current node.
    * Free entries in the first column are set to 1.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxfixes, nrowsused) );
   for (i = 0; i < nrowsused; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxfixes[i], n) ); /*lint !e866*/
   }

   for (i = 0; i < nrowsused; ++i)
   {
      int origrow;

      origrow = roworder[i];

      for (j = 0; j < n; ++j)
      {
         if ( SCIPvarGetUbLocal(vars[origrow][j]) < 0.5 )
            lexmaxfixes[i][j] = 0;
         else if ( SCIPvarGetLbLocal(vars[origrow][j]) > 0.5 || j == 0 )
            lexmaxfixes[i][j] = 1;
         else
            lexmaxfixes[i][j] = 2;
      }
   }

   /* find lexicographically maximal face of hypercube containing lexmax fixes */
   SCIP_CALL( findLexMaxFace(vars, lexmaxfixes, NULL, infeasible, m, n, nrowsused, FALSE) );

   if ( *infeasible )
      goto FREELEXMAX;

   /* Find for each column j the minimal row in which lexminfixes and lexmaxfixes differ. Fix all entries above this
    * row to the corresponding value in lexminfixes (or lexmaxfixes).
    */
   for (j = 0; j < n; ++j)
   {
      for (i = 0; i < nrowsused; ++i)
      {
         int origrow;

         origrow = roworder[i];

         if ( lexminfixes[i][j] != lexmaxfixes[i][j] )
            break;

         if ( SCIPvarGetLbLocal(vars[origrow][j]) < 0.5 && SCIPvarGetUbLocal(vars[origrow][j]) > 0.5 )
         {
            SCIP_Bool success;

            SCIP_CALL( SCIPinferBinvarCons(scip, vars[origrow][j], (SCIP_Bool) lexminfixes[i][j],
                  cons, 0, infeasible, &success) );

            if ( success )
               *nfixedvars += 1;
         }
      }
   }

 FREELEXMAX:
   for (i = 0; i < nrowsused; ++i)
      SCIPfreeBufferArray(scip, &lexmaxfixes[i]);
   SCIPfreeBufferArray(scip, &lexmaxfixes);

 FREELEXMIN:
   for (i = 0; i < nrowsused; ++i)
      SCIPfreeBufferArray(scip, &lexminfixes[i]);
   SCIPfreeBufferArray(scip, &lexminfixes);

   return SCIP_OKAY;
}


/** propagation method for a single orbitope constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be processed */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_ORBITOPETYPE orbitopetype;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   orbitopetype = consdata->orbitopetype;

   if ( orbitopetype == SCIP_ORBITOPETYPE_FULL )
   {
      SCIP_CALL( propagateFullOrbitopeCons(scip, cons, infeasible, nfixedvars,
            consdata->usedynamicprop && !consdata->ismodelcons) );
   }
   else
   {
      assert( orbitopetype == SCIP_ORBITOPETYPE_PACKING || orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING );
      SCIP_CALL( propagatePackingPartitioningCons(scip, cons, infeasible, nfixedvars) );
   }

   return SCIP_OKAY;
}


/** Propagation conflict resolving method of propagator
 *
 *  In this function we use that all variable reductions that can be found by the propagation algorithm
 *  are only due to the fixed variables that are in or above the minimum fixed row of each pair of adjacent
 *  columns of the lexmin and lexmax matrices.
 *
 *  Since the storage of an integer is not enough to store the complete information about the fixing,
 *  we have to use the linear time algorithm for finding the lexmin and lexmax
 *  matrices and determine from this the minimum fixed rows.
 */
static
SCIP_RETCODE resolvePropagationFullOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler of the corresponding constraint */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   int                   inferinfo,          /**< inference information */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int** lexminfixes;
   int** lexmaxfixes;
   int* roworder;
   int* minfixedrowlexmin;
   int* minfixedrowlexmax;
   int i;
   int j;
   int m;
   int n;
   int nrowsused;
   SCIP_Bool dynamic;
   SCIP_Bool terminate;

   *result = SCIP_DIDNOTFIND;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( result != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->orbitopetype == SCIP_ORBITOPETYPE_FULL );

   dynamic = consdata->usedynamicprop && !consdata->ismodelcons;
   m = consdata->nspcons;
   n = consdata->nblocks;
   vars = consdata->vars;

   if ( dynamic )
   {
      assert( consdata->roworder != NULL );
      assert( consdata->nrowsused > 0 );

      nrowsused = consdata->nrowsused;
   }
   else
      nrowsused = m;
   roworder = consdata->roworder;

   assert( inferinfo <= consdata->nspcons );

   /* Initialize lexicographically minimal matrix by fixed entries at the current node.
    * Free entries in the last column are set to 0.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexminfixes, nrowsused) );
   for (i = 0; i < nrowsused; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lexminfixes[i], n) ); /*lint !e866*/
   }

   /* store minimum fixed row for each column */
   SCIP_CALL( SCIPallocBufferArray(scip, &minfixedrowlexmin, n) );
   minfixedrowlexmin[n - 1] = -1;

   for (i = 0; i < nrowsused; ++i)
   {
      int origrow;

      origrow = roworder[i];

      for (j = 0; j < n; ++j)
      {
         if ( SCIPvarGetLbAtIndex(vars[origrow][j], bdchgidx, FALSE) > 0.5 )
            lexminfixes[i][j] = 1;
         else if ( SCIPvarGetUbAtIndex(vars[origrow][j], bdchgidx, FALSE) < 0.5 || j == n - 1 )
            lexminfixes[i][j] = 0;
         else
            lexminfixes[i][j] = 2;
      }
   }

   /* find lexicographically minimal face of hypercube containing lexmin fixes */
   SCIP_CALL( findLexMinFace(vars, lexminfixes, minfixedrowlexmin, &terminate, m, n, nrowsused, TRUE) );

   if ( terminate )
      goto FREELEXMIN;

   /* Initialize lexicographically maximal matrix by fixed entries at the current node.
    * Free entries in the first column are set to 1.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxfixes, nrowsused) );
   for (i = 0; i < nrowsused; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxfixes[i], n) ); /*lint !e866*/
   }

   /* store minimum fixed row for each column */
   SCIP_CALL( SCIPallocBufferArray(scip, &minfixedrowlexmax, n) );
   minfixedrowlexmax[0] = -1;

   for (i = 0; i < nrowsused; ++i)
   {
      int origrow;

      origrow = roworder[i];

      for (j = 0; j < n; ++j)
      {
         if ( SCIPvarGetUbAtIndex(vars[origrow][j], bdchgidx, FALSE) < 0.5 )
            lexmaxfixes[i][j] = 0;
         else if ( SCIPvarGetLbAtIndex(vars[origrow][j], bdchgidx, FALSE) > 0.5 || j == 0 )
            lexmaxfixes[i][j] = 1;
         else
            lexmaxfixes[i][j] = 2;
      }
   }

   /* find lexicographically maximal face of hypercube containing lexmax fixes */
   SCIP_CALL( findLexMaxFace(vars, lexmaxfixes, minfixedrowlexmax, &terminate, m, n, nrowsused, TRUE) );

   if ( terminate )
      goto FREELEXMAX;

   /* Find for each column j the minimal row in which lexminfixes and lexmaxfixes differ. Fix all entries above this
    * row to the corresponding value in lexminfixes (or lexmaxfixes).
    */
   for (j = 0; j < n; ++j)
   {
      int ub = MAX(minfixedrowlexmin[j], minfixedrowlexmax[j]);

      for (i = 0; i <= ub; ++i)
      {
         int origrow;

         origrow = roworder[i];

         if ( SCIPvarGetLbAtIndex(vars[origrow][j], bdchgidx, FALSE) > 0.5 ||
            SCIPvarGetUbAtIndex(vars[origrow][j], bdchgidx, FALSE) < 0.5 )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[origrow][j]) );
            *result = SCIP_SUCCESS;
         }
      }
   }

 FREELEXMAX:
   SCIPfreeBufferArray(scip, &minfixedrowlexmax);
   for (i = 0; i < nrowsused; ++i)
      SCIPfreeBufferArray(scip, &lexmaxfixes[i]);
   SCIPfreeBufferArray(scip, &lexmaxfixes);

 FREELEXMIN:
   SCIPfreeBufferArray(scip, &minfixedrowlexmin);
   for (i = 0; i < nrowsused; ++i)
      SCIPfreeBufferArray(scip, &lexminfixes[i]);
   SCIPfreeBufferArray(scip, &lexminfixes);

   return SCIP_OKAY;
}


/** Propagation conflict resolving method of propagator
 *
 *  In this function we use that the propagation method above implicitly propagates SCIs, i.e., every
 *  fixing can also be gotten via an SCI-fixing.
 *
 *  Since the storage of an integer is not enough to store the complete information about the fixing
 *  nor a complete shifted column, we have to use the linear time algorithm for SCIs.
 *
 *  The inferinfo integer is set as follows:
 *
 *  - If a shifted column is fixed to 0 and the corresponding bar does not necessarily has value 1
 *    then we fix these entries to 0 and inferinfo is i * nblocks + j, where (i,j) is the leader of the
 *    bar. The SCI depends on whether i is in Gamma or not (see Lemma 1 in the paper and the comments
 *    above).
 *
 *  - If a bar has value 1 and the shifted column has one entry that is not fixed, it can be fixed to
 *    1 and inferinfo is (nspcons*nblocks) + i * nblocks + j, where (i,j) is the leader of the bar; see
 *    Proposition 1 (2c).
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   int                   inferinfo,          /**< inference information */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real** vals;
   SCIP_Real** weights;
   SCIP_VAR*** vars;
   SCIP_ORBITOPETYPE orbitopetype;
   int** cases;

   int i;
   int j;
   int nspcons;
   int nblocks;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( result != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->vals != NULL );
   assert( consdata->weights != NULL );
   assert( consdata->cases != NULL );
   assert( consdata->istrianglefixed );

   *result = SCIP_DIDNOTFIND;
   if ( ! consdata->resolveprop )
      return SCIP_OKAY;

   nspcons = consdata->nspcons;
   nblocks = consdata->nblocks;
   vars = consdata->vars;
   vals = consdata->vals;
   weights = consdata->weights;
   orbitopetype = consdata->orbitopetype;
   cases = consdata->cases;

   SCIPdebugMsg(scip, "Propagation resolution method of orbitope constraint using orbitopal fixing\n");

   /* fill table */
   for (i = 0; i < nspcons; ++i)
   {
      int lastcolumn;

      /* last column considered as part of the bar: */
      lastcolumn = nblocks - 1;
      if ( lastcolumn > i )
         lastcolumn = i;
      for (j = 0; j <= lastcolumn; ++j)
      {
         /* if the variable was fixed to zero at conflict time */
         if ( SCIPgetVarUbAtIndex(scip, vars[i][j], bdchgidx, FALSE) < 0.5 )
            vals[i][j] = 0.0;
         else
         {
            /* if the variable was fixed to one at conflict time */
            if ( SCIPgetVarLbAtIndex(scip, vars[i][j], bdchgidx, FALSE) > 0.5 )
               vals[i][j] = 2.0;
            else
               vals[i][j] = 1.0;
         }
      }
   }

#ifdef PRINT_MATRIX
   SCIPdebugMsg(scip, "Matrix:\n");
   printMatrix(scip, consdata);
#endif

   /* computation of table: this now minimizes the value of the shifted column */
   assert( consdata->istrianglefixed );
   computeSCTable(scip, nspcons, nblocks, weights, cases, vals);

   /* if we fixed variables in the bar to zero */
   assert( inferinfo >= 0 && inferinfo < 2 * nspcons * nblocks );
   if ( inferinfo < nspcons * nblocks )
   {
      int p1;
      int p2;
#ifdef SCIP_DEBUG
      char str[SCIP_MAXSTRLEN];
      char tmpstr[SCIP_MAXSTRLEN];
#endif

      i = (int) (inferinfo / nblocks);
      j = inferinfo % nblocks;
      assert( 0 <= i && i < nspcons );
      assert( 0 <= j && j < nblocks );

      /* find SCI with value 0 */
      assert( weights[i-1][j-1] < 0.5 );

      SCIPdebugMsg(scip, " -> reason for x[%d][%d] = ... = x[%d][%d] = 0 was the following SC:\n", i, j, i, MIN(i,nblocks));
#ifdef SCIP_DEBUG
      str[0] = '\0';
#endif

      p1 = i-1;
      p2 = j-1;
      do
      {
         assert( cases[p1][p2] != -1 );
         assert( p1 >= 0 && p1 < i );
         assert( p2 >= 0 && p2 < j );

         /* if case 1 */
         if ( cases[p1][p2] == 1 )
            --p2;   /* decrease column */
         else
         {
            /* case 2 or 3: */
            assert( cases[p1][p2] == 2 || cases[p1][p2] == 3 );
            assert( SCIPgetVarUbAtIndex(scip, vars[p1][p2], bdchgidx, FALSE) < 0.5 );
            SCIP_CALL( SCIPaddConflictUb(scip, vars[p1][p2], bdchgidx) );
            *result = SCIP_SUCCESS;

#ifdef SCIP_DEBUG
            (void) SCIPsnprintf(tmpstr, SCIP_MAXSTRLEN, " (%d,%d)", p1, p2);
            (void) strncat(str, tmpstr, SCIP_MAXSTRLEN);
#endif

            if ( cases[p1][p2] == 3 )
               break;
         }
         --p1;  /* decrease row */
      }
      while ( p1 >= 0 );   /* should always be true, i.e. the break should end the loop */
      assert( cases[p1][p2] == 3 );

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "%s\n", str);
#endif
   }
   else
   {
      int k;
      int p1;
      int p2;
#ifndef NDEBUG
      int pos1;
      int pos2;
#endif
#ifdef SCIP_DEBUG
      char str[SCIP_MAXSTRLEN];
      char tmpstr[SCIP_MAXSTRLEN];
#endif

      /* if we fixed a variable in the SC to 1 */
      inferinfo -= nspcons * nblocks;
      i = (int) inferinfo / nblocks;
      j = inferinfo % nblocks;
      assert( 0 <= i && i < nspcons );
      assert( 0 <= j && j < nblocks );

      /* In rare cases it might happen that we fixed a variable to 1, but the node later becomes infeasible by globally
       * fixing variables to 0. In this case, it might happen that we find a SC with value 0 instead of 1. We then
       * cannot use this SC to repropagate (and do not know how to reconstruct the original reasoning). */
      if ( weights[i-1][j-1] > 0.5 && weights[i-1][j-1] < 1.5 )
      {
         SCIPdebugMsg(scip, " -> reason for x[%d][%d] = 1 was the following SC:\n", i, j);
#ifdef SCIP_DEBUG
         (void) SCIPsnprintf(str, SCIP_MAXSTRLEN, "SC:");
#endif

         p1 = i-1;
         p2 = j-1;
#ifndef NDEBUG
         pos1 = -1;
         pos2 = -1;
#endif
         do
         {
            assert( cases[p1][p2] != -1 );
            assert( p1 >= 0 && p1 < i );
            assert( p2 >= 0 && p2 < j );

            /* if case 1 */
            if ( cases[p1][p2] == 1 )
               --p2;   /* decrease column */
            else
            {
               /* case 2 or 3: reason are formed by variables in SC fixed to 0 */
               assert( cases[p1][p2] == 2 || cases[p1][p2] == 3 );
               if ( SCIPgetVarUbAtIndex(scip, vars[p1][p2], bdchgidx, FALSE) < 0.5 )
               {
                  SCIP_CALL( SCIPaddConflictUb(scip, vars[p1][p2], bdchgidx) );
                  *result = SCIP_SUCCESS;

#ifdef SCIP_DEBUG
                  (void) SCIPsnprintf(tmpstr, SCIP_MAXSTRLEN, " (%d,%d)", p1, p2);
                  (void) strncat(str, tmpstr, SCIP_MAXSTRLEN);
#endif
               }
#ifndef NDEBUG
               else
               {
                  assert( SCIPgetVarLbAtIndex(scip, vars[p1][p2], bdchgidx, FALSE) < 0.5 );
                  assert( pos1 == -1 && pos2 == -1 );
                  pos1 = p1;
                  pos2 = p2;
               }
#endif
               if ( cases[p1][p2] == 3 )
                  break;
            }
            --p1;  /* decrease row */
         }
         while ( p1 >= 0 );   /* should always be true, i.e., the break should end the loop */
         assert( cases[p1][p2] == 3 );
         assert( pos1 >= 0 && pos2 >= 0 );

         /* distinguish partitioning/packing */
         if ( orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
         {
            /* partitioning case */
#ifdef SCIP_DEBUG
            (void) SCIPsnprintf(tmpstr, SCIP_MAXSTRLEN, "  before bar: ");
            (void) strncat(str, tmpstr, SCIP_MAXSTRLEN);
#endif
            /* add variables before the bar in the partitioning case */
            for (k = 0; k < j; ++k)
            {
               assert( SCIPgetVarUbAtIndex(scip, vars[i][k], bdchgidx, FALSE) < 0.5 );
               SCIP_CALL( SCIPaddConflictUb(scip, vars[i][k], bdchgidx) );
               *result = SCIP_SUCCESS;
#ifdef SCIP_DEBUG
               (void) SCIPsnprintf(tmpstr, SCIP_MAXSTRLEN, " (%d,%d)", i, k);
               (void) strncat(str, tmpstr, SCIP_MAXSTRLEN);
#endif
            }

#ifdef SCIP_DEBUG
            SCIPdebugMsg(scip, "%s\n", str);
#endif
         }
         else
         {
            /* packing case */
            int lastcolumn;

            /* last column considered as part of the bar: */
            lastcolumn = nblocks - 1;
            if ( lastcolumn > i )
               lastcolumn = i;

            /* search for variable in the bar that is fixed to 1 in the packing case */
            for (k = j; k <= lastcolumn; ++k)
            {
               if ( SCIPgetVarLbAtIndex(scip, vars[i][k], bdchgidx, FALSE) > 0.5 )
               {
                  SCIP_CALL( SCIPaddConflictLb(scip, vars[i][k], bdchgidx) );
                  *result = SCIP_SUCCESS;
                  SCIPdebugMsg(scip, "   and variable x[%d][%d] fixed to 1.\n", i, k);
                  break;
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** check packing/partitioning orbitope solution for feasibility */
static
SCIP_RETCODE enfopsPackingPartitioningOrbitopeSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pointer to orbitope constraint */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real** weights;
   SCIP_Real** vals;
   int** cases;
   int nspcons;
   int nblocks;
   int i;
   int j;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vals != NULL );
   assert( consdata->weights != NULL );
   assert( consdata->cases != NULL );

   /* check for upper right triangle */
   if ( ! consdata->istrianglefixed )
   {
      SCIP_Bool infeasible = FALSE;
      int nfixedvars = 0;

      SCIP_CALL( fixTriangle(scip, cons, &infeasible, &nfixedvars) );
      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if ( nfixedvars > 0 )
      {
         *result = SCIP_REDUCEDDOM;
         return SCIP_OKAY;
      }
   }

   nspcons = consdata->nspcons;
   nblocks = consdata->nblocks;
   vals = consdata->vals;
   weights = consdata->weights;
   cases = consdata->cases;

   /* get solution */
   copyValues(scip, consdata, NULL);
   SCIPdebugMsg(scip, "Enforcing (pseudo solutions) for orbitope constraint <%s>\n", SCIPconsGetName(cons));

   /* compute table */
   assert( consdata->istrianglefixed );
   computeSCTable(scip, nspcons, nblocks, weights, cases, vals);

   /* loop through rows */
   for (i = 1; i < nspcons; ++i)
   {
      SCIP_Real bar = 0.0;
      int lastcolumn;

      lastcolumn = nblocks - 1;

      /* last column considered as part of the bar: */
      if ( lastcolumn > i )
         lastcolumn = i;

      /* traverse row from right to left */
      for (j = lastcolumn; j > 0; --j)
      {
         bar += vals[i][j];
         assert( SCIPisIntegral(scip, vals[i][j]) );

         /* check whether weights[i-1][j-1] < bar  (<=> bar - weights[i-1][j-1] > 0), i.e. cut is violated) */
         if ( SCIPisGT(scip, bar - weights[i-1][j-1], 0.0) )
         {
            SCIPdebugMsg(scip, "Solution is infeasible.\n");
            *result = SCIP_INFEASIBLE;
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}


/** check packing/partitioning orbitope solution for feasibility */
static
SCIP_RETCODE checkPackingPartitioningOrbitopeSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pointer to orbitope constraint */
   SCIP_SOL*             sol,                /**< solution to be checked */
   SCIP_RESULT*          result,             /**< pointer to store the result of the enforcing call */
   SCIP_Bool             printreason         /**< whether reason for infeasibility should be printed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   SCIP_Real** vals;
   SCIP_Real** weights;
   int** cases;
   int nspcons;
   int nblocks;
   int i;
   int j;

   /* get data of constraint */
   assert( cons != 0 );
   consdata = SCIPconsGetData(cons);

   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->vals != NULL );
   assert( consdata->weights != NULL );
   assert( consdata->cases != NULL );

   nspcons = consdata->nspcons;
   nblocks = consdata->nblocks;
   vars = consdata->vars;
   vals = consdata->vals;
   weights  = consdata->weights;
   cases = consdata->cases;

   /* get solution */
   copyValues(scip, consdata, sol);
   SCIPdebugMsg(scip, "Checking orbitope constraint <%s> ...\n", SCIPconsGetName(cons));

   /* check upper right triangle (if not yet fixed to zero or in debug mode */
#ifdef NDEBUG
   if ( ! consdata->istrianglefixed )
#endif
   {
      int diagsize;

      /* get last row of triangle */
      diagsize = nblocks;
      if ( nspcons < nblocks )
         diagsize = nspcons;

      /* check variables */
      for (i = 0; i < diagsize; ++i)
      {
         for (j = i+1; j < nblocks; ++j)
         {
            if ( ! SCIPisFeasZero(scip, vals[i][j]) )
            {
               if ( printreason )
                  SCIPinfoMessage(scip, NULL, "variable x[%d][%d] = %f on upper right nonzero.\n", i, j, vals[i][j]);
               *result = SCIP_INFEASIBLE;
            }
         }
      }
   }

   /* compute table */
   computeSCTable(scip, nspcons, nblocks, weights, cases, vals);

   /* loop through rows */
   for (i = 1; i < nspcons; ++i)
   {
      SCIP_Real bar;
      int lastcolumn;

      lastcolumn = nblocks - 1;
      bar = 0.0;
      /* last column considered as part of the bar: */
      if ( lastcolumn > i )
         lastcolumn = i;

      /* traverse row from right to left */
      for (j = lastcolumn; j > 0; --j)
      {
         bar += vals[i][j];
         assert( SCIPisFeasIntegral(scip, vals[i][j]) );

         /* check whether weights[i-1][j-1] < bar  (<=> bar - weights[i-1][j-1] > 0), i.e. cut is violated) */
         if ( SCIPisGT(scip, bar - weights[i-1][j-1], 0.0) )
         {
            SCIPdebugMsg(scip, "Solution is infeasible.\n");
            *result = SCIP_INFEASIBLE;

            if ( printreason )
            {
               int l;
               int p1;
               int p2;

               SCIPinfoMessage(scip, NULL, "violated SCI: bar(");

               /* first output bar */
               for (l = j; l < nblocks; ++l)
                  SCIPinfoMessage(scip, NULL, "<%s> (%f)", SCIPvarGetName(vars[i][l]), consdata->vals[i][l]);

               SCIPinfoMessage(scip, NULL, ")  SC(");

               /* output shifted column */
               p1 = i-1;
               p2 = j-1;
               do
               {
                  assert( cases[p1][p2] != -1 );
                  assert( p1 >= 0 && p1 < i );
                  assert( p2 >= 0 && p2 < j );

                  /* if case 1 */
                  if (cases[p1][p2] == 1)
                     --p2;   /* decrease column */
                  else
                  {
                     /* case 2 or 3: */
                     assert( cases[p1][p2] == 2 || cases[p1][p2] == 3 );
                     SCIPinfoMessage(scip, NULL, "<%s> (%f)", SCIPvarGetName(vars[p1][p2]), consdata->vals[p1][p2]);
                     if ( cases[p1][p2] == 3 )
                        break;
                  }
                  --p1;  /* decrease row */
               }
               while ( p1 >= 0 );   /* should always be true, i.e. the break should end the loop */
               assert( cases[p1][p2] == 3 );

               SCIPinfoMessage(scip, NULL, ")");
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** check full orbitope solution for feasibility */
static
SCIP_RETCODE checkFullOrbitopeSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_SOL*             sol,                /**< solution to be checked */
   SCIP_Bool             printreason,        /**< whether reason for infeasibility should be printed */
   SCIP_Bool*            feasible            /**< memory address to store whether solution is feasible */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   int nrows;
   int ncols;
   int j;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( feasible != NULL );

   consdata = SCIPconsGetData(cons);

   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( ! consdata->ismodelcons ); /* non-model constraints are never checked */

   vars = consdata->vars;
   nrows = consdata->nspcons;
   ncols = consdata->nblocks;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars1, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars2, nrows) );

   /* iterate over adjacent columns of orbitope and check whether the first column in this
    * column pair is lexicographically not smaller than the second column in the pair */
   *feasible = TRUE;
   for (j = 1; j < ncols && *feasible; ++j)
   {
      for (i = 0; i < nrows; ++i)
      {
         vars1[i] = vars[i][j - 1];
         vars2[i] = vars[i][j];
      }

      SCIP_CALL( SCIPcheckSolutionOrbisack(scip, sol, vars1, vars2, nrows, printreason, feasible) );
   }

   SCIPfreeBufferArray(scip, &vars2);
   SCIPfreeBufferArray(scip, &vars1);

   return SCIP_OKAY;
}


/** separate orbisack cover inequalities */
static
SCIP_RETCODE separateCoversOrbisack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_SOL*             sol,                /**< solution to separate (NULL for the LP solution) */
   SCIP_Bool             dynamic,            /**< whether we use a dynamic row order */
   int*                  ngen,               /**< pointer to store number of generated cuts */
   SCIP_Bool*            infeasible          /**< pointer to store whether infeasibility has been detected */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int* roworder;
   int nrowsused;
   int nrows;
   int ncols;
   int i;
   int j;
   int origrow;
   SCIP_Real rhs;
   SCIP_Real lhs;
   SCIP_Real* coeffs1;
   SCIP_Real* coeffs2;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( ngen != NULL );
   assert( infeasible != NULL );

   *ngen = 0;
   *infeasible = FALSE;

   /* get basic data */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   vars = consdata->vars;
   nrows = consdata->nspcons;
   ncols = consdata->nblocks;
   nrowsused = dynamic ? consdata->nrowsused : nrows;
   roworder = consdata->roworder;

   /* allocate memory for cover inequalities */
   SCIP_CALL( SCIPallocBufferArray(scip, &coeffs1, nrowsused) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coeffs2, nrowsused) );

   lhs = 0.0;
   rhs = 0.0;

   /* separate orbisack cover inequalities for adjacent columns */
   for (j = 0; j < ncols - 1 && ! *infeasible; ++j)
   {
      SCIP_Real rowval;

      for (i = 0; i < nrowsused; ++i)
      {
         origrow = roworder[i];

         assert( origrow >= 0 );
         assert( origrow < nrows );

         rowval = SCIPgetSolVal(scip, sol, vars[origrow][j + 1]) - SCIPgetSolVal(scip, sol, vars[origrow][j]);

         /* check whether cover inequality is violated */
         if ( SCIPisEfficacious(scip, rowval + lhs - rhs) )
         {
            SCIP_ROW* row;
            int k;

            /* set coefficients for current inequality */
            coeffs1[i] = -1.0;
            coeffs2[i] = 1.0;

            /* add violated orbisack cover inequality */
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "orbisackcover", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            for (k = 0; k <= i; ++k)
            {
               int origrow2;

               origrow2 = roworder[k];

               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[origrow2][j], coeffs1[k]) );
               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[origrow2][j + 1], coeffs2[k]) );
            }
            SCIP_CALL( SCIPflushRowExtensions(scip, row) );

            SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
            SCIP_CALL( SCIPreleaseRow(scip, &row) );

            *ngen += 1;
            if ( *infeasible )
               break;

            /* reset coefficients for next inequality */
            coeffs1[i] = 0.0;
            coeffs2[i] = 0.0;
         }

         /* add argmax( 1 - vals[i][0], vals[i][1] ) as coefficient and ensure that both vars1[0] and vars2[0] are
          * contained in the LIFTED cover inequality */
         rowval = SCIPgetSolVal(scip, sol, vars[origrow][j]) + SCIPgetSolVal(scip, sol, vars[origrow][j + 1]);
         if ( SCIPisEfficacious(scip, 1.0 - rowval) )
         {
            coeffs1[i] = -1.0;
            coeffs2[i] = 0.0;
            lhs -= SCIPgetSolVal(scip, sol, vars[origrow][j]);

            /* apply lifting? */
            if ( i == 0 )
            {
               coeffs2[i] = 1.0;
               lhs += SCIPgetSolVal(scip, sol, vars[origrow][j + 1]);
            }
         }
         else
         {
            coeffs1[i] = 0.0;
            coeffs2[i] = 1.0;
            lhs += SCIPgetSolVal(scip, sol, vars[origrow][j]);
            rhs += 1.0;

            /* apply lifting? */
            if ( i == 0 )
            {
               coeffs1[i] = -1.0;
               lhs -= SCIPgetSolVal(scip, sol, vars[origrow][j]);
               rhs -= 1.0;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &coeffs1);
   SCIPfreeBufferArray(scip, &coeffs2);

   return SCIP_OKAY;
}


/** separate or enforce constraints */
static
SCIP_RETCODE separateConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_SOL*             sol,                /**< solution to separate (NULL for the LP solution) */
   SCIP_RESULT*          result,             /**< pointer to store the result (should be initialized) */
   SCIP_Bool             enforce             /**< whether we enforce orbitope constraints */
   )
{
   SCIP_Bool infeasible = FALSE;
   int nfixedvars = 0;
   int ncuts = 0;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* loop through constraints */
   for (c = 0; c < nconss && ! infeasible; c++)
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_CONSDATA* consdata;
      int nconsfixedvars = 0;
      int nconscuts = 0;
      SCIP_ORBITOPETYPE orbitopetype;

      assert( conss[c] != NULL );

      /* get data of constraint */
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* do not enforce non-model constraints */
      if ( enforce && !consdata->ismodelcons )
         continue;

      /* get solution */
      copyValues(scip, consdata, sol);

      /* separate */
      orbitopetype = consdata->orbitopetype;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      if ( orbitopetype == SCIP_ORBITOPETYPE_PACKING || orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
      {
         SCIP_CALL( separateSCIs(scip, conshdlr, conss[c], consdata, &infeasible, &nconsfixedvars, &nconscuts) );
         nfixedvars += nconsfixedvars;
      }
      else if ( conshdlrdata->sepafullorbitope )
      {
         SCIP_CALL( separateCoversOrbisack(scip, conss[c], sol, consdata->usedynamicprop && !consdata->ismodelcons, &nconscuts, &infeasible) );
      }
      ncuts += nconscuts;

      /* stop after the useful constraints if we found cuts of fixed variables */
      if ( c >= nusefulconss && (ncuts > 0 || nfixedvars > 0) )
         break;
   }

   if ( infeasible )
   {
      SCIPdebugMsg(scip, "Infeasible node.\n");
      *result = SCIP_CUTOFF;
   }
   else if ( nfixedvars > 0 )
   {
      SCIPdebugMsg(scip, "Fixed %d variables.\n", nfixedvars);
      *result = SCIP_REDUCEDDOM;
   }
   else if ( ncuts > 0 )
   {
      SCIPdebugMsg(scip, "Separated %dinequalities.\n", ncuts);
      *result = SCIP_SEPARATED;
   }
   else
   {
      SCIPdebugMsg(scip, "No violated inequality found during separation.\n");
   }

   return SCIP_OKAY;
}


/** check whether all variables in an orbitope constraint are fixed */
static
SCIP_RETCODE checkRedundantCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be processed */
   SCIP_Bool*            redundant           /**< pointer to store whether constraint is redundant (contains no active vars) */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int i;
   int j;
   int nrows;
   int ncols;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( redundant != NULL );

   *redundant = FALSE;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );

   vars = consdata->vars;
   nrows = consdata->nspcons;
   ncols = consdata->nblocks;

   /* check whether there exists an active variable in the orbitope */
   for (i = 0; i < nrows; ++i)
   {
      for (j = 0; j < ncols; ++j)
      {
         if ( SCIPvarIsActive(vars[i][j]) )
            return SCIP_OKAY;
      }
   }

   *redundant = TRUE;

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyOrbitope)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrOrbitope(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** frees constraint handler */
static
SCIP_DECL_CONSFREE(consFreeOrbitope)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteOrbitope)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransOrbitope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create linear constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->vars, sourcedata->nspcons, sourcedata->nblocks,
         sourcedata->orbitopetype, sourcedata->resolveprop, sourcedata->usedynamicprop, sourcedata->ismodelcons,
         sourcedata->mayinteract) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpOrbitope)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation of orbitope constraint handler <%s> for LP solution.\n", SCIPconshdlrGetName(conshdlr));

   *result = SCIP_DIDNOTRUN;

   /* if solution is integer, skip separation */
   if ( SCIPgetNLPBranchCands(scip) <= 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, nusefulconss, NULL, result, FALSE) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolOrbitope)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation of orbitope constraint handler <%s> for primal solution.\n", SCIPconshdlrGetName(conshdlr));

   *result = SCIP_DIDNOTFIND;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, nusefulconss, sol, result, FALSE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpOrbitope)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( result != NULL );

   /* we have a negative priority, so we should come after the integrality conshdlr */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   SCIPdebugMsg(scip, "Enforcement for orbitope constraint handler <%s> for LP solution.\n", SCIPconshdlrGetName(conshdlr));

   *result = SCIP_FEASIBLE;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, nusefulconss, NULL, result, TRUE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxOrbitope)
{  /*lint --e{715}*/
   assert( result != NULL );
   assert( scip != NULL );

   SCIPdebugMsg(scip, "Enforcement for orbitope constraint handler <%s> for relaxation solution.\n", SCIPconshdlrGetName(conshdlr));

   *result = SCIP_FEASIBLE;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, nusefulconss, sol, result, TRUE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsOrbitope)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;
   if ( objinfeasible || solinfeasible )
      return SCIP_OKAY;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;
      SCIP_ORBITOPETYPE orbitopetype;
      SCIP_Bool feasible;

      /* get data of constraint */
      cons = conss[c];
      assert( cons != 0 );
      consdata = SCIPconsGetData(cons);

      assert( consdata != NULL );

      /* do not enforce non-model constraints */
      if ( !consdata->ismodelcons )
         continue;

      orbitopetype = consdata->orbitopetype;

      if ( orbitopetype == SCIP_ORBITOPETYPE_PACKING || orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
      {
         SCIP_CALL( enfopsPackingPartitioningOrbitopeSolution(scip, cons, result) );
      }
      else
      {
         SCIP_CALL( checkFullOrbitopeSolution(scip, cons, NULL, FALSE, &feasible) );

         if ( ! feasible )
            *result = SCIP_INFEASIBLE;
      }

      if ( *result == SCIP_INFEASIBLE )
         break;
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckOrbitope)
{  /*lint --e{715}*/
   int c;
   SCIP_CONSDATA* consdata;
   SCIP_ORBITOPETYPE orbitopetype;
   SCIP_Bool feasible;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   /* loop through constraints */
   for( c = 0; c < nconss && (*result == SCIP_FEASIBLE || completely); ++c )
   {
      assert( conss[c] != 0 );
      consdata = SCIPconsGetData(conss[c]);

      assert( consdata != NULL );

      /* do not check non-model constraints */
      if ( !consdata->ismodelcons )
         continue;

      orbitopetype = consdata->orbitopetype;

      if ( orbitopetype == SCIP_ORBITOPETYPE_PACKING || orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
      {
         SCIP_CALL( checkPackingPartitioningOrbitopeSolution(scip, conss[c], sol, result, printreason) );
      }
      else
      {
         SCIP_CALL( checkFullOrbitopeSolution(scip, conss[c], sol, printreason, &feasible) );

         if ( ! feasible )
            *result = SCIP_INFEASIBLE;
      }
   }
   SCIPdebugMsg(scip, "Solution is feasible.\n");

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropOrbitope)
{  /*lint --e{715}*/
   SCIP_Bool infeasible = FALSE;
   int nfixedvars = 0;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* propagate all useful constraints */
   for (c = 0; c < nusefulconss && !infeasible; ++c)
   {
      assert( conss[c] != 0 );

      SCIPdebugMsg(scip, "Propagation of orbitope constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( propagateCons(scip, conss[c], &infeasible, &nfixedvars) );
   }

   /* return the correct result */
   if ( infeasible )
   {
      *result = SCIP_CUTOFF;
      SCIPdebugMsg(scip, "Propagation via orbitopal fixing proved node to be infeasible.\n");
   }
   else if ( nfixedvars > 0 )
   {
      *result = SCIP_REDUCEDDOM;
      SCIPdebugMsg(scip, "Propagated %d variables via orbitopal fixing.\n", nfixedvars);
   }
   else if ( nusefulconss > 0 )
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMsg(scip, "Propagation via orbitopal fixing did not find anything.\n");
   }

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolOrbitope)
{  /*lint --e{715}*/
   SCIP_Bool infeasible = FALSE;
   int noldfixedvars;
   int c;
   SCIP_Bool redundant;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;
   noldfixedvars = *nfixedvars;

   /* propagate all useful constraints
    *
    * @todo use an event handler to only propagate if a variable in the orbitope has been fixed
    */
   for (c = 0; c < nconss && !infeasible; ++c)
   {
      int nfixed = 0;

      assert( conss[c] != 0 );

      SCIPdebugMsg(scip, "Presolving of orbitope constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      /* first propagate */
      SCIP_CALL( propagateCons(scip, conss[c], &infeasible, &nfixed) );
      *nfixedvars += nfixed;

      if ( ! infeasible )
      {
         SCIP_CALL( checkRedundantCons(scip, conss[c], &redundant) );

         if ( redundant )
         {
            SCIPdebugMsg(scip, "Orbitope constraint <%s> is redundant: it does not contain active variables\n",
               SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelCons(scip, conss[c]) );
            assert( ! SCIPconsIsActive(conss[c]) );
            (*ndelconss)++;
            continue;
         }
      }
   }

   if ( infeasible )
   {
      *result = SCIP_CUTOFF;
      SCIPdebugMsg(scip, "Presolving detected infeasibility.\n");
   }
   else if ( *nfixedvars > noldfixedvars )
   {
      *result = SCIP_SUCCESS;
   }
   else if ( nconss > 0 )
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMsg(scip, "Presolving via orbitopal fixing did not find anything.\n");
   }

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropOrbitope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_ORBITOPETYPE orbitopetype;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   orbitopetype = consdata->orbitopetype;

   /* resolution for full orbitopes not availabe yet */
   if ( orbitopetype == SCIP_ORBITOPETYPE_PACKING || orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
   {
      SCIP_CALL( resolvePropagation(scip, cons, inferinfo, bdchgidx, result) );
   }
   else
   {
      SCIP_CALL( resolvePropagationFullOrbitope(scip, conshdlr, cons, inferinfo, bdchgidx, result) );
   }

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockOrbitope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int i;
   int j;
   int nspcons;
   int nblocks;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( locktype == SCIP_LOCKTYPE_MODEL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vars != NULL );

   SCIPdebugMsg(scip, "Locking method for orbitope constraint handler\n");

   nspcons = consdata->nspcons;
   nblocks = consdata->nblocks;
   vars = consdata->vars;

   /* add up locks and down locks on each variable */
   for (i = 0; i < nspcons; ++i)
   {
      for (j = 0; j < nblocks; ++j)
         SCIP_CALL( SCIPaddVarLocksType(scip, vars[i][j], locktype, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintOrbitope)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int i;
   int j;
   int nspcons;
   int nblocks;
   SCIP_ORBITOPETYPE orbitopetype;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nspcons > 0 );
   assert( consdata->nblocks > 0 );
   assert( consdata->vars != NULL );

   nspcons = consdata->nspcons;
   nblocks = consdata->nblocks;
   vars = consdata->vars;
   orbitopetype = consdata->orbitopetype;

   SCIPdebugMsg(scip, "Printing method for orbitope constraint handler\n");

   switch ( orbitopetype )
   {
   case SCIP_ORBITOPETYPE_PARTITIONING:
      SCIPinfoMessage(scip, file, "partOrbitope(");
      break;
   case SCIP_ORBITOPETYPE_PACKING:
      SCIPinfoMessage(scip, file, "packOrbitope(");
      break;
   case SCIP_ORBITOPETYPE_FULL:
      SCIPinfoMessage(scip, file, "fullOrbitope(");
      break;
   default:
      SCIPABORT();
   }

   for (i = 0; i < nspcons; ++i)
   {
      for (j = 0; j < nblocks; ++j)
      {
         if ( j > 0 )
            SCIPinfoMessage(scip, file, ",");
         SCIP_CALL( SCIPwriteVarName(scip, file, vars[i][j], TRUE) );
      }
      if ( i < nspcons-1 )
         SCIPinfoMessage(scip, file, ".");
   }
   SCIPinfoMessage(scip, file, ")");

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyOrbitope)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_VAR*** sourcevars;
   SCIP_VAR*** vars;
   int nspcons;
   int nblocks;
   int i;
   int k;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( sourcescip != NULL );
   assert( sourceconshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(sourceconshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( varmap != NULL );
   assert( valid != NULL );

   *valid = TRUE;

   SCIPdebugMsg(scip, "Copying method for orbitope constraint handler.\n");

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->nspcons > 0 );
   assert( sourcedata->nblocks > 0 );
   assert( sourcedata->vars != NULL );

   conshdlrdata = SCIPconshdlrGetData(sourceconshdlr);
   assert( conshdlrdata != NULL );

   /* do not copy non-model constraints */
   if ( !sourcedata->ismodelcons && !conshdlrdata->forceconscopy )
   {
      *valid = FALSE;

      return SCIP_OKAY;
   }

   nspcons = sourcedata->nspcons;
   nblocks = sourcedata->nblocks;
   sourcevars = sourcedata->vars;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nspcons) );
   for (i = 0; i < nspcons && *valid; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(vars[i]), nblocks) );  /*lint !e866*/

      for (j = 0; j < nblocks && *valid; ++j)
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[i][j], &(vars[i][j]), varmap, consmap, global, valid) );
         assert( !(*valid) || vars[i][j] != NULL );
      }
   }

   /* only create the target constraint, if all variables could be copied */
   if ( *valid )
   {
      /* create copied constraint */
      if ( name == NULL )
         name = SCIPconsGetName(sourcecons);

      SCIP_CALL( SCIPcreateConsOrbitope(scip, cons, name,
            vars, sourcedata->orbitopetype, nspcons, nblocks, sourcedata->usedynamicprop,
            sourcedata->resolveprop, sourcedata->ismodelcons, sourcedata->mayinteract,
            initial, separate, enforce, check, propagate,
            local, modifiable, dynamic, removable, stickingatnode) );
   }

   /* free space; only up to row i if copying failed */
   assert( 0 <= i && i <= nspcons );
   for (k = i - 1; k >= 0; --k)
      SCIPfreeBufferArray(scip, &vars[k]);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseOrbitope)
{  /*lint --e{715}*/
   const char* s;
   char* endptr;
   SCIP_ORBITOPETYPE orbitopetype;
   SCIP_VAR*** vars;
   SCIP_VAR* var;
   int nspcons;
   int maxnspcons;
   int nblocks;
   int maxnblocks;
   int k;
   int j;

   assert( success != NULL );

   *success = TRUE;
   s = str;

   /* skip white space */
   while ( *s != '\0' && isspace((unsigned char)*s) )
      ++s;

   if ( strncmp(s, "partOrbitope(", 13) == 0 )
      orbitopetype = SCIP_ORBITOPETYPE_PARTITIONING;
   else if ( strncmp(s, "packOrbitope(", 13) == 0 )
      orbitopetype = SCIP_ORBITOPETYPE_PACKING;
   else
   {
      if ( strncmp(s, "fullOrbitope(", 13) != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - expected \"fullOrbitope(\", \"partOrbitope\" or \"packOrbitope\": %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      orbitopetype = SCIP_ORBITOPETYPE_FULL;
   }
   s += 13;

   /* loop through string */
   nspcons = 0;
   nblocks = 0;
   maxnspcons = 10;
   maxnblocks = 10;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, maxnspcons) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(vars[0]), maxnblocks) );

   j = 0;
   do
   {
      /* skip whitespace */
      while ( isspace((int)*s) )
         ++s;

      /* parse variable name */
      SCIP_CALL( SCIPparseVarName(scip, s, &var, &endptr) );
      if ( var == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable name at '%s'\n", str);
         *success = FALSE;
         return SCIP_OKAY;
      }
      vars[nspcons][j++] = var;
      s = endptr;

      if ( j > nblocks )
      {
         if ( nspcons > 0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "variables per row do not match.\n");
            *success = FALSE;
            return SCIP_OKAY;
         }
         nblocks = j;

         if ( nblocks > maxnblocks )
         {
            int newsize;

            newsize = SCIPcalcMemGrowSize(scip, nblocks);
            SCIP_CALL( SCIPreallocBufferArray(scip, &(vars[nspcons]), newsize) );    /*lint !e866*/
            maxnblocks = newsize;
         }
      }
      assert( nblocks <= maxnblocks );

      /* skip white space and ',' */
      while ( *s != '\0' && ( isspace((unsigned char)*s) ||  *s == ',' ) )
         ++s;

      /* begin new row if required */
      if ( *s == '.' )
      {
         ++nspcons;
         ++s;

         if ( nspcons >= maxnspcons )
         {
            int newsize;

            newsize = SCIPcalcMemGrowSize(scip, nspcons+1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &vars, newsize) );
            maxnspcons = newsize;
         }
         assert(nspcons < maxnspcons);

         SCIP_CALL( SCIPallocBufferArray(scip, &(vars[nspcons]), nblocks) );  /*lint !e866*/
         j = 0;
      }
   }
   while ( *s != ')' );
   ++nspcons;

   /* to ensure consistency, we disable dynamic propagation and tell SCIP that the orbitope could potentially
    * interact with other symmetry handling constraints
    */
   SCIP_CALL( SCIPcreateConsOrbitope(scip, cons, name, vars, orbitopetype, nspcons, nblocks, FALSE, TRUE, TRUE, TRUE,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   for (k = nspcons - 1; k >= 0; --k)
      SCIPfreeBufferArray(scip, &vars[k]);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsOrbitope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( success != NULL );
   assert( vars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( varssize < consdata->nblocks * consdata->nspcons )
      (*success) = FALSE;
   else
   {
      int cnt = 0;
      int i;
      int j;

      for (i = 0; i < consdata->nspcons; ++i)
      {
         for (j = 0; j < consdata->nblocks; ++j)
            vars[cnt++] = consdata->vars[i][j];
      }
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsOrbitope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   (*nvars) = consdata->nblocks * consdata->nspcons;
   (*success) = TRUE;

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for orbitope constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrOrbitope(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create orbitope constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpOrbitope, consEnfopsOrbitope, consCheckOrbitope, consLockOrbitope,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyOrbitope, consCopyOrbitope) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeOrbitope) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteOrbitope) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsOrbitope) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsOrbitope) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseOrbitope) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolOrbitope, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintOrbitope) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropOrbitope, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropOrbitope) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpOrbitope, consSepasolOrbitope, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransOrbitope) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxOrbitope) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/checkpporbitope",
         "Strengthen orbitope constraints to packing/partioning orbitopes?",
         &conshdlrdata->checkpporbitope, TRUE, DEFAULT_PPORBITOPE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/sepafullorbitope",
         "Whether we separate inequalities for full orbitopes?",
         &conshdlrdata->sepafullorbitope, TRUE, DEFAULT_SEPAFULLORBITOPE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/forceconscopy",
         "Whether orbitope constraints should be forced to be copied to sub SCIPs.",
         &conshdlrdata->forceconscopy, TRUE, DEFAULT_FORCECONSCOPY, NULL, NULL) );

   return SCIP_OKAY;
}


/** creates and captures a orbitope constraint
 *
 *  @pre If packing/partitioning orbitopes are used, this constraint handler assumes that constraints which enforce
 *  the packing/partitioning constraints are contained in the problem. It does not implement, e.g., separation and
 *  propagation of set packing/partitioning constraints, since this would just copy large parts of the code of the
 *  setppc constraint handler.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR***           vars,               /**< matrix of variables on which the symmetry acts */
   SCIP_ORBITOPETYPE     orbitopetype,       /**< type of orbitope constraint */
   int                   nspcons,            /**< number of set partitioning/packing constraints  <=> p */
   int                   nblocks,            /**< number of symmetric variable blocks             <=> q */
   SCIP_Bool             usedynamicprop,     /**< whether dynamic propagation should be used */
   SCIP_Bool             mayinteract,        /**< whether symmetries corresponding to orbitope might interact
                                              *   with symmetries handled by other routines */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved? */
   SCIP_Bool             ismodelcons,        /**< whether the orbitope is a model constraint */
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
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the orbitope constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("orbitope constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* check for consistency */
   if ( usedynamicprop && mayinteract )
   {
      SCIPwarningMessage(scip, "Dynamic propagation is only possible if orbitope does not interact with \
                          other symmetry handling constraints. Ignore value of <usedynamicprop>.\n");
   }

   assert( nspcons > 0 );
   assert( nblocks > 0 );

   /* run some checks */
#ifndef NDEBUG
   {
      SCIP_Real obj;
      int i;
      int j;
      for (i = 0; i < nspcons; ++i)
      {
         /* initialize obj to infinity */
         obj = SCIPinfinity(scip);
         for (j = 0; j < nblocks; ++j)
         {
            SCIP_Bool fixedZero;
            SCIP_VAR* var;

            var = vars[i][j];
            assert(var != NULL);

            if ( SCIPvarIsNegated(var) )
               var = SCIPvarGetNegatedVar(var);

            /* all variables need to be binary */
            assert( SCIPvarIsBinary(var) );

            /* fixed variables have obj = 0; for variables fixed to 0, we assume that there is no
               problem (but we cannot always check it, e.g., when in the original problem
               variables were fixed and this problem was copied.) */
            fixedZero = ( SCIPisZero(scip, SCIPvarGetLbGlobal(var)) && SCIPisZero(scip, SCIPvarGetUbGlobal(var)) );

            /* @todo adapt correctness of the following check for sub-scips */
            if ( SCIPgetSubscipDepth(scip) == 0 )
            {
               /* check whether all variables in a row have the same objective */
               if ( ! fixedZero && SCIPisInfinity(scip, obj) )
                  obj = SCIPvarGetObj(var);
               else
               {
                  assert( fixedZero || ! SCIPvarIsActive(var) || SCIPisEQ(scip, obj, SCIPvarGetObj(var)) );
               }
            }
         }
      }
   }
#endif

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   if ( conshdlrdata->checkpporbitope && orbitopetype != SCIP_ORBITOPETYPE_PARTITIONING
      && orbitopetype != SCIP_ORBITOPETYPE_PACKING )
   {
      SCIP_CALL( strengthenOrbitopeConstraint(scip, vars, &nspcons, nblocks, &orbitopetype, mayinteract) );
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, vars, nspcons, nblocks, orbitopetype,
         resolveprop, usedynamicprop && ! mayinteract, ismodelcons, mayinteract) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures an orbitope constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR***           vars,               /**< matrix of variables on which the symmetry acts */
   SCIP_ORBITOPETYPE     orbitopetype,       /**< type of orbitope constraint */
   int                   nspcons,            /**< number of set partitioning/packing constraints  <=> p */
   int                   nblocks,            /**< number of symmetric variable blocks             <=> q */
   SCIP_Bool             usedynamicprop,     /**< whether dynamic propagation should be used */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved? */
   SCIP_Bool             ismodelcons,        /**< whether the orbitope is a model constraint */
   SCIP_Bool             mayinteract         /**< whether symmetries corresponding to orbitope might interact
                                              *   with symmetries handled by other routines */
   )
{
   SCIP_CALL( SCIPcreateConsOrbitope(scip, cons, name, vars, orbitopetype, nspcons, nblocks, usedynamicprop,
         resolveprop, ismodelcons, mayinteract, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
