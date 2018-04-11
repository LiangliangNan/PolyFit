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

/**@file   presol_symbreak.c
 * @brief  presolver for adding symmetry breaking constraints
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Christopher Hojny
 *
 * This presolver adds the following symmetry breaking constraints:
 *
 * - minimal cover inequalities for symresacks via a constraint handler
 *
 * @note It is important to control the order of presolvers within SCIP in order to avoid contraditions. First, one needs
 * to take care of presolvers that have an effect on symmetry, e.g., presol_domcol. If symmetry is computed on the
 * original formulation, we perform this presolver at the very beginning. Otherwise, we try to compute symmetry as late
 * as possible and then add constraints based on this information.
 *
 * @note Currently, we only allocate memory for pointers to symresack constraints for group generators. If further
 * constraints are considered, we have to reallocate memory.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include <scip/cons_linear.h>
#include <scip/cons_orbitope.h>
#include <scip/cons_setppc.h>
#include <scip/cons_symresack.h>
#include <scip/misc.h>
#include <scip/prop_probing.h>
#include <scip/presol_symbreak.h>
#include <scip/presol_symmetry.h>
#include <symmetry/type_symmetry.h>

/* presolver properties */
#define PRESOL_NAME            "symbreak"
#define PRESOL_DESC            "presolver for adding symmetry breaking constraints"
#define PRESOL_PRIORITY         -10000000    /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               -1    /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING   SCIP_PRESOLTIMING_EXHAUSTIVE   /**< timing for presolving */

/* default parameters */
#define DEFAULT_CONSSADDLP           TRUE    /**< Should the symmetry breaking constraints be added to the LP? */
#define DEFAULT_ADDSYMRESACKS        TRUE    /**< Add inequalities for symresacks for each generator? */
#define DEFAULT_COMPUTEORBITS       FALSE    /**< Should the orbits of the symmetry group be computed? */
#define DEFAULT_DETECTORBITOPES      TRUE    /**< Should we check whether the components of the symmetry group can be handled by orbitopes? */

/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   int**                 perms;              /**< list of permutations */
   int                   nperms;             /**< number of permutations in perms */
   int                   npermvars;          /**< number of variables affected by permutations */
   SCIP_Real             log10groupsize;     /**< log10 of group size */
   SCIP_VAR**            permvars;           /**< array of variables on which permutations act */
   SCIP_Bool             addedconss;         /**< whether we already added symmetry breaking constraints */
   SCIP_Bool             computedsymmetry;   /**< whether symmetry has been computed already */
   SCIP_Bool             conssaddlp;         /**< Should the symmetry breaking constraints be added to the LP? */
   SCIP_Bool             addsymresacks;      /**< Add symresack constraints for each generator? */
   SCIP_Bool             enabled;            /**< run presolver? */
   SCIP_Bool             early;              /**< run presolver as early as possible if symmetry has been detected in initpre() */
   SCIP_CONS**           genconss;           /**< list of generated constraints */
   int                   ngenconss;          /**< number of generated constraints */
   int                   nsymresacks;        /**< number of symresack constraints */
   SCIP_Bool             detectorbitopes;    /**< Should we check whether the components of the symmetry group can be handled by orbitopes? */
   int                   norbitopes;         /**< number of orbitope constraints */
   int                   norbits;            /**< number of non-trivial orbits of permutation group */
   SCIP_Bool             computeorbits;      /**< whether the orbits of the symmetry group should be computed */
   int*                  orbits;             /**< array containing the indices of variables sorted by orbits */
   int*                  orbitbegins;        /**< array containing in i-th position the first position of orbit i in orbits array */
   int                   ncomponents;        /**< number of components of symmetry group */
   int*                  npermsincomponent;  /**< array containing number of permutations per component */
   int**                 components;         /**< array containing for each components the corresponding permutations */
   SCIP_Shortbool*       componentblocked;   /**< array to store whether a component is blocked to be considered by symmetry handling techniques */
};



/*
 * Local methods
 */


/** compute components of symmetry group */
static
SCIP_RETCODE computeComponents(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOLDATA*      presoldata          /**< data of symmetry breaking presolver */
   )
{
   SCIP_DISJOINTSET* componentstoperm = NULL;
   SCIP_Bool startchanged;
   int** perms;
   int npermsincomponent;
   int curcomponent;
   int ncomponents;
   int componentcnt;
   int npermvars;
   int nperms;
   int newstart;
   int start;
   int i;
   int j;
   int k;

   assert( scip != NULL );
   assert( presoldata != NULL );

   nperms = presoldata->nperms;

   presoldata->ncomponents = -1;

   if ( nperms <= 0 )
      return SCIP_OKAY;

   /* at this point, we have at least one non-trivial permutation */
   npermvars = presoldata->npermvars;
   perms = presoldata->perms;

   /* if there exists only one orbit, we have exactly one component */
   if ( presoldata->norbits == 1 )
      ncomponents = 1;
   else
   {
      /* init array that assigns to each permutation its component of the group */
      SCIP_CALL( SCIPdisjointsetCreate(&componentstoperm, SCIPblkmem(scip), nperms) );
      ncomponents = nperms;

      /* check whether two permutations belong to the same component */
      for (i = 0; i < nperms && ncomponents > 1; ++i)
      {
         for (j = i + 1; j < nperms && ncomponents > 1; ++j)
         {
            int componentI;
            int componentJ;
            int* permI;
            int* permJ;

            componentI = SCIPdisjointsetFind(componentstoperm, i);
            componentJ = SCIPdisjointsetFind(componentstoperm, j);
            if ( componentI == componentJ )
               continue;

            permI = perms[i];
            permJ = perms[j];

            /* Do perms[i] and perms[j] belong to the same component? */
            for (k = 0; k < npermvars; ++k)
            {
               /* both permutations belong to the same component */
               if ( permI[k] != k && permJ[k] != k )
               {
                  /* Keep the smallest identifier to keep track of where a new component starts.
                   * Note that it is necessary to store the smaller identifier to be able to iterate
                   * over all ordered pairs (i,j), i < j, of permutations, instead of all unordered
                   * pairs {i,j}. */
                  if ( componentI < componentJ )
                     SCIPdisjointsetUnion(componentstoperm, i, j, TRUE);
                  else
                     SCIPdisjointsetUnion(componentstoperm, j, i, TRUE);

                  --ncomponents;
                  break;
               }
            }
         }
      }
      assert( ncomponents > 0 );
   }

   /* store data */
   presoldata->ncomponents = ncomponents;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->npermsincomponent, ncomponents) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->components, ncomponents) );

   /* copy componentstoperm adequatly to components and npermsincomponent */
   if ( ncomponents == 1 )
   {
      presoldata->npermsincomponent[0] = nperms;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata->components[0]), nperms) );
      for (i = 0; i < nperms; ++i)
         presoldata->components[0][i] = i;
   }
   else
   {
      componentcnt = 0;
      start = 0;
      newstart = 0;
      startchanged = TRUE;
      while ( startchanged )
      {
         startchanged = FALSE;
         npermsincomponent = 0;
         assert( componentstoperm != NULL );
         curcomponent = SCIPdisjointsetFind(componentstoperm, start);

         /* find number of permutations in current component and detect first perm in another permutation */
         for (i = start; i < nperms; ++i)
         {
            if ( SCIPdisjointsetFind(componentstoperm, i) == curcomponent )
               ++npermsincomponent;
            /* only store other component if it has the same identifier as its node (mark that new component starts) */
            else if ( ! startchanged && SCIPdisjointsetFind(componentstoperm, i) == i )
            {
               newstart = i;
               startchanged = TRUE;
            }
         }

         /* store number of permutations and permutation per component */
         presoldata->npermsincomponent[componentcnt] = npermsincomponent;

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata->components[componentcnt]), npermsincomponent) );

         j = 0;
         for (i = start; i < nperms; ++i)
         {
            if ( SCIPdisjointsetFind(componentstoperm, i) == curcomponent )
               presoldata->components[componentcnt][j++] = i;
         }

         ++componentcnt;
         start = newstart;
      }
   }

   if ( presoldata->norbits != 1 )
      SCIPdisjointsetFree(&componentstoperm, SCIPblkmem(scip));

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata->componentblocked), ncomponents) );
   for (i = 0; i < ncomponents; ++i)
      presoldata->componentblocked[i] = FALSE;

#ifdef SCIP_OUTPUT
   printf("\n\n\nTESTS\n\n");
   printf("Number of components:\t\t%d\n", presoldata->ncomponents);
   for (i = 0; i < presoldata->ncomponents; ++i)
   {
      printf("Component %d: Number of perms: %d\n", i, presoldata->npermsincomponent[i]);
      for (j = 0; j < presoldata->npermsincomponent[i]; ++j)
      {
         printf("%d ", presoldata->components[i][j]);
      }
      printf("\n");
   }

   printf("GENERATORS:");
   for (i = 0; i < presoldata->nperms; ++i)
   {
      printf("generator %d:\n", i);
      for (j = 0; j < presoldata->npermvars; ++j)
      {
         if ( presoldata->perms[i][j] != j )
            printf("%d ", presoldata->perms[i][j]);
      }
      printf("\n");
   }
#endif

   return SCIP_OKAY;
}


/** check whether a permutation is a composition of 2-cycles of binary variables and in this case determines the number of 2-cycles */
static
SCIP_RETCODE getPermProperties(
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< array of variables perm is acting on */
   int                   nvars,              /**< number of variables */
   SCIP_Bool*            iscompoftwocycles,  /**< pointer to store whether permutation is a composition of 2-cycles */
   int*                  ntwocyclesperm,     /**< pointer to store number of 2-cycles */
   SCIP_Bool*            allvarsbinary       /**< pointer to store whether perm is acting on binary variables only */
   )
{
   int ntwocycles = 0;
   int i;

   assert( perm != NULL );
   assert( iscompoftwocycles != NULL );
   assert( ntwocyclesperm != NULL );
   assert( allvarsbinary != NULL );

   *iscompoftwocycles = FALSE;
   *ntwocyclesperm = 0;
   *allvarsbinary = TRUE;
   for (i = 0; i < nvars; ++i)
   {
      /* skip fixed points and avoid treating the same 2-cycle twice */
      if ( perm[i] <= i )
         continue;

      if ( perm[perm[i]] == i )
      {
         if ( SCIPvarIsBinary(vars[i]) && SCIPvarIsBinary(vars[perm[i]]) )
            ++ntwocycles;
         else
         {
            /* at least one variable is not binary */
            *allvarsbinary = FALSE;
            return SCIP_OKAY;
         }
      }
      else
      {
         /* we do not have a 2-cycle */
         return SCIP_OKAY;
      }
   }

   /* at this point the permutation is a composition of 2-cycles on binary variables */
   *iscompoftwocycles = TRUE;
   *ntwocyclesperm = ntwocycles;

   return SCIP_OKAY;
}


/** Given a matrix with nrows and \#perms + 1 columns whose first nfilledcols columns contain entries of variables, this routine
 *  checks whether the 2-cycles of perm intersect each row of column coltoextend in exactly one position. In this case,
 *  we add one column to the suborbitope of the first nfilledcols columns.
 *
 *  @pre Every non-trivial cycle of perm is a 2-cycle.
 *  @pre perm has nrows many 2-cycles
 */
static
SCIP_RETCODE extendSubOrbitope(
   int**                 suborbitope,        /**< matrix containing suborbitope entries */
   int                   nrows,              /**< number of rows of suborbitope */
   int                   nfilledcols,        /**< number of columns of suborbitope which are filled with entries */
   int                   coltoextend,        /**< index of column that should be extended by perm */
   int*                  perm,               /**< permutation */
   SCIP_Bool             leftextension,      /**< whether we extend the suborbitope to the left */
   int**                 nusedelems,         /**< pointer to array storing how often an element was used in the orbitope */
   SCIP_Bool*            success,            /**< pointer to store whether extension was successful */
   SCIP_Bool*            infeasible          /**< pointer to store if the number of intersecting cycles is too small */
   )
{
   int nintersections = 0;
   int row;
   int idx1;
   int idx2;

   assert( suborbitope != NULL );
   assert( nfilledcols > 0 );
   assert( perm != NULL );
   assert( nusedelems != NULL );
   assert( success != NULL );
   assert( infeasible != NULL );

   *success = FALSE;
   *infeasible = FALSE;

   /* if we try to extend the first orbitope generator by another one,
    * we can change the order of entries in each row of suborbitope */
   if ( nfilledcols == 2 )
   {
      /* check whether each cycle of perm intersects with a row of suborbitope */
      for (row = 0; row < nrows; ++row)
      {
         idx1 = suborbitope[row][0];
         idx2 = suborbitope[row][1];

         /* if idx1 or idx2 is affected by perm, we can extend the row of the orbitope */
         if ( idx1 != perm[idx1] )
         {
            /* change order of idx1 and idx2 for right extensions */
            if ( ! leftextension )
            {
               suborbitope[row][0] = idx2;
               suborbitope[row][1] = idx1;
            }
            suborbitope[row][2] = perm[idx1];
            ++nintersections;

            (*nusedelems)[idx1] += 1;
            (*nusedelems)[perm[idx1]] += 1;

            /* if an element appears too often in the orbitope matrix */
            if ( (*nusedelems)[idx1] > 2 || (*nusedelems)[perm[idx1]] > 2 )
            {
               *infeasible = TRUE;
               break;
            }
         }
         else if ( idx2 != perm[idx2] )
         {
            /* change order of idx1 and idx2 for left extensions */
            if ( leftextension )
            {
               suborbitope[row][0] = idx2;
               suborbitope[row][1] = idx1;
            }
            suborbitope[row][2] = perm[idx2];
            ++nintersections;

            (*nusedelems)[idx2] += 1;
            (*nusedelems)[perm[idx2]] += 1;

            /* if an element appears too often in the orbitope matrix */
            if ( (*nusedelems)[idx2] > 2 || (*nusedelems)[perm[idx2]] > 2 )
            {
               *infeasible = TRUE;
               break;
            }
         }
      }
   }
   else
   {
      /* check whether each cycle of perm intersects with a row of suborbitope */
      for (row = 0; row < nrows; ++row)
      {
         idx1 = suborbitope[row][coltoextend];

         /* if idx1 is affected by perm, we can extend the row of the orbitope */
         if ( idx1 != perm[idx1] )
         {
            suborbitope[row][nfilledcols] = perm[idx1];
            ++nintersections;

            (*nusedelems)[idx1] += 1;
            (*nusedelems)[perm[idx1]] += 1;

            /* if an element appears to often in the orbitope matrix */
            if ( (*nusedelems)[idx1] > 2 || (*nusedelems)[perm[idx1]] > 2 )
            {
               *infeasible = TRUE;
               break;
            }

         }
      }
   }

   /* if there are too few intersection, this is not an orbitope */
   if ( nintersections > 0 && nintersections < nrows )
      *infeasible = TRUE;
   else if ( nintersections == nrows )
      *success = TRUE;

   return SCIP_OKAY;
}


/** generate variable matrix for orbitope constraint handler */
static
SCIP_RETCODE generateOrbitopeVarsMatrix(
   SCIP_VAR****          vars,               /**< pointer to matrix of orbitope variables */
   int                   nrows,              /**< number of rows of orbitope */
   int                   ncols,              /**< number of columns of orbitope */
   SCIP_VAR**            permvars,           /**< superset of variables that are contained in orbitope */
   int                   npermvars,          /**< number of variables in permvars array */
   int**                 orbitopevaridx,     /**< permuted index table of variables in permvars that are contained in orbitope */
   int*                  columnorder,        /**< permutation to reorder column of orbitopevaridx */
   int*                  nusedelems,         /**< array storing how often an element was used in the orbitope */
   SCIP_Bool*            infeasible          /**< pointer to store whether the potential orbitope is not an orbitope */
   )
{
   int nfilledcols = 0;
   int curcolumn;
   int i;

   assert( vars != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( orbitopevaridx != NULL );
   assert( columnorder != NULL );
   assert( nusedelems != NULL );
   assert( infeasible != NULL );

   curcolumn = ncols - 1;

   /* start filling vars matrix with the right-most column w.r.t. columnorder */
   while ( curcolumn >= 0 && columnorder[curcolumn] >= 0 )
   {
      for (i = 0; i < nrows; ++i)
      {
         assert( orbitopevaridx[i][curcolumn] < npermvars );
         assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][curcolumn]]) );

         /* elements in first column of orbitope have to appear exactly once in the orbitope */
         if ( nusedelems[orbitopevaridx[i][curcolumn]] > 1 )
         {
            *infeasible = TRUE;
            break;
         }

         (*vars)[i][nfilledcols] = permvars[orbitopevaridx[i][curcolumn]];
      }
      --curcolumn;
      ++nfilledcols;
   }

   /* There are three possibilities for the structure of columnorder:
    * 1)  [0, 1, -1, -1, ..., -1]
    * 2)  [0, 1, 1, 1, ..., 1]
    * 3)  [0, 1, -1, -1, ...., -1, 1, 1, ..., 1]
    */
   /* Either we are in case 1) or case 3), or all columns should have been added to vars in case 2) */
   assert( curcolumn > 1 || (curcolumn < 0 && nfilledcols == ncols) );

   if ( curcolumn > 1 )
   {
      /* add column with columnorder 1 to vars */
      for (i = 0; i < nrows; ++i)
      {
         assert( orbitopevaridx[i][1] < npermvars );
         assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][1]]) );

         (*vars)[i][nfilledcols] = permvars[orbitopevaridx[i][1]];
      }
      ++nfilledcols;

      /* add column with columnorder 0 to vars */
      for (i = 0; i < nrows; ++i)
      {
         assert( orbitopevaridx[i][0] < npermvars );
         assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][0]]) );

         (*vars)[i][nfilledcols] = permvars[orbitopevaridx[i][0]];
      }
      ++nfilledcols;

      /* add columns with a negative column order to vars */
      if ( nfilledcols < ncols )
      {
         assert( ncols > 2 );

         curcolumn = 2;
         while ( nfilledcols < ncols )
         {
            assert( columnorder[curcolumn] < 0 );

            for (i = 0; i < nrows; ++i)
            {
               assert( orbitopevaridx[i][curcolumn] < npermvars );
               assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][curcolumn]]) );

               /* elements in last column of orbitope have to appear exactly once in the orbitope */
               if ( nfilledcols == ncols - 1 && nusedelems[orbitopevaridx[i][curcolumn]] > 1 )
               {
                  *infeasible = TRUE;
                  break;
               }

               (*vars)[i][nfilledcols] = permvars[orbitopevaridx[i][curcolumn]];
            }
            ++curcolumn;
            ++nfilledcols;
         }
      }
   }

   return SCIP_OKAY;
}


/** check whether components of the symmetry group can be completely handled by orbitopes */
static
SCIP_RETCODE detectOrbitopes(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOLDATA*      presoldata          /**< pointer to data of symbreak presolver */
   )
{
   SCIP_VAR** permvars;
   int** components;
   int** perms;
   int* npermsincomponent;
   int ncomponents;
   int npermvars;
   int i;

   assert( scip != NULL );
   assert( presoldata != NULL );

   /* exit if no symmetry is present */
   if ( presoldata->nperms == 0 )
      return SCIP_OKAY;

   /* compute components if not done yet */
   if ( presoldata->ncomponents == -1 )
   {
      SCIP_CALL( computeComponents(scip, presoldata) );
   }
   assert( presoldata->nperms > 0 );
   assert( presoldata->perms != NULL );
   assert( presoldata->npermvars > 0 );
   assert( presoldata->permvars != NULL );
   assert( presoldata->ncomponents > 0 );
   assert( presoldata->components != NULL );
   assert( presoldata->npermsincomponent != NULL );

   perms = presoldata->perms;
   npermvars = presoldata->npermvars;
   permvars = presoldata->permvars;
   ncomponents = presoldata->ncomponents;
   npermsincomponent = presoldata->npermsincomponent;
   components = presoldata->components;

   /* iterate over components */
   for (i = 0; i < ncomponents; ++i)
   {
      SCIP_VAR*** vars;
      SCIP_CONS* cons;
      SCIP_Bool isorbitope = TRUE;
      SCIP_Bool* usedperm;
      int** orbitopevaridx;
      int* columnorder;
      int ntwocyclescomp = -1;
      int nfilledcols;
      int nusedperms;
      int* nusedelems;
      int coltoextend;
      int j;
      int row;
      SCIP_Bool infeasibleorbitope;

      /* get properties of permutations */
      assert( npermsincomponent[i] > 0 );
      for (j = 0; j < npermsincomponent[i]; ++j)
      {
         SCIP_Bool iscompoftwocycles = FALSE;
         SCIP_Bool allvarsbinary = TRUE;
         int ntwocyclesperm = 0;

         SCIP_CALL( getPermProperties(perms[components[i][j]], permvars, npermvars, &iscompoftwocycles, &ntwocyclesperm, &allvarsbinary) );

         /* if we are checking the first permutation */
         if ( ntwocyclescomp == - 1 )
            ntwocyclescomp = ntwocyclesperm;

         /* no or different number of 2-cycles or not all vars binary: permutations cannot generate orbitope */
         if ( ntwocyclescomp == 0 || ntwocyclescomp != ntwocyclesperm || ! allvarsbinary )
         {
            isorbitope = FALSE;
            break;
         }
      }

      /* if no orbitope was detected or orbitope has too few rows */
      if ( ! isorbitope )
         continue;
      assert( ntwocyclescomp > 0 );

      /* iterate over permutations and check whether for each permutation there exists
       * another permutation whose 2-cycles intersect pairwise in exactly one element */

      /* whether a permutation was considered to contribute to orbitope */
      SCIP_CALL( SCIPallocBufferArray(scip, &usedperm, npermsincomponent[i]) );
      for (j = 0; j < npermsincomponent[i]; ++j)
         usedperm[j] = FALSE;
      nusedperms = 0;

      /* orbitope matrix for indices of variables in permvars array */
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx, ntwocyclescomp) );
      for (j = 0; j < ntwocyclescomp; ++j)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx[j], npermsincomponent[i] + 1) ); /*lint !e866*/
      }

      /* order of columns of orbitopevaridx */
      SCIP_CALL( SCIPallocBufferArray(scip, &columnorder, npermsincomponent[i] + 1) );
      for (j = 0; j < npermsincomponent[i] + 1; ++j)
         columnorder[j] = npermsincomponent[i] + 2;

      /* count how often an element was used in the potential orbitope */
      SCIP_CALL( SCIPallocBufferArray(scip, &nusedelems, npermvars) );
      for (j = 0; j < npermvars; ++j)
         nusedelems[j] = 0;

      /* fill first two columns of orbitopevaridx matrix */
      row = 0;
      for (j = 0; j < npermvars; ++j)
      {
         int permidx = components[i][0];

         /* avoid adding the same 2-cycle twice */
         if ( perms[permidx][j] > j )
         {
            orbitopevaridx[row][0] = j;
            orbitopevaridx[row++][1] = perms[permidx][j];
            nusedelems[j] += 1;
            nusedelems[perms[permidx][j]] += 1;
         }

         if ( row == ntwocyclescomp )
            break;
      }
      assert( row == ntwocyclescomp );

      usedperm[0] = TRUE;
      ++nusedperms;
      columnorder[0] = 0;
      columnorder[1] = 1;
      nfilledcols = 2;

      /* extend orbitopevaridx matrix to the left, i.e., iteratively find new permutations that
       * intersect the last added left column in each row in exactly one entry, starting with
       * column 0 */
      coltoextend = 0;
      for (j = 0; j < npermsincomponent[i]; ++j)
      {  /* lint --e{850} */
         SCIP_Bool success = FALSE;
         SCIP_Bool infeasible = FALSE;

         if ( nusedperms == npermsincomponent[i] )
            break;

         if ( usedperm[j] )
            continue;

         SCIP_CALL( extendSubOrbitope(orbitopevaridx, ntwocyclescomp, nfilledcols, coltoextend,
               perms[components[i][j]], TRUE, &nusedelems, &success, &infeasible) );

         if ( infeasible )
         {
            isorbitope = FALSE;
            break;
         }
         else if ( success )
         {
            usedperm[j] = TRUE;
            ++nusedperms;
            coltoextend = nfilledcols;
            columnorder[nfilledcols++] = -1; /* mark column to be filled from the left */
            j = 0; /* reset j since previous permutations can now intersect with the latest added column */
         }
      }

      if ( ! isorbitope ) /*lint !e850*/
         goto FREEDATASTRUCTURES;

      coltoextend = 1;
      for (j = 0; j < npermsincomponent[i]; ++j)
      {  /*lint --e(850)*/
         SCIP_Bool success = FALSE;
         SCIP_Bool infeasible = FALSE;

         if ( nusedperms == npermsincomponent[i] )
            break;

         if ( usedperm[j] )
            continue;

         SCIP_CALL( extendSubOrbitope(orbitopevaridx, ntwocyclescomp, nfilledcols, coltoextend,
               perms[components[i][j]], FALSE, &nusedelems, &success, &infeasible) );

         if ( infeasible )
         {
            isorbitope = FALSE;
            break;
         }
         else if ( success )
         {
            usedperm[j] = TRUE;
            ++nusedperms;
            coltoextend = nfilledcols;
            columnorder[nfilledcols] = 1; /* mark column to be filled from the right */
            ++nfilledcols;
            j = 0; /* reset j since previous permutations can now intersect with the latest added column */
         }
      }

      if ( nusedperms < npermsincomponent[i] ) /*lint !e850*/
         isorbitope = FALSE;

      if ( ! isorbitope )
         goto FREEDATASTRUCTURES;

      /* we have found a potential orbitope, prepare data for orbitope conshdlr */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, ntwocyclescomp) );
      for (j = 0; j < ntwocyclescomp; ++j)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &vars[j], npermsincomponent[i] + 1) ); /*lint !e866*/
      }

      /* prepare variable matrix (reorder columns of orbitopevaridx) */
      infeasibleorbitope = FALSE;
      SCIP_CALL( generateOrbitopeVarsMatrix(&vars, ntwocyclescomp, npermsincomponent[i] + 1, permvars, npermvars,
            orbitopevaridx, columnorder, nusedelems, &infeasibleorbitope) );

      if ( ! infeasibleorbitope )
      {
         SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, "orbitope", vars, SCIP_ORBITOPETYPE_FULL, ntwocyclescomp, npermsincomponent[i] + 1, TRUE,
               presoldata->conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         /* do not release constraint here - will be done later */
         presoldata->genconss[presoldata->ngenconss] = cons;
         ++presoldata->ngenconss;
         ++presoldata->norbitopes;
         presoldata->componentblocked[i] = TRUE;
      }

      /* free data structures */
      for (j = 0; j < ntwocyclescomp; ++j)
         SCIPfreeBufferArray(scip, &vars[j]);
      SCIPfreeBufferArray(scip, &vars);

   FREEDATASTRUCTURES:
      SCIPfreeBufferArray(scip, &nusedelems);
      SCIPfreeBufferArray(scip, &columnorder);
      for (j = 0; j < ntwocyclescomp; ++j)
         SCIPfreeBufferArray(scip, &orbitopevaridx[j]);
      SCIPfreeBufferArray(scip, &orbitopevaridx);
      SCIPfreeBufferArray(scip, &usedperm);
   }

   return SCIP_OKAY;
}


/** add symresack constraints */
static
SCIP_RETCODE addSymresackConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOL*          presol              /**< symmetry breaking presolver */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_VAR** permvars;
   SCIP_Bool conssaddlp;
   int** perms;
   int nperms;
   int nsymresackcons = 0;
   int npermvars;
   int ncomponents;
   int i;
   int p;

   assert( scip != NULL );
   assert( presol != NULL );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   perms = presoldata->perms;
   nperms = presoldata->nperms;
   permvars = presoldata->permvars;
   npermvars = presoldata->npermvars;
   conssaddlp = presoldata->conssaddlp;
   ncomponents = presoldata->ncomponents;

   assert( nperms <= 0 || perms != NULL );
   assert( permvars != NULL );
   assert( npermvars > 0 );

   /* if we use different approaches for components of symmetry group */
   if ( ncomponents > 0 )
   {
      /* loop through components */
      for (i = 0; i < ncomponents; ++i)
      {
         /* skip components that were treated by different symemtry handling techniques */
         if ( presoldata->componentblocked[i] )
            continue;

         /* loop through perms in component i and add symresack constraints */
         for (p = 0; p < presoldata->npermsincomponent[i]; ++p)
         {
            SCIP_CONS* cons;
            int permidx = presoldata->components[i][p];
            char name[SCIP_MAXSTRLEN];

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symbreakcons_component%d_perm%d", i, p);
            SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name, perms[permidx], permvars, npermvars,
                  conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );

            /* do not release constraint here - will be done later */
            presoldata->genconss[presoldata->ngenconss] = cons;
            ++presoldata->ngenconss;
            ++presoldata->nsymresacks;
            ++nsymresackcons;
            SCIPdebugMsg(scip, "Added symresack constraint: %d.\n", nsymresackcons);
         }
      }
   }
   else
   {
      /* loop through perms and add symresack constraints */
      for (p = 0; p < nperms; ++p)
      {
         SCIP_CONS* cons;
         char name[SCIP_MAXSTRLEN];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symbreakcons_perm%d", p);
         SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name, perms[p], permvars, npermvars,
                  conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         /* do not release constraint here - will be done later */
         presoldata->genconss[presoldata->ngenconss] = cons;
         ++presoldata->ngenconss;
         ++presoldata->nsymresacks;
         ++nsymresackcons;
         SCIPdebugMsg(scip, "Added symresack constraint: %d.\n", nsymresackcons);
      }
   }

   return SCIP_OKAY;
}


/** analyze generators and add symmetry breaking constraints */
static
SCIP_RETCODE addSymmetryBreakingConstraints(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   /* exit if no or only trivial symmetry group is available */
   if ( presoldata->nperms < 1 )
      return SCIP_OKAY;

   if ( presoldata->addsymresacks )
   {
      SCIP_CALL( addSymresackConss(scip, presol) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeSymbreak)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   SCIPfreeMemory(scip, &presoldata);

   return SCIP_OKAY;
}


/** initialization method of presolver (called after problem was transformed) */
static
SCIP_DECL_PRESOLINIT(presolInitSymbreak)
{
   SCIP_PRESOLDATA* presoldata;
   int usesymmetry;

   assert( scip != NULL );
   assert( presol != NULL );

   SCIPdebugMsg(scip, "Init method of symmetry breaking presolver ...\n");

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   /* check whether we should run */
   SCIP_CALL( SCIPgetIntParam(scip, "misc/usesymmetry", &usesymmetry) );
   if ( usesymmetry == (int) SYM_HANDLETYPE_SYMBREAK )
      presoldata->enabled = TRUE;
   else
   {
      presoldata->enabled = FALSE;
   }

   if ( presoldata->enabled )
   {
      /* register presolver for symmetry information */
      SCIP_CALL( SCIPregisterSymmetry(scip, SYM_HANDLETYPE_SYMBREAK, SYM_SPEC_BINARY | SYM_SPEC_INTEGER | SYM_SPEC_REAL, 0) );
   }

   return SCIP_OKAY;
}


/** deinitialization method of presolver (called before transformed problem is freed) */
static
SCIP_DECL_PRESOLEXIT(presolExitSymbreak)
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_CONS** genconss;
   int ngenconss;
   int i;

   assert( scip != NULL );
   assert( presol != NULL );

   SCIPdebugMsg(scip, "Exit method of symmetry breaking presolver ...\n");

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   ngenconss = presoldata->ngenconss;

   /* release constraints */
   genconss = presoldata->genconss;
   for (i = 0; i < ngenconss; ++i)
   {
      assert( genconss[i] != NULL );
      SCIP_CALL( SCIPreleaseCons(scip, &genconss[i]) );
   }

   /* reset data (necessary if another problem is read in the SCIP shell) */

   /* free pointers to symmetry group and binary variables */
   SCIPfreeBlockMemoryArrayNull(scip, &presoldata->genconss, presoldata->nperms);

   presoldata->genconss = NULL;
   presoldata->ngenconss = 0;

   /* free orbit structures */
   if ( presoldata->norbits >= 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &presoldata->orbitbegins, presoldata->npermvars);
      SCIPfreeBlockMemoryArray(scip, &presoldata->orbits, presoldata->npermvars);
   }

   presoldata->orbitbegins = NULL;
   presoldata->orbits = NULL;
   presoldata->norbits = -1;

   /* free components */
   if ( presoldata->ncomponents > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &presoldata->componentblocked, presoldata->ncomponents);
      for (i = 0; i < presoldata->ncomponents; ++i)
         SCIPfreeBlockMemoryArray(scip, &presoldata->components[i], presoldata->npermsincomponent[i]);
      SCIPfreeBlockMemoryArray(scip, &presoldata->components, presoldata->ncomponents);
      SCIPfreeBlockMemoryArray(scip, &presoldata->npermsincomponent, presoldata->ncomponents);
   }

   presoldata->componentblocked = NULL;
   presoldata->components = NULL;
   presoldata->npermsincomponent = NULL;
   presoldata->ncomponents = -1;

   /* reset basic data */
   presoldata->addedconss = FALSE;
   presoldata->computedsymmetry = FALSE;
   presoldata->enabled = TRUE;
   presoldata->early = FALSE;
   presoldata->nperms = -1;
   presoldata->log10groupsize = -1.0;
   presoldata->norbitopes = 0;
   presoldata->nsymresacks = 0;

   return SCIP_OKAY;
}


/** presolving initialization method of presolver (called when presolving is about to begin) */
static
SCIP_DECL_PRESOLINITPRE(presolInitpreSymbreak)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   int usesymmetry;

   assert( scip != NULL );
   assert( presol != NULL );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   /* check whether we have to run the presolver at the beginning of presolving */
   SCIP_CALL( SCIPgetTimingSymmetry(scip, &presoldata->early) );
   presoldata->early = ! presoldata->early;

   /* check whether we should run */
   SCIP_CALL( SCIPgetIntParam(scip, "misc/usesymmetry", &usesymmetry) );
   if ( usesymmetry == (int) SYM_HANDLETYPE_SYMBREAK )
      presoldata->enabled = TRUE;
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "presolving/symbreak/maxrounds", 0) );
      presoldata->enabled = FALSE;
   }

   if ( presoldata->enabled )
   {
      /* register presolver for symmetry information */
      SCIP_CALL( SCIPregisterSymmetry(scip, SYM_HANDLETYPE_SYMBREAK, SYM_SPEC_BINARY | SYM_SPEC_INTEGER | SYM_SPEC_REAL, 0) );
   }

   if ( presoldata->early && presoldata->enabled )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Executing presolver <%s> early, since symmetries are computed early.\n\n", SCIPpresolGetName(presol));

      SCIP_CALL( SCIPsetIntParam(scip, "presolving/" PRESOL_NAME "/priority", 90000000) );
   }

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecSymbreak)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   int noldfixedvars = *nfixedvars;
   int noldaggrvars = *naggrvars;
   int noldbdchgs = *nchgbds;
   int noldaddconss = *naddconss;
   int i;

   assert( scip != NULL );
   assert( presol != NULL );
   assert( result != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING );

   *result = SCIP_DIDNOTRUN;

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   /* possibly skip presolver */
   if ( ! presoldata->enabled )
      return SCIP_OKAY;

   /* skip presolving if we are not at the end */
   if ( ! presoldata->early && ! SCIPisPresolveFinished(scip) )
      return SCIP_OKAY;

   /* possibly stop */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Presolving method of symmetry breaking presolver ...\n");

   /* get symmetry information, if not already computed */
   if ( ! presoldata->computedsymmetry )
   {
      SCIPdebugMsg(scip, "Symmetry breaking presolver: computing symmetry ...\n");
      assert( presoldata->nperms < 0 );

      /* get symmetries */
      SCIP_CALL( SCIPgetGeneratorsSymmetry(scip, &(presoldata->npermvars), &(presoldata->permvars),
            &(presoldata->nperms), &(presoldata->perms), &(presoldata->log10groupsize)) );

      presoldata->computedsymmetry = TRUE;

      if ( presoldata->nperms <= 0 )
      {
         SCIPdebugMsg(scip, "Symmetry breaking presolver: no symmetry has been found, turning presolver off.\n");
         presoldata->enabled = FALSE;
         return SCIP_OKAY;
      }
      else
      {
         assert( presoldata->nperms > 0 );

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata->genconss), presoldata->nperms) );

         if ( presoldata->computeorbits )
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->orbits, presoldata->npermvars) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->orbitbegins, presoldata->npermvars) );

            SCIP_CALL( SCIPcomputeGroupOrbitsSymbreak(scip, presoldata->permvars, presoldata->npermvars, presoldata->perms, presoldata->nperms, NULL,
                  presoldata->orbits, presoldata->orbitbegins, &presoldata->norbits) );
         }

         if ( presoldata->detectorbitopes )
         {
            SCIP_CALL( detectOrbitopes(scip, presoldata) );
         }
      }
   }
   else if ( presoldata->nperms <= 0 )
      return SCIP_OKAY;

   /* at this point, the symmetry group should be computed and nontrivial */
   assert( presoldata->nperms > 0 );
   *result = SCIP_DIDNOTFIND;

   /* deactivate presolvers that may conflict with symmetry handling routines */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/domcol/maxrounds", 0) );

   /* possibly stop */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* if not already done, add symmetry breaking constraints */
   if ( ! presoldata->addedconss )
   {
      /* add symmetry breaking constraints */
      int noldngenconns = presoldata->ngenconss;

      SCIP_CALL( addSymmetryBreakingConstraints(scip, presol) );

      presoldata->addedconss = TRUE;
      *naddconss += presoldata->ngenconss - noldngenconns;
      SCIPdebugMsg(scip, "Added symmetry breaking constraints: %d.\n", presoldata->ngenconss - noldngenconns);

      /* if constraints have been added, loop through generated constraints and presolve each */
      for (i = 0; i < presoldata->ngenconss; ++i)
      {
         SCIP_CALL( SCIPpresolCons(scip, presoldata->genconss[i], nrounds, SCIP_PRESOLTIMING_ALWAYS, nnewfixedvars, nnewaggrvars, nnewchgvartypes,
               nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars,
               nchgvartypes, nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

         /* exit if cutoff has been detected */
         if ( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
         {
            SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n", SCIPconsGetName(presoldata->genconss[i]));
            return SCIP_OKAY;
         }
      }
      SCIPdebugMsg(scip, "Presolved %d constraints generated by symbreak.\n", presoldata->ngenconss);
   }

   /* determine success */
   if ( *naddconss + *nchgbds + *naggrvars + *nfixedvars > noldaddconss + noldbdchgs + noldaggrvars + noldfixedvars )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the symbreak presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolSymbreak(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata = NULL;
   SCIP_PRESOL* presol;

   /* create presolver data */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );

   /* we must call the constructor explictly, because memory was m'alloced and not new'ed */
   presoldata->addedconss = FALSE;
   presoldata->computedsymmetry = FALSE;
   presoldata->enabled = TRUE;
   presoldata->early = FALSE;
   presoldata->nsymresacks = 0;
   presoldata->norbitopes = 0;
   presoldata->ngenconss = 0;
   presoldata->genconss = NULL;
   presoldata->nperms = -1;
   presoldata->log10groupsize = -1.0;
   presoldata->norbits = -1;
   presoldata->orbits = NULL;
   presoldata->orbitbegins = NULL;
   presoldata->ncomponents = -1;
   presoldata->npermsincomponent = NULL;
   presoldata->components = NULL;
   presoldata->componentblocked = NULL;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecSymbreak, presoldata) );

   /* set additional callbacks */
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeSymbreak) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitSymbreak) );
   SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitSymbreak) );
   SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreSymbreak) );

   /* add symmetry breaking presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/conssaddlp",
         "Should the symmetry breaking constraints be added to the LP?",
         &presoldata->conssaddlp, TRUE, DEFAULT_CONSSADDLP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/addsymresacks",
         "Add inequalities for symresacks for each generator?",
         &presoldata->addsymresacks, TRUE, DEFAULT_ADDSYMRESACKS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/computeorbits",
         "Should the orbits of the symmetry group be computed?",
         &presoldata->computeorbits, TRUE, DEFAULT_COMPUTEORBITS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/detectorbitopes",
         "Should we check whether the components of the symmetry group can be handled by orbitopes?",
         &presoldata->detectorbitopes, TRUE, DEFAULT_DETECTORBITOPES, NULL, NULL) );

   return SCIP_OKAY;
}


/** compute non-trivial orbits of symmetry group
 *
 *  The non-tivial orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same orbit appear
 *  consecutively in the orbits array. The variables of the i-th orbit have indices
 *  orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1].
 *  Note that the description of the orbits ends at orbitbegins[norbits] - 1.
 */
SCIP_RETCODE SCIPcomputeGroupOrbitsSymbreak(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            permvars,           /**< variables considered by symbreak presolver */
   int                   npermvars,          /**< length of a permutation array */
   int**                 perms,              /**< matrix containing in each row a permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   SCIP_Shortbool*       activeperms,        /**< array for marking active permutations (or NULL) */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits             /**< pointer to number of orbits currently stored in orbits */
   )
{
   SCIP_Shortbool* varadded;
   int* curorbit;
   int i;
   int curorbitsize;
   int beginneworbit = 0;

   assert( scip != NULL );
   assert( permvars != NULL );
   assert( perms != NULL );
   assert( nperms > 0 );
   assert( npermvars > 0 );
   assert( orbits != NULL );
   assert( norbits != NULL );

   /* init data structures*/
   SCIP_CALL( SCIPallocBufferArray(scip, &curorbit, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varadded, npermvars) );

   /* initially, every variable is contained in no orbit */
   for (i = 0; i < npermvars; ++i)
      varadded[i] = FALSE;

   /* find variable orbits */
   *norbits = 0;
   for (i = 0; i < npermvars; ++i)
   {
      /* if variable is not contained in an orbit of a previous variable */
      if ( ! varadded[i] )
      {
         int j;
         int p;
         int curelem;
         int image;

         /* compute and store orbit if it is non-trivial */
         curorbit[0] = i;
         curorbitsize = 1;
         varadded[i] = TRUE;

         /* iterate over variables in curorbit and compute their images */
         for (j = 0; j < curorbitsize; ++j)
         {
            curelem = curorbit[j];

            for (p = 0; p < nperms; ++p)
            {
               if ( activeperms == NULL || activeperms[p] )
               {
                  image = perms[p][curelem];

                  /* found new element of the orbit of i */
                  if ( ! varadded[image] )
                  {
                     curorbit[curorbitsize++] = image;
                     varadded[image] = TRUE;
                  }
               }
            }
         }

         if ( curorbitsize > 1 )
         {
            orbitbegins[*norbits] = beginneworbit;
            for (j = 0; j < curorbitsize; ++j)
               orbits[beginneworbit++] = curorbit[j];

            ++(*norbits);
         }
      }
   }

   /* store end in "last" orbitbegins entry */
   assert( *norbits < npermvars );
   orbitbegins[*norbits] = beginneworbit;

#ifdef SCIP_OUTPUT
   printf("Orbits (total number: %d):\n", *norbits);
   for (i = 0; i < *norbits; ++i)
   {
      int j;

      printf("%d: ", i);
      for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
         printf("%d ", orbits[j]);
      printf("\n");
   }
#endif

   /* free memory */
   SCIPfreeBufferArray(scip, &varadded);
   SCIPfreeBufferArray(scip, &curorbit);

   return SCIP_OKAY;
}
