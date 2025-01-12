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

/**@file   symmetry.h
 * @ingroup PUBLICCOREAPI
 * @brief  methods for handling symmetries
 * @author Christopher Hojny
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SYMMETRY_H__
#define __SCIP_SYMMETRY_H__

#include "scip/def.h"
#include "scip/pub_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include <symmetry/type_symmetry.h>

#ifdef __cplusplus
extern "C" {
#endif


/**@addtogroup PublicSymmetryMethods
 *
 * @{
 */


/** compute non-trivial orbits of symmetry group
 *
 *  The non-tivial orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same orbit appear
 *  consecutively in the orbits array. The variables of the i-th orbit have indices
 *  orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1].
 *  Note that the description of the orbits ends at orbitbegins[norbits] - 1.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeOrbitsSym(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            permvars,           /**< variables considered in a permutation array */
   int                   npermvars,          /**< length of a permutation array */
   int**                 perms,              /**< matrix containing in each row a permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits             /**< pointer to number of orbits currently stored in orbits */
   );


/** compute non-trivial orbits of symmetry group using filtered generators
 *
 *  The non-trivial orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same orbit appear
 *  consecutively in the orbits array. The variables of the i-th orbit have indices
 *  orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1].
 *  Note that the description of the orbits ends at orbitbegins[norbits] - 1.
 *
 *  Only permutations that are not inactive (as marked by @p inactiveperms) are used. Thus, one can use this array to
 *  filter out permutations.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeOrbitsFilterSym(
   SCIP*                 scip,               /**< SCIP instance */
   int                   npermvars,          /**< length of a permutation array */
   int**                 permstrans,         /**< transposed matrix containing in each column a
                                              *   permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   SCIP_Shortbool*       inactiveperms,      /**< array to store whether permutations are inactive */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits,            /**< pointer to number of orbits currently stored in orbits */
   int*                  components,         /**< array containing the indices of permutations sorted by components */
   int*                  componentbegins,    /**< array containing in i-th position the first position of
                                              *   component i in components array */
   int*                  vartocomponent,     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   unsigned*             componentblocked,   /**< array to store which symmetry methods have been used on a component
                                              *   using the same bitset information as for misc/usesymmetry */
   int                   ncomponents,        /**< number of components of symmetry group */
   int                   nmovedpermvars      /**< number of variables moved by any permutation in a symmetry component
                                              *   that is handled by orbital fixing */
   );

/** compute non-trivial orbits of symmetry group
 *
 *  The non-tivial orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same orbit appear
 *  consecutively in the orbits array. The variables of the i-th orbit have indices
 *  orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1].
 *  Note that the description of the orbits ends at orbitbegins[norbits] - 1.
 *
 *  This function is adapted from SCIPcomputeOrbitsFilterSym().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeOrbitsComponentsSym(
   SCIP*                 scip,               /**< SCIP instance */
   int                   npermvars,          /**< length of a permutation array */
   int**                 permstrans,         /**< transposed matrix containing in each column a permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   int*                  components,         /**< array containing the indices of permutations sorted by components */
   int*                  componentbegins,    /**< array containing in i-th position the first position of component i in components array */
   int*                  vartocomponent,     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   int                   ncomponents,        /**< number of components of symmetry group */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits,            /**< pointer to number of orbits currently stored in orbits */
   int*                  varorbitmap         /**< array for storing the orbits for each variable */
   );

/** Compute orbit of a given variable and store it in @p orbit. The first entry of the orbit will
 *  be the given variable index and the rest is filled with the remaining variables excluding
 *  the ones specified in @p ignoredvars.
 *
 *  @pre orbit is an initialized array of size propdata->npermvars
 *  @pre at least one of @p perms and @p permstrans should not be NULL
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeOrbitVar(
   SCIP*                 scip,               /**< SCIP instance */
   int                   npermvars,          /**< number of variables in permvars */
   int**                 perms,              /**< the generators of the permutation group (or NULL) */
   int**                 permstrans,         /**< the transposed matrix of generators (or NULL) */
   int*                  components,         /**< the components of the permutation group */
   int*                  componentbegins,    /**< array containing the starting index of each component */
   SCIP_Shortbool*       ignoredvars,        /**< array indicating which variables should be ignored */
   SCIP_Shortbool*       varfound,           /**< bitmap to mark which variables have been added (or NULL) */
   int                   varidx,             /**< index of variable for which the orbit is requested */
   int                   component,          /**< component that var is in */
   int *                 orbit,              /**< array in which the orbit should be stored */
   int*                  orbitsize           /**< buffer to store the size of the orbit */
   );

/** Checks whether a permutation is a composition of 2-cycles and in this case determine the number of overall
 *  2-cycles and binary 2-cycles. It is a composition of 2-cycles iff @p ntwocyclesperm > 0 upon termination.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPisInvolutionPerm(
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< array of variables perm is acting on */
   int                   nvars,              /**< number of variables */
   int*                  ntwocyclesperm,     /**< pointer to store number of 2-cycles */
   int*                  nbincyclesperm,     /**< pointer to store number of binary cycles */
   SCIP_Bool             earlytermination    /**< whether we terminate early if not all affected variables are binary */
   );

/** determine number of variables affected by symmetry group */
SCIP_EXPORT
SCIP_RETCODE SCIPdetermineNVarsAffectedSym(
   SCIP*                 scip,               /**< SCIP instance */
   int**                 perms,              /**< permutations */
   int                   nperms,             /**< number of permutations in perms */
   SCIP_VAR**            permvars,           /**< variables corresponding to permutations */
   int                   npermvars,          /**< number of permvars in perms */
   int*                  nvarsaffected       /**< pointer to store number of all affected variables */
   );

/** compute components of symmetry group */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeComponentsSym(
   SCIP*                 scip,               /**< SCIP instance */
   int**                 perms,              /**< permutation generators as
                                              *   (either nperms x npermvars or npermvars x nperms) matrix */
   int                   nperms,             /**< number of permutations */
   SCIP_VAR**            permvars,           /**< variables on which permutations act */
   int                   npermvars,          /**< number of variables for permutations */
   SCIP_Bool             transposed,         /**< transposed permutation generators as (npermvars x nperms) matrix */
   int**                 components,         /**< array containing the indices of permutations sorted by components */
   int**                 componentbegins,    /**< array containing in i-th position the first position of
                                              *   component i in components array */
   int**                 vartocomponent,     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   unsigned**            componentblocked,   /**< array to store which symmetry methods have been used on a component
                                              *   using the same bitset information as for misc/usesymmetry */
   int*                  ncomponents         /**< pointer to store number of components of symmetry group */
   );

/** Given a matrix with nrows and \#perms + 1 columns whose first nfilledcols columns contain entries of variables, this routine
 *  checks whether the 2-cycles of perm intersect each row of column coltoextend in exactly one position. In this case,
 *  we add one column to the suborbitope of the first nfilledcols columns.
 *
 *  @pre Every non-trivial cycle of perm is a 2-cycle.
 *  @pre perm has nrows many 2-cycles
 */
SCIP_EXPORT
SCIP_RETCODE SCIPextendSubOrbitope(
   int**                 suborbitope,        /**< matrix containing suborbitope entries */
   int                   nrows,              /**< number of rows of suborbitope */
   int                   nfilledcols,        /**< number of columns of suborbitope which are filled with entries */
   int                   coltoextend,        /**< index of column that should be extended by perm */
   int*                  perm,               /**< permutation */
   SCIP_Bool             leftextension,      /**< whether we extend the suborbitope to the left */
   int**                 nusedelems,         /**< pointer to array storing how often an element was used in the orbitope */
   SCIP_VAR**            permvars,           /**< permutation vars array */
   SCIP_Shortbool*       rowisbinary,        /**< array encoding whether variables in an orbitope row are binary */
   SCIP_Bool*            success,            /**< pointer to store whether extension was successful */
   SCIP_Bool*            infeasible          /**< pointer to store if the number of intersecting cycles is too small */
   );

/** generate variable matrix for orbitope constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPgenerateOrbitopeVarsMatrix(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR****          vars,               /**< pointer to matrix of orbitope variables */
   int                   nrows,              /**< number of rows of orbitope */
   int                   ncols,              /**< number of columns of orbitope */
   SCIP_VAR**            permvars,           /**< superset of variables that are contained in orbitope */
   int                   npermvars,          /**< number of variables in permvars array */
   int**                 orbitopevaridx,     /**< permuted index table of variables in permvars that are contained in orbitope */
   int*                  columnorder,        /**< permutation to reorder column of orbitopevaridx */
   int*                  nusedelems,         /**< array storing how often an element was used in the orbitope */
   SCIP_Shortbool*       rowisbinary,        /**< array encoding whether a row contains only binary variables */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the potential orbitope is not an orbitope */
   SCIP_Bool             storelexorder,      /**< whether the lexicographic order induced by the orbitope shall be stored */
   int**                 lexorder,           /**< pointer to array storing the lexorder (or NULL) */
   int*                  nvarsorder,         /**< pointer to store number of variables in lexorder (or NULL) */
   int*                  maxnvarsorder       /**< pointer to store maximum number of variables in lexorder (or NULL) */
   );

/** checks whether an orbitope is a packing or partitioning orbitope */
SCIP_RETCODE SCIPisPackingPartitioningOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< variable matrix of orbitope constraint */
   int                   nrows,              /**< pointer to number of rows of variable matrix */
   int                   ncols,              /**< number of columns of variable matrix */
   SCIP_Bool**           pprows,             /**< pointer to store which rows are are contained in
                                              *   packing/partitioning constraints or NULL if not needed */
   int*                  npprows,            /**< pointer to store how many rows are contained
                                              *   in packing/partitioning constraints or NULL if not needed */
   SCIP_ORBITOPETYPE*    type                /**< pointer to store type of orbitope constraint after strengthening */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
