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

/**@file   presol_symmetry.h
 * @ingroup PRESOLVERS
 * @brief  presolver for storing symmetry information about current problem
 * @author Marc Pfetsch
 * @author Thomas Rehn
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_SYMMETRY_H_
#define __SCIP_PRESOL_SYMMETRY_H_

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <symmetry/type_symmetry.h>

/** include symmetry presolver */
EXTERN
SCIP_RETCODE SCIPincludePresolSymmetry(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return symmetry group generators */
EXTERN
SCIP_RETCODE SCIPgetGeneratorsSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   int*                  nperms,             /**< pointer to store number of permutations */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize      /**< pointer to store log10 of group size (or NULL) */
   );

/** return objective coefficients of permuted variables at time of symmetry computation */
EXTERN
SCIP_RETCODE SCIPgetPermvarsObjSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           permvarsobj         /**< pointer to store objective coefficients of permuted variables (NULL if not available) */
   );

/** register that a specific symmetry is needed */
EXTERN
SCIP_RETCODE SCIPregisterSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_HANDLETYPE        symtype,            /**< type of symmetry handling of callee */
   SYM_SPEC              type,               /**< variable types the callee is interested in */
   SYM_SPEC              fixedtype           /**< variable types that callee wants to have fixed */
   );

/** return at what time symmetry is computed (before or after presolving) */
EXTERN
SCIP_RETCODE SCIPgetTimingSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            afterpresolve       /**< pointer to store whether symmetry is computed in stage initpre or exitpre */
   );

#ifdef __cplusplus
}
#endif

#endif
