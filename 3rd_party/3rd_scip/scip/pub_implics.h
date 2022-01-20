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

/**@file   pub_implics.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for implications, variable bounds, and cliques
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_IMPLICS_H__
#define __SCIP_PUB_IMPLICS_H__


#include "scip/def.h"
#include "scip/type_var.h"
#include "scip/type_implics.h"

#ifdef NDEBUG
#include "scip/struct_implics.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * methods for cliques
 */

/** returns the position of the given variable/value pair in the clique; returns -1 if variable/value pair is not member
 *  of the clique
 */
EXTERN
int SCIPcliqueSearchVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   SCIP_VAR*             var,                /**< variable to search for */
   SCIP_Bool             value               /**< value of the variable in the clique */
   );

/** returns whether the given variable/value pair is member of the given clique */
EXTERN
SCIP_Bool SCIPcliqueHasVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   SCIP_VAR*             var,                /**< variable to remove from the clique */
   SCIP_Bool             value               /**< value of the variable in the clique */
   );

/** gets number of variables in the cliques */
EXTERN
int SCIPcliqueGetNVars(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   );

/** gets array of active problem variables in the cliques */
EXTERN
SCIP_VAR** SCIPcliqueGetVars(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   );

/** gets array of values of active problem variables in the cliques, i.e. whether the variable is fixed to FALSE or
 *  to TRUE in the clique
 */
EXTERN
SCIP_Bool* SCIPcliqueGetValues(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   );

/** gets unique identifier of the clique */
EXTERN
unsigned int SCIPcliqueGetId(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   );

/** gets index of the clique in the clique table */
EXTERN
int SCIPcliqueGetIndex(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   );

/** returns whether the given clique is cleaned up */
EXTERN
SCIP_Bool SCIPcliqueIsCleanedUp(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   );

/** return whether the given clique is an equation */
EXTERN
SCIP_Bool SCIPcliqueIsEquation(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   );


#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPcliqueGetNVars(clique)                   ((clique)->nvars)
#define SCIPcliqueGetVars(clique)                    ((clique)->vars)
#define SCIPcliqueGetValues(clique)                  ((clique)->values)
#define SCIPcliqueGetId(clique)                      ((clique)->id)
#define SCIPcliqueGetIndex(clique)                      ((clique)->index)
#define SCIPcliqueIsCleanedUp(clique)                ((clique)->startcleanup == -1)
#define SCIPcliqueIsEquation(clique)                 ((SCIP_Bool)(clique)->equation)

#endif

#ifdef __cplusplus
}
#endif

#endif
