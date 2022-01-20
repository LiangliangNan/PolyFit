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

/**@file   struct_sol.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_STRUCT_SOL_H__
#define __SCIP_STRUCT_SOL_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_sol.h"
#include "scip/type_heur.h"


#ifdef __cplusplus
extern "C" {
#endif

/** maximum violations of problem constraints */
struct SCIP_Viol
{
   SCIP_Real             absviolbounds;       /**< maximum absolute violation of variable bounds */
   SCIP_Real             relviolbounds;       /**< maximum relative violation of variable bounds */
   SCIP_Real             absviollprows;       /**< maximum absolute violation of LP rows */
   SCIP_Real             relviollprows;       /**< maximum relative violation of LP rows */
   SCIP_Real             absviolintegrality;  /**< maximum absolute violation integrality */
   SCIP_Real             absviolcons;         /**< maximum absolute violation of constraints */
   SCIP_Real             relviolcons;         /**< maximum relative violation of constraints */
};

/** primal CIP solution
 *
 *  For reasons of efficiency, a working solution only stores values that have been accessed at least once,
 *  or that have been changed from the value in the solution's source.
 *  The user has to call SCIPsolUnlink() in order to retrieve all non-cached elements from the solution's source
 *  and to store the values in the solution's own array. This changes the solution's origin to SCIP_SOLORIGIN_ZERO.
 *  A linked solution with origin SCIP_SOLORIGIN_LPSOL or SCIP_SOLORIGIN_PSEUDOSOL becomes invalid after the
 *  next node is focused (i.e. the LP and pseudo solutions changed) and cannot be accessed anymore.
 *
 *  Solutions with origin ORIGINAL contain the values for original variables. The stored objective value also
 *  corresponds to the original problem.
 */
struct SCIP_Sol
{
   SCIP_Real             obj;                /**< objective value of solution */
   SCIP_Real             time;               /**< clock time, when the solution was discovered */
   SCIP_Longint          nodenum;            /**< last node number of current run, where this solution was modified */
#ifndef NDEBUG
   SCIP_Longint          lpcount;            /**< number of LPs solved when this solution was created, needed for debug checks
                                              *   concerning solutions linked to the LP solution
                                              */
#endif
   SCIP_REALARRAY*       vals;               /**< solution values for variables */
   SCIP_BOOLARRAY*       valid;              /**< is value in vals array valid? otherwise it has to be retrieved from
                                              *   origin */
   SCIP_HEUR*            heur;               /**< heuristic that found the solution (or NULL if it's an LP solution) */
   SCIP_VIOL             viol;               /**< maximum violations of problem constraints */
   int                   runnum;             /**< branch and bound run number in which the solution was found */
   int                   depth;              /**< depth at which the solution was found */
   int                   primalindex;        /**< index of solution in array of existing solutions of primal data */
   int                   index;              /**< consecutively numbered unique index of all created solutions */
   SCIP_SOLORIGIN        solorigin;          /**< origin of solution: where to retrieve uncached elements */
   SCIP_Bool             hasinfval;          /**< does the solution (potentially) contain an infinite value? Note: this
                                              * could also be implemented as a counter for the number of infinite
                                              * values, to avoid redundant checks when resetting inf. solution values
                                              */
};

#ifdef __cplusplus
}
#endif

#endif
