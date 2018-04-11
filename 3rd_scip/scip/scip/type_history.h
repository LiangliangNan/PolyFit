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

/**@file   type_history.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for branching and inference history
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_HISTORY_H__
#define __SCIP_TYPE_HISTORY_H__

#ifdef __cplusplus
extern "C" {
#endif

/** branching direction for branching on variables */
enum SCIP_BranchDir
{
   SCIP_BRANCHDIR_DOWNWARDS = 0,        /**< downwards branching: decreasing upper bound  */
   SCIP_BRANCHDIR_UPWARDS   = 1,        /**< upwards branching: increasing lower bound    */
   SCIP_BRANCHDIR_FIXED     = 2,        /**< fixed branching: both bounds changed         */
   SCIP_BRANCHDIR_AUTO      = 3         /**< automatic setting for choosing bound changes */
};
typedef enum SCIP_BranchDir SCIP_BRANCHDIR;       /**< branching direction for branching on variables */

typedef struct SCIP_History SCIP_HISTORY;         /**< branching and inference history information for single variable */

/** Value history data structure
 *
 *  branching and inference history informations for single variable dependent on the domain value
 *
 *  - \ref SCIP_VALUEHISTORY "List of all available methods"
 */
typedef struct SCIP_ValueHistory SCIP_VALUEHISTORY;


#ifdef __cplusplus
}
#endif

#endif
