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

/**@file   pub_history.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for branching and inference history structure
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_HISTORY_H__
#define __SCIP_PUB_HISTORY_H__

#ifdef NDEBUG
#include "scip/struct_history.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NDEBUG

/** gets the conflict score of the history entry */
EXTERN
SCIP_Real SCIPhistoryGetVSIDS(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction */
   );

/** get number of cutoffs counter */
EXTERN
SCIP_Real SCIPhistoryGetCutoffSum(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** return the number of (domain) values for which a history exists */
EXTERN
int SCIPvaluehistoryGetNValues(
   SCIP_VALUEHISTORY*    valuehistory        /**< value based history */
   );

/** return the array containing the histories for the individual (domain) values */
EXTERN
SCIP_HISTORY** SCIPvaluehistoryGetHistories(
   SCIP_VALUEHISTORY*    valuehistory        /**< value based history */
   );

/** return the array containing the (domain) values for which a history exists */
EXTERN
SCIP_Real* SCIPvaluehistoryGetValues(
   SCIP_VALUEHISTORY*    valuehistory        /**< value based history */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPhistoryGetVSIDS(history,dir)   ((history)->vsids[dir])

#define SCIPvaluehistoryGetNValues(valuehistory)     (valuehistory)->nvalues
#define SCIPvaluehistoryGetHistories(valuehistory)      (valuehistory)->histories
#define SCIPvaluehistoryGetValues(valuehistory)      (valuehistory)->values

#endif


#ifdef __cplusplus
}
#endif

#endif
