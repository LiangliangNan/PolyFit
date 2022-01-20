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

/**@file   type_stat.h
 * @brief  type definitions for problem statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_STAT_H__
#define __SCIP_TYPE_STAT_H__

#ifdef __cplusplus
extern "C" {
#endif

/** SCIP solving status */
enum SCIP_Status
{
   SCIP_STATUS_UNKNOWN        =  0,     /**< the solving status is not yet known */
   SCIP_STATUS_USERINTERRUPT  =  1,     /**< the user interrupted the solving process (by pressing CTRL-C) */
   SCIP_STATUS_NODELIMIT      =  2,     /**< the solving process was interrupted because the node limit was reached */
   SCIP_STATUS_TOTALNODELIMIT =  3,     /**< the solving process was interrupted because the total node limit was
                                         *   reached (incl. restarts)
                                         */
   SCIP_STATUS_STALLNODELIMIT =  4,     /**< the solving process was interrupted because the stalling node limit was
                                         *   reached (no inprovement w.r.t. primal bound)
                                         */
   SCIP_STATUS_TIMELIMIT      =  5,     /**< the solving process was interrupted because the time limit was reached */
   SCIP_STATUS_MEMLIMIT       =  6,     /**< the solving process was interrupted because the memory limit was reached */
   SCIP_STATUS_GAPLIMIT       =  7,     /**< the solving process was interrupted because the gap limit was reached */
   SCIP_STATUS_SOLLIMIT       =  8,     /**< the solving process was interrupted because the solution limit was
                                          *  reached
                                          */
   SCIP_STATUS_BESTSOLLIMIT   =  9,    /**< the solving process was interrupted because the solution improvement limit
                                          *   was reached
                                          */
   SCIP_STATUS_RESTARTLIMIT   = 10,     /**< the solving process was interrupted because the restart limit was reached */
   SCIP_STATUS_OPTIMAL        = 11,     /**< the problem was solved to optimality, an optimal solution is available */
   SCIP_STATUS_INFEASIBLE     = 12,     /**< the problem was proven to be infeasible */
   SCIP_STATUS_UNBOUNDED      = 13,     /**< the problem was proven to be unbounded */
   SCIP_STATUS_INFORUNBD      = 14      /**< the problem was proven to be either infeasible or unbounded */
};
typedef enum SCIP_Status SCIP_STATUS;

typedef struct SCIP_Stat SCIP_STAT;               /**< problem and runtime specific statistics */

#ifdef __cplusplus
}
#endif

#endif
