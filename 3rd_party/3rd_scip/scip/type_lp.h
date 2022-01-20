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

/**@file   type_lp.h
 * @brief  type definitions for LP management
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_LP_H__
#define __SCIP_TYPE_LP_H__

#ifdef __cplusplus
extern "C" {
#endif

/** solution status after solving LP */
enum SCIP_LPSolStat
{
   SCIP_LPSOLSTAT_NOTSOLVED    = 0,     /**< LP was not solved, no solution exists */
   SCIP_LPSOLSTAT_OPTIMAL      = 1,     /**< LP was solved to optimality */
   SCIP_LPSOLSTAT_INFEASIBLE   = 2,     /**< LP is primal infeasible */
   SCIP_LPSOLSTAT_UNBOUNDEDRAY = 3,     /**< LP has a primal unbounded ray */
   SCIP_LPSOLSTAT_OBJLIMIT     = 4,     /**< objective limit was reached during optimization */
   SCIP_LPSOLSTAT_ITERLIMIT    = 5,     /**< iteration limit was reached during optimization */
   SCIP_LPSOLSTAT_TIMELIMIT    = 6,     /**< time limit was reached during optimization */
   SCIP_LPSOLSTAT_ERROR        = 7      /**< an error occured during optimization */
};
typedef enum SCIP_LPSolStat SCIP_LPSOLSTAT;

/** type of variable bound: lower or upper bound */
enum SCIP_BoundType
{
   SCIP_BOUNDTYPE_LOWER = 0,            /**< lower bound */
   SCIP_BOUNDTYPE_UPPER = 1             /**< upper bound */
};
typedef enum SCIP_BoundType SCIP_BOUNDTYPE;

/** type of row side: left hand or right hand side */
enum SCIP_SideType
{
   SCIP_SIDETYPE_LEFT  = 0,             /**< left hand side */
   SCIP_SIDETYPE_RIGHT = 1              /**< right hand side */
};
typedef enum SCIP_SideType SCIP_SIDETYPE;

/** type of origin of row */
enum SCIP_RowOriginType
{
   SCIP_ROWORIGINTYPE_UNSPEC = 0,       /**< unspecified origin of row */
   SCIP_ROWORIGINTYPE_CONS   = 1,       /**< row created by constraint handler */
   SCIP_ROWORIGINTYPE_SEPA   = 2,       /**< row created by separator */
   SCIP_ROWORIGINTYPE_REOPT  = 3        /**< row created by reoptimization */
};
typedef enum SCIP_RowOriginType SCIP_ROWORIGINTYPE;

/** type of LP algorithm */
enum SCIP_LPAlgo
{
   SCIP_LPALGO_PRIMALSIMPLEX    = 0,    /**< primal simplex */
   SCIP_LPALGO_DUALSIMPLEX      = 1,    /**< dual simplex */
   SCIP_LPALGO_BARRIER          = 2,    /**< barrier algorithm */
   SCIP_LPALGO_BARRIERCROSSOVER = 3     /**< barrier algorithm with crossover */
};
typedef enum SCIP_LPAlgo SCIP_LPALGO;

typedef struct SCIP_ColSolVals SCIP_COLSOLVALS;   /**< collected values of a column which depend on the LP solution */
typedef struct SCIP_RowSolVals SCIP_ROWSOLVALS;   /**< collected values of a row which depend on the LP solution */
typedef struct SCIP_LpSolVals SCIP_LPSOLVALS;     /**< collected values of the LP data which depend on the LP solution */

/** column of an LP
 *
 *  - \ref PublicColumnMethods "List of all available methods"
 */
typedef struct SCIP_Col SCIP_COL;

/** row of an LP
 *
 *  - \ref PublicRowMethods "List of all available methods"
 */
typedef struct SCIP_Row SCIP_ROW;

/** LP structure
 *
 *  - \ref PublicLPMethods "List of all available methods"
 */
typedef struct SCIP_Lp SCIP_LP;

#ifdef __cplusplus
}
#endif

#endif
