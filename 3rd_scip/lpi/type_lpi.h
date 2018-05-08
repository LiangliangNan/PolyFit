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

/**@file   type_lpi.h
 * @brief  type definitions for specific LP solvers interface
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_LPI_H__
#define __SCIP_TYPE_LPI_H__

#ifdef __cplusplus
extern "C" {
#endif

/** objective sense */
enum SCIP_ObjSen
{
   SCIP_OBJSEN_MAXIMIZE = -1,           /**< maximize objective function */
   SCIP_OBJSEN_MINIMIZE = +1            /**< minimize objective function */
};
typedef enum SCIP_ObjSen SCIP_OBJSEN;

/** LP solver parameters */
enum SCIP_LPParam
{
   SCIP_LPPAR_FROMSCRATCH    =  0,      /**< solver should start from scratch at next call? */
   SCIP_LPPAR_FASTMIP        =  1,      /**< fast mip setting of LP solver */
   SCIP_LPPAR_SCALING        =  2,      /**< should LP solver use scaling? */
   SCIP_LPPAR_PRESOLVING     =  3,      /**< should LP solver use presolving? */
   SCIP_LPPAR_PRICING        =  4,      /**< pricing strategy */
   SCIP_LPPAR_LPINFO         =  5,      /**< should LP solver output information to the screen? */
   SCIP_LPPAR_FEASTOL        =  6,      /**< feasibility tolerance for primal variables and slacks */
   SCIP_LPPAR_DUALFEASTOL    =  7,      /**< feasibility tolerance for dual variables and reduced costs */
   SCIP_LPPAR_BARRIERCONVTOL =  8,      /**< convergence tolerance used in barrier algorithm */
   SCIP_LPPAR_OBJLIM         =  9,      /**< objective limit (stop if objective is known be larger/smaller than limit for min/max-imization) */
   SCIP_LPPAR_LPITLIM        = 10,      /**< LP iteration limit */
   SCIP_LPPAR_LPTILIM        = 11,      /**< LP time limit */
   SCIP_LPPAR_MARKOWITZ      = 12,      /**< Markowitz tolerance */
   SCIP_LPPAR_ROWREPSWITCH   = 13,      /**< simplex algorithm shall use row representation of the basis
                                         *   if number of rows divided by number of columns exceeds this value */
   SCIP_LPPAR_THREADS        = 14,      /**< number of threads used to solve the LP */
   SCIP_LPPAR_CONDITIONLIMIT = 15,      /**< maximum condition number of LP basis counted as stable */
   SCIP_LPPAR_TIMING         = 16,      /**< type of timer (1 - cpu, 2 - wallclock, 0 - off) */
   SCIP_LPPAR_RANDOMSEED     = 17,      /**< inital random seed, e.g. for perturbations in the simplex (0: LP default) */
   SCIP_LPPAR_POLISHING      = 18,      /**< set solution polishing (0 - disable, 1 - enable) */
   SCIP_LPPAR_REFACTOR       = 19       /**< set refactorization interval (0 - automatic) */
};
typedef enum SCIP_LPParam SCIP_LPPARAM;

/** LP pricing strategy */
enum SCIP_Pricing
{
   SCIP_PRICING_LPIDEFAULT  = 0,        /**< the SCIP/LP interface should use its preferred strategy */
   SCIP_PRICING_AUTO        = 1,        /**< the LP solver should use its preferred strategy */
   SCIP_PRICING_FULL        = 2,        /**< full pricing */
   SCIP_PRICING_PARTIAL     = 3,        /**< partial pricing */
   SCIP_PRICING_STEEP       = 4,        /**< steepest edge pricing */
   SCIP_PRICING_STEEPQSTART = 5,        /**< steepest edge pricing without initial dual norms */
   SCIP_PRICING_DEVEX       = 6         /**< devex pricing */
};
typedef enum SCIP_Pricing SCIP_PRICING;

/** basis status for columns and rows */
enum SCIP_BaseStat
{
   SCIP_BASESTAT_LOWER = 0,             /**< (slack) variable is at its lower bound */
   SCIP_BASESTAT_BASIC = 1,             /**< (slack) variable is basic */
   SCIP_BASESTAT_UPPER = 2,             /**< (slack) variable is at its upper bound */
   SCIP_BASESTAT_ZERO  = 3              /**< free variable is non-basic and set to zero */
};
typedef enum SCIP_BaseStat SCIP_BASESTAT;

/** LP solution quality quantities */
enum SCIP_LPSolQuality
{
   SCIP_LPSOLQUALITY_ESTIMCONDITION = 0,    /**< estimated condition number of (scaled) basis matrix (SCIP_Real) */
   SCIP_LPSOLQUALITY_EXACTCONDITION = 1     /**< exact condition number of (scaled) basis matrix (SCIP_Real) */
};
typedef enum SCIP_LPSolQuality SCIP_LPSOLQUALITY;

typedef struct SCIP_LPi SCIP_LPI;                 /**< solver dependent LP interface */
typedef struct SCIP_LPiState SCIP_LPISTATE;       /**< complete LP state (i.e. basis information) */
typedef struct SCIP_LPiNorms SCIP_LPINORMS;       /**< LP pricing norms information */

#ifdef __cplusplus
}
#endif

#endif
