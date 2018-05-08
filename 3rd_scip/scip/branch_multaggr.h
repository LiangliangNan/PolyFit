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

/**@file   branch_multaggr.h
 * @ingroup BRANCHINGRULES
 * @brief  fullstrong branching on fractional and multi-aggregated variables
 * @author Anna Melchiori
 * @author Gerald Gamrath
 *
 * This branching rule uses all fractional binary and integer variables as candidates,
 * as well as fractional multiaggregated binary and integer variables. Although not
 * directly contained in the presolved problem anymore, the multi-aggregation provides
 * an affine linear sum of integer variables, on which branching can be performed.
 *
 * For more details, see
 * G.Gamrath, A.Melchiori, T.Berthold, A.M.Gleixner, D.Salvagnin: Branching on Multi-aggregated Variables
 * (http://dx.doi.org/10.1007/978-3-319-18008-3_10)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_MULTAGGR_H__
#define __SCIP_BRANCH_MULTAGGR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the multi-aggregated branching rule and includes it in SCIP
 *
 *  @ingroup BranchingRuleIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleMultAggr(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
