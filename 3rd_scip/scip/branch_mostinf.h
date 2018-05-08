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

/**@file   branch_mostinf.h
 * @ingroup BRANCHINGRULES
 * @brief  most infeasible LP branching rule
 * @author Tobias Achterberg
 *
 * The most infeasible branching rule selects a candidate variable $j$ with fractional solution value \f$ \hat{x}_j\f$
 * which maximizes
 * \f[
 *      \min \left\{ \lceil \hat{x}_j \rceil - \hat{x}_j, \hat{x}_j - \lfloor \hat{x}_j \rfloor  \right\}.
 * \f]
 * i. e., a variable which still is farthest from taking an integer value among all branching candidates.
 *
 * The most infeasible branching rule and many other branching rules of SCIP are explained and compared in
 *
 * @par
 * Tobias Achterberg@n
 * Constraint Integer Programming@n
 * PhD Thesis, Technische Universität Berlin, 2007@n
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_MOSTINF_H__
#define __SCIP_BRANCH_MOSTINF_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the most infeasible LP branching rule and includes it in SCIP
 *
 *  @ingroup BranchingRuleIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleMostinf(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
