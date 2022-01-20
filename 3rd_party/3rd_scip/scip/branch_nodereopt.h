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

/**@file   branch_nodereopt.h
 * @ingroup BRANCHINGRULES
 * @brief  nodereopt branching rule
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_NODEREOPT_H__
#define __SCIP_BRANCH_NODEREOPT_H__


#include "scip/scip.h"
#include "scip/type_branch.h"
#include "scip/type_reopt.h"
#include "scip/type_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the nodereopt branching rule and includes it in SCIP
 *
 *  @ingroup BranchingRuleIncludes
 */
extern
SCIP_RETCODE SCIPincludeBranchruleNodereopt(
   SCIP*                 scip                /**< SCIP data structure */
   );


#ifdef __cplusplus
}
#endif

#endif
