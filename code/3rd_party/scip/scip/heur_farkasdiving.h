/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_farkasdiving.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LP diving heuristic that tries to construct a Farkas-proof
 * @author Jakob Witzig
 *
 * The heuristic dives into the direction of the pseudosolution, i.e., variables get rounded
 * towards their best bound w.r.t there objective coefficient. This strategy is twofold, if
 * a feasible solution is found the solution has potentially a very good objective value; on the other
 * hand, the left-hand side of a potential Farkas-proof \f$y^Tb - y^TA{l',u'} > 0\f$ (i.e., infeasibility proof)
 * gets increased, where \f$l',u'\f$ are the local bounds. The contribution of each variable \f$x_i\f$ to the
 * Farkas-proof can be approximated by \f$c_i = y^TA_i\f$ because we only dive on basic variables with
 * reduced costs \f$c_i - y^TA_i = 0\f$.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_FARKASDIVING_H__
#define __SCIP_HEUR_FARKASDIVING_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the farkasdiving heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurFarkasdiving(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
