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

/**@file   sepa_mcf.h
 * @ingroup SEPARATORS
 * @brief  multi-commodity-flow network cut separator
 * @author Tobias Achterberg
 * @author Christian Raack
 *
 * We try to identify a multi-commodity flow structure in the LP relaxation of the
 * following type:
 *
 *  (1)  sum_{a in delta^+(v)} f_a^k  - sum_{a in delta^-(v)} f_a^k  <=  -d_v^k   for all v in V and k in K
 *  (2)  sum_{k in K} f_a^k - c_a x_a                                <=  0        for all a in A
 *
 * Constraints (1) are flow conservation constraints, which say that for each commodity k and node v the
 * outflow (delta^+(v)) minus the inflow (delta^-(v)) of a node v must not exceed the negative of the demand of
 * node v in commodity k. To say it the other way around, inflow minus outflow must be at least equal to the demand.
 * Constraints (2) are the arc capacity constraints, which say that the sum of all flow over an arc a must not
 * exceed its capacity c_a x_a, with x being a binary or integer variable.
 * c_a x_a does not need to be a single product of a capacity and an integer variable; we also accept general scalar
 * products.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_MCF_H__
#define __SCIP_SEPA_MCF_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the mcf separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeSepaMcf(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
