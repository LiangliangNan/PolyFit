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

/**@file   sepa_aggregation.h
 * @ingroup SEPARATORS
 * @brief  flow cover and complemented mixed integer rounding cuts separator (Marchand's version)
 * @author Robert Lion Gottwald
 * @author Kati Wolter
 * @author Tobias Achterberg
 *
 * For an overview see:
 *
 * Marchand, H., & Wolsey, L. A. (2001).@n
 * Aggregation and mixed integer rounding to solve MIPs.@n
 * Operations research, 49(3), 363-371.
 *
 * Some remarks:
 * - In general, continuous variables are less prefered than integer variables, since their cut
 *   coefficient is worse.
 * - We seek for aggregations that project out continuous variables that are far away from their bound,
 *   since if it is at its bound then it doesn't contribute to the violation
 * - These aggregations are also useful for the flowcover separation, so after building an aggregation
 *   we try to generate a MIR cut and a flowcover cut.
 * - We only keep the best cut.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_CMIR_H__
#define __SCIP_SEPA_CMIR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the aggregation separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeSepaAggregation(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
