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

/**@file   sepa_disjunctive.h
 * @ingroup SEPARATORS
 * @brief  disjunctive cut separator
 * @author Tobias Fischer
 * @author Marc Pfetsch
 *
 * We separate disjunctive cuts for two term disjunctions of the form \f$x_1 = 0 \vee x_2 = 0\f$. They can be generated
 * directly from the simplex tableau. For further information, we refer to@n
 * "A complementarity-based partitioning and disjunctive cut algorithm for mathematical programming problems with
 * equilibrium constraints"@n
 * Júdice, J.J., Sherali, H.D., Ribeiro, I.M., Faustino, A.M., Journal of Global Optimization 36(1), 89–114 (2006)
 *
 * Cut coefficients belonging to integer variables can be strengthened by the 'monoidal cut strengthening' procedure, see@n
 * "Strengthening cuts for mixed integer programs"@n
 * Egon Balas, Robert G. Jeroslow, European Journal of Operational Research, Volume 4, Issue 4, 1980, Pages 224-234
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_DISJUNCTIVE_H__
#define __SCIP_SEPA_DISJUNCTIVE_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the disjunctive cut separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeSepaDisjunctive(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
