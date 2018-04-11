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

/**@file   presol_implics.h
 * @ingroup PRESOLVERS
 * @brief  implication graph presolver which checks for aggregations
 * @author Tobias Achterberg
 *
 * This presolver looks for variable implications in \f$x == 0\f$ and \f$x == 1\f$ with the same implied variable.
 * There are four possible cases:
 * \f[
 *  x = 0 \Rightarrow y = lb,\; \mathrm{and}\; x = 1 \Rightarrow y = lb:\; \mathrm{fix}\; y\; \mathrm{to}\; lb
 * \f]
 * \f[
 *  x = 0 \Rightarrow y = lb,\; \mathrm{and}\; x = 1 \Rightarrow y = ub:\; \mathrm{aggregate}\; y == lb + (ub-lb)x
 * \f]
 * \f[
 *  x = 0 \Rightarrow y = ub,\; \mathrm{and}\; x = 1 \Rightarrow y = lb:\; \mathrm{aggregate}\;  y == ub - (ub-lb)x
 * \f]
 * \f[
 *  x = 0 \Rightarrow y = ub,\; \mathrm{and}\; x = 1 \Rightarrow y = ub:\; \mathrm{fix}\; y\; \mathrm{to}\; ub
 * \f]
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_IMPLICS_H__
#define __SCIP_PRESOL_IMPLICS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the implics presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludePresolImplics(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
