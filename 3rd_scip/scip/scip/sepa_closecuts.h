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

/**@file   sepa_closecuts.h
 * @ingroup SEPARATORS
 * @brief  closecuts meta separator
 * @author Marc Pfetsch
 *
 * This separator generates a convex combination of the current LP solution and either the best
 * primal feasible solution or an interior point of the LP relaxation. If the convex combination is
 * proper, the new point is closer to the convex hull of the feasible points. The separator then
 * calls all other separators to separate this point. The idea is that in this way possibly "deeper"
 * cuts are generated. Note, however, that the new point is not a basic solution, i.e., separators
 * relying basis information, e.g., Gomory cuts, will not work.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_CLOSECUTS_H__
#define __SCIP_SEPA_CLOSECUTS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the closecuts separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeSepaClosecuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup SEPARATORS
 *
 * @{
 */

/** sets point to be used as base point for computing the point to be separated
 *
 *  The point is only stored if separation of relative interior points is used. The solution is copied.
 */
EXTERN
SCIP_RETCODE SCIPsetBasePointClosecuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< base point solution */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
