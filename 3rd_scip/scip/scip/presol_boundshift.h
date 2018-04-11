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

/**@file   presol_boundshift.h
 * @ingroup PRESOLVERS
 * @brief  presolver that converts integer variables with domain [a,b] to integer variables with domain [0,b-a]
 * @author Tobias Achterberg
 * @author Michael Winkler
 *
 * This presolver converts all integer variables with domain \f$[a,b]\f$ to integer variables with domain
 * \f$[0,b-a]\f$. This is done by creating a new integer variable \f$y\f$ which will be aggregated to the old variable
 * \f$x\f$ such that
 * \f[
 * x = y + a
 * \f]
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_BOUNDSHIFT_H__
#define __SCIP_PRESOL_BOUNDSHIFT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the boundshift presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludePresolBoundshift(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
