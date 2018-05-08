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

/**@file   reader_mps.h
 * @ingroup FILEREADERS
 * @brief  (extended) MPS file reader
 * @author Thorsten Koch
 * @author Tobias Achterberg
 *
 * This reader allows to parse and write MPS files with linear and quadratic constraints and objective,
 * special ordered sets of type 1 and 2, indicators on linear constraints, and semicontinuous variables.
 * For writing, linear (general and specialized), indicator, quadratic, second order cone, and
 * special ordered set constraints are supported.
 *
 * See http://en.wikipedia.org/wiki/MPS_%28format%29 for a description.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_MPS_H__
#define __SCIP_READER_MPS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the mps file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeReaderMps(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
