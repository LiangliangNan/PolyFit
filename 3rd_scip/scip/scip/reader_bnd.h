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

/**@file   reader_bnd.h
 * @ingroup FILEREADERS
 * @brief  file reader for variable bounds
 * @author Ambros Gleixner
 *
 * This reader allows to read a file containing new bounds for variables of the current problem.  Each line of the file
 * should have format
 *
 *    \<variable name\> \<lower bound\> \<upper bound\>
 *
 * where infinite bounds can be written as inf, +inf or -inf.  Note that only a subset of the variables may appear in
 * the file.  Lines with unknown variable names are ignored.  The writing functionality is currently not supported.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_BND_H__
#define __SCIP_READER_BND_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the bnd file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeReaderBnd(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
