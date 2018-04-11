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

/**@file   reader_osil.h
 * @ingroup FILEREADERS
 * @brief  OS instance language (OSiL) format file reader
 * @author Stefan Vigerske
 *
 * This reader allows to parse OSiL files with linear and nonlinear constraints and objective.
 * Writing is not implemented yet.
 *
 * The OSiL format is an XML based format to represent a broad class of mathematical programming instances, see http://www.coin-or.org/OS/OSiL.html .
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_OSIL_H__
#define __SCIP_READER_OSIL_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the osil file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeReaderOsil(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
