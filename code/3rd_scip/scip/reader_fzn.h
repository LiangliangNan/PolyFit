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

/**@file   reader_fzn.h
 * @ingroup FILEREADERS
 * @brief  FlatZinc file reader
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * FlatZinc is a low-level solver input language that is the target language for MiniZinc. It is designed to be easy to
 * translate into the form required by a solver. For more details see http://www.g12.cs.mu.oz.au/minizinc/ .
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __READER_FZN_H__
#define __READER_FZN_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the FlatZinc file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeReaderFzn(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup FILEREADERS
 *
 * @{
 */

/** print given solution in Flatzinc format w.r.t. the output annotation */
EXTERN
SCIP_RETCODE SCIPprintSolReaderFzn(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
