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

/**@file   reader_mst.h
 * @ingroup FILEREADERS
 * @brief  file reader for partial primal solutions
 * @author Jakob Witzig
 *
 * This reader handles solutions in two formats:
 *
 * - <b>SCIP raw format</b>@n
 *   The format is as follows:@n@n
 *   line 1: "solution status: <status>"@n
 *   line 2: "objective value: <value>"@n
 *   line 3+i: \<variable name\> \<value\> (obj: \<objective coefficient of variable\>)
 *   @n@n
 *   Only known values need to be listed.
 *   @par
 *   Example:
 *   @code
 *     solution status: optimal
 *     objective value: 1
 *     x1  1 (obj:1)
 *     x2  1 (obj:0)
 *   @endcode
 * - <b>XML format</b>@n
 *   This format is used by CPLEX, for example. For reading we require a section of @p
 *   \<variables\>. Each entry in this section consists of@n
 *   \<variable name="<name>" index="<number>" value="<value>"/>
 *   @par
 *   Example:
 *   @code
 *   <?xml version = "1.0" standalone="yes"?>
 *   <variables>
 *      <variable name="x1" index="1" value="1"/>
 *      <variable name="x2" index="2" value="1"/>
 *   </variables>
 *   </xml>
 *   @endcode
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_MST_H__
#define __SCIP_READER_MST_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the mst file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeReaderMst(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
