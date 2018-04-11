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

/**@file    nlpi_all.h
 * @brief   ALL NLP interface
 * @brief   NLP interface that uses all available NLP interfaces
 * @author  Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_ALL_H__
#define __SCIP_NLPI_ALL_H__

#include "nlpi/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup NLPIS
 *
 * @{
 */

/** create solver interface for All solver */
extern
SCIP_RETCODE SCIPcreateNlpSolverAll(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi,               /**< pointer to buffer for nlpi address */
   SCIP_NLPI**           nlpis,              /**< array containing existing nlpis */
   int                   nnlpis              /**< total number of nlpis */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_ALL_H__ */
