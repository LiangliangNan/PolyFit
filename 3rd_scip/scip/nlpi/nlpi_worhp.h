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

/**@file    nlpi_worhp.h
 * @brief   Worhp NLP interface
 * @ingroup NLPIS
 * @author  Benjamin Mueller
 * @author  Renke Kuhlmann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_WORHP_H__
#define __SCIP_NLPI_WORHP_H__

#include "nlpi/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup NLPIS
 *
 * @{
 */

/** create solver interface for Worhp solver */
extern
SCIP_RETCODE SCIPcreateNlpSolverWorhp(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi,               /**< pointer to buffer for nlpi address */
   SCIP_Bool             useip               /**< TRUE for using Interior Point, FALSE for SQP */
   );

/** gets string that identifies Worhp (version number) */
extern
const char* SCIPgetSolverNameWorhp(
   void
   );

/** gets string that describes Worhp (version number) */
extern
const char* SCIPgetSolverDescWorhp(
   void
   );

/** returns whether Worhp is available, i.e., whether it has been linked in */
extern
SCIP_Bool SCIPisWorhpAvailableWorhp(
   void
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_WORHP_H__ */
