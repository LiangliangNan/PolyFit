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

/**@file    nlpi_worhp_dummy.c
 * @ingroup NLPIS
 * @brief   dummy WORHP NLP interface
 * @author  Benjamin Mueller
 */

#include "nlpi/nlpi_worhp.h"

/** create solver interface for Worhp solver */
SCIP_RETCODE SCIPcreateNlpSolverWorhp(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi,               /**< pointer to buffer for nlpi address */
   SCIP_Bool             useip               /**< TRUE for using Interior Point, FALSE for SQP */
   )
{
   *nlpi = NULL;

   return SCIP_OKAY;
} /*lint !e715*/

/** gets string that identifies Worhp (version number) */
const char* SCIPgetSolverNameWorhp(void)
{
   return "WORHP";
}

/** gets string that describes Worhp (version number) */
const char* SCIPgetSolverDescWorhp(void)
{
   return "this is WORHP";
}

/** returns whether Worhp is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisWorhpAvailableWorhp(void)
{
   return FALSE;
}
