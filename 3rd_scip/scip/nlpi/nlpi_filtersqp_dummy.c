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

/**@file    nlpi_filtersqp_dummy.c
 * @ingroup NLPIS
 * @brief   dummy filterSQP NLP interface for the case that FilterSQP is not available
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "nlpi/nlpi_filtersqp.h"

/** create solver interface for filterSQP solver */
SCIP_RETCODE SCIPcreateNlpSolverFilterSQP(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi                /**< pointer to buffer for nlpi address */
   )
{
   assert(nlpi != NULL);

   *nlpi = NULL;

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets string that identifies filterSQP */
const char* SCIPgetSolverNameFilterSQP(void)
{
   return "";
}

/** gets string that describes filterSQP */
const char* SCIPgetSolverDescFilterSQP(void)
{
   return "";
}

/** returns whether filterSQP is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisFilterSQPAvailableFilterSQP(void)
{
   return FALSE;
}
