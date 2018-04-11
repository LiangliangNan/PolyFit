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

/**@file   tpi_none.c
 * @ingroup TASKINTERFACE
 * @brief  the interface functions for dummy tpi
 * @author Stephen J. Maher
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "tpi/tpi.h"

/** creates a job for parallel processing */
SCIP_RETCODE SCIPtpiCreateJob(
   SCIP_JOB**            job,                /**< pointer to the job that will be created */
   int                   jobid,              /**< the id for the current job */
   SCIP_RETCODE          (*jobfunc)(void* args),/**< pointer to the job function */
   void*                 jobarg              /**< the job's argument */
   )
{
   SCIP_UNUSED( job );
   SCIP_UNUSED( jobid );
   SCIP_UNUSED( jobfunc );
   SCIP_UNUSED( jobarg );

   return SCIP_ERROR;
}

/** submit a job for parallel processing */
/* the return is a globally defined status */
SCIP_RETCODE SCIPtpiSumbitJob(
   SCIP_JOB*             job,                /**< pointer to the job to be submitted */
   SCIP_SUBMITSTATUS*    status              /**< pointer to store the job's submit status */
   )
{
   SCIP_UNUSED( job );
   SCIP_UNUSED( status );

   return SCIP_ERROR;
}

/** Blocks until all jobs with the given jobid have finished
 * and then returns the smallest SCIP_RETCODE of all the jobs */
SCIP_RETCODE SCIPtpiCollectJobs(
   int                   jobid               /**< the id of the jobs to collect */
   )
{
   SCIP_UNUSED( jobid );

   return SCIP_ERROR;
}

/** initializes tpi */
SCIP_RETCODE SCIPtpiInit(
   int         nthreads,                     /**< the number of threads to be used */
   int         queuesize,                    /**< the size of the queue */
   SCIP_Bool   blockwhenfull                 /**< should the queue block when full */
   )
{
   SCIP_UNUSED( nthreads );
   SCIP_UNUSED( queuesize );
   SCIP_UNUSED( blockwhenfull );

   return SCIP_ERROR;
}

/** deinitializes the tpi */
SCIP_RETCODE SCIPtpiExit(
   void
   )
{
   return SCIP_ERROR;
}
