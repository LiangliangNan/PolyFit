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

/**@file   tpi_openmp.c
 * @ingroup TASKINTERFACE
 * @brief  the interface functions for openmp
 * @author Stephen J. Maher
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "tpi/tpi.h"
#include "blockmemshell/memory.h"

/** A job added to the queue */
struct SCIP_Job
{
   int                   jobid;              /**< id to identify jobs from a common process */
   struct                SCIP_Job* nextjob;  /**< pointer to the next job in the queue */
   SCIP_RETCODE          (*jobfunc)(void* args);/**< pointer to the job function */
   void*                 args;               /**< pointer to the function arguements */
   SCIP_RETCODE          retcode;            /**< return code of the job */
};

/** the thread pool job queue */
struct SCIP_JobQueue
{
   SCIP_JOB*             firstjob;           /**< pointer to the first job in the queue */
   SCIP_JOB*             lastjob;            /**< pointer to the last job in the queue */
   int                   njobs;              /**< number of jobs in the queue */
};
typedef struct SCIP_JobQueue SCIP_JOBQUEUE;

struct SCIP_JobQueues
{
   SCIP_JOBQUEUE         jobqueue;           /**< queue of unprocessed jobs */
   SCIP_JOB**            currentjobs;        /**< array with slot for each thread to store the currently running job */
   int                   ncurrentjobs;       /**< number of currently running jobs */
   int                   nthreads;           /**< number of threads */
   SCIP_JOBQUEUE         finishedjobs;       /**< jobqueue containing the finished jobs */
   SCIP_LOCK             lock;               /**< lock to protect this stucture from concurrent access */
   SCIP_CONDITION        jobfinished;        /**< condition to signal if a job was finished */
};
typedef struct SCIP_JobQueues SCIP_JOBQUEUES;

static SCIP_JOBQUEUES* _jobqueues = NULL;



static
SCIP_RETCODE createJobQueue(
   int                   nthreads,           /**< the number of threads */
   int                   qsize,              /**< the queue size */
   SCIP_Bool             blockwhenfull       /**< should the queue be blocked from new jobs when full */
   )
{
   int i;

   assert(nthreads >= 0);
   assert(qsize >= 0);
   SCIP_UNUSED( blockwhenfull );

   /* allocting memory for the job queue */
   SCIP_ALLOC( BMSallocMemory(&_jobqueues) );
   _jobqueues->jobqueue.firstjob = NULL;
   _jobqueues->jobqueue.lastjob = NULL;
   _jobqueues->jobqueue.njobs = 0;
   _jobqueues->finishedjobs.firstjob = NULL;
   _jobqueues->finishedjobs.lastjob = NULL;
   _jobqueues->finishedjobs.njobs = 0;
   _jobqueues->ncurrentjobs = 0;

   _jobqueues->nthreads = nthreads;
   SCIP_ALLOC( BMSallocMemoryArray(&_jobqueues->currentjobs, nthreads) );

   for( i = 0; i < nthreads; ++i )
      _jobqueues->currentjobs[i] = NULL;

   SCIP_CALL( SCIPtpiInitLock(&_jobqueues->lock) );
   SCIP_CALL( SCIPtpiInitCondition(&_jobqueues->jobfinished) );

   return SCIP_OKAY;
}



static
SCIP_RETCODE freeJobQueue(
   void
   )
{
   assert(_jobqueues != NULL);

   SCIPtpiDestroyLock(&_jobqueues->lock);
   SCIPtpiDestroyCondition(&_jobqueues->jobfinished);
   BMSfreeMemoryArray(&_jobqueues->currentjobs);

   BMSfreeMemory(&_jobqueues);

   return SCIP_OKAY;
}


static
void executeJob(
   SCIP_JOB*             job                 /**< the job to be executed in parallel */
   )
{
   int threadnum;

   threadnum = SCIPtpiGetThreadNum();

   SCIP_CALL_ABORT( SCIPtpiAcquireLock(&_jobqueues->lock) );
   _jobqueues->currentjobs[threadnum] = job;
   SCIP_CALL_ABORT( SCIPtpiReleaseLock(&_jobqueues->lock) );

   job->retcode = (*(job->jobfunc))(job->args);

   SCIP_CALL_ABORT( SCIPtpiAcquireLock(&_jobqueues->lock) );
   _jobqueues->ncurrentjobs--;
   _jobqueues->currentjobs[threadnum] = NULL;

   /* insert job into finished jobs */
   if( _jobqueues->finishedjobs.njobs == 0 )
   {
      _jobqueues->finishedjobs.firstjob = job;
      _jobqueues->finishedjobs.lastjob = job;
   }
   else
   {
      _jobqueues->finishedjobs.lastjob->nextjob = job;
      _jobqueues->finishedjobs.lastjob = job;
   }

   ++_jobqueues->finishedjobs.njobs;

   SCIP_CALL_ABORT( SCIPtpiBroadcastCondition(&_jobqueues->jobfinished) );

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(&_jobqueues->lock) );

}




/** this is a job that will be executed on to process the job queue */
/* the job will only be added when the number of active jobs is equal to the number of threads.
 * As such, there will always be number of threads + 1 tasks available for the scheduler to run. */
static
void jobQueueProcessJob(
   void
   )
{
   SCIP_JOB* job;
   SCIP_CALL_ABORT( SCIPtpiAcquireLock(&_jobqueues->lock) );

   while( _jobqueues->ncurrentjobs == SCIPtpiGetNumThreads() )
   {
      SCIP_CALL_ABORT( SCIPtpiWaitCondition(&_jobqueues->jobfinished, &_jobqueues->lock) );
   }

   if( _jobqueues->jobqueue.njobs == 1 )
   {
      job = _jobqueues->jobqueue.firstjob;
      _jobqueues->jobqueue.firstjob = NULL;
      _jobqueues->jobqueue.lastjob = NULL;
      --_jobqueues->jobqueue.njobs;
   }
   else if( _jobqueues->jobqueue.njobs > 1 )
   {
      job = _jobqueues->jobqueue.firstjob;
      _jobqueues->jobqueue.firstjob = job->nextjob;
      --_jobqueues->jobqueue.njobs;
   }
   else
   {
      job = NULL;
   }

   ++_jobqueues->ncurrentjobs;
   SCIP_CALL_ABORT( SCIPtpiReleaseLock(&_jobqueues->lock) );

   if( job )
   {
      executeJob(job);
   }
}




/** adding a job to the job queue.
 * This gives some more flexibility in the handling of new jobs.
 * IMPORTANT: This function MUST be called from within a mutex. */
static
SCIP_RETCODE jobQueueAddJob(
   SCIP_JOB*             newjob
   )
{
   /* @todo we want to work out what to do with a full job queue. Is there a problem if the limit is hit? */
   /* @note it is important to have a queuesize. This will stop the code submitting infinitely many jobs. */
   assert(newjob != NULL);

   newjob->nextjob = NULL;

   /* this function queries the current job list. This could change by other threads writing to the list. So a lock is
    * required to ensure that the current joblist remains static. */
   SCIP_CALL( SCIPtpiAcquireLock(&_jobqueues->lock) );

   /* checking the status of the job queue */
   if( _jobqueues->ncurrentjobs == SCIPtpiGetNumThreads() )
   {
      if( _jobqueues->jobqueue.njobs == 0 )
      {
         _jobqueues->jobqueue.firstjob = newjob;
         _jobqueues->jobqueue.lastjob = newjob;
      }
      else /* it is assumed that the jobqueue is not full */
      {
         _jobqueues->jobqueue.lastjob->nextjob = newjob;
         _jobqueues->jobqueue.lastjob = newjob;
      }

      _jobqueues->jobqueue.njobs++;

      SCIP_CALL( SCIPtpiReleaseLock(&_jobqueues->lock) );

      #pragma omp task
      jobQueueProcessJob();
   }
   else
   {
      assert(_jobqueues->ncurrentjobs < SCIPtpiGetNumThreads());

      _jobqueues->ncurrentjobs++;

      SCIP_CALL( SCIPtpiReleaseLock(&_jobqueues->lock) );
      /* running the new job */
      #pragma omp task firstprivate(newjob)
      executeJob(newjob);
   }

   return SCIP_OKAY;
}


SCIP_RETCODE SCIPtpiSignalCondition(
   SCIP_CONDITION* condition
   )
{
   SCIP_CALL( SCIPtpiAcquireLock(&condition->_lock) );

   if( condition->_waitnum > condition->_signals )
      ++condition->_signals;

   SCIP_CALL( SCIPtpiReleaseLock(&condition->_lock) );

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPtpiBroadcastCondition(
   SCIP_CONDITION* condition
   )
{

   SCIP_CALL( SCIPtpiAcquireLock(&condition->_lock) );
   condition->_signals = condition->_waitnum;
   SCIP_CALL( SCIPtpiReleaseLock(&condition->_lock) );

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPtpiWaitCondition(
   SCIP_CONDITION* condition,
   SCIP_LOCK*      lock
   )
{
   int waitnum;

   SCIP_CALL( SCIPtpiReleaseLock(lock) );

   SCIP_CALL( SCIPtpiAcquireLock(&condition->_lock) );
   waitnum = ++condition->_waitnum;

   ++condition->_waiters;

   do
   {
      SCIP_CALL( SCIPtpiReleaseLock(&condition->_lock) );
      #pragma omp taskyield
      SCIP_CALL( SCIPtpiAcquireLock(&condition->_lock) );
   }
   while( condition->_signals < waitnum );

   --condition->_waiters;

   if( condition->_waiters == 0 )
   {
      condition->_signals = 0;
      condition->_waitnum = 0;
   }

   SCIP_CALL( SCIPtpiReleaseLock(&condition->_lock) );

   SCIP_CALL( SCIPtpiAcquireLock(lock) );

   return SCIP_OKAY;
}

/** Returns the number of threads */
int SCIPtpiGetNumThreads(
   )
{
   return omp_get_num_threads();
}

/** Returns the thread number */
int SCIPtpiGetThreadNum(
   )
{
   return omp_get_thread_num();
}

/** creates a job for parallel processing*/
SCIP_RETCODE SCIPtpiCreateJob(
   SCIP_JOB**            job,                /**< pointer to the job that will be created */
   int                   jobid,              /**< the id for the current job */
   SCIP_RETCODE          (*jobfunc)(void* args),/**< pointer to the job function */
   void*                 jobarg              /**< the job's argument */
   )
{
   SCIP_ALLOC( BMSallocMemory(job) );

   (*job)->jobid = jobid;
   (*job)->jobfunc = jobfunc;
   (*job)->args = jobarg;
   (*job)->nextjob = NULL;

   return SCIP_OKAY;
}

/** get a new job id for the new set of submitted jobs */
int SCIPtpiGetNewJobID(
   void
   )
{
   static int currentjobid = 0;
   int jobid;

   #pragma omp atomic capture
   jobid = ++currentjobid;

   return jobid;
}

/** submit a job for parallel processing */
/* the return is a globally defined status */
SCIP_RETCODE SCIPtpiSumbitJob(
   SCIP_JOB*             job,                /**< pointer to the job to be submitted */
   SCIP_SUBMITSTATUS*    status              /**< pointer to store the submit status */

   )
{
   assert(_jobqueues != NULL);

   *status = SCIP_SUBMIT_SUCCESS;
   SCIP_CALL( jobQueueAddJob(job) );

   return SCIP_OKAY;
}

static
SCIP_Bool isJobRunning(
   int                   jobid
   )
{
   int i;

   if( _jobqueues->ncurrentjobs > 0 )
   {
      for( i = 0; i < _jobqueues->nthreads; ++i )
      {
         if( _jobqueues->currentjobs[i] != NULL && _jobqueues->currentjobs[i]->jobid == jobid )
            return TRUE;
      }
   }

   return FALSE;
}

static
SCIP_Bool isJobWaiting(
   int                   jobid
   )
{

   if( _jobqueues->jobqueue.njobs > 0 )
   {
      SCIP_JOB* currjob;
      currjob = _jobqueues->jobqueue.firstjob;

      do
      {
         if( currjob->jobid == jobid )
            return TRUE;

         if( currjob == _jobqueues->jobqueue.lastjob )
            break;

         currjob = currjob->nextjob;
      }
      while( TRUE ); /*lint !e506*/
   }

   return FALSE;
}


/** Blocks until all jobs of the given jobid have finished
 * and then returns the smallest SCIP_RETCODE of all the jobs */
SCIP_RETCODE SCIPtpiCollectJobs(
   int                   jobid
   )
{
   SCIP_RETCODE retcode;

   retcode = SCIP_OKAY;
   SCIP_CALL( SCIPtpiAcquireLock(&_jobqueues->lock) );

   while( isJobRunning(jobid) || isJobWaiting(jobid) )
   {
      SCIP_CALL( SCIPtpiWaitCondition(&_jobqueues->jobfinished, &_jobqueues->lock) );
   }

   if( _jobqueues->finishedjobs.njobs > 0 )
   {
      SCIP_JOB* currjob = _jobqueues->finishedjobs.firstjob;
      SCIP_JOB* prevjob = NULL;

      /* finding the location of the processed job in the currentjobs queue */
      do
      {
         if( currjob->jobid == jobid )
         {
            SCIP_JOB* nextjob;

            /** if the job has the right jobid collect its retcode,
             *  remove it from the finished job list, and free it */
            retcode = MIN(retcode, currjob->retcode);

            /* removing the finished job from finished jobs list */
            if( currjob == _jobqueues->finishedjobs.firstjob )
               _jobqueues->finishedjobs.firstjob = currjob->nextjob;
            else
               prevjob->nextjob = currjob->nextjob; /*lint !e613*/

            if( currjob == _jobqueues->finishedjobs.lastjob )
               _jobqueues->finishedjobs.lastjob = prevjob;

            _jobqueues->finishedjobs.njobs--;

            /* update currjob and free finished job; prevjob stays the same */
            nextjob = currjob->nextjob;
            BMSfreeMemory(&currjob);
            currjob = nextjob;
         }
         else
         {
            prevjob = currjob;
            currjob = prevjob->nextjob;
         }
      }
      while( prevjob != _jobqueues->finishedjobs.lastjob );

   }
   else
   {
      /* given jobid was not submitted */
      printf("err1");
      retcode = SCIP_ERROR;
   }

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(&_jobqueues->lock) );

   return retcode;
}

/** initializes tpi */
SCIP_RETCODE SCIPtpiInit(
   int         nthreads,
   int         queuesize,
   SCIP_Bool   blockwhenfull
   )
{
   omp_set_num_threads(nthreads);
   assert(_jobqueues == NULL);

   SCIP_CALL( createJobQueue(nthreads, queuesize, blockwhenfull) );
   return SCIP_OKAY;
}

/** deinitializes tpi */
SCIP_RETCODE SCIPtpiExit(
   void
   )
{
   assert(_jobqueues != NULL);
   assert(_jobqueues->finishedjobs.njobs == 0);
   assert(_jobqueues->jobqueue.njobs == 0);
   assert(_jobqueues->ncurrentjobs == 0);

   SCIP_CALL( freeJobQueue() );
   return SCIP_OKAY;
}
