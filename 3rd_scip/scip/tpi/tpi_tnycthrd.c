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

/**@file   tpi_tnycthrd.c
 * @ingroup TASKINTERFACE
 * @brief  a TPI implementation using tinycthreads
 * @author Stephen J. Maher
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "tpi/tpi.h"
#include "blockmemshell/memory.h"

typedef struct SCIP_ThreadPool SCIP_THREADPOOL;
static SCIP_THREADPOOL* _threadpool = NULL;
_Thread_local int _threadnumber; /*lint !e129*/

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

/** The thread pool */
struct SCIP_ThreadPool
{
   /* Pool Characteristics */
   int                   nthreads;           /**< number of threads in the pool */
   int                   queuesize;          /**< the total number of items to enter the queue */

   /* Current pool state */
   thrd_t*               threads;            /**< the threads included in the pool */
   SCIP_JOBQUEUE*        jobqueue;           /**< the job queue */
   SCIP_JOBQUEUE*        currentjobs;        /**< the jobs currently being processed on a thread.
                                                  Only a single job is allowed per thread. */
   SCIP_JOBQUEUE*        finishedjobs;       /**< finished jobs that are not yet collected */
   int                   currworkingthreads; /**< the threads currently processing jobs */
   SCIP_Bool             blockwhenfull;      /**< indicates that the queue can only be as large as nthreads */
   int                   currentid;          /**< current job id */

   /* Control indicators */
   SCIP_Bool             shutdown;           /**< indicates whether the pool needs to be shutdown */
   SCIP_Bool             queueopen;          /**< indicates whether the queue is open */

   /* mutex and locks for the thread pool */
   SCIP_LOCK             poollock;           /**< mutex to allow read and write of the pool features */
   SCIP_CONDITION        queuenotempty;      /**< condition to broadcast the queue has jobs */
   SCIP_CONDITION        queuenotfull;       /**< condition to broadcast the queue is not full */
   SCIP_CONDITION        queueempty;         /**< condition to broadcast that the queue is empty */
   SCIP_CONDITION        jobfinished;        /**< condition to broadcast that a job has been finished */
};

/** this function controls the execution of each of the threads */
static
SCIP_RETCODE threadPoolThreadRetcode(
   void*                 threadnum           /**< thread number is passed in as argument stored inside a void pointer */
   )
{
   SCIP_JOB* newjob;
   SCIP_JOB* prevjob;
   SCIP_JOB* currjob;

   _threadnumber = (int)(uintptr_t) threadnum;

   /* Increase the number of active threads */
   SCIP_CALL( SCIPtpiAcquireLock(&(_threadpool->poollock)) );
   _threadpool->currworkingthreads += 1;
   SCIP_CALL( SCIPtpiReleaseLock(&(_threadpool->poollock)) );

   /* this is an endless loop that runs until the thrd_exit function is called */
   while( TRUE ) /*lint !e716*/
   {
      SCIP_CALL( SCIPtpiAcquireLock(&(_threadpool->poollock)) );

      /* the queue is empty but the shutdown command has not been given */
      while( _threadpool->jobqueue->njobs == 0 && !_threadpool->shutdown )
      {
         SCIP_CALL( SCIPtpiWaitCondition(&(_threadpool->queuenotempty), &(_threadpool->poollock)) );
      }

      /* if the shutdown command has been given, the exit the thread */
      if( _threadpool->shutdown )
      {
         /* Decrease the thread count when execution of job queue has completed */
         _threadpool->currworkingthreads -= 1;
         SCIP_CALL( SCIPtpiReleaseLock(&(_threadpool->poollock)) );

         thrd_exit((int)SCIP_OKAY);
      }

      /* getting the next job in the queue */
      newjob = _threadpool->jobqueue->firstjob;
      _threadpool->jobqueue->njobs--;            /* decreasing the number of jobs in the queue */

      if( _threadpool->jobqueue->njobs == 0 )
      {
         _threadpool->jobqueue->firstjob = NULL;
         _threadpool->jobqueue->lastjob = NULL;
      }
      else
         _threadpool->jobqueue->firstjob = newjob->nextjob;     /* updating the queue */

      /* if we want to wait when the queue is full, then we broadcast that the queue can now take new jobs */
      if( _threadpool->blockwhenfull &&
          _threadpool->jobqueue->njobs == _threadpool->queuesize - 1 )
      {
         SCIP_CALL( SCIPtpiBroadcastCondition(&(_threadpool->queuenotfull)) );
      }

      /* indicating that the queue is empty */
      if( _threadpool->jobqueue->njobs == 0 )
      {
         SCIP_CALL( SCIPtpiBroadcastCondition(&(_threadpool->queueempty)) );
      }

      /* updating the current job list */
      if( _threadpool->currentjobs->njobs == 0 )
      {
         _threadpool->currentjobs->firstjob = newjob;
         _threadpool->currentjobs->lastjob = newjob;
      }
      else
      {
         _threadpool->currentjobs->lastjob->nextjob = newjob;
         _threadpool->currentjobs->lastjob = newjob;
      }

      _threadpool->currentjobs->njobs++;

      SCIP_CALL( SCIPtpiReleaseLock(&(_threadpool->poollock)) );

      /* setting the job to run on this thread */
      newjob->retcode = (*(newjob->jobfunc))(newjob->args);

      /* setting the current job on this thread to NULL */
      SCIP_CALL( SCIPtpiAcquireLock(&(_threadpool->poollock)) );

      /* finding the location of the processed job in the currentjobs queue */
      currjob = _threadpool->currentjobs->firstjob;
      prevjob = NULL;

      while( currjob != newjob )
      {
         prevjob = currjob;
         currjob = prevjob->nextjob;
      }

      /* removing the processed job from current jobs list */
      if( currjob == _threadpool->currentjobs->firstjob )
         _threadpool->currentjobs->firstjob = currjob->nextjob;
      else
         prevjob->nextjob = currjob->nextjob; /*lint !e794*/

      if( currjob == _threadpool->currentjobs->lastjob )
         _threadpool->currentjobs->lastjob = prevjob;

      _threadpool->currentjobs->njobs--;

      /* updating the finished job list */
      if( _threadpool->finishedjobs->njobs == 0 )
      {
         _threadpool->finishedjobs->firstjob = newjob;
         _threadpool->finishedjobs->lastjob = newjob;
      }
      else
      {
         _threadpool->finishedjobs->lastjob->nextjob = newjob;
         _threadpool->finishedjobs->lastjob = newjob;
      }

      _threadpool->finishedjobs->njobs++;

      /* signalling that a job has been finished */
      SCIP_CALL( SCIPtpiBroadcastCondition(&(_threadpool)->jobfinished) );

      SCIP_CALL( SCIPtpiReleaseLock(&(_threadpool->poollock)) );
   }
}

/** this function controls the execution of each of the threads */
static
int threadPoolThread(
   void*                 threadnum           /**< thread number is passed in as argument stored inside a void pointer */
   )
{
   return (int) threadPoolThreadRetcode(threadnum);
}

/** creates a threadpool */
static
SCIP_RETCODE createThreadPool(
   SCIP_THREADPOOL**     thrdpool,           /**< pointer to store threadpool */
   int                   nthreads,           /**< number of threads in the threadpool */
   int                   qsize,              /**< maximum size of the jobqueue */
   SCIP_Bool             blockwhenfull       /**< should the jobqueue block if it is full */
   )
{
   uintptr_t i;

   assert(nthreads >= 0);
   assert(qsize >= 0);

   /* @todo think about the correct memory here */
   SCIP_ALLOC( BMSallocMemory(thrdpool) );
   (*thrdpool)->currentid = 0;
   (*thrdpool)->queuesize = qsize;
   (*thrdpool)->nthreads = nthreads;
   (*thrdpool)->blockwhenfull = blockwhenfull;
   (*thrdpool)->shutdown = FALSE;
   (*thrdpool)->queueopen = TRUE;

   /* allocating memory for the job queue */
   SCIP_ALLOC( BMSallocMemory(&(*thrdpool)->jobqueue) );
   (*thrdpool)->jobqueue->firstjob = NULL;
   (*thrdpool)->jobqueue->lastjob = NULL;
   (*thrdpool)->jobqueue->njobs = 0;

   /* allocating memory for the job queue */
   SCIP_ALLOC( BMSallocMemory(&(*thrdpool)->currentjobs) );
   (*thrdpool)->currentjobs->firstjob = NULL;
   (*thrdpool)->currentjobs->lastjob = NULL;
   (*thrdpool)->currentjobs->njobs = 0;

   /* allocating memory for the job queue */
   SCIP_ALLOC( BMSallocMemory(&(*thrdpool)->finishedjobs) );
   (*thrdpool)->finishedjobs->firstjob = NULL;
   (*thrdpool)->finishedjobs->lastjob = NULL;
   (*thrdpool)->finishedjobs->njobs = 0;

   /* initialising the mutex */
   SCIP_CALL( SCIPtpiInitLock(&(*thrdpool)->poollock) );

   /* initialising the conditions */
   SCIP_CALL( SCIPtpiInitCondition(&(*thrdpool)->queuenotempty) );
   SCIP_CALL( SCIPtpiInitCondition(&(*thrdpool)->queuenotfull) );
   SCIP_CALL( SCIPtpiInitCondition(&(*thrdpool)->queueempty) );
   SCIP_CALL( SCIPtpiInitCondition(&(*thrdpool)->jobfinished) );

   /* creating the threads */
   (*thrdpool)->currworkingthreads = 0;

   /* allocating memory for the threads */
   SCIP_ALLOC( BMSallocMemoryArray(&((*thrdpool)->threads), nthreads) );

   /* create the threads */
   for( i = 0; i < (unsigned)nthreads; i++ )
   {
      if( thrd_create(&((*thrdpool)->threads[i]), threadPoolThread, (void*)i) != thrd_success )
         return SCIP_ERROR;
   }

   _threadnumber = nthreads;
   /* halt while all threads are not active TODO: is synchronization required here ? */
   /*TODO: this caused a deadlock, is it important to wait for all threads to start?
    * while( (*thrdpool)->currworkingthreads != nthreads )
   {}*/

   return SCIP_OKAY;
}

/** adding a job to the job queue.
 * This gives some more flexibility in the handling of new jobs.
 * This function needs to be called from within a mutex. */
static
void jobQueueAddJob(
   SCIP_THREADPOOL*      threadpool,
   SCIP_JOB*             newjob
   )
{
   /* @todo we want to work out what to do with a full job queue. Is there a problem if the limit is hit? */
   /* @note it is important to have a queuesize. This will stop the code submitting infinitely many jobs. */
   assert(threadpool->jobqueue->njobs < threadpool->queuesize);

   newjob->nextjob = NULL;

   /* checking the status of the job queue */
   if( threadpool->jobqueue->njobs == 0 )
   {
      threadpool->jobqueue->firstjob = newjob;
      threadpool->jobqueue->lastjob = newjob;
   }
   else /* it is assumed that the jobqueue is not full */
   {
      threadpool->jobqueue->lastjob->nextjob = newjob;
      threadpool->jobqueue->lastjob = newjob;
   }

   SCIP_CALL_ABORT( SCIPtpiSignalCondition(&(threadpool->queuenotempty)) );    /* signalling to all threads that the queue has jobs
                                                                       * using the signal instead of broadcast because only one
                                                                       * thread should be awakened */

   threadpool->jobqueue->njobs++;
}

/** adds a job to the threadpool */
static
SCIP_RETCODE threadPoolAddWork(
   SCIP_JOB*             newjob,             /**< job to add to threadpool */
   SCIP_SUBMITSTATUS*    status              /**< pointer to store the job's submit status */
   )
{
   assert(newjob != NULL);
   assert(_threadpool != NULL);

   SCIP_CALL( SCIPtpiAcquireLock(&(_threadpool->poollock)) );

   /* if the queue is full and we are blocking, then return an error. */
   if( _threadpool->jobqueue->njobs == _threadpool->queuesize && _threadpool->blockwhenfull )
   {
      SCIP_CALL( SCIPtpiReleaseLock(&(_threadpool->poollock)) );
      *status = SCIP_SUBMIT_QUEUEFULL;
      return SCIP_OKAY;
   }

   /* Wait until the job queue is not full. If the queue is closed or the thread pool is shutdown, then stop waiting */
   /* @todo this needs to be checked. It is possible that a job can be submitted and then the queue is closed or the
    * thread pool is shutdown. Need to work out the best way to handle this. */
   while( _threadpool->jobqueue->njobs == _threadpool->queuesize && !(_threadpool->shutdown || !_threadpool->queueopen) )
   {
      SCIP_CALL( SCIPtpiWaitCondition(&(_threadpool->queuenotfull), &(_threadpool->poollock)) );
   }

   /* if the thread pool is shutdown or the queue is closed, then we need to leave the job submission */
   if( !_threadpool->queueopen )
   {
      SCIP_CALL( SCIPtpiReleaseLock(&(_threadpool->poollock)) );
      *status = SCIP_SUBMIT_QUEUECLOSED;
      return SCIP_OKAY;
   }
   else if( _threadpool->shutdown )
   {
      SCIP_CALL( SCIPtpiReleaseLock(&(_threadpool->poollock)) );
      *status =  SCIP_SUBMIT_SHUTDOWN;
      return SCIP_OKAY;
   }

   /* creating the job for submission */
   newjob->nextjob = NULL;

   /* adding the job to the queue */
   /* this can only happen if the queue is not full
    */
   assert(_threadpool->jobqueue->njobs != _threadpool->queuesize);
   jobQueueAddJob(_threadpool, newjob);

   SCIP_CALL( SCIPtpiReleaseLock(&(_threadpool->poollock)) );

   *status = SCIP_SUBMIT_SUCCESS;
   return SCIP_OKAY;
}

/** frees the jobqueue of the threadpool */
static
void freeJobQueue(
   SCIP_THREADPOOL*      thrdpool
   )
{
   SCIP_JOB* currjob;

   assert(!thrdpool->queueopen);
   assert(thrdpool->shutdown);

   /* iterating through all jobs until all have been freed. */
   while( thrdpool->jobqueue->firstjob != NULL )
   {
      currjob = thrdpool->jobqueue->firstjob->nextjob;
      thrdpool->jobqueue->firstjob = thrdpool->jobqueue->firstjob->nextjob;
      BMSfreeMemory(&currjob);
   }

   assert(thrdpool->jobqueue->firstjob == NULL);
   assert(thrdpool->jobqueue->lastjob == NULL);

   BMSfreeMemory(&thrdpool->jobqueue);
}




static
SCIP_RETCODE freeThreadPool(
   SCIP_THREADPOOL**     thrdpool,
   SCIP_Bool             finishjobs,
   SCIP_Bool             completequeue
   )
{
   int          i;
   SCIP_RETCODE retcode;

   /*TODO remove argument? */
   SCIP_UNUSED( finishjobs );

   SCIP_CALL( SCIPtpiAcquireLock(&((*thrdpool)->poollock)) );

   /* if the shutdown is already in progress, then we don't need to completed this function */
   if( !(*thrdpool)->queueopen || (*thrdpool)->shutdown )
   {
      SCIP_CALL( SCIPtpiReleaseLock(&((*thrdpool)->poollock)) );

      return SCIP_OKAY;
   }

   /* indicating that the job queue is now closed for new jobs */
   (*thrdpool)->queueopen = FALSE;

   /* if the jobs in the queue should be completed, then we wait until the queueempty condition is set */
   if( completequeue )
   {
      while( (*thrdpool)->jobqueue->njobs > 0 )
      {
         SCIP_CALL( SCIPtpiWaitCondition(&((*thrdpool)->queueempty), &((*thrdpool)->poollock)) );
      }
   }

   /* indicating that the tpi has commenced the shutdown process */
   (*thrdpool)->shutdown = TRUE;

   SCIP_CALL( SCIPtpiReleaseLock(&((*thrdpool)->poollock)) );

   /* waking up all threads so that they can check the shutdown condition
    * this requires that the conditions queuenotempty and queuenotfull is broadcast
    */
   SCIP_CALL( SCIPtpiBroadcastCondition(&((*thrdpool)->queuenotempty)) );
   SCIP_CALL( SCIPtpiBroadcastCondition(&((*thrdpool)->queuenotfull)) );

   retcode = SCIP_OKAY;

   /* calling a join to ensure that all worker finish before the thread pool is closed */
   for( i = 0; i < (*thrdpool)->nthreads; i++ )
   {
      int thrdretcode;

      if( thrd_join((*thrdpool)->threads[i], &thrdretcode) != thrd_success )
         retcode = (SCIP_RETCODE) MIN((int)SCIP_ERROR, (int)retcode);
      else
         retcode = (SCIP_RETCODE) MIN(thrdretcode, (int)retcode);
   }

   /* freeing memory and data structures */
   BMSfreeMemoryArray(&(*thrdpool)->threads);

   /* freeing the current jobs list. This assumes that all jobs complete before the tpi is closed. */
   assert((*thrdpool)->currentjobs->njobs == 0);
   BMSfreeMemory(&(*thrdpool)->currentjobs);
   assert((*thrdpool)->finishedjobs->njobs == 0);
   BMSfreeMemory(&(*thrdpool)->finishedjobs);

   freeJobQueue(*thrdpool);

   /* destroying the conditions */
   SCIPtpiDestroyCondition(&(*thrdpool)->jobfinished);
   SCIPtpiDestroyCondition(&(*thrdpool)->queueempty);
   SCIPtpiDestroyCondition(&(*thrdpool)->queuenotfull);
   SCIPtpiDestroyCondition(&(*thrdpool)->queuenotempty);

   /* destroying the mutex */
   SCIPtpiDestroyLock(&(*thrdpool)->poollock);

   BMSfreeMemory(thrdpool);

   return retcode;
}


/* checking a job queue */
static
SCIP_JOBSTATUS checkJobQueue(
   SCIP_JOBQUEUE*        jobqueue,
   int                   jobid
   )
{
   SCIP_JOB* currjob = jobqueue->firstjob;

   /* checking the job ids */
   if( currjob != NULL )
   {
      while( currjob != jobqueue->lastjob )
      {
         if( currjob->jobid == jobid )
            return SCIP_JOB_INQUEUE;

         currjob = currjob->nextjob;
      }

      if( currjob->jobid == jobid )
         return SCIP_JOB_INQUEUE;
   }

   return SCIP_JOB_DOESNOTEXIST;
}



/* returns whether the job id is running */
static
SCIP_Bool isJobRunning(
   SCIP_JOBQUEUE*        currentjobs,
   int                   jobid
   )
{
   if( checkJobQueue(currentjobs, jobid) == SCIP_JOB_INQUEUE )
      return TRUE;
   else
      return FALSE;
}

/** returns the number of threads */
int SCIPtpiGetNumThreads(
   void
   )
{
   return _threadpool->nthreads;
}

/** initializes tpi */
SCIP_RETCODE SCIPtpiInit(
   int                   nthreads,
   int                   queuesize,
   SCIP_Bool             blockwhenfull
   )
{
   assert(_threadpool == NULL);
   SCIP_CALL( createThreadPool(&_threadpool, nthreads, queuesize, blockwhenfull) );
   return SCIP_OKAY;
}

/** deinitializes tpi */
SCIP_RETCODE SCIPtpiExit(
   void
   )
{
   assert(_threadpool != NULL);

   SCIP_CALL( freeThreadPool(&_threadpool, TRUE, TRUE) );

   return SCIP_OKAY;
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
   int id;
   assert(_threadpool != NULL);

   SCIP_CALL_ABORT( SCIPtpiAcquireLock(&_threadpool->poollock) );
   id = ++_threadpool->currentid;
   SCIP_CALL_ABORT( SCIPtpiReleaseLock(&_threadpool->poollock) );

   return id;
}

/** submit a job for parallel processing
 *  the return is a globally defined status */
SCIP_RETCODE SCIPtpiSumbitJob(
   SCIP_JOB*             job,                /**< pointer to the job to be submitted */
   SCIP_SUBMITSTATUS*    status              /**< pointer to store the job's submit status */
   )
{
   assert(job != NULL);

   assert(job->jobid == _threadpool->currentid);    /* the job id must be set before submitting the job. The submitter controls
                                                       whether a new id is required. */

   SCIP_CALL( threadPoolAddWork(job, status) );

   return SCIP_OKAY;
}


/** blocks until all jobs of the given jobid have finished
 *  and then returns the smallest SCIP_RETCODE of all the jobs
 */
SCIP_RETCODE SCIPtpiCollectJobs(
   int                   jobid               /**< the jobid of the jobs to wait for */
   )
{
   SCIP_RETCODE retcode;
   SCIP_JOB* currjob;
   SCIP_JOB* prevjob;

   SCIP_CALL( SCIPtpiAcquireLock(&(_threadpool->poollock)) );

   while( isJobRunning(_threadpool->currentjobs, jobid) ||
          isJobRunning(_threadpool->jobqueue, jobid) )
   {
      SCIP_CALL( SCIPtpiWaitCondition(&_threadpool->jobfinished, &_threadpool->poollock) );
   }

   /* finding the location of the processed job in the currentjobs queue */
   retcode = SCIP_OKAY;
   currjob = _threadpool->finishedjobs->firstjob;
   prevjob = NULL;

   while( currjob )
   {
      if( currjob->jobid == jobid )
      {
         SCIP_JOB* nextjob;

         /* if the job has the right jobid collect its retcode,
          * remove it from the finished job list, and free it
          */
         retcode = MIN(retcode, currjob->retcode);

         /* removing the finished job from finished jobs list */
         if( currjob == _threadpool->finishedjobs->firstjob )
         {
            _threadpool->finishedjobs->firstjob = currjob->nextjob;
         }
         else
         {
            assert(prevjob != NULL);
            prevjob->nextjob = currjob->nextjob;
         }

         if( currjob == _threadpool->finishedjobs->lastjob )
            _threadpool->finishedjobs->lastjob = prevjob;

         _threadpool->finishedjobs->njobs--;

         /* update currjob and free finished job; prevjob stays the same */
         nextjob = currjob->nextjob;
         BMSfreeMemory(&currjob);
         currjob = nextjob;
      }
      else
      {
         /* otherwise leave job untouched */
         prevjob = currjob;
         currjob = prevjob->nextjob;
      }
   }

   SCIP_CALL( SCIPtpiReleaseLock(&_threadpool->poollock) );

   return retcode;
}
