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

/**@file   clock.c
 * @brief  methods for clocks and timing issues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#else
#include <sys/times.h>
#include <sys/time.h>
#include <unistd.h>
#endif
#include <time.h>

#include "scip/def.h"
#include "scip/pub_message.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/clock.h"

#include "scip/struct_clock.h"

/** converts CPU clock ticks into seconds */
static
SCIP_Real cputime2sec(
   clock_t               cputime             /**< clock ticks for CPU time */
   )
{
   clock_t clocks_per_second;

#if defined(_WIN32) || defined(_WIN64)
   clocks_per_second = 100;
#else
#ifndef CLK_TCK
   clocks_per_second = sysconf(_SC_CLK_TCK);
#else
   clocks_per_second = CLK_TCK;
#endif
#endif

   return (SCIP_Real)cputime / (SCIP_Real)clocks_per_second;
}

/*lint -esym(*,timeval)*/
/*lint -esym(*,gettimeofday)*/

/** converts wall clock time into seconds */
static
SCIP_Real walltime2sec(
   long                  sec,                /**< seconds counter */
   long                  usec                /**< microseconds counter */
   )
{
   return (SCIP_Real)sec + 0.000001 * (SCIP_Real)usec;
}

/** converts seconds into CPU clock ticks */
static
void sec2cputime(
   SCIP_Real             sec,                /**< seconds */
   clock_t*              cputime             /**< pointer to store clock ticks for CPU time */
   )
{
   clock_t clocks_per_second;

   assert(cputime != NULL);

#if defined(_WIN32) || defined(_WIN64)
   clocks_per_second = 100;
#else
#ifndef CLK_TCK
   clocks_per_second = sysconf(_SC_CLK_TCK);
#else
   clocks_per_second = CLK_TCK;
#endif
#endif
   *cputime = (clock_t)(sec * clocks_per_second);
}

/** converts wall clock time into seconds */
static
void sec2walltime(
   SCIP_Real             sec,                /**< seconds */
   long*                 wallsec,            /**< pointer to store seconds counter */
   long*                 wallusec            /**< pointer to store microseconds counter */
   )
{
   assert(wallsec != NULL);
   assert(wallusec != NULL);

   *wallsec = (long)sec;
   *wallusec = (long)((sec  - *wallsec) * 1000000.0);
}


/** sets the clock's type and converts the clock timer accordingly */
static
void clockSetType(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_CLOCKTYPE        newtype             /**< new clock type */
   )
{
   assert(clck != NULL);
   assert(newtype != SCIP_CLOCKTYPE_DEFAULT);

   if( clck->clocktype != newtype )
   {
      if( clck->clocktype == SCIP_CLOCKTYPE_DEFAULT )
      {
         assert(clck->nruns == 0);
         clck->clocktype = newtype;
         SCIPclockReset(clck);
         SCIPdebugMessage("switched clock type to %d\n", newtype);
      }
      else
      {
         SCIP_Real sec;

         sec = SCIPclockGetTime(clck);
         clck->clocktype = newtype;
         SCIPclockSetTime(clck, sec);
         SCIPdebugMessage("switched clock type to %d (%g seconds -> %g seconds)\n", newtype, sec, SCIPclockGetTime(clck));
      }
   }
}

/** if the clock uses the default clock type and the default changed, converts the clock timer to the new type */
static
void clockUpdateDefaultType(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_CLOCKTYPE        defaultclocktype    /**< default type of clock to use */
   )
{
   assert(clck != NULL);
   assert(defaultclocktype != SCIP_CLOCKTYPE_DEFAULT);

   if( clck->usedefault && clck->clocktype != defaultclocktype )
      clockSetType(clck, defaultclocktype);
}

/** creates a clock and initializes it */
SCIP_RETCODE SCIPclockCreate(
   SCIP_CLOCK**          clck,               /**< pointer to clock timer */
   SCIP_CLOCKTYPE        clocktype           /**< type of clock */
   )
{
   assert(clck != NULL);

   SCIP_ALLOC( BMSallocMemory(clck) );

   SCIPclockInit(*clck, clocktype);

   return SCIP_OKAY;
}

/** frees a clock */
void SCIPclockFree(
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   )
{
   assert(clck != NULL);

   BMSfreeMemory(clck);
}

/** initializes and resets a clock */
void SCIPclockInit(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_CLOCKTYPE        clocktype           /**< type of clock */
   )
{
   assert(clck != NULL);

   SCIPdebugMessage("initializing clock %p of type %d\n", (void*)clck, clocktype);
   clck->enabled = TRUE;
   clck->lasttime = 0.0;
   SCIPclockSetType(clck, clocktype);
}

/** completely stop the clock and reset the clock's counter to zero */
void SCIPclockReset(
   SCIP_CLOCK*           clck                /**< clock timer */
   )
{
   assert(clck != NULL);

   SCIPdebugMessage("resetting clock %p of type %d (usedefault=%u)\n", (void*)clck, clck->clocktype, clck->usedefault);
   switch( clck->clocktype )
   {
   case SCIP_CLOCKTYPE_DEFAULT:
      break;
   case SCIP_CLOCKTYPE_CPU:
      clck->data.cpuclock.user = 0;
      break;
   case SCIP_CLOCKTYPE_WALL:
      clck->data.wallclock.sec = 0;
      clck->data.wallclock.usec = 0;
      break;
   default:
      SCIPerrorMessage("invalid clock type\n");
      SCIPABORT();
   }
   clck->nruns = 0;
}

/** enables the clock */
void SCIPclockEnable(
   SCIP_CLOCK*           clck                /**< clock timer */
   )
{
   assert(clck != NULL);

   SCIPdebugMessage("enabling clock %p of type %d (usedefault=%u)\n", (void*)clck, clck->clocktype, clck->usedefault);

   clck->enabled = TRUE;
}

/** disables and resets the clock */
void SCIPclockDisable(
   SCIP_CLOCK*           clck                /**< clock timer */
   )
{
   assert(clck != NULL);

   SCIPdebugMessage("disabling clock %p of type %d (usedefault=%u)\n", (void*)clck, clck->clocktype, clck->usedefault);

   clck->enabled = FALSE;
   SCIPclockReset(clck);
}

/** enables or disables \p clck, depending on the value of the flag */
void SCIPclockEnableOrDisable(
   SCIP_CLOCK*           clck,               /**< the clock to be disabled/enabled */
   SCIP_Bool             enable              /**< should the clock be enabled? */
   )
{
   assert(clck != NULL);

   if( enable )
      SCIPclockEnable(clck);
   else
      SCIPclockDisable(clck);
}

/** sets the type of the clock, overriding the default clock type, and resets the clock */
void SCIPclockSetType(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_CLOCKTYPE        clocktype           /**< type of clock */
   )
{
   assert(clck != NULL);

   SCIPdebugMessage("setting type of clock %p (type %d, usedefault=%u) to %d\n", 
      (void*)clck, clck->clocktype, clck->usedefault, clocktype);

   clck->clocktype = clocktype;
   clck->usedefault = (clocktype == SCIP_CLOCKTYPE_DEFAULT);
   SCIPclockReset(clck);
}

/** starts measurement of time in the given clock */
void SCIPclockStart(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(clck != NULL);
   assert(set != NULL);

   if( set->time_enabled && clck->enabled )
   {
      clockUpdateDefaultType(clck, set->time_clocktype);

      if( clck->nruns == 0 )
      {
#if defined(_WIN32) || defined(_WIN64)
         FILETIME creationtime;
         FILETIME exittime;
         FILETIME kerneltime;
         FILETIME usertime;
#else
         struct timeval tp; /*lint !e86*/
         struct tms now;
#endif

         SCIPdebugMessage("starting clock %p (type %d, usedefault=%u)\n", (void*)clck, clck->clocktype, clck->usedefault);

         switch( clck->clocktype )
         {
         case SCIP_CLOCKTYPE_CPU:
#if defined(_WIN32) || defined(_WIN64)
            GetProcessTimes(GetCurrentProcess(), &creationtime, &exittime, &kerneltime, &usertime);
            clck->data.cpuclock.user -= usertime.dwHighDateTime * 42950 + usertime.dwLowDateTime / 100000L;
#else
            (void)times(&now);
            clck->data.cpuclock.user -= now.tms_utime;
#endif
            clck->lasttime = cputime2sec(clck->data.cpuclock.user);
            break;

         case SCIP_CLOCKTYPE_WALL:
#if defined(_WIN32) || defined(_WIN64)
            clck->data.wallclock.sec -= time(NULL);
#else
            gettimeofday(&tp, NULL);
            if( tp.tv_usec > clck->data.wallclock.usec ) /*lint !e115 !e40*/
            {
               clck->data.wallclock.sec -= (tp.tv_sec + 1); /*lint !e115 !e40*/
               clck->data.wallclock.usec += (1000000 - tp.tv_usec); /*lint !e115 !e40*/
            }
            else
            {
               clck->data.wallclock.sec -= tp.tv_sec; /*lint !e115 !e40*/
               clck->data.wallclock.usec -= tp.tv_usec; /*lint !e115 !e40*/
            }
#endif
            clck->lasttime = walltime2sec(clck->data.wallclock.sec, clck->data.wallclock.usec);
            break;

         case SCIP_CLOCKTYPE_DEFAULT:
         default:
            SCIPerrorMessage("invalid clock type\n");
            SCIPABORT();
         }
      }

      clck->nruns++;
   }
}

/** stops measurement of time in the given clock */
void SCIPclockStop(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(clck != NULL);
   assert(set != NULL);

   if( set->time_enabled && clck->enabled )
   {
      assert(clck->nruns >= 1);

      clck->nruns--;
      if( clck->nruns == 0 )
      {
#if defined(_WIN32) || defined(_WIN64)
         FILETIME creationtime;
         FILETIME exittime;
         FILETIME kerneltime;
         FILETIME usertime;
#else
         struct timeval tp; /*lint !e86*/
         struct tms now;
#endif

         SCIPdebugMessage("stopping clock %p (type %d, usedefault=%u)\n", (void*)clck, clck->clocktype, clck->usedefault);

         switch( clck->clocktype )
         {
         case SCIP_CLOCKTYPE_CPU:
#if defined(_WIN32) || defined(_WIN64)
            GetProcessTimes(GetCurrentProcess(), &creationtime, &exittime, &kerneltime, &usertime);
            clck->data.cpuclock.user += usertime.dwHighDateTime * 42950 + usertime.dwLowDateTime / 100000L;
#else
            (void)times(&now);
            clck->data.cpuclock.user += now.tms_utime;
#endif
            break;

         case SCIP_CLOCKTYPE_WALL:
#if defined(_WIN32) || defined(_WIN64)
            clck->data.wallclock.sec += time(NULL);
#else
            gettimeofday(&tp, NULL);
            if( tp.tv_usec + clck->data.wallclock.usec > 1000000 ) /*lint !e115 !e40*/
            {
               clck->data.wallclock.sec += (tp.tv_sec + 1); /*lint !e115 !e40*/
               clck->data.wallclock.usec -= (1000000 - tp.tv_usec); /*lint !e115 !e40*/
            }
            else
            {
               clck->data.wallclock.sec += tp.tv_sec; /*lint !e115 !e40*/
               clck->data.wallclock.usec += tp.tv_usec; /*lint !e115 !e40*/
            }
#endif
            break;

         case SCIP_CLOCKTYPE_DEFAULT:
         default:
            SCIPerrorMessage("invalid clock type\n");
            SCIPABORT();
         }
      }
   }
}

/** returns whether the clock is currently running */
SCIP_Bool SCIPclockIsRunning(
   SCIP_CLOCK*           clck                /**< clock timer */
   )
{
   assert(clck != NULL);

   return (clck->nruns > 0);
}


/** gets the used time of this clock in seconds */
SCIP_Real SCIPclockGetTime(
   SCIP_CLOCK*           clck                /**< clock timer */
   )
{
   SCIP_Real result;
   assert(clck != NULL);
   result = 0.0;

   SCIPdebugMessage("getting time of clock %p (type %d, usedefault=%u, nruns=%d)\n",
      (void*)clck, clck->clocktype, clck->usedefault, clck->nruns);

   if( clck->nruns == 0 )
   {
      /* the clock is not running: convert the clocks timer into seconds */
      switch( clck->clocktype )
      {
      case SCIP_CLOCKTYPE_DEFAULT:
         break;
      case SCIP_CLOCKTYPE_CPU:
         result = cputime2sec(clck->data.cpuclock.user);
         break;
      case SCIP_CLOCKTYPE_WALL:
         result = walltime2sec(clck->data.wallclock.sec, clck->data.wallclock.usec);
         break;
      default:
         SCIPerrorMessage("invalid clock type\n");
         SCIPABORT();
         result = 0.0; /*lint !e527*/
      }
   }
   else
   {
#if defined(_WIN32) || defined(_WIN64)
      FILETIME creationtime;
      FILETIME exittime;
      FILETIME kerneltime;
      FILETIME usertime;
#else
      struct timeval tp; /*lint !e86*/
      struct tms now;
#endif

      /* the clock is currently running: we have to add the current time to the clocks timer */
      switch( clck->clocktype )
      {
      case SCIP_CLOCKTYPE_CPU:
#if defined(_WIN32) || defined(_WIN64)
          GetProcessTimes(GetCurrentProcess(), &creationtime, &exittime, &kerneltime, &usertime);
          result = cputime2sec(clck->data.cpuclock.user + usertime.dwHighDateTime * 42950 + usertime.dwLowDateTime / 100000L);
#else
         (void)times(&now);
         result = cputime2sec(clck->data.cpuclock.user + now.tms_utime);
#endif
         break;
      case SCIP_CLOCKTYPE_WALL:
#if defined(_WIN32) || defined(_WIN64)
         result = walltime2sec(clck->data.wallclock.sec + time(NULL), 0);
#else
         gettimeofday(&tp, NULL);
         if( tp.tv_usec + clck->data.wallclock.usec > 1000000 ) /*lint !e115 !e40*/
            result = walltime2sec(clck->data.wallclock.sec + tp.tv_sec + 1, /*lint !e115 !e40*/
               (clck->data.wallclock.usec - 1000000) + tp.tv_usec); /*lint !e115 !e40*/
         else
            result = walltime2sec(clck->data.wallclock.sec + tp.tv_sec, /*lint !e115 !e40*/
               clck->data.wallclock.usec + tp.tv_usec); /*lint !e115 !e40*/
#endif
         break;
      case SCIP_CLOCKTYPE_DEFAULT:
      default:
         SCIPerrorMessage("invalid clock type\n");
         SCIPABORT();
         result = 0.0; /*lint !e527*/
      }
   }

   clck->lasttime = result;
   return result;
}

/** gets the last validated time of this clock in seconds */
SCIP_Real SCIPclockGetLastTime(
   SCIP_CLOCK*           clck                /**< clock timer */
   )
{
   assert(clck != NULL);

   return clck->lasttime;
}

/** sets the used time of this clock in seconds */
void SCIPclockSetTime(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_Real             sec                 /**< time in seconds to set the clock's timer to */
   )
{
   assert(clck != NULL);

   SCIPdebugMessage("setting time of clock %p (type %d, usedefault=%u, nruns=%d) to %g\n",
      (void*)clck, clck->clocktype, clck->usedefault, clck->nruns, sec);

   /* if the clock type is not yet set, set it to an arbitrary value to be able to store the number */
   if( clck->clocktype == SCIP_CLOCKTYPE_DEFAULT )
      clockSetType(clck, SCIP_CLOCKTYPE_WALL);

   switch( clck->clocktype )
   {
   case SCIP_CLOCKTYPE_CPU:
      sec2cputime(sec, &clck->data.cpuclock.user);
      break;

   case SCIP_CLOCKTYPE_WALL:
      sec2walltime(sec, &clck->data.wallclock.sec, &clck->data.wallclock.usec);
      break;

   case SCIP_CLOCKTYPE_DEFAULT:
   default:
      SCIPerrorMessage("invalid clock type\n");
      SCIPABORT();
   }

   if( clck->nruns >= 1 )
   {
#if defined(_WIN32) || defined(_WIN64)
      FILETIME creationtime;
      FILETIME exittime;
      FILETIME kerneltime;
      FILETIME usertime;
#else
      struct timeval tp; /*lint !e86*/
      struct tms now;
#endif

      /* the clock is currently running: we have to subtract the current time from the new timer value */
      switch( clck->clocktype )
      {
      case SCIP_CLOCKTYPE_CPU:
#if defined(_WIN32) || defined(_WIN64)
         GetProcessTimes(GetCurrentProcess(), &creationtime, &exittime, &kerneltime, &usertime);
         clck->data.cpuclock.user -= usertime.dwHighDateTime * 42950 + usertime.dwLowDateTime / 100000L;
#else
         (void)times(&now);
         clck->data.cpuclock.user -= now.tms_utime;
#endif
         break;

      case SCIP_CLOCKTYPE_WALL:
#if defined(_WIN32) || defined(_WIN64)
         clck->data.wallclock.sec -= time(NULL);
#else
         gettimeofday(&tp, NULL);
         if( tp.tv_usec > clck->data.wallclock.usec ) /*lint !e115 !e40*/
         {
            clck->data.wallclock.sec -= (tp.tv_sec + 1); /*lint !e115 !e40*/
            clck->data.wallclock.usec += (1000000 - tp.tv_usec); /*lint !e115 !e40*/
         }
         else
         {
            clck->data.wallclock.sec -= tp.tv_sec; /*lint !e115 !e40*/
            clck->data.wallclock.usec -= tp.tv_usec; /*lint !e115 !e40*/
         }
#endif
         break;

      case SCIP_CLOCKTYPE_DEFAULT:
      default:
         SCIPerrorMessage("invalid clock type\n");
         SCIPABORT();
      }
   }
}

/** gets current time of day in seconds (standard time zone) */
SCIP_Real SCIPclockGetTimeOfDay(
   void
   )
{
#if defined(_WIN32) || defined(_WIN64)
   time_t now;
   now = time(NULL);
   return (SCIP_Real)(now % (24*3600));
#else
   struct timeval tp; /*lint !e86*/

   gettimeofday(&tp, NULL);

   return (SCIP_Real)(tp.tv_sec % (24*3600)) + (SCIP_Real)tp.tv_usec / 1e+6; /*lint !e40 !e115*/
#endif
}
