/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <assert.h>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#else
#include <sys/times.h>
#include <sys/time.h>
#endif
#include <time.h>

#include "spxdefines.h"
#include "wallclocktimer.h"

namespace soplex
{

// start timer, resume accounting user, system and real time.
void WallclockTimer::start()
{
   // ignore start request if timer is running
   if (status != RUNNING)
   {
#if !defined(_WIN32) && !defined(_WIN64)
      struct timeval tp; /*lint !e86*/
#endif
#if defined(_WIN32) || defined(_WIN64)
      sec = -::time(NULL);
#else
      gettimeofday(&tp, NULL);
      if( tp.tv_usec > usec ) /*lint !e115 !e40*/
      {
         sec = -(tp.tv_sec + 1); /*lint !e115 !e40*/
         usec = (1000000 - tp.tv_usec); /*lint !e115 !e40*/
      }
      else
      {
         sec = -tp.tv_sec; /*lint !e115 !e40*/
         usec = -tp.tv_usec; /*lint !e115 !e40*/
      }
#endif
      status = RUNNING;
   }
   lasttime = 0.0;
}

// stop timer, return accounted wallclock time.
Real WallclockTimer::stop()
{
   // status remains unchanged if timer is not running
   if (status == RUNNING)
   {
#if !defined(_WIN32) && !defined(_WIN64)
      struct timeval tp; /*lint !e86*/
#endif

#if defined(_WIN32) || defined(_WIN64)
      // we need the blank specifier to distiguish this method from WallclockTimer::time
      sec += ::time(NULL);
#else
      gettimeofday(&tp, NULL);
      if( tp.tv_usec + usec > 1000000 ) /*lint !e115 !e40*/
      {
         sec += (tp.tv_sec + 1); /*lint !e115 !e40*/
         usec -= (1000000 - tp.tv_usec); /*lint !e115 !e40*/
      }
      else
      {
         sec += tp.tv_sec; /*lint !e115 !e40*/
         usec += tp.tv_usec; /*lint !e115 !e40*/
      }
#endif
      status = STOPPED;
      lasttime = wall2sec(sec, usec);
   }
   return lasttime;
}


Real WallclockTimer::time() const
{
#if !defined(_WIN32) && !defined(_WIN64)
   struct timeval tp; /*lint !e86*/
#endif
   // only update times if timer is still running
   if( status == RUNNING )
   {
#if defined(_WIN32) || defined(_WIN64)
      // we need the blank specifier to distiguish this method from WallclockTimer::time
      lasttime = wall2sec(sec + ::time(NULL), 0);
#else
      gettimeofday(&tp, NULL);
      // check whether the microseconds add up to more than a second
      if( tp.tv_usec + usec > 1000000 ) /*lint !e115 !e40*/
         lasttime = wall2sec(sec + tp.tv_sec + 1, /*lint !e115 !e40*/
                         (usec - 1000000) + tp.tv_usec); /*lint !e115 !e40*/
      else
         lasttime = wall2sec(sec + tp.tv_sec, /*lint !e115 !e40*/
                         usec + tp.tv_usec); /*lint !e115 !e40*/
#endif
   }
   return lasttime;
}

Real WallclockTimer::lastTime() const
{
   return lasttime;
}

} // namespace soplex
