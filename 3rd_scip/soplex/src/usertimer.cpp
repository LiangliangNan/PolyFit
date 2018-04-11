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

#include <time.h>

#else   // !(_WIN32 || _WIN64)

#include <sys/types.h>
#include <sys/times.h>
//#include <sys/param.h>
#include <unistd.h>

#endif  // !(_WIN32 || _WIN64)

#include "spxdefines.h"
#include "usertimer.h"

namespace soplex
{
/* determine TIMES_TICKS_PER_SEC for clock ticks delivered by times().
 * (don't use CLOCKS_PER_SEC since this is related to clock() only).
 */
#if defined(CLK_TCK)
#define TIMES_TICKS_PER_SEC CLK_TCK
#elif defined(_SC_CLK_TCK)
#define TIMES_TICKS_PER_SEC sysconf(_SC_CLK_TCK)
#elif defined(HZ)
#define TIMES_TICKS_PER_SEC HZ
#else // !CLK_TCK && !_SC_CLK_TCK && !HZ
#define TIMES_TICKS_PER_SEC 60
#endif // !CLK_TCK && !_SC_CLK_TCK && !HZ

const long UserTimer::ticks_per_sec = long(TIMES_TICKS_PER_SEC);

// get actual user, system and real time from system
void UserTimer::updateTicks() const
{
#if defined(_WIN32) || defined(_WIN64)

   uTicks = clock();

#else   /* !(_WIN32 || _WIN64) */

   struct tms now;
   clock_t    ret = times(&now);

   if (int(ret) == -1)
      now.tms_utime = now.tms_stime = ret = 0;

   uTicks = long(now.tms_utime);

#endif  /* !(_WIN32 || _WIN64) */
}

// start timer, resume accounting user, system and real time.
void UserTimer::start()
{
   // ignore start request if timer is running
   if (status != RUNNING)
   {
      updateTicks();

      uAccount -= uTicks;
      status    = RUNNING;
   }
   lasttime = 0;
}

// stop timer, return accounted user time.
Real UserTimer::stop()
{
   // status remains unchanged if timer is not running
   if (status == RUNNING)
   {
      updateTicks();

      uAccount += uTicks;
      status    = STOPPED;
   }
   return ticks2sec(uAccount);
}

// get accounted user time.
Real UserTimer::time() const
{
   if (status == RUNNING)
   {
      updateTicks();
      lasttime = ticks2sec(uTicks + uAccount);
   }
   else
   {
      lasttime = ticks2sec(uAccount);
   }
   return lasttime;
}

Real UserTimer::lastTime() const
{
   return lasttime;
}

} // namespace soplex
