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

/**@file  timer.h
 * @brief Timer class.
 */

#ifndef _TIMER_H_
#define _TIMER_H_

#include "spxdefines.h"

namespace soplex
{
/**@class   Timer
   @ingroup Elementary

   @brief Wrapper for the system time query methods.

    In C or C++ programs, the usual way to measure time intervals,
    e.g., running times of some complex computations, is to call one
    of the provided system functions like %clock(), %time(), %times(),
    %gettimeofday(), %getrusage() etc.  By these functions one can
    gather information about the process' user and system time and the
    system clock (real time).

    Unfortunately, these functions are rather clumsy.  The programmer
    determines computation times by querying a (virtual) clock value
    at the beginning and another one at the end of some computation
    and converting the difference of these values into seconds.  Some
    functions impose restrictions; for instance, the values of
    the ANSI C function %clock() are of high resolution but will wrap
    around after about 36 minutes (cpu time).  Most timing functions
    take some data structure as argument that has to be allocated
    before the call and from which the user has to pick up the
    information of interest after the call.  Problems can arise when
    porting programs to other operating systems that do not support
    standards like POSIX etc.

    In order to simplify measuring computation times and to hide the
    system-dependencies involved, a concept of \em timers accounting the
    process' system and real time is implemented.  C and C++ interfaces
    are provided as a set of functions operating on timers and a timer class
    respectively.

    Look into the file timerfactory.h to see how to switch between different
    timing types or to disable timing altogether.

    The idea is to provide a type Timer for objects that act like a stopwatch.
    Operations on such an objects include: start accounting time, stop
    accounting, read the actual time account and reset the objects time account
    to zero.

    After initialization, accounting for the time can be
    started by calling a function start(). Accounting is suspended by calling
    a function stop() and can be resumed at any time by calling start()
    again.

    For convenience, the actually accounted user time is returned by stop()
    too.  Function reset() re-initializes a timer clearing all time
    accounts.

*/
class Timer
{
protected:

   //------------------------------------
   /**@name Types */
   //@{
   /// status of the timer
   enum
   {
      RESET,                   ///< reset
      STOPPED,                 ///< stopped
      RUNNING                  ///< running
   } status;                   ///< timer status

   //@}

public:

   //------------------------------------
   /**@name Timers */
   //@{
   /// types of timers
   typedef enum
   {
      OFF = 0,
      USER_TIME = 1,
      WALLCLOCK_TIME = 2
   } TYPE;
   //@}

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   Timer()
      : status(RESET)
   {}
   /// copy constructor
   Timer(const Timer& old)
      : status(old.status)
   {}
   /// assignment operator
   Timer& operator=(const Timer& old)
   {
      status = old.status;
      return *this;
   }
   virtual ~Timer()
   {}
   //@}

   //------------------------------------
   /**@name Control */
   //@{
   /// initialize timer, set timing accounts to zero.
   virtual void reset() = 0;

   /// start timer, resume accounting user, system and real time.
   virtual void start() = 0;

   /// stop timer, return accounted user time.
   virtual Real stop() = 0;

   /// return type of timer
   virtual TYPE type() = 0;
   //@}

   //------------------------------------
   /**@name Access */
   //@{
   /// return accounted time.
   /// get accounted user, system, or real time when ticks were updated last
   void getLastTimes(Real* userTime, Real* systemTime, Real* realTime) const;


   virtual Real time() const = 0;

   /// return last accounted time without rechecking the clock
   virtual Real lastTime() const = 0;


   /// return accounted real time without rechecking the clock
   Real realTimeLast() const;

   //@}
};
} // namespace soplex
#endif // _TIMER_H_
