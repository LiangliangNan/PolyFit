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

/**@file  wallclocktimer.h
 * @brief WallclockTimer class.
 */

#ifndef _WALLCLOCK_TIMER_H_
#define _WALLCLOCK_TIMER_H_

#include "spxdefines.h"
#include "timer.h"

namespace soplex
{

class WallclockTimer : public Timer
{
private:

   //------------------------------------
   /**@name Data */
   //@{
   mutable long sec;           ///< seconds
   mutable long usec;          ///< microseconds

   mutable Real lasttime;
   //@}

   //------------------------------------
   /**@name Internal helpers */
   //@{
   /// convert wallclock time to secounds.
   Real wall2sec(long s, long us) const
   {
      return (Real)s+ 0.000001 * (Real)us;
   }

   //@}

public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   WallclockTimer()
      : Timer(), sec(0), usec(0), lasttime(0.0)
   {}
   /// copy constructor
   WallclockTimer(const WallclockTimer& old)
      : Timer(), sec(old.sec), usec(old.usec), lasttime(old.lasttime)
   {}
   /// assignment operator
   WallclockTimer& operator=(const WallclockTimer& old)
   {
      sec = old.sec;
      usec = old.usec;
      lasttime = old.lasttime;
      return *this;
   }

   virtual ~WallclockTimer()
   {}
   //@}

   //------------------------------------
   /**@name Control */
   //@{
   /// initialize timer, set timing accounts to zero.
   virtual void reset()
   {
      status   = RESET;
      sec = usec = 0;
      lasttime = 0.0;
   }

   /// start timer, resume accounting user, system and real time.
   virtual void start();

   /// stop timer, return accounted user time.
   virtual Real stop();

   /// return type of timer
   virtual TYPE type()
   {
      return WALLCLOCK_TIME;
   }
   //@}

   //------------------------------------
   /**@name Access */
   //@{
   virtual Real time() const;

   virtual Real lastTime() const;

   //@}
};
} // namespace soplex
#endif // _WALLCLOCK_TIMER_H_
