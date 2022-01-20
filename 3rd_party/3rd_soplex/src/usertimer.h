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

/**@file  usertimer.h
 * @brief UserTimer class.
 */

#ifndef _USER_TIMER_H_
#define _USER_TIMER_H_

#include "spxdefines.h"
#include "timer.h"

namespace soplex
{

class UserTimer : public Timer
{
private:

   //------------------------------------
   /**@name number of ticks per second */
   //@{
   static const long ticks_per_sec;  ///< ticks per secound, should be constant
   //@}

   //------------------------------------
   /**@name Data */
   //@{
   mutable long uAccount;      ///< user time
   mutable long uTicks;        ///< user ticks

   mutable Real lasttime;
   //@}

   //------------------------------------
   /**@name Internal helpers */
   //@{
   /// convert ticks to secounds.
   Real ticks2sec(long ticks) const
   {
      return (Real(ticks) * 1000.0 / Real(ticks_per_sec)) / 1000.0;
   }

   /// get actual user ticks from the system.
   void updateTicks() const;

   //@}

public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   UserTimer()
      : Timer(), uAccount(0), uTicks(0), lasttime(0.0)
   {
      assert(ticks_per_sec > 0);
   }
   /// copy constructor
   UserTimer(const UserTimer& old)
      : Timer(), uAccount(old.uAccount), uTicks(old.uTicks), lasttime(old.lasttime)
   {
      assert(ticks_per_sec > 0);
   }
   /// assignment operator
   UserTimer& operator=(const UserTimer& old)
   {
      assert(ticks_per_sec > 0);
      uAccount = old.uAccount;
      uTicks = old.uTicks;
      lasttime = old.lasttime;
      return *this;
   }

   virtual ~UserTimer()
   {}
   //@}

   //------------------------------------
   /**@name Control */
   //@{
   /// initialize timer, set timing accounts to zero.
   virtual void reset()
   {
      status   = RESET;
      uAccount = 0;
      lasttime = 0.0;
   }

   /// start timer, resume accounting user, system and real time.
   virtual void start();

   /// stop timer, return accounted user time.
   virtual Real stop();

   /// return type of timer
   virtual TYPE type()
   {
      return USER_TIME;
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
#endif // _USER_TIMER_H_
