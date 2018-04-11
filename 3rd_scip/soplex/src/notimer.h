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

/**@file  notimer.h
 * @brief NoTimer class.
 */

#ifndef _NO_TIMER_H_
#define _NO_TIMER_H_

#include "spxdefines.h"
#include "timer.h"

namespace soplex
{

class NoTimer : public Timer
{

public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   NoTimer()
      : Timer()
   {}
   /// copy constructor
   NoTimer(const NoTimer&)
      : Timer()
   {}
   /// assignment operator
   NoTimer& operator=(const NoTimer&)
   {
      return *this;
   }

   virtual ~NoTimer()
   {}
   //@}

   //------------------------------------
   /**@name Control */
   //@{
   /// initialize timer
   virtual void reset()
   {}

   /// start timer
   virtual void start()
   {}

   /// stop timer
   virtual Real stop()
   {
      return 0.0;
   }

   /// return type of timer
   virtual TYPE type()
   {
      return OFF;
   }
   //@}

   //------------------------------------
   /**@name Access */
   //@{
   virtual Real time() const
   {
      return 0.0;
   }

   virtual Real lastTime() const
   {
      return 0.0;
   }

   //@}
};
} // namespace soplex
#endif // _NO_TIMER_H_
