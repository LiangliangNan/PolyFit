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

/**@file  timerfactory.h
 * @brief TimerFactory class.
 */

#ifndef _TIMERFACTORY_H_
#define _TIMERFACTORY_H_

#include "spxdefines.h"
#include "spxalloc.h"
#include "notimer.h"
#include "usertimer.h"
#include "wallclocktimer.h"

namespace soplex
{
/**@class   TimerFactory
   @ingroup Elementary

   @brief Class to create new timers and to switch types of exiting ones
   */

class TimerFactory
{

public:

   /// create timers and allocate memory for them
   static Timer* createTimer(Timer::TYPE ttype)
   {
      Timer* timer = 0;
      switch( ttype )
      {
      case Timer::OFF:
         spx_alloc(timer, sizeof(NoTimer));
         timer = new (timer) NoTimer();
         break;
      case Timer::USER_TIME:
         spx_alloc(timer, sizeof(UserTimer));
         timer = new (timer) UserTimer();
         break;
      case Timer::WALLCLOCK_TIME:
         spx_alloc(timer, sizeof(WallclockTimer));
         timer = new (timer) WallclockTimer();
         break;
      default:
         MSG_ERROR( std::cerr << "wrong timer specified" << std::endl; )
      }
      return timer;
   }

   static Timer* switchTimer(Timer* timer, Timer::TYPE ttype)
   {
      // check whether the type is different from the current one
      if( ttype != timer->type() )
      {
         // @todo transfer the old times
         spx_free(timer);
         timer = createTimer(ttype);
      }
      return timer;
   }

};
} // namespace soplex
#endif // _TIMERFACTORY_H_
