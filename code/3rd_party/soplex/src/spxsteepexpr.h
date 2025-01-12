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


/**@file  spxsteepexpr.h
 * @brief Steepest edge pricer with exact initialization of weights.
 */
#ifndef _SPXSTEEPEXPR_H_
#define _SPXSTEEPEXPR_H_


#include <assert.h>

#include "spxdefines.h"
#include "spxsteeppr.h"

namespace soplex
{

/**@brief   Steepest edge pricer.
   @ingroup Algo

   Class SPxSteepExPR implements a steepest edge pricer to be used with
   SoPlex. Exact initialization of weights is used.

   See SPxPricer for a class documentation.
*/
class SPxSteepExPR : public SPxSteepPR
{

public:

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   ///
   SPxSteepExPR()
      : SPxSteepPR("SteepEx", EXACT)
   {
      assert(isConsistent());
   }
   /// copy constructor
   SPxSteepExPR( const SPxSteepExPR& old)
      : SPxSteepPR(old)
   {
      assert(isConsistent());
   }
   /// assignment operator
   SPxSteepExPR& operator=( const SPxSteepExPR& rhs)
   {
      if(this != &rhs)
      {
         SPxSteepPR::operator=(rhs);

         assert(isConsistent());
      }

      return *this;
   }
   /// destructor
   virtual ~SPxSteepExPR()
   {}
   /// clone function for polymorphism
   inline virtual SPxSteepPR* clone()  const
   {
      return new SPxSteepExPR(*this);
   }
   //@}
};

} // namespace soplex
#endif // _SPXSTEEPPR_H_
