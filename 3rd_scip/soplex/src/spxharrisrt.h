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

/**@file  spxharrisrt.h
 * @brief Harris pricing with shifting.
 */
#ifndef _SPXHARRISRT_H_
#define _SPXHARRISRT_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxratiotester.h"

namespace soplex
{

/**@brief   Harris pricing with shifting.
   @ingroup Algo
   
   Class SPxHarrisRT is a stable implementation of a SPxRatioTester class
   along the lines of Harris' two phase algorithm. Additionally it uses
   shifting of bounds in order to avoid cycling.

   See SPxRatioTester for a class documentation.
*/
/**@todo HarrisRT leads to cycling in dcmulti.sub.lp */
class SPxHarrisRT : public SPxRatioTester
{
private:

   //-------------------------------------
   /**@name Private helpers */
   //@{
   ///
   Real degenerateEps() const;

   ///
   int maxDelta(
      Real* /*max*/,        ///< max abs value in \p upd
      Real* val,            ///< initial and chosen value
      int num,              ///< number of indices in \p idx
      const int* idx,       ///< nonzero indices in \p upd
      const Real* upd,      ///< update vector for \p vec
      const Real* vec,      ///< current vector
      const Real* low,      ///< lower bounds for \p vec
      const Real* up,       ///< upper bounds for \p vec
      Real epsilon          ///< what is 0?
      ) const;

   ///
   int minDelta(
      Real* /*max*/,        ///< max abs value in \p upd
      Real* val,            ///< initial and chosen value
      int num,              ///< of indices in \p idx
      const int* idx,       ///< nonzero indices in \p upd
      const Real* upd,      ///< update vector for \p vec
      const Real* vec,      ///< current vector
      const Real* low,      ///< lower bounds for \p vec
      const Real* up,       ///< upper bounds for \p vec
      Real epsilon          ///< what is 0?
      ) const;
   //@}

public:

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   SPxHarrisRT() 
      : SPxRatioTester("Harris")
   {}
   /// copy constructor
   SPxHarrisRT(const SPxHarrisRT& old) 
      : SPxRatioTester(old)
   {}
   /// assignment operator
   SPxHarrisRT& operator=( const SPxHarrisRT& rhs)
   {
      if(this != &rhs)
      {
         SPxRatioTester::operator=(rhs);
      }

      return *this;
   }
   /// destructor
   virtual ~SPxHarrisRT()
   {}
   /// clone function for polymorphism
   inline virtual SPxRatioTester* clone() const
   {
      return new SPxHarrisRT(*this);
   }
   //@}

   //-------------------------------------
   /**@name Leave / enter */
   //@{
   ///
   virtual int selectLeave(Real& val, Real, bool);
   ///
   virtual SPxId selectEnter(Real& val, int, bool);
   //@}

};

} // namespace soplex
#endif // _SPXHARRISRT_H_
