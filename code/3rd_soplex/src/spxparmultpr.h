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

/**@file  spxparmultpr.h
 * @brief Partial multiple pricing.
 */
#ifndef _SPXPARMULTPR_H_
#define _SPXPARMULTPR_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxpricer.h"
#include "dataarray.h"
#include "array.h"
#include "ssvector.h"

namespace soplex
{

/**@brief   Partial multiple pricing.
   @ingroup Algo

   Class SPxParMultPr is an implementation class for SPxPricer implementing
   Dantzig's default pricing strategy with partial multiple pricing.
   Partial multiple pricing applies to the entering Simplex only. A set of
   #partialSize eligible pivot indices is selected (partial pricing). In the
   following Simplex iterations pricing is restricted to these indices
   (multiple pricing) until no more eliiable pivots are available. Partial
   multiple pricing significantly reduces the computation time for computing
   the matrix-vector-product in the Simplex algorithm.

   See SPxPricer for a class documentation.
*/
class SPxParMultPR : public SPxPricer
{
private:

   //-------------------------------------
   /**@name Private types */
   //@{
   /// Helper structure.
   struct SPxParMultPr_Tmp
   {
      ///
      SPxId id;
      ///
      Real test;
   };
   //@}

   //-------------------------------------
   /**@name Helper data */
   //@{
   ///
   DataArray < SPxParMultPr_Tmp > pricSet;
   ///
   int multiParts;
   ///
   int used;
   ///
   int min;
   ///
   int last;
   /// Set size for partial pricing.
   int partialSize;
   //@}

public:

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   SPxParMultPR() 
      : SPxPricer("ParMult")
      , multiParts(0)
      , used(0)
      , min(0)
      , last(0)
      , partialSize(17)
   {}
   /// copy constructor
   SPxParMultPR( const SPxParMultPR& old)
      : SPxPricer(old)
      , pricSet(old.pricSet)
      , multiParts(old.multiParts)
      , used(old.used)
      , min(old.min)
      , last(old.last)
      , partialSize(old.partialSize)
   {}
   /// assignment operator
   SPxParMultPR& operator=( const SPxParMultPR& rhs)
   {
      if(this != &rhs)
      {
         SPxPricer::operator=(rhs);
         pricSet = rhs.pricSet;
         multiParts = rhs.multiParts;
         used = rhs.used;
         min = rhs.min;
         last = rhs.last;
         partialSize = rhs.partialSize;
      }

      return *this;
   }
   /// destructor
   virtual ~SPxParMultPR()
   {}
   /// clone function for polymorphism
   inline virtual SPxPricer* clone()  const
   {
      return new SPxParMultPR(*this);
   }
   //@}

   //-------------------------------------
   /**@name Interface */
   //@{
   /// set the solver
   virtual void load(SPxSolver* solver);
   /// set entering or leaving algorithm
   virtual void setType(SPxSolver::Type tp);
   /// 
   virtual int selectLeave();
   ///
   virtual SPxId selectEnter();
   //@}

};

} // namespace soplex
#endif // _SPXPARMULTPRR_H_
