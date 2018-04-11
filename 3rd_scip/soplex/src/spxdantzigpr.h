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

/**@file  spxdantzigpr.h
 * @brief Dantzig pricer.
 */
#ifndef _SPXDEFAULTPR_H_
#define _SPXDEFAULTPR_H_

#include <assert.h>

#include "spxpricer.h"

namespace soplex
{

/**@brief   Dantzig pricer.
   @ingroup Algo

   Class SPxDantzigPR is an implementation class of an SPxPricer implementing
   Dantzig's default pricing strategy, i.e., maximal/minimal reduced cost or
   maximally violated constraint.

   See SPxPricer for a class documentation.
*/
class SPxDantzigPR : public SPxPricer
{
private:
   int                   selectLeaveSparse();/**< sparse pricing method for leaving Simplex */

   SPxId                 selectEnterX();                                /**< choose the best entering index among columns and rows but prefer sparsity */
   SPxId                 selectEnterSparseDim(Real& best,SPxId& id);    /**< sparse pricing method for entering Simplex (slack variables)*/
   SPxId                 selectEnterSparseCoDim(Real& best, SPxId& id); /**< sparse pricing method for entering Simplex */
   SPxId                 selectEnterDenseDim(Real& best, SPxId& id);    /**< selectEnter() in dense case (slack variables) */
   SPxId                 selectEnterDenseCoDim(Real& best,SPxId& id);   /**< selectEnter() in dense case */
public:

   //-------------------------------------
   /**@name Constructors / destructors */
   //@{
   /// default constructor
   SPxDantzigPR() 
      : SPxPricer("Dantzig")
   {} 
   /// copy constructor
   SPxDantzigPR(const SPxDantzigPR& old ) 
      : SPxPricer(old)
   {}
   /// assignment operator
   SPxDantzigPR& operator=( const SPxDantzigPR& rhs)
   {
      if(this != &rhs)
      {
         SPxPricer::operator=(rhs);
      }

      return *this;
   }   
   /// destructor
   virtual ~SPxDantzigPR()
   {}
   /// clone function for polymorphism
   inline virtual SPxPricer* clone() const
   {
      return new SPxDantzigPR(*this);
   }
   //@}


   //-------------------------------------
   /**@name Select enter/leave */
   //@{
   ///
   virtual int selectLeave();
   ///
   virtual SPxId selectEnter();
   //@}
};
} // namespace soplex
#endif // _SPXDEFAULTPRR_H_
