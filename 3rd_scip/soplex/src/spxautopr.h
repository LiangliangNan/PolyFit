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

/**@file  spxautopr.h
 * @brief Auto pricer.
 */
#ifndef _SPXAUTOPR_H_
#define _SPXAUTOPR_H_

#include <assert.h>

#include "spxpricer.h"
#include "spxdevexpr.h"
#include "spxsteeppr.h"
#include "spxsteepexpr.h"


namespace soplex
{

/**@brief   Auto pricer.
   @ingroup Algo

   This pricer switches between Devex and Steepest edge pricer based on the difficulty of the problem
   which is determined by the number of iterations.

   See SPxPricer for a class documentation.
*/
class SPxAutoPR : public SPxPricer
{
private:

   int            switchIters;   ///< number of iterations before switching pricers
   SPxPricer*     activepricer;  ///< pointer to currently selected pricer
   SPxDevexPR     devex;         ///< internal Devex pricer
   SPxSteepExPR     steep;         ///< internal Steepest edge pricer

   bool setActivePricer(SPxSolver::Type type);          ///< switches active pricing method

public:

   //-------------------------------------
   /**@name Constructors / destructors */
   //@{
   /// default constructor
   SPxAutoPR()
      : SPxPricer("Auto")
      , switchIters(10000)
      , activepricer(&devex)
      , devex()
      , steep()
   {}
   /// copy constructor
   SPxAutoPR(const SPxAutoPR& old )
      : SPxPricer(old)
      , switchIters(old.switchIters)
      , devex(old.devex)
      , steep(old.steep)
   {
      assert(old.activepricer == &old.devex || old.activepricer == &old.steep);
      if( old.activepricer == &old.devex )
         activepricer = &devex;
      else
         activepricer = &steep;
   }
   /// assignment operator
   SPxAutoPR& operator=( const SPxAutoPR& rhs)
   {
      if(this != &rhs)
      {
         SPxPricer::operator=(rhs);
         switchIters = rhs.switchIters;
         devex = rhs.devex;
         steep = rhs.steep;

         assert(rhs.activepricer == &rhs.devex || rhs.activepricer == &rhs.steep);
         if( rhs.activepricer == &rhs.devex )
            activepricer = &devex;
         else
            activepricer = &steep;
      }

      return *this;
   }
   /// destructor
   virtual ~SPxAutoPR()
   {}
   /// clone function for polymorphism
   inline virtual SPxPricer* clone() const
   {
      return new SPxAutoPR(*this);
   }
   //@}

   //-------------------------------------
   /**@name Access / modification */
   //@{
   /// set max number of iterations before switching pricers
   void setSwitchIters(int iters);
   /// clear the data
   void clear();
   /// set epsilon of internal pricers
   void setEpsilon(Real eps);
   /// set the solver
   virtual void load(SPxSolver* base);
   /// set entering/leaving algorithm
   virtual void setType(SPxSolver::Type);
   /// set row/column representation
   virtual void setRep(SPxSolver::Representation);
   ///
   virtual int selectLeave();
   ///
   virtual SPxId selectEnter();
   ///
   virtual void left4(int n, SPxId id);
   ///
   virtual void entered4(SPxId id, int n);
   //@}
};
} // namespace soplex
#endif // _SPXAUTOPR_H_
