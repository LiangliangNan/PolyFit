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


/**@file  spxstarter.h
 * @brief SoPlex start basis generation base class.
 */
#ifndef _SPXDSTARTER_H_
#define _SPXDSTARTER_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxsolver.h"

namespace soplex
{

/**@brief   SoPlex start basis generation base class.
   @ingroup Algo
   
   SPxStarter is the virtual base class for classes generating a starter basis
   for the Simplex solver SoPlex. When a SPxStarter object has been loaded
   to a SoPlex solver, the latter will call method #generate() in order to
   have a start basis generated. Implementations of method #generate() must
   terminate by \ref soplex::SPxSolver::load() "loading" the generated basis to 
   SoPlex. Loaded bases must be nonsingular.
*/
class SPxStarter
{
protected:

   //-------------------------------------
   /**@name Data */
   //@{
   /// name of the starter
   const char* m_name;
   //@}

public:

   //-------------------------------------
   /**@name Data */
   //@{
   /// constructor
   explicit SPxStarter(const char* name)
      : m_name(name)
   {}
   /// copy constructor
   SPxStarter( const SPxStarter& old)
      : m_name(old.m_name)
   {}
   /// assignment operator
   SPxStarter& operator=( const SPxStarter& rhs)
   {
      if(this != &rhs)
      {
         m_name = rhs.m_name;
      }

      return *this;
   }
   /// destructor.
   virtual ~SPxStarter()
   {
      m_name = 0;
   }
   /// clone function for polymorphism
   virtual SPxStarter* clone()const = 0;
   //@}

   //-------------------------------------
   /**@name Access */
   //@{
   /// get name of starter.
   virtual const char* getName() const
   {
      return m_name;
   }
   //@}

   //-------------------------------------
   /**@name Starting */
   //@{
   /// generates start basis for loaded basis.
   virtual void generate(SPxSolver& base) = 0;
   //@}

   //-------------------------------------
   /**@name Misc */
   //@{
   /// checks consistency.
   virtual bool isConsistent() const;
   //@}

private:

   //------------------------------------
   /**@name Blocked */
   //@{
   /// we have no default constructor.
   SPxStarter();
   //@}

};
} // namespace soplex
#endif // _SPXDSTARTER_H_
