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

/**@file  spxvectorst.h
 * @brief Solution vector based start basis.
 */
#ifndef _SPXVECTORST_H_
#define _SPXVECTORST_H_

#include <assert.h>

#include "spxweightst.h"
#include "vector.h"

namespace soplex
{

/**@brief   Solution vector based start basis.
   @ingroup Algo

   This version of SPxWeightST can be used to construct a starting basis for
   an LP to be solved with SoPlex if an approximate solution vector or dual
   vector (possibly optained by a heuristic) is available. This is done by
   setting up weights for the SPxWeightST it is derived from.
   
   The primal vector to be used is loaded by calling method #primal() while
   #dual() setups for the dual vector. Methods #primal() or #dual() must be
   called \em before #generate() is called by SoPlex to set up a
   starting basis. If more than one call of method #primal() or #dual()
   occurred only the most recent one is valid for generating the starting base.
*/
class SPxVectorST : public SPxWeightST
{
private:

   //-------------------------------------
   /**@name Types */
   //@{
   /// specifies whether to work on the primal, the dual, or not at all.
   enum { NONE, PVEC, DVEC } state;
   //@}

   //-------------------------------------
   /**@name Data */
   //@{
   /// the current (approximate) primal or dual vector
   DVector vec;
   //@}

protected:

   //-------------------------------------
   /**@name Protected helpers */
   //@{
   /// sets up variable weights.
   void setupWeights(SPxSolver& base);
   //@}

public:

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor.
   SPxVectorST() 
      : state(NONE)
   {
      m_name = "Vector";
   }
   /// copy constructor
   SPxVectorST( const SPxVectorST& old)
      : SPxWeightST(old)
      , state(old.state)
      , vec(old.vec)
   {
      assert(isConsistent());
   }
   /// assignment operator
   SPxVectorST& operator=( const SPxVectorST& rhs)
   {
      if(this != &rhs)
      {
         SPxWeightST::operator=(rhs);
         state = rhs.state;
         vec = rhs.vec;

         assert(isConsistent());
      }

      return *this;
   }
   /// destructor.
   virtual ~SPxVectorST()
   {}  
   /// clone function for polymorphism
   inline virtual SPxStarter* clone() const
   {
      return new SPxVectorST(*this);
   }
   //@}

   //-------------------------------------
   /**@name Modification */
   //@{
   /// sets up primal solution vector.
   void primal(const Vector& v)
   {
      vec = v;
      state = PVEC;
   }
   /// sets up primal solution vector.
   void dual(const Vector& v)
   {
      vec = v;
      state = DVEC;
   }
   //@}

};

} // namespace soplex
#endif // _SPXVECTORST_H_
