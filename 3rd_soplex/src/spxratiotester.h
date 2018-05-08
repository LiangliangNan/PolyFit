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

/**@file  spxratiotester.h
 * @brief Abstract ratio test base class.
 */
#ifndef _SPXRATIOTESTER_H_
#define _SPXRATIOTESTER_H_


#include <assert.h>

#include "spxdefines.h"
#include "spxsolver.h"

namespace soplex
{

/**@brief Abstract ratio test base class.
   @ingroup Algo

   Class SPxRatioTester is the virtual base class for computing the ratio
   test within the Simplex algorithm driven by SoPlex. After a SoPlex
   solver has been #load()%ed to an SPxRatioTester, the solver calls
   #selectLeave() for computing the ratio test for the entering simplex and
   #selectEnter() for computing the ratio test in leaving simplex.
*/
class SPxRatioTester
{
protected:

   //-------------------------------------
   /**@name Data */
   //@{
   /// the solver
   SPxSolver*  thesolver;
   /// name of the ratio tester
   const char* m_name;
   /// internal storage of type
   SPxSolver::Type m_type;
   /// allowed bound violation
   Real delta;
   //@}

public:

   //-------------------------------------
   /**@name Access / modification */
   //@{
   /// get name of ratio tester.
   virtual const char* getName() const
   {
      return m_name;
   }
   /// loads LP.
   /** Load the solver and LP for which pricing steps are to be performed.
    */
   virtual void load(SPxSolver* p_solver)
   {
      thesolver = p_solver;
   }

   /// unloads LP.
   virtual void clear()
   {
      thesolver = 0;
   }

   /// returns loaded LP solver.
   virtual SPxSolver* solver() const
   {
      return thesolver;
   }

   /// set allowed bound violation
   virtual void setDelta( Real newDelta )
   {
      if( newDelta <= DEFAULT_EPS_ZERO )
         delta = DEFAULT_EPS_ZERO;
      else
         delta = newDelta;
   }

   /// get allowed bound violation
   virtual Real getDelta()
   {
      return delta;
   }
   //@}

   //-------------------------------------
   /**@name Entering / leaving */
   //@{
   /// selects index to leave the basis.
   /** Method #selectLeave() is called by the loaded SoPlex solver when
       computing the entering simplex algorithm. Its task is to select and
       return the index of the basis variable that is to leave the basis.
       When being called, 
       \ref SPxSolver::fVec() "fVec()" fullfills the basic bounds 
       \ref SPxSolver::lbBound() "lbBound()" and
       \ref SPxSolver::ubBound() "ubBound()" within 
       \ref SPxSolver::entertol() "entertol()".
       fVec().delta() is the vector by
       which fVec() will be updated in this simplex step. Its nonzero
       indices are stored in sorted order in fVec().idx().
       
       If \p val > 0, \p val is the maximum allowed update value for fVec(),
       otherwise the minimum. Method #selectLeave() must chose \p val of the
       same sign as passed, such that updating fVec() by \p val yields a
       new vector that satisfies all basic bounds (within entertol). The
       returned index, must be the index of an element of fVec(), that
       reaches one of its bounds with this update.
   */
   virtual int selectLeave(Real& val, Real enterTest, bool polish = false) = 0;

   /// selects variable Id to enter the basis.
   /** Method #selectEnter() is called by the loaded SoPlex solver, when
       computing the leaving simplex algorithm. It's task is to select and
       return the Id of the basis variable that is to enter the basis.
       When being called, 
       \ref SPxSolver::pVec() "pVec()" fullfills the bounds 
       \ref SPxSolver::lbBound() "lbBound()" and
       \ref SPxSolver::ubBound() "ubBound()" within 
       \ref SPxSolver::leavetol() "leavetol()".
       Similarly, 
       \ref SPxSolver::coPvec() "coPvec()" fulfills the bounds
       \ref SPxSolver::lbBound() "lbBound()" and
       \ref SPxSolver::ubBound() "ubBound()" within
       \ref SPxSolver::leavetol() "leavetol()".
       pVec().delta() and coPvec().delta() are
       the vectors by which pVec() and coPvec() will be updated in this
       simplex step. Their nonzero indices are stored in sorted order in
       pVec().idx() and coPvec().idx().
       
       If \p val > 0, \p val is the maximum allowed update value for pVec()
       and coPvec(), otherwise the minimum. Method #selectEnter() must
       chose \p val of the same sign as passed, such that updating pVec()
       and coPvec() by \p val yields a new vector that satisfies all basic
       bounds (within leavetol). The returned Id must be the Id of an
       element of pVec() or coPvec(), that reaches one of its bounds
       with this update.
    */
   virtual SPxId selectEnter(Real& val, int leaveIdx, bool polish = false) = 0;

   /// sets Simplex type.
   /** Informs pricer about (a change of) the loaded SoPlex's Type. In
       the sequel, only the corresponding select methods may be called.
   */
   virtual void setType(SPxSolver::Type)
   {}
   //@}

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   explicit SPxRatioTester(const char* name)
      : thesolver(0)
      , m_name(name)
      , m_type(SPxSolver::LEAVE)
      , delta(1e-6)
   {}
   /// copy constructor
   SPxRatioTester( const SPxRatioTester& old)
      : thesolver(old.thesolver)
      , m_name(old.m_name)
      , m_type(old.m_type)
      , delta(old.delta)
   {}
   /// assignment operator
   SPxRatioTester& operator=( const SPxRatioTester& rhs)
   {
      if(this != &rhs)
      {
         m_name = rhs.m_name;
         thesolver = rhs.thesolver;
         m_type = rhs.m_type;
         delta = rhs.delta;
      }

      return *this;
   }
   /// destructor.
   virtual ~SPxRatioTester()
   {
      thesolver = 0;
      m_name    = 0;
   }
   /// clone function for polymorphism
   virtual SPxRatioTester* clone() const = 0;
   //@}

};


} // namespace soplex
#endif // _SPXRATIOTESTER_H_
