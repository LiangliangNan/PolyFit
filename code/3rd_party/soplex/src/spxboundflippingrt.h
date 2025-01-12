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

/**@file   spxboundflippingrt.h
 * @brief  Bound flipping ratio test (long step dual) for SoPlex.
 * @author Matthias Miltenberger
 * @author Eva Ramlow
 */
#ifndef _SPXBOUNDFLIPPINGRT_H_
#define _SPXBOUNDFLIPPINGRT_H_


#include <assert.h>
#include "spxdefines.h"
#include "spxratiotester.h"
#include "spxfastrt.h"

namespace soplex
{
struct Compare;

/**@brief   Bound flipping ratio test ("long step dual") for SoPlex.
   @ingroup Algo

   Class SPxBoundFlippingRT provides an implementation of the bound flipping
   ratio test as a derived class of SPxRatioTester. Dual step length is
   increased beyond some breakpoints and corresponding primal nonbasic
   variables are set to their other bound to handle the resulting dual infeasibility.

   The implementation mostly follows the papers
   - I. Maros, "A generalized dual phase-2 simplex algorithm",
     European Journal of Operational Research Vol 149, Issue 1, pp. 1-16, 2003
   - A. Koberstein, "Progress in the dual simplex algorithm for solving large scale LP problems:
     techniques for a fast and stable implementation",
     Computational Optimization and Applications Vol 41, Nr 2, pp. 185-204, 2008

   See SPxRatioTester for a class documentation.
*/
class SPxBoundFlippingRT : public SPxFastRT
{
private:
   /**@name substructures */
   //@{
   /** enumerator to remember which vector we have been searching to find a breakpoint
    */
   enum BreakpointSource
   {
      FVEC               = -1,
      PVEC               = 0,
      COPVEC             = 1
   };

   /** breakpoint structure
    */
   struct Breakpoint
   {
      Real               val;                /**< breakpoint value (step length) */
      int                idx;                /**< index of corresponding row/column */
      BreakpointSource   src;                /**< origin of breakpoint, i.e. vector which was searched */
   };

   /** Compare class for breakpoints
    */
   struct BreakpointCompare
   {
   public:
      /** constructor
       */
      BreakpointCompare()
         : entry(0)
      {
      }

      const Breakpoint*  entry;

      Real operator() (
         Breakpoint      i,
         Breakpoint      j
         ) const
      {
         return i.val - j.val;
      }
   };
   //@}

   /**@name Data
    */
   //@{
   bool                  enableBoundFlips;   /**< enable or disable long steps in BoundFlippingRT */
   bool                  enableRowBoundFlips;/**< enable bound flips also for row representation */
   Real                  flipPotential;      /**< tracks bound flip history and decides which ratio test to use */
   int                   relax_count;        /**< count rounds of ratio test */
   DataArray<Breakpoint> breakpoints;        /**< array of breakpoints */
   SSVector              updPrimRhs;         /**< right hand side vector of additional system to be solved after the ratio test */
   SSVector              updPrimVec;         /**< allocation of memory for additional solution vector */
   //@}

   /** store all available pivots/breakpoints in an array (positive pivot search direction) */
   void collectBreakpointsMax(
      int&               nBp,                /**< number of found breakpoints so far */
      int&               minIdx,             /**< index to current minimal breakpoint */
      const int*         idx,                /**< pointer to indices of current vector */
      int                nnz,                /**< number of nonzeros in current vector */
      const Real*        upd,                /**< pointer to update values of current vector */
      const Real*        vec,                /**< pointer to values of current vector */
      const Real*        upp,                /**< pointer to upper bound/rhs of current vector */
      const Real*        low,                /**< pointer to lower bound/lhs of current vector */
      BreakpointSource   src                 /**< type of vector (pVec or coPvec)*/
      );

   /** store all available pivots/breakpoints in an array (negative pivot search direction) */
   void collectBreakpointsMin(
      int&               nBp,                /**< number of found breakpoints so far */
      int&               minIdx,             /**< index to current minimal breakpoint */
      const int*         idx,                /**< pointer to indices of current vector */
      int                nnz,                /**< number of nonzeros in current vector */
      const Real*        upd,                /**< pointer to update values of current vector */
      const Real*        vec,                /**< pointer to values of current vector */
      const Real*        upp,                /**< pointer to upper bound/rhs of current vector */
      const Real*        low,                /**< pointer to lower bound/lhs of current vector */
      BreakpointSource   src                 /**< type of vector (pVec or coPvec)*/
      );

   /** get values for entering index and perform shifts if necessary */
   bool getData(
      Real&              val,
      SPxId&             enterId,
      int                idx,
      Real               stab,
      Real               degeneps,
      const Real*        upd,
      const Real*        vec,
      const Real*        low,
      const Real*        upp,
      BreakpointSource   src,
      Real               max
      );

   /** get values for leaving index and perform shifts if necessary */
   bool getData(
      Real&              val,
      int&             leaveIdx,
      int                idx,
      Real               stab,
      Real               degeneps,
      const Real*        upd,
      const Real*        vec,
      const Real*        low,
      const Real*        upp,
      BreakpointSource   src,
      Real               max
      );

   /** perform necessary bound flips to restore dual feasibility */
   void flipAndUpdate(
      int&               usedBp              /**< number of bounds that should be flipped */
      );

   /** comparison method for breakpoints */
   static bool isSmaller(
      Breakpoint         x,
      Breakpoint         y
      )
   {
      return (spxAbs(x.val) < spxAbs(y.val));
   };

public:

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   SPxBoundFlippingRT()
      : SPxFastRT("Bound Flipping")
      , enableBoundFlips(true)
      , enableRowBoundFlips(false)
      , flipPotential(1)
      , relax_count(0)
      , breakpoints(10)
      , updPrimRhs(0)
      , updPrimVec(0)
   {}
   /// copy constructor
   SPxBoundFlippingRT(const SPxBoundFlippingRT& old)
      : SPxFastRT(old)
      , enableBoundFlips(old.enableBoundFlips)
      , enableRowBoundFlips(old.enableRowBoundFlips)
      , flipPotential(1)
      , relax_count(0)
      , breakpoints(10)
      , updPrimRhs(0)
      , updPrimVec(0)
   {}
   /// assignment operator
   SPxBoundFlippingRT& operator=( const SPxBoundFlippingRT& rhs)
   {
      if(this != &rhs)
      {
         SPxFastRT::operator=(rhs);
      }

      enableBoundFlips = rhs.enableBoundFlips;
      enableRowBoundFlips = rhs.enableRowBoundFlips;
      flipPotential = rhs.flipPotential;

      return *this;
   }
   /// destructor
   virtual ~SPxBoundFlippingRT()
   {}
   /// clone function for polymorphism
   inline virtual SPxRatioTester* clone() const
   {
      return new SPxBoundFlippingRT(*this);
   }
   //@}

   //-------------------------------------
   /**@name Select enter/leave */
   //@{
   ///
   virtual int selectLeave(
      Real&              val,
      Real               enterTest,
      bool               polish = false
      );
   ///
   virtual SPxId selectEnter(
      Real&              val,
      int                leaveIdx,
      bool               polish = false
      );

   void useBoundFlips(bool bf)
   {
      enableBoundFlips = bf;
   }

   void useBoundFlipsRow(bool bf)
   {
      enableRowBoundFlips = bf;
   }
   //@}
};

} // namespace soplex
#endif // _SPXBOUNDFLIPPINGRT_H_
