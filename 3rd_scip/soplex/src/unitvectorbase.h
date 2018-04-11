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

/**@file  unitvector.h
 * @brief Sparse vector \f$e_i\f$.
 */

#ifndef _UNITVECTORBASE_H_
#define _UNITVECTORBASE_H_

#include <assert.h>
#include "spxdefines.h"
#include "svectorbase.h"

namespace soplex
{


/**@brief   Sparse vector \f$e_i\f$.
   @ingroup Algebra

   A UnitVectorBase is an SVectorBase that can take only one nonzero value with
   value 1 but arbitrary index.

   \todo Several SVectorBase modification methods are still accessible for UnitVector.
   They might be used to change the vector.

   \todo UnitVectorBase memory management must be changed when SVectorBase is redesigned.
*/

template < class R >
class UnitVectorBase : public SVectorBase<R>
{
private:

   //------------------------------------
   /**@name Data */
   //@{
   typename SVectorBase<R>::Element themem;  ///< memory for sparse vector entry
   //@}

   using SVectorBase<R>::mem;

   using SVectorBase<R>::size;

   using SVectorBase<R>::max;

public:

   //------------------------------------
   /**@name Access */
   //@{
   /// returns value = 1
   /**\pre \c n must be 0.
    */
   /* ARGSUSED n */
   R value(int n) const
   {
      assert( n == 0 );
      return 1;
   }
   //@}

   //------------------------------------
   /**@name Constructors / destructors */
   //@{
   /// construct \c i 'th unit vector.
   explicit
   UnitVectorBase<R>(int i = 0)
      : SVectorBase<R>(1, &themem)
   {
      SVectorBase<R>::add(i, 1.0);

      assert(isConsistent());
   }
   ///  copy constructor
   UnitVectorBase<R>(const UnitVectorBase<R>& rhs)
      : SVectorBase<R>(1, &themem)
   {
      themem = rhs.themem;

      assert(isConsistent());
   }
   /// assignment
   UnitVectorBase<R>& operator=(const UnitVectorBase<R>& rhs)
   {
      if ( this != &rhs )
      {
         themem = rhs.themem;

         assert(isConsistent());
      }
      return *this;
   }
   /// destructor
   ~UnitVectorBase<R>()
   {}
   //@}

   //------------------------------------
   /**@name Miscellaneous */
   //@{
   /// consistency check
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      if (mem() != &themem)
         return MSGinconsistent("UnitVectorBase");
      if (size() != 1)
         return MSGinconsistent("UnitVectorBase");
      if (max() != 1)
         return MSGinconsistent("UnitVectorBase");

      return SVectorBase<R>::isConsistent();
#else
      return true;
#endif
   }
   //@}
};


} // namespace soplex
#endif // _UNITVECTORBASE_H_
