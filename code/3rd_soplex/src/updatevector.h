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

/**@file  updatevector.h
 * @brief Dense vector with semi-sparse vector for updates
 */

#ifndef _UPDATEVECTOR_H_
#define _UPDATEVECTOR_H_

#include <assert.h>


#include "spxdefines.h"
#include "dvector.h"
#include "ssvector.h"

namespace soplex
{

/**@brief   Dense vector with semi-sparse vector for updates
   @ingroup Algebra

    In many algorithms vectors are updated in every iteration, by
    adding a multiple of another vector to it, i.e., given a vector \c
    x, a scalar \f$\alpha\f$ and another vector \f$\delta\f$, the
    update to \c x constists of substituting it by \f$x \leftarrow x +
    \alpha\cdot\delta\f$.
 
    While the update itself can easily be expressed with methods of
    the class Vector, it is often desirable to save the last update
    vector \f$\delta\f$ and value \f$\alpha\f$. This is provided by
    class UpdateVector.
 
    UpdateVectors are derived from DVector and provide additional
    methods for saving and setting the multiplicator \f$\alpha\f$ and
    the update vector \f$\delta\f$. Further, it allows for efficient
    sparse updates, by providing an IdxSet idx() containing the
    nonzero indices of \f$\delta\f$.  
*/
class UpdateVector : public DVector
{
private:

   //------------------------------------
   /**@name Data */
   //@{
   Real     theval;      ///< update multiplicator 
   SSVector thedelta;    ///< update vector
   //@}

public:

   //------------------------------------
   /**@name Constructors / destructors */
   //@{
   /// default constructor.
   explicit
   UpdateVector(int p_dim /*=0*/, Real p_eps /*=1e-16*/)
      : DVector (p_dim)
      , theval (0)
      , thedelta(p_dim, p_eps)
   {
      assert(isConsistent());
   }
   ///
   ~UpdateVector()
   {}
   /// copy constructor
   UpdateVector( const UpdateVector& );
   /// assignment from DVector
   UpdateVector& operator=(const DVector& rhs)
   {
      if ( this != & rhs )
         DVector::operator=(rhs);

      assert(isConsistent());

      return *this;
   }
   /// assignment from Vector
   UpdateVector& operator=(const Vector& rhs)
   {
      if ( this != & rhs )
         DVector::operator=(rhs);

      assert(isConsistent());

      return *this;
   }
   /// assignment
   UpdateVector& operator=(const UpdateVector& rhs);
   //@}

   //------------------------------------
   /**@name Access */
   //@{
   /// update multiplicator \f$\alpha\f$, writeable
   Real& value()
   {
      return theval;
   }
   /// update multiplicator \f$\alpha\f$
   Real value() const
   {
      return theval;
   }

   /// update vector \f$\delta\f$, writeable
   SSVector& delta()
   {
      return thedelta;
   }
   /// update vector \f$\delta\f$
   const SSVector& delta() const
   {
      return thedelta;
   }

   /// nonzero indices of update vector \f$\delta\f$
   const IdxSet& idx() const
   {
      return thedelta.indices();
   }
   //@}

   //------------------------------------
   /**@name Modification */
   //@{
   /// Perform the update
   /**  Add \c value() * \c delta() to the UpdateVector. Only the indices 
    *  in idx() are affected. For all other indices, delta() is asumed
    *  to be 0.
    */
   void update()
   {
      multAdd(theval, thedelta);
   }

   /// clear vector and update vector
   void clear()
   {
      DVector::clear();
      clearUpdate();
   }

   /// clear \f$\delta\f$, \f$\alpha\f$
   void clearUpdate()
   {
      thedelta.clear();
      theval = 0;
   }

   /// reset dimension
   void reDim(int newdim)
   {
      DVector::reDim(newdim);
      thedelta.reDim(newdim);
   }
   //@}

   //------------------------------------
   /**@name Consistency check */
   //@{
   /// 
   bool isConsistent() const;
   //@}
};


} // namespace soplex
#endif // _UPDATEVECTOR_H_
