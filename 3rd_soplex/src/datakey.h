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

/**@file  datakey.h
 * @brief Entry identifier class for items of a DataSet.
 */
#ifndef _DATAKEY_H_
#define _DATAKEY_H_

#include <assert.h>

namespace soplex
{
/**@brief   Entry identifier class for items of a DataSet.
   @ingroup Elementary

   Every item in a DataSet is assigned a DataKey by which it can be
   accessed (using DataSet::operator[]()). A DataKey consists of an integer
   member #idx, which is a positive number for any valid DataKey. No
   #idx of an element in a DataSet may exceed the sets max().
   This property may be used to build arrays with additional information to
   the elements of a DataSet.

   In addition, #DataKey%s provide a member #info which can be used to store 
   further information.
   
   Each DataKey is unique for one DataSet but different DataSets may (and
   generally will) manage the same #DataKey%s. When an element is removed from
   a DataSet its DataKey may (and generally will) be reused for other
   elements added to the DataSet later on.

   @todo data members should be private.
*/
class DataKey
{
public:

   //-------------------------------------
   /**@name Data */
   //@{
   /* This was originally implemented as bitfield "signed int info: 2; signed int idx: (8 * sizeof(int) - 2);",
      however, this seems to trigger a bug with old versions of GCC/glibc on 32bit machines. */
   int info;                                  ///< user information to store values -1, 0, +1
   int idx;                                   ///< (locally) unique key index
   //@}

public:

   //-------------------------------------
   /**@name Constructors / destructors */
   //@{
   /// Default constructor. Constructs an invalid DataKey.
   DataKey() 
      : info(0), idx(-1) 
   {}
   // Full constructor
   DataKey(int p_info, int p_idx)
      : info(p_info)
      , idx(p_idx)
   {
      assert( p_info <= 1 && p_info >= -1 );
   }
   /// Assignment operator.
   DataKey& operator=(const DataKey& rhs)
   {
      if ( this != &rhs ) {
         info = rhs.info;
         idx  = rhs.idx;
      }
      
      return *this;
   }
   /// Copy constructor.
   DataKey(const DataKey& old) 
      : info(old.info) 
      , idx(old.idx)
   {}
   //@}

   //-------------------------------------
   /**@name Access / modification */
   //@{
   /// gets the index number (\ref soplex::DataKey::idx "idx") of the DataKey.
   inline int getIdx() const
   {
      return idx;
   }
   /// sets the index number (\ref soplex::DataKey::idx "idx") of the DataKey.
   inline void setIdx(int p_idx) 
   {
      idx = p_idx;
   }   
   /// returns TRUE, iff the DataKey is valid.
   inline bool isValid() const
   {
      return idx >= 0;
   }
   /// makes the DataKey invalid and clears the \ref soplex::DataKey::info "info" field.
   inline void inValidate()
   {
      idx  = -1;
      info = 0;
   }
   //@}

};

} // namespace soplex
#endif // _DATAKEY_H_
