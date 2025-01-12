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

/**@file  nameset.h
 * @brief Set of strings.
 */
#ifndef _NAMESET_H_
#define _NAMESET_H_

#include <assert.h>

#include "spxdefines.h"
#include "dataset.h"
#include "datahashtable.h"
#include "datakey.h"

namespace soplex
{
/**@class NameSet
   @brief   Set of strings.
   @ingroup Elementary

   Class NameSet implements a symbol or name table. It allows to store or
   remove names (i.e., \c char*), but does not provide means for manipulating
   stored names.
   
   Names in a NameSet may be accessed via numbers from 0 through num()-1
   and via \ref soplex::DataKey "DataKeys". See DataSet for a description of 
   these concepts.
   
   At a time a NameSet can hold a maximum of max() entries. This can be
   reset with method reMax(). If more than max() names are added to a
   NameSet, it adjusts itself automatically to the required size.  This
   implies, that references to names within a NameSet may become invalid if
   the NameSet is expanded.
   
   All names (i.e., the actual char strings) in a NameSet are stored in one
   continuous memory block of size memMax(). At one time memSize() bytes of
   it are used for actually saving names; the remaining memory is free to hold
   additional names. memRemax() can be used to reset memMax() but not lower
   than to memSize(). Method memPack() performs a garbage collection to
   gain free memory resulting from removed names.

   @warning Since the keys the NameSet uses to reference the strings are
            generated internally, it is extremly important that the calls
            to DataSet from within NameSet are synchronous to any calls
            outside to DataSet, such as in row or column adding.
*/
class NameSet
{
public:

   /**@brief   Handles of names in a NameSet.
    * @ingroup Elementary
    *
    *  Class Name provides the handles (i.e., char*s) of names in a
    *  NameSet.
    */   
   class Name
   {
   private:

      //------------------------------
      /**@name Private Data */
      //@{
      static const char deflt;     ///< default zero string.
      //@}

   public:

      //------------------------------
      /**@name Public Data */
      //@{
      const char *name;      ///< pointer to the name string.
      //@}

      //------------------------------
      /**@name Friends */
      //@{
      /// equality operator.
      friend int operator==(const Name& n1, const Name& n2)
      {
         return (strcmp (n1.name, n2.name) == 0);
      }
      /// output operator.
      friend std::ostream& operator<<(std::ostream& out, const Name& n)
      {
         return out << n.name;
      }
      //@}

      //------------------------------
      /**@name Debugging */
      //@{
      /// consistency check.
      bool isConsistent () const
      {
         return (name != 0);
      }
      //@}

      //------------------------------------
      /**@name Constructors / destructors */
      //@{
      /// default constructor.
      Name() 
         : name(&deflt)
      {}
      /// copy constructor.
      /** Only the pointer to the name is copied, but not the name itself.
       */
      Name (const Name& str) 
         : name(str.name)
      {}
      /// implictly constructs a Name out of a C style character string.
      Name (const char* str) 
         : name(str)
      {}
      //@}
   };

private:

   //------------------------------
   /**@name Data */
   //@{
   DataSet < int > set;  ///< name set.
   char* mem;            ///< string memory
   int memmax;           ///< size of string memory
   int memused;          ///< size of used string memory
   /** Every name in a NameSet is assigned a DataKey by which it can be
       accessed (see NameSet::operator[]()). See DataKey for a more
       detailed description of the concept of Keys.
   */
   DataHashTable < Name, DataKey > hashtab;  ///< hashtable for names
   //@}

public:

   //------------------------------
   /**@name Inquiry */
   //@{
   /// returns \p num 'th name of NameSet.
   const char* operator[](int pnum) const
   {
      return &mem[set[pnum]];
   }

   /// returns name for DataKey \p pkey of NameSet.
   const char* operator[](const DataKey& pkey) const
   {
      return &mem[set[pkey]];
   }

   /// returns nr. of names in NameSet.
   int num() const
   {
      return set.num();
   }

   /// returns maximum nr. of names that fit into NameSet.
   int max() const
   {
      return set.max();
   }

   /// returns maximum DataKey::idx used in NameSet.
   int size() const
   {
      return set.size();
   }

   /// returns maximum length of string memory.
   int memMax() const
   {
      return memmax;
   }

   /// returns used length of string memory.
   int memSize() const
   {
      return memused;
   }

   /// returns DataKey of the \p pnum 'th name in NameSet.
   DataKey key(int pnum) const
   {
      return set.key(pnum);
   }

   /// returns DataKey of name \p str in NameSet.
   DataKey key(const char* str) const
   {
      const Name nam(str);
      const DataKey* result = hashtab.get(nam);
      return result == 0 ? DataKey() : *hashtab.get(nam);
   }

   /// returns number of name with DataKey \p pkey in NameSet.
   int number(const DataKey& pkey) const
   {
      return set.number(pkey);
   }

   /// returns number of name \p str in NameSet.
   int number(const char *str) const
   {
      const Name nam(str);
      if (hashtab.has(nam))
      {
         assert(hashtab.get(nam) != 0);
         return number(*hashtab.get(nam));
      }
      else
         return -1;
   }

   /// does NameSet has a name with number \p pnum?
   bool has(int pnum) const
   {
      return set.has(pnum);
   }

   /// does NameSet has a name \p str?
   bool has(const char* str) const
   {
      const Name nam(str);
      return hashtab.has(nam);
   }

   /// does NameSet has a name with DataKey \p pkey?
   bool has(const DataKey& pkey) const
   {
      return set.has(pkey);
   }
   //@}

   //----------------------------
   /**@name Extension */
   //@{
   ///
   void add(const char* str);
   /// adds name \p str to NameSet.
   void add(DataKey& key, const char* str);

   ///
   void add(const NameSet& set);
   /// adds all names in \p set to NameSet.
   void add(DataKey key[], const NameSet& nset);
   //@}


   //----------------------------
   /**@name Shrinking */
   //@{
   /// removes name with DataKey \p key from NameSet.
   void remove(const DataKey& key);

   /// removes \p pnum 'th name from NameSet.
   void remove(int pnum)
   {
      remove(key(pnum));
   }

   /// removes name \p str from NameSet.
   void remove(const char* str);

   /// removes \p n names with DataKeys \p keys from NameSet.
   void remove(const DataKey keys[], int n);

   /// removes \p n names with numbers \p nums from NameSet.
   void remove(const int nums[], int n);

   /// remove all entries where \p dstat is less than zero.
   void remove(int dstat[]);

   /// removes all names from NameSet.
   void clear();
   //@}


   //----------------------------
   /**@name Memory Control */
   //@{
   /// resets max() to \p newmax.
   void reMax(int newmax = 0);

   /// resets memMax() to \p newmax.
   void memRemax(int newmax = 0);

   /// garbage collection.
   void memPack();
   //@}


   //----------------------------
   /**@name Control Parameters */
   //@{
   /// memory extension factor for entries.
   /** When more than max() names are added to a NameSet, it is
       automatically resized to fit the additional names. Parameter
       \p factor is the factor by which the element memory is extended to do
       so.
    */
   Real factor;

   /// memory extension factor for names.
   /** When the names added to a NameSet do no longer fit into the name
       memory (i.e. the memory for saving the strings), it is automatically
       resized to fit the additional names. Parameter \p memFactor is the
       factor by which this memory is extended to do so.
    */
   Real memFactor;
   //@}

   //----------------------------
   /**@name Miscellaneous */
   //@{
   /// consistency check.
   bool isConsistent() const;
   //@}

   //--------------------------------------
   /**@name Constructors / Destructors */
   //@{
   /// default constructor.
   /** @param      max     start value for max()
    *  @param      mmax    start value for memMax()
    *  @param      fac     start value for #factor
    *  @param      memFac  start value for #memFactor
    */
   explicit
   NameSet(int max = 10000,
           int mmax = -1,
           Real fac = 2,
           Real memFac = 2);

   /// destructor.
   ~NameSet();
   //@}

private:

   //--------------------------------------
   /**@name Blocked */
   //@{
   /// copy constructor.
   NameSet(const NameSet& old);

   /// assignment operator.
   NameSet& operator=(const NameSet& rhs);
   //@}
};

extern std::ostream& operator<<(std::ostream& s, const NameSet& nset);

} // namespace soplex
#endif // _NAMESET_H_
