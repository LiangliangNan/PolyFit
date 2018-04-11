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

/**@file  rational.cpp
 * @brief Wrapper for GMP types.
 */

#ifndef SOPLEX_LEGACY
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <limits.h>


#include "rational.h"
#include "spxalloc.h"
#include "spxdefines.h"

#ifdef SOPLEX_WITH_GMP
#include "gmp.h"
#endif

/// use mpq_sgn instead of mpq_equal to compare rationals with zero
#define SOPLEX_PERFALT_1

/// turn on checking rationals for zero in arithmetic plus and minus operators
#define SOPLEX_PERFALT_2a

/// turn on checking rationals for zero, posone, negone in arithmetic multiplication and division operators
#define SOPLEX_PERFALT_2b

/// turn on checking rationals for zero, posone, negone before conversion to double and long double
#define SOPLEX_PERFALT_3

/// turn on checking for equality before assigning in operator=()
#define SOPLEX_PERFALT_4

/// turn on checking whether assignment in assignment operators is redundant for 0
#define SOPLEX_PERFALT_5a

/// turn on checking whether assignment in assignment operators is redundant for +1 and -1
#define SOPLEX_PERFALT_5b

namespace soplex
{

#ifdef SOPLEX_WITH_GMP

/// rational zero
const Rational Rational::ZERO(0, true);

/// rational plus one
const Rational Rational::POSONE(1, true);

/// rational minus one
const Rational Rational::NEGONE(-1, true);

/// list of unused Private objects; note that this cannot be used if SOPLEX_WITH_GMP is not defined, since then the
/// Private class has no member next() and prev()
/// should list memory be used?
#ifdef SOPLEX_NOLISTMEM
THREADLOCAL bool Rational::useListMem = false;
#else
THREADLOCAL bool Rational::useListMem = true;
#endif




/// list of unused Private objects

#ifdef SOPLEX_WITH_GMP
   THREADLOCAL IdList< Rational::Private > Rational::unusedPrivateList(0, 0, true);
#endif

/// Defines the "Pimpl"-class Private
class Rational::Private
{
public:

   mpq_t privatevalue;  ///< actual value of the Rational object
   Private* theprev;    ///< pointer to the previous element in the list
   Private* thenext;    ///< pointer to the next element in the list

   /// default constructor
   Private()
      : theprev(0)
      , thenext(0)
   {
      mpq_init(privatevalue);
   }

   /// copy constructor
   Private(const Private& p)
      : theprev(0)
      , thenext(0)
   {
      // a newly constructed element is not in any list, even if the original element (p) is; hence we initialize
      // theprev and thenext to zero
      mpq_init(privatevalue);
      mpq_set(this->privatevalue, p.privatevalue);
   }

   /// constructor from long double
   Private(const long double& r)
      : theprev(0)
      , thenext(0)
   {
      mpq_init(privatevalue);
      if( r == (long double)(1.0) )
         mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      else if( r == (long double)(-1.0) )
         mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      else if( r == (long double)(0.0) )
      {
         assert(mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) != 0);
      }
      else
         mpq_set_d(privatevalue, double(r));
   }

   /// constructor from double
   Private(const double& r)
      : theprev(0)
      , thenext(0)
   {
      mpq_init(privatevalue);
      if( r == 1.0 )
         mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      else if( r == -1.0 )
         mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      else if( r == 0.0 )
      {
         assert(mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) != 0);
      }
      else
         mpq_set_d(privatevalue, r);
   }

   /// constructor from int
   Private(const int& i)
      : theprev(0)
      , thenext(0)
   {
      mpq_init(privatevalue);
      if( i == 1 )
         mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      else if( i == -1 )
         mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      else if( i == 0 )
      {
         assert(mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) != 0);
      }
      else
         mpq_set_si(privatevalue, i, 1);
   }

   /// constructor from mpq_t
   Private(const mpq_t& q)
      : theprev(0)
      , thenext(0)
   {
      mpq_init(privatevalue);
      mpq_set(privatevalue, q);
   }

   /// destructor
   ~Private()
   {
      mpq_clear(privatevalue);
   }

   /// assignment operator
   Private& operator=(const Private& p)
   {
#ifdef SOPLEX_PERFALT_4
      if( mpq_equal(this->privatevalue, p.privatevalue) != 0 )
         return *this;
#endif

      // we only assign the value; the position in the list, i.e., theprev and thenext, must not be modified
      mpq_set(this->privatevalue, p.privatevalue);
      return *this;
   }

   /// assignment operator from long double
   Private& operator=(const long double& r)
   {
      // we only assign the value; the position in the list, i.e., theprev and thenext, must not be modified
      if( r == (long double)(0.0) )
      {
#ifdef SOPLEX_PERFALT_5a
#ifdef SOPLEX_PERFALT_1
         if( mpq_sgn(privatevalue) != 0 )
#else
         if( mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) == 0 )
#endif
#endif
            mpq_set(privatevalue, Rational::ZERO.dpointer->privatevalue);
      }
      else if( r == (long double)(1.0) )
      {
#ifdef SOPLEX_PERFALT_5b
         if( mpq_equal(privatevalue, Rational::POSONE.dpointer->privatevalue) == 0 )
#endif
            mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      }
      else if( r == (long double)(-1.0) )
      {
#ifdef SOPLEX_PERFALT_5b
         if( mpq_equal(privatevalue, Rational::NEGONE.dpointer->privatevalue) == 0 )
#endif
            mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      }
      else
         mpq_set_d(this->privatevalue, double(r));

      return *this;
   }

   /// assignment operator from double
   Private& operator=(const double& r)
   {
      // we only assign the value; the position in the list, i.e., theprev and thenext, must not be modified
      if( r == 0.0 )
      {
#ifdef SOPLEX_PERFALT_5a
#ifdef SOPLEX_PERFALT_1
         if( mpq_sgn(privatevalue) != 0 )
#else
         if( mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) == 0 )
#endif
#endif
            mpq_set(privatevalue, Rational::ZERO.dpointer->privatevalue);
      }
      else if( r == 1.0 )
      {
#ifdef SOPLEX_PERFALT_5b
         if( mpq_equal(privatevalue, Rational::POSONE.dpointer->privatevalue) == 0 )
#endif
            mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      }
      else if( r == -1.0 )
      {
#ifdef SOPLEX_PERFALT_5b
         if( mpq_equal(privatevalue, Rational::NEGONE.dpointer->privatevalue) == 0 )
#endif
            mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      }
      else
         mpq_set_d(privatevalue, r);

      return *this;
   }

   /// assignment operator from int
   Private& operator=(const int& i)
   {
      // we only assign the value; the position in the list, i.e., theprev and thenext, must not be modified
      if( i == 0 )
      {
#ifdef SOPLEX_PERFALT_5a
#ifdef SOPLEX_PERFALT_1
         if( mpq_sgn(privatevalue) != 0 )
#else
         if( mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) == 0 )
#endif
#endif
            mpq_set(privatevalue, Rational::ZERO.dpointer->privatevalue);
      }
      else if( i == 1 )
      {
#ifdef SOPLEX_PERFALT_5b
         if( mpq_equal(privatevalue, Rational::POSONE.dpointer->privatevalue) == 0 )
#endif
            mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      }
      else if( i == -1 )
      {
#ifdef SOPLEX_PERFALT_5b
         if( mpq_equal(privatevalue, Rational::NEGONE.dpointer->privatevalue) == 0 )
#endif
            mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      }
      else
         mpq_set_si(privatevalue, i, 1);

      return *this;
   }

   /// assignment operator from mpq_t
   Private& operator=(const mpq_t& q)
   {
#ifdef SOPLEX_PERFALT_4
   if( mpq_equal(this->privatevalue, q) != 0 )
      return *this;
#endif

      // we only assign the value; the position in the list, i.e., theprev and thenext, must not be modified
      mpq_set(this->privatevalue, q);
      return *this;
   }

   /// previous Private element
   Private*& prev()
   {
      return theprev;
   }

   /// previous Private element
   Private* const& prev() const
   {
      return theprev;
   }

   /// next Private element
   Private*& next()
   {
      return thenext;
   }

   /// next Private element
   Private* const& next() const
   {
      return thenext;
   }
};



/// special constructor only for initializing static rational variables; this is necessary since we need a constructor
/// for Rational::{ZERO, POSONE, NEGONE} that does not use these numbers
Rational::Rational(const int& i, const bool& dummy)
{
   dpointer = 0;
   spx_alloc(dpointer);
   new (dpointer) Private();
   mpq_set_si(dpointer->privatevalue, i, 1);

   assert(dpointer != 0);
}



/// default constructor
Rational::Rational()
{
   if( Rational::useListMem )
   {
      dpointer = unusedPrivateList.last();

      if( dpointer != 0 )
      {
         assert(unusedPrivateList.first() != 0);
         unusedPrivateList.remove(dpointer);
      }
      else
      {
         assert(unusedPrivateList.first() == 0);
         spx_alloc(dpointer);
         new (dpointer) Private();
      }
   }
   else
   {
      assert(unusedPrivateList.length() == 0);
      dpointer = 0;
      spx_alloc(dpointer);
      new (dpointer) Private();
   }

   assert(dpointer != 0);
}



/// copy constructor
Rational::Rational(const Rational& r)
{
   if( Rational::useListMem )
   {
      dpointer = unusedPrivateList.last();

      if( dpointer != 0 )
      {
         assert(unusedPrivateList.first() != 0);
         unusedPrivateList.remove(dpointer);
         *dpointer = *(r.dpointer);
      }
      else
      {
         assert(unusedPrivateList.first() == 0);
         spx_alloc(dpointer);
         new (dpointer) Private(*(r.dpointer));
      }
   }
   else
   {
      assert(unusedPrivateList.length() == 0);
      dpointer = 0;
      spx_alloc(dpointer);
      new (dpointer) Private(*(r.dpointer));
   }

   assert(dpointer != 0);
}



/// constructor from long double
Rational::Rational(const long double& r)
{
   if( Rational::useListMem )
   {
      dpointer = unusedPrivateList.last();

      if( dpointer != 0 )
      {
         assert(unusedPrivateList.first() != 0);
         unusedPrivateList.remove(dpointer);
         *dpointer = r;
      }
      else
      {
         assert(unusedPrivateList.first() == 0);
         spx_alloc(dpointer);
         new (dpointer) Private(r);
      }
   }
   else
   {
      assert(unusedPrivateList.length() == 0);
      dpointer = 0;
      spx_alloc(dpointer);
      new (dpointer) Private(r);
   }

   assert(dpointer != 0);
}



/// constructor from double
Rational::Rational(const double& r)
{
   if( Rational::useListMem )
   {
      dpointer = unusedPrivateList.last();

      if( dpointer != 0 )
      {
         assert(unusedPrivateList.first() != 0);
         unusedPrivateList.remove(dpointer);
         *dpointer = r;
      }
      else
      {
         assert(unusedPrivateList.first() == 0);
         spx_alloc(dpointer);
         new (dpointer) Private(r);
      }
   }
   else
   {
      assert(unusedPrivateList.length() == 0);
      dpointer = 0;
      spx_alloc(dpointer);
      new (dpointer) Private(r);
   }

   assert(dpointer != 0);
}



/// constructor from int
Rational::Rational(const int& i)
{
   if( Rational::useListMem )
   {
      dpointer = unusedPrivateList.last();

      if( dpointer != 0 )
      {
         assert(unusedPrivateList.first() != 0);
         unusedPrivateList.remove(dpointer);
         *dpointer = i;
      }
      else
      {
         assert(unusedPrivateList.first() == 0);
         spx_alloc(dpointer);
         new (dpointer) Private(i);
      }
   }
   else
   {
      assert(unusedPrivateList.length() == 0);
      dpointer = 0;
      spx_alloc(dpointer);
      new (dpointer) Private(i);
   }

   assert(dpointer != 0);
}



/// constructor from mpq_t
Rational::Rational(const mpq_t& q)
{
   if( Rational::useListMem )
   {
      dpointer = unusedPrivateList.last();

      if( dpointer != 0 )
      {
         assert(unusedPrivateList.first() != 0);
         unusedPrivateList.remove(dpointer);
         *dpointer = q;
      }
      else
      {
         assert(unusedPrivateList.first() == 0);
         spx_alloc(dpointer);
         new (dpointer) Private(q);
      }
   }
   else
   {
      assert(unusedPrivateList.length() == 0);
      dpointer = 0;
      spx_alloc(dpointer);
      new (dpointer) Private(q);
   }

   assert(dpointer != 0);
}



/// destructor
Rational::~Rational()
{
   assert(Rational::useListMem || unusedPrivateList.length() == 0);

   if( !Rational::useListMem || this == &Rational::ZERO || this == &Rational::POSONE || this == &Rational::NEGONE )
   {
      dpointer->~Private();
      spx_free(dpointer);
   }
   else
   {
      // for memory efficiency, we could free the Private object (or even more Private objects from the list of unused
      // elements) if there are much more unused than used Private objects; this requires counting the used Private
      // objects, though; we do not implement this currently, because we have not encountered memory problems, so far, and
      // because freeing costs time
      unusedPrivateList.append(dpointer);
   }
}



/// enables list memory
void Rational::enableListMem()
{
   assert(Rational::useListMem || unusedPrivateList.length() == 0);
   Rational::useListMem = true;
}



/// frees the unused rational elements in the memory list
/// frees the unused rational elements in the memory list
/** this can be useful when you want to save memory or needed when working with a GMP memory manager like the one
 *  in EGlib that frees GMP memory before the destructor of the static memory list is called; in most cases this
 *  method is optional; note that this does not free the Rational elements that are currently in use
 */
void Rational::freeListMem()
{
   unusedPrivateList.clear(true);
   assert(unusedPrivateList.length() == 0);
}



/// disables list memory
void Rational::disableListMem()
{
   Rational::freeListMem();
   Rational::useListMem = false;
}



/// assignment operator
Rational& Rational::operator=(const Rational &r)
{
   *(this->dpointer) = *(r.dpointer);
   return *this;
}



/// assignment operator from long double
Rational& Rational::operator=(const long double &r)
{
   *(this->dpointer) = r;
   return *this;
}



/// assignment operator from double
Rational& Rational::operator=(const double &r)
{
   *(this->dpointer) = r;
   return *this;
}



/// assignment operator from int
Rational& Rational::operator=(const int &i)
{
   *(this->dpointer) = i;
   return *this;
}



/// assignment operator from mpq_t
Rational& Rational::operator=(const mpq_t &q)
{
   *(this->dpointer) = q;
   return *this;
}



/// typecasts Rational to double (allows only explicit typecast)
Rational::operator double() const
{
#ifdef SOPLEX_PERFALT_3
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(this->dpointer->privatevalue) == 0 )
#else
   if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
#endif
      return 0.0;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      return 1.0;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
      return -1.0;
#endif

   return mpq_get_d(this->dpointer->privatevalue);
}



/// typecasts Rational to long double (allows only explicit typecast)
Rational::operator long double() const
{
#ifdef SOPLEX_PERFALT_3
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(this->dpointer->privatevalue) == 0 )
#else
   if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
#endif
      return 0.0;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      return 1.0;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
      return -1.0;
#endif

   return (long double)mpq_get_d(this->dpointer->privatevalue);
}



/// provides read-only access to underlying mpq_t
const mpq_t* Rational::getMpqPtr() const
{
   return &(this->dpointer->privatevalue);
}



/// provides read-only access to underlying mpq_t
const mpq_t& Rational::getMpqRef() const
{
   return this->dpointer->privatevalue;
}



/// provides write access to underlying mpq_t; use with care
mpq_t* Rational::getMpqPtr_w() const
{
   return &(this->dpointer->privatevalue);
}



/// provides write access to underlying mpq_t; use with care
mpq_t& Rational::getMpqRef_w() const
{
   return this->dpointer->privatevalue;
}



/// addition operator
Rational Rational::operator+(const Rational& r) const
{
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(r.dpointer->privatevalue) == 0 )
      return *this;
   else if( mpq_sgn(this->dpointer->privatevalue) == 0 )
      return r;
#else
   if( mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return r;
#endif
#endif

   Rational retval;
   mpq_add(retval.dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// addition assignment operator
Rational& Rational::operator+=(const Rational& r)
{
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(r.dpointer->privatevalue) == 0 )
      return *this;
   else if( mpq_sgn(this->dpointer->privatevalue) == 0 )
      return (*this = r);
#else
   if( mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return (*this = r);
#endif
#endif

   mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return *this;
}



/// addition operator for doubles
Rational Rational::operator+(const double& d) const
{
   if( d == 0.0 )
      return *this;
   else
   {
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return d;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return d;
#endif
#endif

      Rational retval(d);
      mpq_add(retval.dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
      return retval;
   }
}



/// addition assignment operator for doubles
Rational& Rational::operator+=(const double& d)
{
   if( d == 1.0 )
      return (*this += Rational::POSONE);
   else if( d == -1.0 )
      return (*this += Rational::NEGONE);
   else if( d != 0.0 )
   {
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return (*this = d);
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return (*this = d);
#endif
#endif

      Rational retval(d);
      mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
   }

   return *this;
}



/// addition operator for ints
Rational Rational::operator+(const int& d) const
{
   if( d == 0 )
      return *this;
   else
   {
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return d;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return d;
#endif
#endif

      Rational retval(d);
      mpq_add(retval.dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
      return retval;
   }
}



/// addition assignment operator for ints
Rational& Rational::operator+=(const int& d)
{
   if( d == 1 )
      return (*this += Rational::POSONE);
   else if( d == -1 )
      return (*this += Rational::NEGONE);
   else if( d != 0 )
   {
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return (*this = d);
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return (*this = d);
#endif
#endif

      Rational retval(d);
      mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
   }

   return *this;
}



/// subtraction operator
Rational Rational::operator-(const Rational& r) const
{
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(r.dpointer->privatevalue) == 0 )
      return *this;
   else if( mpq_sgn(this->dpointer->privatevalue) == 0 )
      return -r;
#else
   if( mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return -r;
#endif
#endif

   Rational retval;
   mpq_sub(retval.dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// subtraction assignment operator
Rational& Rational::operator-=(const Rational& r)
{
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(r.dpointer->privatevalue) == 0 )
      return *this;
   else if( mpq_sgn(this->dpointer->privatevalue) == 0 )
   {
      *this = r;
      *this *= -1;
      return *this;
   }
#else
   if( mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
   {
      *this = r;
      *this *= -1;
      return *this;
   }
#endif
#endif

   mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return *this;
}



/// subtraction operator for doubles
Rational Rational::operator-(const double& d) const
{
   if( d == 0.0 )
      return *this;
   else
   {
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return -d;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return -d;
#endif
#endif

      Rational retval(d);
      mpq_sub(retval.dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
      return retval;
   }
}



/// subtraction assignment operator for doubles
Rational& Rational::operator-=(const double& d)
{

   if( d == 1.0 )
      return (*this -= Rational::POSONE);
   else if( d == -1.0 )
      return (*this -= Rational::NEGONE);
   else if( d != 0.0 )
   {
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return (*this = -d);
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return (*this = -d);
#endif
#endif

      Rational retval(d);
      mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
   }

   return *this;
}



/// subtraction operator for ints
Rational Rational::operator-(const int& d) const
{
   if( d == 0 )
      return *this;
   else
   {
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return -d;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return -d;
#endif
#endif

      Rational retval(d);
      mpq_sub(retval.dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
      return retval;
   }
}



/// subtraction assignment operator for ints
Rational& Rational::operator-=(const int& d)
{

   if( d == 1 )
      return (*this -= Rational::POSONE);
   else if( d == -1 )
      return (*this -= Rational::NEGONE);
   else if( d != 0 )
   {
#ifdef SOPLEX_PERFALT_2a
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return (*this = -d);
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return (*this = -d);
#endif
#endif

      Rational retval(d);
      mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
   }

   return *this;
}



/// multiplication operator
Rational Rational::operator*(const Rational& r) const
{
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(r.dpointer->privatevalue) == 0 )
      return Rational::ZERO;
   else if( mpq_sgn(this->dpointer->privatevalue) == 0 )
      return Rational::ZERO;
#else
   if( mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return Rational::ZERO;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return Rational::ZERO;
#endif
   else if( mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      return *this;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      return r;
   else if( mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
      return -*this;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
      return -r;
#endif

   Rational retval;
   mpq_mul(retval.dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// multiplication assignment operator
Rational& Rational::operator*=(const Rational& r)
{
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(r.dpointer->privatevalue) == 0 )
      return (*this = Rational::ZERO);
   else if( mpq_sgn(this->dpointer->privatevalue) == 0 )
      return *this;
#else
   if( mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return (*this = Rational::ZERO);
   else if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
#endif
   else if( mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      return *this;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      return (*this = r);
   else if( mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
   {
      mpq_neg(this->dpointer->privatevalue, this->dpointer->privatevalue);
      return *this;
   }
   else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
   {
      mpq_neg(this->dpointer->privatevalue, r.dpointer->privatevalue);
      return *this;
   }
#endif

   mpq_mul(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return *this;
}



/// multiplication operator for doubles
Rational Rational::operator*(const double& d) const
{
   if( d == 0.0 )
      return Rational::ZERO;
   else if( d == 1.0 )
      return *this;
   else if( d == -1.0 )
   {
      Rational retval;
      mpq_neg(retval.dpointer->privatevalue, this->dpointer->privatevalue);
      return retval;
   }
   else
   {
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return Rational::ZERO;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return Rational::ZERO;
#endif
      else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
         return d;
      else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
         return -d;
#endif

      Rational retval(d);
      mpq_mul(retval.dpointer->privatevalue, retval.dpointer->privatevalue, this->dpointer->privatevalue);
      return retval;
   }
}



/// multiplication assignment operator for doubles
Rational& Rational::operator*=(const double& d)
{
   if( d == 0.0 )
      return (*this = Rational::ZERO);
   else if( d == 1.0 )
      return *this;
   else if( d == -1.0 )
   {
      mpq_neg(this->dpointer->privatevalue, this->dpointer->privatevalue);
      return *this;
   }
   else
   {
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return *this;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return *this;
#endif
      else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
         return (*this = d);
      else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
      {
         *this = d;
         mpq_neg(this->dpointer->privatevalue, this->dpointer->privatevalue);
         return *this;
      }
#endif

      Rational retval(d);
      mpq_mul(this->dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
      return *this;
   }
}



/// multiplication operator for ints
Rational Rational::operator*(const int& d) const
{
   if( d == 0 )
      return Rational::ZERO;
   else if( d == 1 )
      return *this;
   else if( d == -1 )
   {
      Rational retval;
      mpq_neg(retval.dpointer->privatevalue, this->dpointer->privatevalue);
      return retval;
   }
   else
   {
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return Rational::ZERO;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return Rational::ZERO;
#endif
      else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
         return d;
      else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
         return -d;
#endif

      Rational retval(d);
      mpq_mul(retval.dpointer->privatevalue, retval.dpointer->privatevalue, this->dpointer->privatevalue);
      return retval;
   }
}



/// multiplication assignment operator for ints
Rational& Rational::operator*=(const int& d)
{
   if( d == 0 )
      return (*this = Rational::ZERO);
   else if( d == 1 )
      return *this;
   else if( d == -1 )
   {
      mpq_neg(this->dpointer->privatevalue, this->dpointer->privatevalue);
      return *this;
   }
   else
   {
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return *this;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return *this;
#endif
      else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
         return (*this = d);
      else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
      {
         *this = d;
         mpq_neg(this->dpointer->privatevalue, this->dpointer->privatevalue);
         return *this;
      }
#endif

      Rational retval(d);
      mpq_mul(this->dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
      return *this;
   }
}



/// division operator
Rational Rational::operator/(const Rational& r) const
{
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(this->dpointer->privatevalue) == 0 )
      return Rational::ZERO;
#else
   if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return Rational::ZERO;
#endif
   else if( mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      return *this;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
   {
      Rational retval(r);
      retval.invert();
      return retval;
   }
   else if( mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
      return -*this;
   else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
   {
      Rational retval(r);
      retval.invert();
      mpq_neg(retval.dpointer->privatevalue, retval.dpointer->privatevalue);
      return retval;
   }
#endif

   Rational retval;
   mpq_div(retval.dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// division assignment operator
Rational& Rational::operator/=(const Rational& r)
{
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(this->dpointer->privatevalue) == 0 )
      return *this;
#else
   if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
#endif
   else if( mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      return *this;
   else if( mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
   {
      mpq_neg(this->dpointer->privatevalue, this->dpointer->privatevalue);
      return *this;
   }
   else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
   {
      mpq_inv(this->dpointer->privatevalue, r.dpointer->privatevalue);
      return *this;
   }
   else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
   {
      mpq_inv(this->dpointer->privatevalue, r.dpointer->privatevalue);
      mpq_neg(this->dpointer->privatevalue, this->dpointer->privatevalue);
      return *this;
   }
#endif

   mpq_div(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return *this;
}



/// division operator for doubles
Rational Rational::operator/(const double& d) const
{
   if( d == 1.0 )
      return *this;
   else if( d == -1.0 )
      return -(*this);
   else
   {
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return Rational::ZERO;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return Rational::ZERO;
#endif
      else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      {
         Rational retval(d);
         retval.invert();
         return retval;
      }
      else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
      {
         Rational retval(d);
         retval.invert();
         mpq_neg(retval.dpointer->privatevalue, retval.dpointer->privatevalue);
         return retval;
      }
#endif

      Rational retval(d);
      mpq_div(retval.dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
      return retval;
   }
}



/// division assignment operator for doubles
Rational& Rational::operator/=(const double& d)
{
   if( d == 1.0 )
      return *this;
   else if( d == -1.0 )
   {
      mpq_neg(this->dpointer->privatevalue, this->dpointer->privatevalue);
      return *this;
   }
   else
   {
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return *this;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return *this;
#endif
      else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      {
         *this = d;
         this->invert();
         return *this;
      }
      else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
      {
         *this = -d;
         this->invert();
         return *this;
      }
#endif

      Rational retval(d);
      mpq_div(this->dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
      return *this;
   }
}



/// division operator for ints
Rational Rational::operator/(const int& d) const
{
   if( d == 1 )
      return *this;
   else if( d == -1 )
      return -(*this);
   else
   {
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return Rational::ZERO;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return Rational::ZERO;
#endif
      else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      {
         Rational retval(d);
         retval.invert();
         return retval;
      }
      else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
      {
         Rational retval(d);
         retval.invert();
         mpq_neg(retval.dpointer->privatevalue, retval.dpointer->privatevalue);
         return retval;
      }
#endif

      Rational retval(d);
      mpq_div(retval.dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
      return retval;
   }
}



/// division assignment operator for ints
Rational& Rational::operator/=(const int& d)
{
   if( d == 1 )
      return *this;
   else if( d == -1 )
   {
      mpq_neg(this->dpointer->privatevalue, this->dpointer->privatevalue);
      return *this;
   }
   else
   {
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
      if( mpq_sgn(this->dpointer->privatevalue) == 0 )
         return *this;
#else
      if( mpq_equal(this->dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
         return *this;
#endif
      else if( mpq_equal(this->dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
      {
         *this = d;
         this->invert();
         return *this;
      }
      else if( mpq_equal(this->dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
      {
         *this = -d;
         this->invert();
         return *this;
      }
#endif

      Rational retval(d);
      mpq_div(this->dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
      return *this;
   }
}



/// add product of two rationals
Rational& Rational::addProduct(const Rational& r, const Rational& s)
{
#ifdef SOPLEX_PERFALT_2b
   if( mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
   {
      mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, s.dpointer->privatevalue);
      return *this;
   }
   else if( mpq_equal(s.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
   {
      mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
      return *this;
   }
   else if( mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
   {
      mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, s.dpointer->privatevalue);
      return *this;
   }
   else if( mpq_equal(s.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
   {
      mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
      return *this;
   }
#if 0 // currently, SoPlex calls this method only with nonzero r and s, hence we do not check this case
#ifdef SOPLEX_PERFALT_1
   else if( mpq_sgn(r.dpointer->privatevalue) == 0 )
      return *this;
   else if( mpq_sgn(s.dpointer->privatevalue) == 0 )
      return *this;
#else
   else if( mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
   else if( mpq_equal(s.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
#endif
#endif
#endif

   Rational product(r);
   mpq_mul(product.dpointer->privatevalue, product.dpointer->privatevalue, s.dpointer->privatevalue);
   mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, product.dpointer->privatevalue);
   return *this;
}



/// subtract product of two rationals
Rational& Rational::subProduct(const Rational& r, const Rational& s)
{
#ifdef SOPLEX_PERFALT_2b
   if( mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
   {
      mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, s.dpointer->privatevalue);
      return *this;
   }
   else if( mpq_equal(s.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
   {
      mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
      return *this;
   }
   else if( mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
   {
      mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, s.dpointer->privatevalue);
      return *this;
   }
   else if( mpq_equal(s.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
   {
      mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
      return *this;
   }
#if 0 // currently, SoPlex calls this method only with nonzero r and s, hence we do not check this case
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(r.dpointer->privatevalue) == 0 )
      return *this;
   else if( mpq_sgn(s.dpointer->privatevalue) == 0 )
      return *this;
#else
   if( mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
   else if( mpq_equal(s.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
#endif
#endif
#endif

   Rational product(r);
   mpq_mul(product.dpointer->privatevalue, product.dpointer->privatevalue, s.dpointer->privatevalue);
   mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, product.dpointer->privatevalue);
   return *this;
}



/// add quotient of two rationals, r divided by s
Rational& Rational::addQuotient(const Rational& r, const Rational& s)
{
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(r.dpointer->privatevalue) == 0 )
      return *this;
#else
   if( mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
#endif
   else if( mpq_equal(s.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
   {
      mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
      return *this;
   }
   else if( mpq_equal(s.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
   {
      mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
      return *this;
   }
#endif

   Rational quotient(r);
   mpq_div(quotient.dpointer->privatevalue, quotient.dpointer->privatevalue, s.dpointer->privatevalue);
   mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, quotient.dpointer->privatevalue);
   return *this;
}



/// subtract quotient of two rationals, r divided by s
Rational& Rational::subQuotient(const Rational& r, const Rational& s)
{
#ifdef SOPLEX_PERFALT_2b
#ifdef SOPLEX_PERFALT_1
   if( mpq_sgn(r.dpointer->privatevalue) == 0 )
      return *this;
#else
   if( mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0 )
      return *this;
#endif
   else if( mpq_equal(s.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0 )
   {
      mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
      return *this;
   }
   else if( mpq_equal(s.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0 )
   {
      mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
      return *this;
   }
#endif

   Rational quotient(r);
   mpq_div(quotient.dpointer->privatevalue, quotient.dpointer->privatevalue, s.dpointer->privatevalue);
   mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, quotient.dpointer->privatevalue);
   return *this;
}



/// inversion
Rational& Rational::invert()
{
   mpq_inv(this->dpointer->privatevalue, this->dpointer->privatevalue);
   return *this;
}



/// round up to next power of two
Rational& Rational::powRound()
{
   mpz_t roundval;

   MSG_DEBUG( std::cout << "rounding " << rationalToString(this->dpointer->privatevalue) << " to power of two" << "\n" );

   mpz_init(roundval);
   mpz_cdiv_q(roundval, mpq_numref(this->dpointer->privatevalue), mpq_denref(this->dpointer->privatevalue));
   mpz_sub_ui(roundval, roundval, 1);

   MSG_DEBUG( std::cout << "   --> " << mpz_get_str(0, 10, roundval) << "\n" );

   size_t binlog = mpz_sizeinbase(roundval, 2);

   MSG_DEBUG( std::cout << "   --> 2^" << binlog << "\n" );

   mpz_ui_pow_ui(roundval, 2, binlog);

   MSG_DEBUG( std::cout << "   --> " << mpz_get_str(0, 10, roundval) << "\n" );

   mpq_set_z(this->dpointer->privatevalue, roundval);
   mpz_clear(roundval);

   MSG_DEBUG( std::cout << "   --> " << rationalToString(this->dpointer->privatevalue) << "\n" );

   return *this;
}



/// checks if d is the closest possible double
bool Rational::isNextTo(const double& d)
{
   // get intervall [a,b] of doubles that the Rational is in
   double x = mpq_get_d(this->dpointer->privatevalue);
   double a;
   double b;

   if( Rational(x) < *this )
   {
      a = x;
      b = (double)spxNextafter(a, infinity);
   }
   else
   {
      b = x;
      a = (double)spxNextafter(b, -infinity);
   }

   // check if d equals the closer end of the intervall
   bool result = (spxAbs(*this - a) < spxAbs(*this - b))
      ? (d == a)
      : (d == b);

   return result;
}



/// checks if d is exactly equal to the Rational and if not, if it is one of the two adjacent doubles
bool Rational::isAdjacentTo(const double& d) const
{
   double x = mpq_get_d(this->dpointer->privatevalue);
   double a;
   double b;
   mpq_t tmp;

   mpq_init(tmp);
   mpq_set_d(tmp, x);

   int cmp = mpq_cmp(tmp, this->dpointer->privatevalue);
   mpq_clear(tmp);

   // the rounded value is smaller than the rational value
   if( cmp < 0 )
   {
      a = x;
      b = (double)spxNextafter(a, infinity);
   }
   // the rounded value is larger than the rational value
   else if( cmp > 0 )
   {
      b = x;
      a = (double)spxNextafter(b, -infinity);
   }
   // the rational value is representable in double precision
   else
      return (x == d);

   return ((a == d) || (b == d));
}



/// Size in specified base (bit size for base 2)
int Rational::sizeInBase(const int base) const
{
   return (int)mpz_sizeinbase(mpq_numref(this->dpointer->privatevalue), base)
      + (int)mpz_sizeinbase(mpq_denref(this->dpointer->privatevalue), base);
}



/// returns precision of Rational implementation, i.e., number of bits used to store Rational numbers (INT_MAX if exact)
int Rational::precision()
{
   return INT_MAX;
}



#define MAX_STR_LEN 10000
/// read Rational from string
bool Rational::readString(const char* s)
{
   assert(s != 0);
   assert(strlen(s) <= MAX_STR_LEN);
   Rational value;
   const char* pos;

   // if there is a slash or there is no dot and exponent (i.e. we
   // have an integer), we may simply call GMP's string reader
   if( strchr(s, '/') != 0 || strpbrk(s, ".eE") == 0 )
   {
      pos = (*s == '+') ? s + 1 : s;
      if( mpq_set_str(value.dpointer->privatevalue, pos, 10) == 0 )
      {
         mpq_canonicalize(value.dpointer->privatevalue);
         mpq_set(this->dpointer->privatevalue, value.dpointer->privatevalue);
         return true;
      }
      else
         return false;
   }

   // otherwise we analyze the string
#ifndef NDEBUG
   bool has_exponent = false;
   bool has_dot = false;
#endif
   bool has_digits = false;
   bool has_emptyexponent = false;
   long int exponent = 0;
   long int decshift = 0;
   mpz_t shiftpower;
   mpz_init(shiftpower);
   mpq_t shiftpowerRational;
   mpq_init(shiftpowerRational);
   char* t;
   char tmp[MAX_STR_LEN];

   pos = s;

   // 1. sign
   if( (*pos == '+') || (*pos == '-') )
      pos++;

   // 2. Digits before the decimal dot
   while( (*pos >= '0') && (*pos <= '9') )
   {
      has_digits = true;
      pos++;
   }

   // 3. Decimal dot
   if( *pos == '.' )
   {
#ifndef NDEBUG
      has_dot = true;
#endif
      pos++;

      // 4. If there was a dot, possible digit behind it
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_digits = true;
         pos++;
      }
   }

   // 5. Exponent
   if( tolower(*pos) == 'e' )
   {
#ifndef NDEBUG
      has_exponent = true;
#endif
      has_emptyexponent = true;
      pos++;

      // 6. Exponent sign
      if( (*pos == '+') || (*pos == '-') )
         pos++;

      // 7. Exponent digits
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_emptyexponent = false;
         pos++;
      }
   }

   if( has_emptyexponent || !has_digits )
      return false;

   assert(has_exponent || has_dot);

   //read up to dot recording digits
   t = tmp;
   pos = s;

   if( *pos == '+' )
      pos++;

   while( ((*pos >= '0') && (*pos <= '9') ) || *pos == '+' || *pos == '-'  )
   {
      *t++ = *pos;
      pos++;
   }
   //record digits after dot, recording positions
   decshift = 0;
   if( *pos == '.' )
   {
      assert(has_dot);
      pos++;
      while( (*pos >= '0') && (*pos <= '9') )
      {
         *t++ = *pos;
         decshift++;
         pos++;
      }
   }
   *t = '\0';

   if( mpq_set_str(value.dpointer->privatevalue, tmp, 10) != 0)
      return false;
   mpq_canonicalize(value.dpointer->privatevalue);

   //record exponent and update final result
   exponent = -decshift;
   if( tolower(*pos) == 'e' )
   {
      pos++;
      assert(has_exponent);
      for( t = tmp; *pos != '\0'; pos++ )
         *t++ = *pos;
      *t = '\0';
      exponent += atol(tmp);
   }
   if( exponent > 0 )
   {
      mpz_ui_pow_ui(shiftpower, 10, exponent);
      mpq_set_z(shiftpowerRational, shiftpower);
      mpq_mul(value.dpointer->privatevalue, value.dpointer->privatevalue, shiftpowerRational);
   }
   else if( exponent < 0 )
   {
      mpz_ui_pow_ui(shiftpower, 10, -exponent);
      mpq_set_z(shiftpowerRational, shiftpower);
      mpq_div(value.dpointer->privatevalue, value.dpointer->privatevalue, shiftpowerRational);
   }

   mpq_canonicalize(value.dpointer->privatevalue);
   mpq_set(this->dpointer->privatevalue, value.dpointer->privatevalue);
   mpz_clear(shiftpower);
   mpq_clear(shiftpowerRational);
   return true;

}



/// convert rational number to string
std::string rationalToString(const Rational& r, const int precision)
{

#if defined(_WIN32) || defined(_WIN64) || defined(__APPLE__)
  std::stringstream sstream;
  sstream << r;
  return sstream.str();
#else
   if( precision <= 0 )
   {
      std::stringstream sstream;
      sstream << r;
      return sstream.str();
   }
   else
   {
      mpf_t tmpFloat;
      char tmpString[64];
      FILE* tmpStream;

      tmpStream = fmemopen(tmpString, 63, "w");
      mpf_init2(tmpFloat, 256);
      mpf_set_q(tmpFloat, r.dpointer->privatevalue);
      mpf_out_str(tmpStream, 10, precision, tmpFloat);
      mpf_clear(tmpFloat);

      fflush(tmpStream);
      std::string retString = std::string(tmpString);
      fclose(tmpStream);
      return retString;
   }
#endif
}



/// read Rational from string
bool readStringRational(const char* s, Rational& value)
{
   assert(s != 0);
   assert(strlen(s) <= MAX_STR_LEN);
   const char* pos;

   // if there is a slash or there is no dot and exponent (i.e. we
   // have an integer), we may simply call GMP's string reader
   if( strchr(s, '/') != 0 || strpbrk(s, ".eE") == 0 )
   {
      pos = (*s == '+') ? s + 1 : s;
      if( mpq_set_str(value.dpointer->privatevalue, pos, 10) == 0 )
      {
         mpq_canonicalize(value.dpointer->privatevalue);
         return true;
      }
      else
         return false;
   }

   // otherwise we analyze the string
#ifndef NDEBUG
   bool has_exponent = false;
   bool has_dot = false;
#endif
   bool has_digits = false;
   bool has_emptyexponent = false;
   long int exponent = 0;
   long int decshift = 0;
   mpz_t shiftpower;
   mpz_init(shiftpower);
   mpq_t shiftpowerRational;
   mpq_init(shiftpowerRational);
   char* t;
   char tmp[MAX_STR_LEN];

   pos = s;

   // 1. sign
   if( (*pos == '+') || (*pos == '-') )
      pos++;

   // 2. Digits before the decimal dot
   while( (*pos >= '0') && (*pos <= '9') )
   {
      has_digits = true;
      pos++;
   }

   // 3. Decimal dot
   if( *pos == '.' )
   {
#ifndef NDEBUG
      has_dot = true;
#endif
      pos++;

      // 4. If there was a dot, possible digit behind it
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_digits = true;
         pos++;
      }
   }

   // 5. Exponent
   if( tolower(*pos) == 'e' )
   {
#ifndef NDEBUG
      has_exponent = true;
#endif
      has_emptyexponent = true;
      pos++;

      // 6. Exponent sign
      if( (*pos == '+') || (*pos == '-') )
         pos++;

      // 7. Exponent digits
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_emptyexponent = false;
         pos++;
      }
   }

   if( has_emptyexponent || !has_digits )
      return false;

   assert(has_exponent || has_dot);

   // read up to dot recording digits
   t = tmp;
   pos = s;

   if( *pos == '+' )
      pos++;

   while( ((*pos >= '0') && (*pos <= '9') ) || *pos == '+' || *pos == '-'  )
   {
      *t++ = *pos;
      pos++;
   }
   // record digits after dot, recording positions
   decshift = 0;
   if( *pos == '.' )
   {
      assert(has_dot);
      pos++;
      while( (*pos >= '0') && (*pos <= '9') )
      {
         *t++ = *pos;
         decshift++;
         pos++;
      }
   }
   *t = '\0';

   if( mpq_set_str(value.dpointer->privatevalue, tmp, 10) != 0)
      return false;

   mpq_canonicalize(value.dpointer->privatevalue);

   // record exponent and update final result
   exponent = -decshift;
   if( tolower(*pos) == 'e' )
   {
      pos++;
      assert(has_exponent);
      for( t = tmp; *pos != '\0'; pos++ )
         *t++ = *pos;
      *t = '\0';
      exponent += atol(tmp);
   }
   if( exponent > 0 )
   {
      mpz_ui_pow_ui(shiftpower, 10, exponent);
      mpq_set_z(shiftpowerRational, shiftpower);
      mpq_mul(value.dpointer->privatevalue, value.dpointer->privatevalue, shiftpowerRational);
   }
   else if( exponent < 0 )
   {
      mpz_ui_pow_ui(shiftpower, 10, -exponent);
      mpq_set_z(shiftpowerRational, shiftpower);
      mpq_div(value.dpointer->privatevalue, value.dpointer->privatevalue, shiftpowerRational);
   }

   mpq_canonicalize(value.dpointer->privatevalue);
   mpz_clear(shiftpower);
   mpq_clear(shiftpowerRational);

   return true;
}



/// print Rational
std::ostream& operator<<(std::ostream& os, const Rational& r)
{
   char* buffer;
   buffer = (char*) malloc (mpz_sizeinbase(mpq_numref(r.dpointer->privatevalue), 10) + mpz_sizeinbase(mpq_denref(r.dpointer->privatevalue), 10) + 3);
   os << mpq_get_str(buffer, 10, r.dpointer->privatevalue);
   free(buffer);
   return os;
}



/// comparison operator returning a positive value if r > s, zero if r = s, and a negative value if r < s
int compareRational(const Rational& r, const Rational& s)
{
   return mpq_cmp(r.dpointer->privatevalue, s.dpointer->privatevalue);
}



/// equality operator
bool operator==(const Rational& r, const Rational& s)
{
   return (mpq_equal(r.dpointer->privatevalue, s.dpointer->privatevalue) != 0);
}



/// inequality operator
bool operator!=(const Rational& r, const Rational& s)
{
   return (mpq_equal(r.dpointer->privatevalue, s.dpointer->privatevalue) == 0);
}



/// less than operator
bool operator<(const Rational& r, const Rational& s)
{
   return (mpq_cmp(r.dpointer->privatevalue, s.dpointer->privatevalue) < 0);
}



/// less than or equal to operator
bool operator<=(const Rational& r, const Rational& s)
{
   return (mpq_cmp(r.dpointer->privatevalue, s.dpointer->privatevalue) <= 0);
}



/// greater than operator
bool operator>(const Rational& r, const Rational& s)
{
   return (mpq_cmp(r.dpointer->privatevalue, s.dpointer->privatevalue) > 0);
}



/// greater than or equal to operator
bool operator>=(const Rational& r, const Rational& s)
{
   return (mpq_cmp(r.dpointer->privatevalue, s.dpointer->privatevalue) >= 0);
}



/// equality operator for Rational and double
bool operator==(const Rational& r, const double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) == 0);
#else
      return (mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0);
#endif
   else if( s == 1.0 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0);
   else if( s == -1.0 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0);
   else
      return (r == Rational(s));
}



/// inequality operator for Rational and double
bool operator!=(const Rational& r, const double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) != 0);
#else
      return (mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) == 0);
#endif
   else if( s == 1.0 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) == 0);
   else if( s == -1.0 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) == 0);
   else
      return (r != Rational(s));
}



/// less than operator for Rational and double
bool operator<(const Rational& r, const double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) == -1);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) < 0);
#endif
   else if( s == 1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) < 0);
   else if( s == -1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) < 0);
   else
      return (r < Rational(s));
}



/// less than or equal to operator for Rational and double
bool operator<=(const Rational& r, const double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) <= 0);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) <= 0);
#endif
   else if( s == 1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) <= 0);
   else if( s == -1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) <= 0);
   else
      return (r <= Rational(s));
}



/// greater than operator for Rational and double
bool operator>(const Rational& r, const double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) == 1);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) > 0);
#endif
   else if( s == 1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) > 0);
   else if( s == -1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) > 0);
   else
      return (r > Rational(s));
}



/// greater than or equal to operator for Rational and double
bool operator>=(const Rational& r, const double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) >= 0);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) >= 0);
#endif
   else if( s == 1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) >= 0);
   else if( s == -1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) >= 0);
   else
      return (r >= Rational(s));
}



/// equality operator for double and Rational
bool operator==(const double& r, const Rational& s)
{
   return (s == r);
}



/// inequality operator double and Rational
bool operator!=(const double& r, const Rational& s)
{
   return (s != r);
}



/// less than operator double and Rational
bool operator<(const double& r, const Rational& s)
{
   return (s > r);
}



/// less than or equal to operator double and Rational
bool operator<=(const double& r, const Rational& s)
{
   return (s >= r);
}



/// greater than operator double and Rational
bool operator>(const double& r, const Rational& s)
{
   return (s < r);
}



/// greater than or equal to operator double and Rational
bool operator>=(const double& r, const Rational& s)
{
   return (s <= r);
}



/// equality operator for Rational and long double
bool operator==(const Rational& r, const long double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) == 0);
#else
      return (mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0);
#endif
   else if( s == 1.0 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0);
   else if( s == -1.0 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0);
   else
      return (r == Rational(s));
}



/// inequality operator for Rational and long double
bool operator!=(const Rational& r, const long double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) != 0);
#else
      return (mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) == 0);
#endif
   else if( s == 1.0 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) == 0);
   else if( s == -1.0 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) == 0);
   else
      return (r != Rational(s));
}



/// less than operator for Rational and long double
bool operator<(const Rational& r, const long double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) == -1);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) < 0);
#endif
   else if( s == 1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) < 0);
   else if( s == -1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) < 0);
   else
      return (r < Rational(s));
}



/// less than or equal to operator for Rational and long double
bool operator<=(const Rational& r, const long double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) <= 0);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) <= 0);
#endif
   else if( s == 1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) <= 0);
   else if( s == -1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) <= 0);
   else
      return (r <= Rational(s));
}



/// greater than operator for Rational and long double
bool operator>(const Rational& r, const long double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) == 1);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) > 0);
#endif
   else if( s == 1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) > 0);
   else if( s == -1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) > 0);
   else
      return (r > Rational(s));
}



/// greater than or equal to operator for Rational and long double
bool operator>=(const Rational& r, const long double& s)
{
   if( s == 0.0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) >= 0);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) >= 0);
#endif
   else if( s == 1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) >= 0);
   else if( s == -1.0 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) >= 0);
   else
      return (r >= Rational(s));
}



/// equality operator for long double and Rational
bool operator==(const long double& r, const Rational& s)
{
   return (s == r);
}



/// inequality operator for long double and Rational
bool operator!=(const long double& r, const Rational& s)
{
   return (s != r);
}



/// less than operator for long double and Rational
bool operator<(const long double& r, const Rational& s)
{
   return (s > r);
}



/// less than or equal to operator for long double and Rational
bool operator<=(const long double& r, const Rational& s)
{
   return (s >= r);
}



/// greater than operator for long double and Rational
bool operator>(const long double& r, const Rational& s)
{
   return (s < r);
}



/// greater than or equal to operator for long double and Rational
bool operator>=(const long double& r, const Rational& s)
{
   return (s <= r);
}



/// addition operator for double and Rational
Rational operator+(const double& d, const Rational& r)
{
   return (r + d);
}



/// subtraction operator for double and Rational
Rational operator-(const double& d, const Rational& r)
{
   Rational res(r);
   mpq_neg(res.dpointer->privatevalue, res.dpointer->privatevalue);
   res += d;
   return res;
}



/// multiplication operator for double and Rational
Rational operator*(const double& d, const Rational& r)
{
   if( d == 0.0 )
      return Rational::ZERO;
   else if( d == 1.0 )
      return r;
   else if( d == -1.0 )
      return -r;
   else
   {
      Rational retval(d);
      mpq_mul(retval.dpointer->privatevalue, retval.dpointer->privatevalue, r.dpointer->privatevalue);
      return retval;
   }
}



/// division operator for double and Rational
Rational operator/(const double& d, const Rational& r)
{
   if( d == 0.0 )
      return Rational::ZERO;
   else if( d == 1.0 )
   {
      Rational retval(r);
      retval.invert();
      return retval;
   }
   else if( d == -1.0 )
   {
      Rational retval(r);
      retval.invert();
      mpq_neg(retval.dpointer->privatevalue, retval.dpointer->privatevalue);
      return retval;
   }
   else
   {
      Rational retval(d);
      mpq_div(retval.dpointer->privatevalue, retval.dpointer->privatevalue, r.dpointer->privatevalue);
      return retval;
   }
}



/// equality operator for Rational and int
bool operator==(const Rational& r, const int& s)
{
   if( s == 0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) == 0);
#else
      return (mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) != 0);
#endif
   else if( s == 1 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) != 0);
   else if( s == -1 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) != 0);
   else
      return (r == Rational(s));
}



/// inequality operator for Rational and int
bool operator!=(const Rational& r, const int& s)
{
   if( s == 0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) != 0);
#else
      return (mpq_equal(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) == 0);
#endif
   else if( s == 1 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) == 0);
   else if( s == -1 )
      return (mpq_equal(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) == 0);
   else
      return (r != Rational(s));
}



/// less than operator for Rational and int
bool operator<(const Rational& r, const int& s)
{
   if( s == 0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) == -1);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) < 0);
#endif
   else if( s == 1 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) < 0);
   else if( s == -1 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) < 0);
   else
      return (r < Rational(s));
}



/// less than or equal to operator for Rational and int
bool operator<=(const Rational& r, const int& s)
{
   if( s == 0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) <= 0);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) <= 0);
#endif
   else if( s == 1 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) <= 0);
   else if( s == -1 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) <= 0);
   else
      return (r <= Rational(s));
}



/// greater than operator for Rational and int
bool operator>(const Rational& r, const int& s)
{
   if( s == 0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) == 1);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) > 0);
#endif
   else if( s == 1 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) > 0);
   else if( s == -1 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) > 0);
   else
      return (r > Rational(s));
}



/// greater than or equal to operator for Rational and int
bool operator>=(const Rational& r, const int& s)
{
   if( s == 0 )
#ifdef SOPLEX_PERFALT_1
      return (mpq_sgn(r.dpointer->privatevalue) >= 0);
#else
      return (mpq_cmp(r.dpointer->privatevalue, Rational::ZERO.dpointer->privatevalue) >= 0);
#endif
   else if( s == 1 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::POSONE.dpointer->privatevalue) >= 0);
   else if( s == -1 )
      return (mpq_cmp(r.dpointer->privatevalue, Rational::NEGONE.dpointer->privatevalue) >= 0);
   else
      return (r >= Rational(s));
}



/// equality operator for int and Rational
bool operator==(const int& r, const Rational& s)
{
   return (s == r);
}



/// inequality operator int and Rational
bool operator!=(const int& r, const Rational& s)
{
   return (s != r);
}



/// less than operator int and Rational
bool operator<(const int& r, const Rational& s)
{
   return (s > r);
}



/// less than or equal to operator int and Rational
bool operator<=(const int& r, const Rational& s)
{
   return (s >= r);
}



/// greater than operator int and Rational
bool operator>(const int& r, const Rational& s)
{
   return (s < r);
}



/// greater than or equal to operator int and Rational
bool operator>=(const int& r, const Rational& s)
{
   return (s <= r);
}



/// addition operator for int and Rational
Rational operator+(const int& d, const Rational& r)
{
   return (r + d);
}



/// subtraction operator for int and Rational
Rational operator-(const int& d, const Rational& r)
{
   Rational res(r);
   mpq_neg(res.dpointer->privatevalue, res.dpointer->privatevalue);
   res += d;
   return res;
}



/// multiplication operator for int and Rational
Rational operator*(const int& d, const Rational& r)
{
   if( d == 0 )
      return Rational::ZERO;
   else if( d == 1 )
      return r;
   else if( d == -1 )
      return -r;
   else
   {
      Rational retval(d);
      mpq_mul(retval.dpointer->privatevalue, retval.dpointer->privatevalue, r.dpointer->privatevalue);
      return retval;
   }
}



/// division operator for int and Rational
Rational operator/(const int& d, const Rational& r)
{
   if( d == 0 )
      return Rational::ZERO;
   else if( d == 1 )
   {
      Rational retval(r);
      retval.invert();
      return retval;
   }
   else if( d == -1 )
   {
      Rational retval(r);
      retval.invert();
      mpq_neg(retval.dpointer->privatevalue, retval.dpointer->privatevalue);
      return retval;
   }
   else
   {
      Rational retval(d);
      mpq_div(retval.dpointer->privatevalue, retval.dpointer->privatevalue, r.dpointer->privatevalue);
      return retval;
   }
}



/// Absolute.
Rational spxAbs(const Rational& r)
{
   Rational res;
   mpq_abs(res.dpointer->privatevalue, r.dpointer->privatevalue);
   return res;
}



/// Sign function; returns 1 if r > 0, 0 if r = 0, and -1 if r < 0.
int sign(const Rational& r)
{
      return mpq_sgn(r.dpointer->privatevalue);
}



/// Negation.
Rational operator-(const Rational& r)
{
   Rational res;
   mpq_neg(res.dpointer->privatevalue, r.dpointer->privatevalue);
   return res;
}



/// Total size of rational vector.
int totalSizeRational(const Rational* vector, const int length, const int base)
{
   assert(vector != 0);
   assert(length >= 0);
   assert(base >= 0);

   int size = 0;

   for( int i = 0; i < length; i++ )
      size += vector[i].sizeInBase(base);

   return size;
}



/// Size of least common multiple of denominators in rational vector.
int dlcmSizeRational(const Rational* vector, const int length, const int base)
{
   assert(vector != 0);
   assert(length >= 0);
   assert(base >= 0);

   mpz_t lcm;

   mpz_init_set_ui(lcm, 1);

   for( int i = 0; i < length; i++ )
      mpz_lcm(lcm, lcm, mpq_denref(vector[i].getMpqRef()));

   int size = (int)mpz_sizeinbase(lcm, base);

   mpz_clear(lcm);

   return size;
}



/// Size of largest denominator in rational vector.
int dmaxSizeRational(const Rational* vector, const int length, const int base)
{
   assert(vector != 0);
   assert(length >= 0);
   assert(base >= 0);

   size_t dmax = 0;

   for( int i = 0; i < length; i++ )
   {
      size_t dsize = mpz_sizeinbase(mpq_denref(vector[i].getMpqRef()), base);
      if( dsize > dmax )
         dmax = dsize;
   }

   return (int)dmax;
}



#else



/// Defines the "Pimpl"-class Private
class Rational::Private
{

public:

   /// value
   long double privatevalue;

   /// default constructor
   Private()
   {
      privatevalue = 0;
   }

   /// copy constructor
   Private(const Private& p)
   {
      *this = p;
   }

   /// constructor from long double
   Private(const long double& r)
   {
      privatevalue = r;
   }

   /// constructor from double
   Private(const double& r)
   {
      privatevalue = r;
   }

   /// constructor from int
   Private(const int& i)
   {
      privatevalue = i;
   }

   /// assignment operator
   Private& operator=(const Private& p)
   {
      this->privatevalue = p.privatevalue;
      return *this;
   }

   /// assignment operator from long double
   Private& operator=(const long double& r)
   {
      this->privatevalue = r;
      return *this;
   }

   /// assignment operator from double
   Private& operator=(const double& r)
   {
      this->privatevalue = (long double)(r);
      return *this;
   }

   /// assignment operator from int
   Private& operator=(const int& i)
   {
      this->privatevalue = (long double)(i);
      return *this;
   }
};



/// default constructor
Rational::Rational()
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private();
}



/// copy constructor
Rational::Rational(const Rational& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(*(r.dpointer));
}



/// constructor from long double
Rational::Rational(const long double& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(r);
}



/// constructor from double
Rational::Rational(const double& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(r);
}



/// constructor from int
Rational::Rational(const int& i)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(i);
}



/// destructor
Rational::~Rational()
{
   spx_free(dpointer);
}



/// enables list memory
void Rational::enableListMem()
{
   // because list memory is not used when SOPLEX_WITH_GMP is not defined, there is nothing to do here
}



/// frees the unused rational elements in the memory list
/** this can be useful when you want to save memory or needed when working with a GMP memory manager like the one
 *  in EGlib that frees GMP memory before the destructor of the static memory list is called; in most cases this
 *  method is optional; note that this does not free the Rational elements that are currently in use
 */
void Rational::freeListMem()
{
   // because list memory is not used when SOPLEX_WITH_GMP is not defined, there is nothing to do here
}



/// disables list memory
void Rational::disableListMem()
{
   // because list memory is not used when SOPLEX_WITH_GMP is not defined, there is nothing to do here
}



/// assignment operator
Rational& Rational::operator=(const Rational &r)
{
   *dpointer = *(r.dpointer);
   return *this;
}



/// assignment operator from long double
Rational& Rational::operator=(const long double &r)
{
   *dpointer = r;
   return *this;
}



/// assignment operator from double
Rational& Rational::operator=(const double &r)
{
   *dpointer = r;
   return *this;
}



/// assignment operator from int
Rational& Rational::operator=(const int &i)
{
   *dpointer = i;
   return *this;
}



/// typecasts Rational to double (allows only explicit typecast)
Rational::operator double() const
{
   return (double)this->dpointer->privatevalue;
}



/// typecasts Rational to long double (allows only explicit typecast)
Rational::operator long double() const

{
   return this->dpointer->privatevalue;
}



/// addition operator
Rational Rational::operator+(const Rational& r) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue += r.dpointer->privatevalue;
   return retval;
}



/// addition assignment operator
Rational& Rational::operator+=(const Rational& r)
{
   this->dpointer->privatevalue += r.dpointer->privatevalue;
   return *this;
}



/// addition operator for doubles
Rational Rational::operator+(const double& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue += d;
   return retval;
}



/// addition assignment operator for doubles
Rational& Rational::operator+=(const double& d)
{
   this->dpointer->privatevalue += d;
   return *this;
}



/// addition operator for ints
Rational Rational::operator+(const int& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue += d;
   return retval;
}



/// addition assignment operator for ints
Rational& Rational::operator+=(const int& d)
{
   this->dpointer->privatevalue += d;
   return *this;
}



/// subtraction operator
Rational Rational::operator-(const Rational& r) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue -= r.dpointer->privatevalue;
   return retval;
}



/// subtraction assignment operator
Rational& Rational::operator-=(const Rational& r)
{
   this->dpointer->privatevalue -= r.dpointer->privatevalue;
   return *this;
}



/// subtraction operator for doubles
Rational Rational::operator-(const double& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue -= d;
   return retval;
}



/// subtraction assignment operator for doubles
Rational& Rational::operator-=(const double& d)
{
   this->dpointer->privatevalue -= d;
   return *this;
}



/// subtraction operator for ints
Rational Rational::operator-(const int& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue -= d;
   return retval;
}



/// subtraction assignment operator for ints
Rational& Rational::operator-=(const int& d)
{
   this->dpointer->privatevalue -= d;
   return *this;
}



/// multiplication operator
Rational Rational::operator*(const Rational& r) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue *= r.dpointer->privatevalue;
   return retval;
}



/// multiplication assignment operator
Rational& Rational::operator*=(const Rational& r)
{
   this->dpointer->privatevalue *= r.dpointer->privatevalue;
   return *this;
}



/// multiplication operator for doubles
Rational Rational::operator*(const double& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue *= d;
   return retval;
}



/// multiplication assignment operator for doubles
Rational& Rational::operator*=(const double& d)
{
   this->dpointer->privatevalue *= d;
   return *this;
}



/// multiplication operator for ints
Rational Rational::operator*(const int& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue *= d;
   return retval;
}



/// multiplication assignment operator for ints
Rational& Rational::operator*=(const int& d)
{
   this->dpointer->privatevalue *= d;
   return *this;
}



/// division operator
Rational Rational::operator/(const Rational& r) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue /= r.dpointer->privatevalue;
   return retval;
}



/// division assignment operator
Rational& Rational::operator/=(const Rational& r)
{
   this->dpointer->privatevalue /= r.dpointer->privatevalue;
   return *this;
}



/// division operator for doubles
Rational Rational::operator/(const double& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue /= d;
   return retval;
}



/// division assignment operator for doubles
Rational& Rational::operator/=(const double& d)
{
   this->dpointer->privatevalue /= d;
   return *this;
}



/// division operator for ints
Rational Rational::operator/(const int& r) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue /= r;
   return retval;
}



/// division assignment operator for ints
Rational& Rational::operator/=(const int& r)
{
   this->dpointer->privatevalue /= r;
   return *this;
}



/// add product of two rationals
Rational& Rational::addProduct(const Rational& r, const Rational& s)
{
   this->dpointer->privatevalue += r.dpointer->privatevalue * s.dpointer->privatevalue;
   return *this;
}



/// subtract product of two rationals
Rational& Rational::subProduct(const Rational& r, const Rational& s)
{
   this->dpointer->privatevalue -= r.dpointer->privatevalue * s.dpointer->privatevalue;
   return *this;
}



/// add quotient of two rationals, r divided by s
Rational& Rational::addQuotient(const Rational& r, const Rational& s)
{
   this->dpointer->privatevalue += r.dpointer->privatevalue / s.dpointer->privatevalue;
   return *this;
}



/// subtract quotient of two rationals, r divided by s
Rational& Rational::subQuotient(const Rational& r, const Rational& s)
{
   this->dpointer->privatevalue -= r.dpointer->privatevalue / s.dpointer->privatevalue;
   return *this;
}



/// inversion
Rational& Rational::invert()
{
   this->dpointer->privatevalue = 1.0 / this->dpointer->privatevalue;
   return *this;
}



/// round up to next power of two
Rational& Rational::powRound()
{
   ///@todo implement
   return *this;
}



/// checks if d is the closest possible double
bool Rational::isNextTo(const double& d)
{
   return (this->dpointer->privatevalue == d);
}



/// checks if d is exactly equal to the Rational and if not, if it is one of the two adjacent doubles
bool Rational::isAdjacentTo(const double& d) const
{
   return (double(this->dpointer->privatevalue) == d);
}



/// Size in specified base (bit size for base 2)
int Rational::sizeInBase(const int base) const
{
   ///@todo this is only correct for base 2
   return precision();
}



/// returns precision of Rational implementation, i.e., number of bits used to store Rational numbers (INT_MAX if exact)
int Rational::precision()
{
   return sizeof(long double);
}



/// read Rational from string
bool Rational::readString(const char* s)
{
   return (sscanf(s, "%Lf", &this->dpointer->privatevalue) == 1 );
}




/// convert rational number to string
std::string rationalToString(const Rational& r, const int precision)
{
   std::stringstream sstream;
   sstream << std::setprecision(precision <= 0 ? 16 : precision) << r;
   return sstream.str();
}



/// read Rational from string
bool readStringRational(const char* s, Rational& value)
{
   return (sscanf(s, "%Lf", &value.dpointer->privatevalue) == 1);
}



/// print Rational
std::ostream& operator<<(std::ostream& os, const Rational& r)
{
   os << r.dpointer->privatevalue;
   return os;
}



/// comparison operator returning a positive value if r > s, zero if r = s, and a negative value if r < s
int compareRational(const Rational& r, const Rational& s)
{
   if( r.dpointer->privatevalue > s.dpointer->privatevalue)
      return 1;
   else if( r.dpointer->privatevalue < s.dpointer->privatevalue)
      return -1;
   else
      return 0;
}



/// equality operator
bool operator==(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue == s.dpointer->privatevalue);
}



/// inequality operator
bool operator!=(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue != s.dpointer->privatevalue);
}



/// less than operator
bool operator<(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue < s.dpointer->privatevalue);
}



/// less than or equal to operator
bool operator<=(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue <= s.dpointer->privatevalue);
}



/// greater than operator
bool operator>(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue > s.dpointer->privatevalue);
}



/// greater than or equal to operator
bool operator>=(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue >= s.dpointer->privatevalue);
}



/// equality operator for Rational and double
bool operator==(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue > s - DEFAULT_EPS_ZERO)
      && (r.dpointer->privatevalue < s + DEFAULT_EPS_ZERO);
}



/// inequality operator for Rational and double
bool operator!=(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue <= s - DEFAULT_EPS_ZERO)
      || (r.dpointer->privatevalue >= s + DEFAULT_EPS_ZERO);
}



/// less than operator for Rational and double
bool operator<(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue < s);
}



/// less than or equal to operator for Rational and double
bool operator<=(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue <= s);
}



/// greater than operator for Rational and double
bool operator>(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue > s);
}



/// greater than or equal to operator for Rational and double
bool operator>=(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue >= s);
}



/// equality operator for double and Rational
bool operator==(const double& r, const Rational& s)
{
   return (s.dpointer->privatevalue > r - DEFAULT_EPS_ZERO)
      && (s.dpointer->privatevalue < r + DEFAULT_EPS_ZERO);
}



/// inequality operator double and Rational
bool operator!=(const double& r, const Rational& s)
{
   return (s.dpointer->privatevalue <= r - DEFAULT_EPS_ZERO)
      || (s.dpointer->privatevalue >= r + DEFAULT_EPS_ZERO);
}



/// less than operator double and Rational
bool operator<(const double& r, const Rational& s)
{
   return (r < s.dpointer->privatevalue);
}



/// less than or equal to operator double and Rational
bool operator<=(const double& r, const Rational& s)
{
   return (r <= s.dpointer->privatevalue);
}



/// greater than operator double and Rational
bool operator>(const double& r, const Rational& s)
{
   return (r > s.dpointer->privatevalue);
}



/// greater than or equal to operator double and Rational
bool operator>=(const double& r, const Rational& s)
{
   return (r >= s.dpointer->privatevalue);
}




/// equality operator for Rational and long double
bool operator==(const Rational& r, const long double& s)
{
   return (r.dpointer->privatevalue > s - DEFAULT_EPS_ZERO)
      && (r.dpointer->privatevalue < s + DEFAULT_EPS_ZERO);
}



/// inequality operator for Rational and long double
bool operator!=(const Rational& r, const long double& s)
{
   return (r.dpointer->privatevalue <= s - DEFAULT_EPS_ZERO)
      || (r.dpointer->privatevalue >= s + DEFAULT_EPS_ZERO);
}



/// less than operator for Rational and long double
bool operator<(const Rational& r, const long double& s)
{
   return (r.dpointer->privatevalue < s);
}



/// less than or equal to operator for Rational and long double
bool operator<=(const Rational& r, const long double& s)
{
   return (r.dpointer->privatevalue <= s);
}



/// greater than operator for Rational and long double
bool operator>(const Rational& r, const long double& s)
{
   return (r.dpointer->privatevalue > s);
}



/// greater than or equal to operator for Rational and long double
bool operator>=(const Rational& r, const long double& s)
{
   return (r.dpointer->privatevalue >= s);
}



/// equality operator for long double and Rational
bool operator==(const long double& r, const Rational& s)
{
   return (s.dpointer->privatevalue > r - DEFAULT_EPS_ZERO)
      && (s.dpointer->privatevalue < r + DEFAULT_EPS_ZERO);
}



/// inequality operator long double and Rational
bool operator!=(const long double& r, const Rational& s)
{
   return (s.dpointer->privatevalue <= r - DEFAULT_EPS_ZERO)
      || (s.dpointer->privatevalue >= r + DEFAULT_EPS_ZERO);
}



/// less than operator long double and Rational
bool operator<(const long double& r, const Rational& s)
{
   return (r < s.dpointer->privatevalue);
}



/// less than or equal to operator long double and Rational
bool operator<=(const long double& r, const Rational& s)
{
   return (r <= s.dpointer->privatevalue);
}



/// greater than operator long double and Rational
bool operator>(const long double& r, const Rational& s)
{
   return (r > s.dpointer->privatevalue);
}



/// greater than or equal to operator long double and Rational
bool operator>=(const long double& r, const Rational& s)
{
   return (r >= s.dpointer->privatevalue);
}



/// equality operator for Rational and int
bool operator==(const Rational& r, const int& s)
{
   return r.dpointer->privatevalue == s;
}



/// inequality operator for Rational and int
bool operator!=(const Rational& r, const int& s)
{
   return r.dpointer->privatevalue != s;
}



/// less than operator for Rational and int
bool operator<(const Rational& r, const int& s)
{
   return r.dpointer->privatevalue < s;
}



/// less than or equal to operator for Rational and int
bool operator<=(const Rational& r, const int& s)
{
   return r.dpointer->privatevalue <= s;
}



/// greater than operator for Rational and int
bool operator>(const Rational& r, const int& s)
{
   return r.dpointer->privatevalue > s;
}



/// greater than or equal to operator for Rational and int
bool operator>=(const Rational& r, const int& s)
{
   return r.dpointer->privatevalue >= s;
}



/// equality operator for int and Rational
bool operator==(const int& r, const Rational& s)
{
   return r == s.dpointer->privatevalue;
}



/// inequality operator for int and Rational
bool operator!=(const int& r, const Rational& s)
{
   return r != s.dpointer->privatevalue;
}



/// less than operator for int and Rational
bool operator<(const int& r, const Rational& s)
{
   return r < s.dpointer->privatevalue;
}



/// less than or equal to operator for int and Rational
bool operator<=(const int& r, const Rational& s)
{
   return r <= s.dpointer->privatevalue;
}



/// greater than operator for int and Rational
bool operator>(const int& r, const Rational& s)
{
   return r > s.dpointer->privatevalue;
}



/// greater than or equal to operator for int and Rational
bool operator>=(const int& r, const Rational& s)
{
   return r >= s.dpointer->privatevalue;
}



/// addition operator for double and Rational
Rational operator+(const double& d, const Rational& r)
{
   Rational retval(d);
   retval.dpointer->privatevalue += r.dpointer->privatevalue;
   return retval;
}



/// subtraction operator for double and Rational
Rational operator-(const double& d, const Rational& r)
{
   Rational retval(d);
   retval.dpointer->privatevalue -= r.dpointer->privatevalue;
   return retval;
}



/// multiplication operator for double and Rational
Rational operator*(const double& d, const Rational& r)
{
   Rational retval(d);
   retval.dpointer->privatevalue *= r.dpointer->privatevalue;
   return retval;
}



/// division operator for double and Rational
Rational operator/(const double& d, const Rational& r)
{
   Rational retval(d);
   retval.dpointer->privatevalue /= r.dpointer->privatevalue;
   return retval;
}



/// addition operator for int and Rational
Rational operator+(const int& d, const Rational& r)
{
   return (r + d);
}



/// subtraction operator for int and Rational
Rational operator-(const int& d, const Rational& r)
{
   return -(r - d);
}



/// multiplication operator for int and Rational
Rational operator*(const int& d, const Rational& r)
{
   return (r * d);
}



/// division operator for int and Rational
Rational operator/(const int& d, const Rational& r)
{
   return 1.0 / r;
}



/// Absolute.
Rational spxAbs(const Rational& r)
{
   Rational res = r;

   if( res.dpointer->privatevalue < 0 )
      res.dpointer->privatevalue *= -1;

   return res;
}



/// Sign function; returns 1 if r > 0, 0 if r = 0, and -1 if r < 0.
int sign(const Rational& r)
{
      return (r.dpointer->privatevalue > 0) - (r.dpointer->privatevalue < 0);
}



/// Negation.
Rational operator-(const Rational& r)
{
   Rational res = r;
   res.dpointer->privatevalue *= -1;
   return res;
}



/// Total size of rational vector.
int totalSizeRational(const Rational* vector, const int length, const int base)
{
   assert(vector != 0);
   assert(length >= 0);
   assert(base >= 0);

   int size = 0;

   for( int i = 0; i < length; i++ )
      size += vector[i].sizeInBase(base);

   return size;
}



/// Size of least common multiple of denominators in rational vector.
int dlcmSizeRational(const Rational* vector, const int length, const int base)
{
   assert(vector != 0);
   assert(length >= 0);
   assert(base >= 0);

   return 0;
}



/// Size of largest denominator in rational vector.
int dmaxSizeRational(const Rational* vector, const int length, const int base)
{
   assert(vector != 0);
   assert(length >= 0);
   assert(base >= 0);

   return 0;
}



#endif // SOPLEX_WITH_GMP
} // namespace soplex
#endif
