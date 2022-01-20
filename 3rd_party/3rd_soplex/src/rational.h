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

/**@file  rational.h
 * @brief Wrapper for GMP types.
 */
#ifndef _RATIONAL_H_
#define _RATIONAL_H_

#ifndef SOPLEX_LEGACY
#include <math.h>
#include <assert.h>
#include <string.h>
#include <iostream>

#include "spxdefines.h"
#include "idlist.h"

#ifdef SOPLEX_WITH_GMP
#include "gmp.h"
#endif


namespace soplex
{
/**@brief   Wrapper for GMP type mpq_class.
 * @ingroup Algebra
 *
 * We wrap mpq_class so that we can replace it by a double type if GMP is not available.
 */

/// If compiled with GMP support, Rational is defined as mpq_class.
   class Rational
   {
   private:
      class Private;
      Private* dpointer;

#ifdef SOPLEX_WITH_GMP
      THREADLOCAL static IdList< Private > unusedPrivateList;
      THREADLOCAL static bool useListMem;


      /// special constructor only for initializing static rational variables; this is necessary since we need a
      /// constructor for Rational::{ZERO, POSONE, NEGONE} that does not use these numbers
      Rational(const int& i, const bool& dummy);

      //**@name Static variables for special rational values */
      //@{

      static const Rational ZERO;
      static const Rational POSONE;
      static const Rational NEGONE;

      //@}
#endif

   public:

      //**@name Construction and destruction */
      //@{

      /// default constructor
      Rational();

      /// copy constructor
      Rational(const Rational& r);

      /// constructor from long double
      Rational(const long double& r);

      /// constructor from double
      Rational(const double& r);

      ///constructor from int
      Rational(const int& i);

#ifdef SOPLEX_WITH_GMP
      /// constructor from mpq_t
      Rational(const mpq_t& q);
#endif

      /// destructor
      ~Rational();

      /// enables list memory
      static void enableListMem();

      /// frees the unused rational elements in the memory list
      /** this can be useful when you want to save memory or needed when working with a GMP memory manager like the one
       *  in EGlib that frees GMP memory before the destructor of the static memory list is called; in most cases this
       *  method is optional; note that this does not free the Rational elements that are currently in use
       */
      static void freeListMem();

      /// disables list memory
      static void disableListMem();

      /// assignment operator
      Rational& operator=(const Rational&);

      /// assignment operator from long double
      Rational& operator=(const long double &r);

      /// assignment operator from double
      Rational& operator=(const double &r);

      /// assignment operator from int
      Rational& operator=(const int &i);

#ifdef SOPLEX_WITH_GMP
      /// assignment operator from mpq_t
      Rational& operator=(const mpq_t &q);
#endif

      //@}


      //**@name Typecasts */
      //@{

      operator double() const;
      operator long double() const;

#ifdef SOPLEX_WITH_GMP
      /// provides read-only access to underlying mpq_t
      const mpq_t* getMpqPtr() const;

      /// provides read-only access to underlying mpq_t
      const mpq_t& getMpqRef() const;

      /// provides write access to underlying mpq_t; use with care
      mpq_t* getMpqPtr_w() const;

      /// provides write access to underlying mpq_t; use with care
      mpq_t& getMpqRef_w() const;
#endif
      //@}


      //**@name Arithmetic operators */
      //@{

      /// addition operator
      Rational operator+(const Rational& r) const;

      /// addition assignment operator
      Rational& operator+=(const Rational& r);

      /// addition operator for doubles
      Rational operator+(const double& r) const;

      /// addition assignment operator  for doubles
      Rational& operator+=(const double& r);

      /// addition operator for ints
      Rational operator+(const int& r) const;

      /// addition assignment operator  for ints
      Rational& operator+=(const int& r);

      /// subtraction operator
      Rational operator-(const Rational& r) const;

      /// subtraction assignment operator
      Rational& operator-=(const Rational& r);

      /// subtraction operator for doubles
      Rational operator-(const double& r) const;

      /// subtraction assignment operator for doubles
      Rational& operator-=(const double& r);

      /// subtraction operator for ints
      Rational operator-(const int& r) const;

      /// subtraction assignment operator for ints
      Rational& operator-=(const int& r);

      /// multiplication operator
      Rational operator*(const Rational& r) const;

      /// multiplication assignment operator operator
      Rational& operator*=(const Rational& r);

      /// multiplication operator for doubles
      Rational operator*(const double& r) const;

      /// multiplication assignment operator for doubles
      Rational& operator*=(const double& r);

      /// multiplication operator for ints
      Rational operator*(const int& r) const;

      /// multiplication assignment operator for ints
      Rational& operator*=(const int& r);

      /// division operator
      Rational operator/(const Rational& r) const;

      /// division assignment operator
      Rational& operator/=(const Rational& r);

      /// division operator for doubles
      Rational operator/(const double& r) const;

      /// division assignment operator for doubles
      Rational& operator/=(const double& r);

      /// division operator for ints
      Rational operator/(const int& r) const;

      /// division assignment operator for ints
      Rational& operator/=(const int& r);

      /// add product of two rationals
      Rational& addProduct(const Rational& r, const Rational& s);

      /// subtract product of two rationals
      Rational& subProduct(const Rational& r, const Rational& s);

      /// add quotient of two rationals, r divided by s
      Rational& addQuotient(const Rational& r, const Rational& s);

      /// subtract quotient of two rationals, r divided by s
      Rational& subQuotient(const Rational& r, const Rational& s);

      /// inversion
      Rational& invert();

      /// round up to next power of two
      Rational& powRound();

      //@}


      //**@name Methods for checking exactness of doubles  */
      //@{

      /// checks if \p d is the closest number that can be represented by double
      bool isNextTo(const double& d);

      /// checks if \p d is exactly equal to the Rational and if not, if it is one of the two adjacent doubles
      bool isAdjacentTo(const double& d) const;

      //@}


      //**@name Methods for querying size */
      //@{

      /// Size in specified base (bit size for base 2)
      int sizeInBase(const int base = 2) const;

      //@}


      //**@name Static methods  */
      //@{

      /// returns precision of Rational implementation, i.e., number of bits used to store Rational numbers (INT_MAX if exact)
      static int precision();

      //@}


      //**@name Conversion from and to String */
      //@{

      /// read Rational from string
      bool readString(const char* s);

      friend std::string rationalToString(const Rational& r, const int precision);
      friend bool readStringRational(const char* s, Rational& value);
      friend std::ostream& operator<<(std::ostream& os, const Rational& q);

      //@}


      //**@name Friends */
      //@{

      friend int compareRational(const Rational& r, const Rational& s);
      friend bool operator!=(const Rational& r, const Rational& s);
      friend bool operator==(const Rational& r, const Rational& s);
      friend bool operator<(const Rational& r, const Rational& s);
      friend bool operator<=(const Rational& r, const Rational& s);
      friend bool operator>(const Rational& r, const Rational& s);
      friend bool operator>=(const Rational& r, const Rational& s);

      friend bool operator!=(const Rational& r, const double& s);
      friend bool operator==(const Rational& r, const double& s);
      friend bool operator<(const Rational& r, const double& s);
      friend bool operator<=(const Rational& r, const double& s);
      friend bool operator>(const Rational& r, const double& s);
      friend bool operator>=(const Rational& r, const double& s);

      friend bool operator!=(const double& r, const Rational& s);
      friend bool operator==(const double& r, const Rational& s);
      friend bool operator<(const double& r, const Rational& s);
      friend bool operator<=(const double& r, const Rational& s);
      friend bool operator>(const double& r, const Rational& s);
      friend bool operator>=(const double& r, const Rational& s);

      friend bool operator!=(const Rational& r, const long double& s);
      friend bool operator==(const Rational& r, const long double& s);
      friend bool operator<(const Rational& r, const long double& s);
      friend bool operator<=(const Rational& r, const long double& s);
      friend bool operator>(const Rational& r, const long double& s);
      friend bool operator>=(const Rational& r, const long double& s);

      friend bool operator!=(const long double& r, const Rational& s);
      friend bool operator==(const long double& r, const Rational& s);
      friend bool operator<(const long double& r, const Rational& s);
      friend bool operator<=(const long double& r, const Rational& s);
      friend bool operator>(const long double& r, const Rational& s);
      friend bool operator>=(const long double& r, const Rational& s);

      friend Rational operator+(const double& d, const Rational& r);
      friend Rational operator-(const double& d, const Rational& r);
      friend Rational operator*(const double& d, const Rational& r);
      friend Rational operator/(const double& d, const Rational& r);

      friend bool operator!=(const Rational& r, const int& s);
      friend bool operator==(const Rational& r, const int& s);
      friend bool operator<(const Rational& r, const int& s);
      friend bool operator<=(const Rational& r, const int& s);
      friend bool operator>(const Rational& r, const int& s);
      friend bool operator>=(const Rational& r, const int& s);

      friend bool operator!=(const int& r, const Rational& s);
      friend bool operator==(const int& r, const Rational& s);
      friend bool operator<(const int& r, const Rational& s);
      friend bool operator<=(const int& r, const Rational& s);
      friend bool operator>(const int& r, const Rational& s);
      friend bool operator>=(const int& r, const Rational& s);

      friend Rational operator+(const int& d, const Rational& r);
      friend Rational operator-(const int& d, const Rational& r);
      friend Rational operator*(const int& d, const Rational& r);
      friend Rational operator/(const int& d, const Rational& r);

      friend Rational spxAbs(const Rational& r);
      friend int sign(const Rational& r);
      friend Rational operator-(const Rational& q);

      //@}
   };


   //**@name Parsing and printing */
   //@{

   /// convert rational number to string
   std::string rationalToString(const Rational& r, const int precision = 32);

   /// read Rational from string
   bool readStringRational(const char* s, Rational& value);

   /// print Rational
   std::ostream& operator<<(std::ostream& os, const Rational& r);

   //@}


   //**@name Relational operators */
   //@{

   /// comparison operator returning a positive value if r > s, zero if r = s, and a negative value if r < s
   int compareRational(const Rational& r, const Rational& s);

   /// equality operator
   bool operator==(const Rational& r, const Rational& s);

   /// inequality operator
   bool operator!=(const Rational& r, const Rational& s);

   /// less than operator
   bool operator<(const Rational& r, const Rational& s);

   /// less than or equal to operator
   bool operator<=(const Rational& r, const Rational& s);

   /// greater than operator
   bool operator>(const Rational& r, const Rational& s);

   /// greater than or equal to operator
   bool operator>=(const Rational& r, const Rational& s);

   /// equality operator for Rational and double
   bool operator==(const Rational& r, const double& s);

   /// inequality operator for Rational and double
   bool operator!=(const Rational& r, const double& s);

   /// less than operator for Rational and double
   bool operator<(const Rational& r, const double& s);

   /// less than or equal to operator for Rational and double
   bool operator<=(const Rational& r, const double& s);

   /// greater than operator for Rational and double
   bool operator>(const Rational& r, const double& s);

   /// greater than or equal to operator for Rational and double
   bool operator>=(const Rational& r, const double& s);

   /// equality operator for double and Rational
   bool operator==(const double& r, const Rational& s);

   /// inequality operator for double and Rational
   bool operator!=(const double& r, const Rational& s);

   /// less than operator for double and Rational
   bool operator<(const double& r, const Rational& s);

   /// less than or equal to operator for double and Rational
   bool operator<=(const double& r, const Rational& s);

   /// greater than operator for double and Rational
   bool operator>(const double& r, const Rational& s);

   /// greater than or equal to operator for double and Rational
   bool operator>=(const double& r, const Rational& s);

   /// equality operator for Rational and long double
   bool operator==(const Rational& r, const long double& s);

   /// inequality operator for Rational and long double
   bool operator!=(const Rational& r, const long double& s);

   /// less than operator for Rational and long double
   bool operator<(const Rational& r, const long double& s);

   /// less than or equal to operator for Rational and long double
   bool operator<=(const Rational& r, const long double& s);

   /// greater than operator for Rational and long double
   bool operator>(const Rational& r, const long double& s);

   /// greater than or equal to operator for Rational and long double
   bool operator>=(const Rational& r, const long double& s);

   /// equality operator for long double and Rational
   bool operator==(const long double& r, const Rational& s);

   /// inequality operator for long double and Rational
   bool operator!=(const long double& r, const Rational& s);

   /// less than operator for long double and Rational
   bool operator<(const long double& r, const Rational& s);

   /// less than or equal to operator for long double and Rational
   bool operator<=(const long double& r, const Rational& s);

   /// greater than operator for long double and Rational
   bool operator>(const long double& r, const Rational& s);

   /// greater than or equal to operator for long double and Rational
   bool operator>=(const long double& r, const Rational& s);

   /// equality operator for Rational and int
   bool operator==(const Rational& r, const int& s);

   /// inequality operator for Rational and int
   bool operator!=(const Rational& r, const int& s);

   /// less than operator for Rational and int
   bool operator<(const Rational& r, const int& s);

   /// less than or equal to operator for Rational and int
   bool operator<=(const Rational& r, const int& s);

   /// greater than operator for Rational and int
   bool operator>(const Rational& r, const int& s);

   /// greater than or equal to operator for Rational and int
   bool operator>=(const Rational& r, const int& s);

   /// equality operator for int and Rational
   bool operator==(const int& r, const Rational& s);

   /// inequality operator for int and Rational
   bool operator!=(const int& r, const Rational& s);

   /// less than operator for int and Rational
   bool operator<(const int& r, const Rational& s);

   /// less than or equal to operator for int and Rational
   bool operator<=(const int& r, const Rational& s);

   /// greater than operator for int and Rational
   bool operator>(const int& r, const Rational& s);

   /// greater than or equal to operator for int and Rational
   bool operator>=(const int& r, const Rational& s);

   //@}


   //**@name Non-member arithmetic operators and functions */
   //@{

   /// addition operator for double and Rational
   Rational operator+(const double& d, const Rational& r);

   /// addition operator for double and Rational
   Rational operator+(const double& d, const Rational& r);

   /// addition operator for double and Rational
   Rational operator+(const double& d, const Rational& r);

   /// addition operator for double and Rational
   Rational operator+(const double& d, const Rational& r);

   /// addition operator for int and Rational
   Rational operator+(const int& d, const Rational& r);

   /// addition operator for int and Rational
   Rational operator+(const int& d, const Rational& r);

   /// addition operator for int and Rational
   Rational operator+(const int& d, const Rational& r);

   /// addition operator for int and Rational
   Rational operator+(const int& d, const Rational& r);

   /// absolute function
   Rational spxAbs(const Rational& r);

   /// Sign function; returns 1 if r > 0, 0 if r = 0, and -1 if r < 0.
   int sign(const Rational& r);

   /// Negation.
   Rational operator-(const Rational& r);

   /// Total size of rational vector.
   int totalSizeRational(const Rational* vector, const int length, const int base = 2);

   /// Size of least common multiple of denominators in rational vector.
   int dlcmSizeRational(const Rational* vector, const int length, const int base = 2);

   /// Size of largest denominator in rational vector.
   int dmaxSizeRational(const Rational* vector, const int length, const int base = 2);

   //@}

} // namespace soplex
#else
namespace soplex
{
   typedef Real Rational;
}
#endif
#endif // _RATIONAL_H_
