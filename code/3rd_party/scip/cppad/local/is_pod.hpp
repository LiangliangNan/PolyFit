# ifndef CPPAD_LOCAL_IS_POD_HPP
# define CPPAD_LOCAL_IS_POD_HPP
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-17 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
// make sure size_t is defined because autotools version of
// is_pod_specialize_98 uses it
# include <cstddef>

/*!
\file is_pod.hpp
File that defines is_pod<Type>(void)
*/
namespace CppAD { namespace local { // BEGIN_CPPAD_LOCAL_NAMESPACE
/*!
Is this type plain old data; i.e., its constructor need not be called.

The default definition is false. This include file defines it as true
for all the fundamental types except for void and nullptr_t.
*/
template <class T> inline bool is_pod(void) { return false; }
// The following command suppresses doxygen processing for the code below
/// \cond
// C++98 Fundamental types
template <> inline bool is_pod<size_t>(void) { return true; }

# if CPPAD_USE_CPLUSPLUS_2011
// C++11 Fundamental types

# endif

/// \endcond
} } // END_CPPAD_LOCAL_NAMESPACE

# endif
