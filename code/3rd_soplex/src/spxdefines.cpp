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

/**@file  spxdefines.cpp
 * @brief Debugging, floating point type and parameter definitions.
 */
#include "spxdefines.h"
#include "assert.h"
#include "spxout.h"

namespace soplex
{

THREADLOCAL const Real infinity                 = DEFAULT_INFINITY;

THREADLOCAL Real Param::s_epsilon               = DEFAULT_EPS_ZERO;
THREADLOCAL Real Param::s_epsilon_factorization = DEFAULT_EPS_FACTOR;
THREADLOCAL Real Param::s_epsilon_update        = DEFAULT_EPS_UPDATE;
THREADLOCAL Real Param::s_epsilon_pivot         = DEFAULT_EPS_PIVOT;

bool msginconsistent(const char* name, const char* file, int line)
{
   assert(name != 0);
   assert(file != 0);
   assert(line >= 0);

   MSG_ERROR( std::cerr << file << "(" << line << ") "
   << "Inconsistency detected in " << name << std::endl; )

   return 0;
}

Real Param::epsilon()
{
   return s_epsilon;
}

void Param::setEpsilon(Real eps)
{
   s_epsilon = eps;
}

Real Param::epsilonFactorization()
{
   return s_epsilon_factorization;
}

void Param::setEpsilonFactorization(Real eps)
{
   s_epsilon_factorization = eps;
}

Real Param::epsilonUpdate()
{
   return s_epsilon_update;
}

void Param::setEpsilonUpdate(Real eps)
{
   s_epsilon_update = eps;
}

Real Param::epsilonPivot()
{
   return s_epsilon_pivot;
}

void Param::setEpsilonPivot(Real eps)
{
   s_epsilon_pivot = eps;
}

} // namespace soplex
