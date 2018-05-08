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

/**@file  vector.h
 * @brief Dense vector for linear algebra.
 */
#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <math.h>

#include "spxdefines.h"
#include "basevectors.h"

namespace soplex
{
typedef VectorBase< Real > Vector;
typedef VectorBase< Real > VectorReal;
typedef VectorBase< Rational > VectorRational;
} // namespace soplex
#endif // _VECTOR_H_
