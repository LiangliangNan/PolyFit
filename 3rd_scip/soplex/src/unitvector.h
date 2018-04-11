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

#ifndef _UNITVECTOR_H_
#define _UNITVECTOR_H_

#include <assert.h>
#include "spxdefines.h"
#include "basevectors.h"

namespace soplex
{
typedef UnitVectorBase< Real > UnitVector;
typedef UnitVectorBase< Real > UnitVectorReal;
typedef UnitVectorBase< Rational > UnitVectorRational;
} // namespace soplex
#endif // _UNITVECTOR_H_
