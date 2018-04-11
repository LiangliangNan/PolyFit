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

/**@file  dvector.h
 * @brief Dynamic vectors.
 */
#ifndef _DVECTOR_H_
#define _DVECTOR_H_

#include "spxdefines.h"
#include "basevectors.h"
#include "vector.h" // for compatibility

namespace soplex
{
typedef DVectorBase< Real > DVector;
typedef DVectorBase< Real > DVectorReal;
typedef DVectorBase< Rational > DVectorRational;
} // namespace soplex
#endif // _DVECTOR_H_
