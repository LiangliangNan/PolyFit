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


/**@file  ssvector.h
 * @brief Semi sparse vector.
 */
#ifndef _SSVECTOR_H_
#define _SSVECTOR_H_

#include "spxdefines.h"
#include "basevectors.h"

namespace soplex
{
typedef SSVectorBase<Real> SSVector;
typedef SSVectorBase<Real> SSVectorReal;
typedef SSVectorBase<Rational> SSVectorRational;
} // namespace soplex
#endif // _SSVECTOR_H_
