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

/**@file  dsvector.h
 * @brief Dynamic sparse vectors.
 */
#ifndef _DSVECTOR_H_
#define _DSVECTOR_H_

#include "spxdefines.h"
#include "basevectors.h"
#include "svector.h" // for compatibility

namespace soplex
{
typedef DSVectorBase< Real > DSVector;
typedef DSVectorBase< Real > DSVectorReal;
typedef DSVectorBase< Rational > DSVectorRational;
} // namespace soplex
#endif // _DSVECTOR_H_
