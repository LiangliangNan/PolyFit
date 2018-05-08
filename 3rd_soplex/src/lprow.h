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

/**@file  lprow.h
 * @brief (In)equality for LPs.
 */
#ifndef _LPROW_H_
#define _LPROW_H_

#include "spxdefines.h"
#include "lprowbase.h"

namespace soplex
{
typedef LPRowBase< Real > LPRow;
typedef LPRowBase< Real > LPRowReal;
typedef LPRowBase< Rational > LPRowRational;
} // namespace soplex
#endif // _LPROW_H_
