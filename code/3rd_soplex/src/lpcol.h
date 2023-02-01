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

/**@file  lpcol.h
 * @brief LP column.
 */
#ifndef _LPCOL_H_
#define _LPCOL_H_

#include <assert.h>

#include "spxdefines.h"
#include "lpcolbase.h"

namespace soplex
{
typedef LPColBase< Real > LPCol;
typedef LPColBase< Real > LPColReal;
typedef LPColBase< Rational > LPColRational;
} // namespace soplex
#endif // _LPCOL_H_
