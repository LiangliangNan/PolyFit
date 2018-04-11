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

/**@file  lprowset.h
 * @brief Set of LP columns.
 */
#ifndef _LPROWSET_H_
#define _LPROWSET_H_

#include "spxdefines.h"
#include "lprowsetbase.h"

namespace soplex
{
typedef LPRowSetBase< Real > LPRowSet;
typedef LPRowSetBase< Real > LPRowSetReal;
typedef LPRowSetBase< Rational > LPRowSetRational;
} // namespace soplex
#endif // _LPROWSET_H_
