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

/**@file  spxlp.h
 * @brief Saving LPs in a form suitable for SoPlex.
 */

#ifndef _SPXLP_H_
#define _SPXLP_H_

#include "spxdefines.h"
#include "spxlpbase.h"
#include "vector.h" // for compatibility
#include "svector.h" // for compatibility
#include "svset.h" // for compatibility
#include "lprowset.h" // for compatibility
#include "lpcolset.h" // for compatibility
#include "lprow.h" // for compatibility
#include "lpcol.h" // for compatibility

namespace soplex
{
typedef SPxLPBase< Real > SPxLP;
typedef SPxLPBase< Real > SPxLPReal;
typedef SPxLPBase< Rational > SPxLPRational;
} // namespace soplex
#endif // _SPXLP_H_
