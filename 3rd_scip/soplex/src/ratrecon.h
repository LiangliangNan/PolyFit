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

/**@file  ratrecon.h
 * @brief Rational reconstruction of solution vector
 */


#ifndef _RATRECON_H_
#define _RATRECON_H_

#ifndef SOPLEX_LEGACY

#include "spxdefines.h"
#include "rational.h"
#include "sol.h"
#include "basevectors.h"
#include "didxset.h"

namespace soplex
{
   /** reconstruct a rational vector */
   bool reconstructVector(VectorRational& input, const Rational& denomBoundSquared, const DIdxSet* indexSet = 0);

   /** reconstruct a rational solution */
   bool reconstructSol(SolRational& solution);
} // namespace soplex
#endif
#endif // _RATRECON_H_
