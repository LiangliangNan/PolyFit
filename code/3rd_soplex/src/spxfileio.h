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

/**@file  spxfileio.h
 * @brief declaration of types for file output
 *
 * This is to make the use of compressed input files transparent
 * when programming.
 *
 * @todo maybe rename this file (it is unrelated to spxfileio.cpp)
 */
#ifndef _SPXFILEIO_H_
#define _SPXFILEIO_H_

#include <iostream>
#include <fstream>

/*-----------------------------------------------------------------------------
 * compressed file support
 *-----------------------------------------------------------------------------
 */
#ifdef SOPLEX_WITH_ZLIB
#include "gzstream.h"
#endif // WITH_GSZSTREAM

namespace soplex
{
#ifdef SOPLEX_WITH_ZLIB
   typedef gzstream::igzstream spxifstream;
#else
   typedef std::ifstream spxifstream;
#endif // SOPLEX_WITH_ZLIB

} // namespace soplex
#endif // _SPXFILEIO_H_
