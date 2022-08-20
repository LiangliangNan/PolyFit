/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_cuts.h
 * @ingroup PUBLICCOREAPI
 * @brief  struct definitions for cuts
 * @author Robert Lion Gottwald
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CUTS_H__
#define __SCIP_STRUCT_CUTS_H__

#include "scip/def.h"
#include "scip/dbldblarith.h"

struct SCIP_AggrRow
{
   SCIP_Real*            vals;               /**< non-zero coefficients of the cut row */
   int*                  inds;               /**< problem indices of variables with a non-zero coefficient in the cut row */
   int*                  rowsinds;           /**< lpposition of rows that have been added to the cutrow */
   int*                  slacksign;          /**< slacksign of rows that have been added to the cutrow */
   SCIP_Real*            rowweights;         /**< weights of rows that have been added to the cutrow */
   QUAD_MEMBER(SCIP_Real rhs);               /**< right hand side of the cut row */
   int                   nnz;                /**< number of non-zeros in the cut row */
   int                   nrows;              /**< number of rows that have been added to the cutrow */
   int                   rowssize;           /**< size of the row and slacksign array */
   int                   rank;               /**< rank of the cut row */
   SCIP_Bool             local;              /**< is the cut row only valid locally? */
};

#endif
