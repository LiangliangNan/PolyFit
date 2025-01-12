/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_cuts.h
 * @ingroup PUBLICCOREAPI
 * @brief  struct definitions for cuts
 * @author Leona Gottwald
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
