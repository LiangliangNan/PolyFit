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

/**@file   bitencode.h
 * @brief  packing single and dual bit values
 * @author Thorsten Koch
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BITENCODE_H__
#define __SCIP_BITENCODE_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned int SCIP_SINGLEPACKET;                /**< storing single bits in packed format */
#define SCIP_SINGLEPACKETSIZE (sizeof(SCIP_SINGLEPACKET)*8) /**< each entry needs one bit of information */
typedef unsigned int SCIP_DUALPACKET;                  /**< storing bit pairs in packed format */
#define SCIP_DUALPACKETSIZE   (sizeof(SCIP_DUALPACKET)*4)   /**< each entry needs two bits of information */


/** encode a single bit vector into packed format */
void SCIPencodeSingleBit(
   const int*            inp,                /**< unpacked input vector */
   SCIP_SINGLEPACKET*    out,                /**< buffer to store the packed vector */
   int                   count               /**< number of elements */
   );

/** decode a packed single bit vector into unpacked format */
void SCIPdecodeSingleBit(
   const SCIP_SINGLEPACKET* inp,             /**< packed input vector */
   int*                  out,                /**< buffer to store unpacked vector */
   int                   count               /**< number of elements */
   );

/** encode a dual bit vector into packed format */
void SCIPencodeDualBit(
   const int*            inp,                /**< unpacked input vector */
   SCIP_DUALPACKET*      out,                /**< buffer to store the packed vector */
   int                   count               /**< number of elements */
   );

/** decode a packed dual bit vector into unpacked format */
void SCIPdecodeDualBit(
   const SCIP_DUALPACKET* inp,               /**< packed input vector */
   int*                  out,                /**< buffer to store unpacked vector */
   int                   count               /**< number of elements */
   );

#ifdef __cplusplus
}
#endif

#endif
