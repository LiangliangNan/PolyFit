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
extern
void SCIPencodeSingleBit(
   const int*            inp,                /**< unpacked input vector */
   SCIP_SINGLEPACKET*    out,                /**< buffer to store the packed vector */
   int                   count               /**< number of elements */
   );

/** decode a packed single bit vector into unpacked format */
extern
void SCIPdecodeSingleBit(
   const SCIP_SINGLEPACKET* inp,             /**< packed input vector */
   int*                  out,                /**< buffer to store unpacked vector */
   int                   count               /**< number of elements */
   );

/** encode a dual bit vector into packed format */
extern
void SCIPencodeDualBit(
   const int*            inp,                /**< unpacked input vector */
   SCIP_DUALPACKET*      out,                /**< buffer to store the packed vector */
   int                   count               /**< number of elements */
   );

/** decode a packed dual bit vector into unpacked format */
extern
void SCIPdecodeDualBit(
   const SCIP_DUALPACKET* inp,               /**< packed input vector */
   int*                  out,                /**< buffer to store unpacked vector */
   int                   count               /**< number of elements */
   );

#ifdef __cplusplus
}
#endif

#endif
