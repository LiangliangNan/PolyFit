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

/**@file   bitencode.c
 * @ingroup OTHER_CFILES
 * @brief  packing single and dual bit values
 * @author Thorsten Koch
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/bitencode.h"


/** encode a single bit vector into packed format */
void SCIPencodeSingleBit(
   const int*            inp,                /**< unpacked input vector */
   SCIP_SINGLEPACKET*    out,                /**< buffer to store the packed vector */
   int                   count               /**< number of elements */
   )
{
   static const SCIP_SINGLEPACKET mask[SCIP_SINGLEPACKETSIZE][2] = {   /* if the packet size changes, the mask has to be updated */
      {0x00000000, 0x00000001},
      {0x00000000, 0x00000002},
      {0x00000000, 0x00000004},
      {0x00000000, 0x00000008},
      {0x00000000, 0x00000010},
      {0x00000000, 0x00000020},
      {0x00000000, 0x00000040},
      {0x00000000, 0x00000080},
      {0x00000000, 0x00000100},
      {0x00000000, 0x00000200},
      {0x00000000, 0x00000400},
      {0x00000000, 0x00000800},
      {0x00000000, 0x00001000},
      {0x00000000, 0x00002000},
      {0x00000000, 0x00004000},
      {0x00000000, 0x00008000},
      {0x00000000, 0x00010000},
      {0x00000000, 0x00020000},
      {0x00000000, 0x00040000},
      {0x00000000, 0x00080000},
      {0x00000000, 0x00100000},
      {0x00000000, 0x00200000},
      {0x00000000, 0x00400000},
      {0x00000000, 0x00800000},
      {0x00000000, 0x01000000},
      {0x00000000, 0x02000000},
      {0x00000000, 0x04000000},
      {0x00000000, 0x08000000},
      {0x00000000, 0x10000000},
      {0x00000000, 0x20000000},
      {0x00000000, 0x40000000},
      {0x00000000, 0x80000000}
   };
   int rest;
   int nfull;
   const int packetsize = (int) SCIP_SINGLEPACKETSIZE;

   assert(inp != NULL || count == 0);
   assert(out != NULL || count == 0);
   assert(count >= 0);
   assert(packetsize == 32);

   rest = count % packetsize;
   nfull = count - rest;

   for( int i = 0; i < nfull; i += packetsize )
   {
      assert(inp != NULL);
      assert(out != NULL);

#ifndef NDEBUG
      {
         for( int j = 0; j < packetsize; ++j )
            assert(0 <= inp[j] && inp[j] <= 1);
      }
#endif
      *out++ =
         mask[0][inp[0]] | mask[1][inp[1]] | mask[2][inp[2]] | mask[3][inp[3]]
         | mask[4][inp[4]] | mask[5][inp[5]] | mask[6][inp[6]] | mask[7][inp[7]]
         | mask[8][inp[8]] | mask[9][inp[9]] | mask[10][inp[10]] | mask[11][inp[11]]
         | mask[12][inp[12]] | mask[13][inp[13]] | mask[14][inp[14]] | mask[15][inp[15]]
         | mask[16][inp[16]] | mask[17][inp[17]] | mask[18][inp[18]] | mask[19][inp[19]]
         | mask[20][inp[20]] | mask[21][inp[21]] | mask[22][inp[22]] | mask[23][inp[23]]
         | mask[24][inp[24]] | mask[25][inp[25]] | mask[26][inp[26]] | mask[27][inp[27]]
         | mask[28][inp[28]] | mask[29][inp[29]] | mask[30][inp[30]] | mask[31][inp[31]];
      inp += packetsize;
   }

   if( rest > 0 )
   {
      SCIP_SINGLEPACKET m = (SCIP_SINGLEPACKET) 0u;

      assert(inp != NULL);
      assert(out != NULL);
      assert(rest <= (int) SCIP_SINGLEPACKETSIZE);

      for( int i = 0; i < rest; i++ )
         m |= mask[i][inp[i]];
      *out = m;
   }
}

/** decode a packed single bit vector into unpacked format */
void SCIPdecodeSingleBit(
   const SCIP_SINGLEPACKET* inp,             /**< packed input vector */
   int*                  out,                /**< buffer to store unpacked vector */
   int                   count               /**< number of elements */
   )
{
   SCIP_SINGLEPACKET m;
   int rest;
   int nfull;
   int i;

   assert(inp != NULL || count == 0);
   assert(out != NULL || count == 0);
   assert(count >= 0);
   assert(SCIP_SINGLEPACKETSIZE == 32); /*lint !e506*/

   rest = count % (int)SCIP_SINGLEPACKETSIZE;
   nfull = count - rest;

   for( i = 0; i < nfull; i += (int)SCIP_SINGLEPACKETSIZE )
   {
      assert(inp != NULL);
      assert(out != NULL);

      m = *inp++;

      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      m >>= 1;
      *out++ = (int)m & 1;
      assert(m >> 1 == 0);
   }

   if( rest > 0 )
   {
      assert(inp != NULL);
      assert(out != NULL);

      m = *inp;
      for( i = 0; i < rest; i++ )
      {
         *out++ = (int)m & 1;
         m >>= 1;
      }
   }
}

/** encode a dual bit vector into packed format */
void SCIPencodeDualBit(
   const int*            inp,                /**< unpacked input vector */
   SCIP_DUALPACKET*      out,                /**< buffer to store the packed vector */
   int                   count               /**< number of elements */
   )
{
   static const SCIP_DUALPACKET mask[SCIP_DUALPACKETSIZE][4] = {   /* if the packet size changes, the mask has to be updated */
      {0x00000000, 0x00000001, 0x00000002, 0x00000003},
      {0x00000000, 0x00000004, 0x00000008, 0x0000000C},
      {0x00000000, 0x00000010, 0x00000020, 0x00000030},
      {0x00000000, 0x00000040, 0x00000080, 0x000000C0},
      {0x00000000, 0x00000100, 0x00000200, 0x00000300},
      {0x00000000, 0x00000400, 0x00000800, 0x00000C00},
      {0x00000000, 0x00001000, 0x00002000, 0x00003000},
      {0x00000000, 0x00004000, 0x00008000, 0x0000C000},
      {0x00000000, 0x00010000, 0x00020000, 0x00030000},
      {0x00000000, 0x00040000, 0x00080000, 0x000C0000},
      {0x00000000, 0x00100000, 0x00200000, 0x00300000},
      {0x00000000, 0x00400000, 0x00800000, 0x00C00000},
      {0x00000000, 0x01000000, 0x02000000, 0x03000000},
      {0x00000000, 0x04000000, 0x08000000, 0x0C000000},
      {0x00000000, 0x10000000, 0x20000000, 0x30000000},
      {0x00000000, 0x40000000, 0x80000000, 0xC0000000}
   };
   int rest;
   int nfull;
   const int dualpacketsize = (int) SCIP_DUALPACKETSIZE;

   assert(inp != NULL || count == 0);
   assert(out != NULL || count == 0);
   assert(count >= 0);
   assert(dualpacketsize == 16);

   rest = count % dualpacketsize;
   nfull = count - rest;

   for( int i = 0; i < nfull; i += dualpacketsize, inp += dualpacketsize ) /*lint !e679*/
   {
      assert(inp != NULL);
      assert(out != NULL);

#ifndef NDEBUG
      {
         for( int j = 0; j < dualpacketsize; ++j )
            assert(0 <= inp[j] && inp[j] <= 3);
      }
#endif
      *out++ =
         mask[0][inp[0]] | mask[1][inp[1]] | mask[2][inp[2]] | mask[3][inp[3]]
         | mask[4][inp[4]] | mask[5][inp[5]] | mask[6][inp[6]]
         | mask[7][inp[7]] | mask[8][inp[8]] | mask[9][inp[9]] 
         | mask[10][inp[10]] | mask[11][inp[11]] | mask[12][inp[12]] 
         | mask[13][inp[13]] | mask[14][inp[14]] | mask[15][inp[15]];
   }

   if( rest > 0 )
   {
      SCIP_DUALPACKET m = (SCIP_DUALPACKET) 0u;

      assert(inp != NULL);
      assert(out != NULL);
      assert(rest <= (int) SCIP_DUALPACKETSIZE);

      for( int i = 0; i < rest; i++ )
         m |= mask[i][inp[i]];
      *out = m;
   }
}

/** decode a packed dual bit vector into unpacked format */
void SCIPdecodeDualBit(
   const SCIP_DUALPACKET* inp,               /**< packed input vector */
   int*                  out,                /**< buffer to store unpacked vector */
   int                   count               /**< number of elements */
   )
{
   SCIP_DUALPACKET m;
   int rest;
   int nfull;
   int i;

   assert(inp != NULL || count == 0);
   assert(out != NULL || count == 0);
   assert(count >= 0);
   assert(SCIP_DUALPACKETSIZE == 16); /*lint !e506*/

   rest = count % (int)SCIP_DUALPACKETSIZE;
   nfull = count - rest;

   for( i = 0; i < nfull; i += (int)SCIP_DUALPACKETSIZE )
   {
      assert(inp != NULL);
      assert(out != NULL);

      m = *inp++;

      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      m >>= 2;
      *out++ = (int)m & 3;
      assert(m >> 2 == 0);
   }

   if( rest > 0 )
   {
      assert(inp != NULL);
      assert(out != NULL);

      m = *inp;
      for( i = 0; i < rest; i++ )
      {
         *out++ = (int)m & 3;
         m >>= 2;
      }
   }
}

