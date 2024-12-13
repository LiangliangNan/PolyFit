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

/**@file   scip_mem.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for memory management
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/mem.h"
#include "scip/pub_message.h"
#include "scip/scip_mem.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"

/** returns block memory to use at the current time
 *
 *  @return the block memory to use at the current time.
 */
BMS_BLKMEM* SCIPblkmem(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);

   return scip->mem->probmem;
}

/** returns buffer memory for short living temporary objects
 *
 *  @return the buffer memory for short living temporary objects
 */
BMS_BUFMEM* SCIPbuffer(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   return scip->mem->buffer;
}

/** returns clean buffer memory for short living temporary objects initialized to all zero
 *
 *  @return the buffer memory for short living temporary objects initialized to all zero
 */
BMS_BUFMEM* SCIPcleanbuffer(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   return scip->mem->cleanbuffer;
}

/** returns the total number of bytes used in block and buffer memory
 *
 *  @return the total number of bytes used in block and buffer memory.
 */
SCIP_Longint SCIPgetMemUsed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return SCIPmemGetUsed(scip->mem);
}

/** returns the total number of bytes in block and buffer memory
 *
 *  @return the total number of bytes in block and buffer memory.
 */
SCIP_Longint SCIPgetMemTotal(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return SCIPmemGetTotal(scip->mem);
}

/** returns the estimated number of bytes used by external software, e.g., the LP solver
 *
 *  @return the estimated number of bytes used by external software, e.g., the LP solver.
 */
SCIP_Longint SCIPgetMemExternEstim(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return SCIPstatGetMemExternEstim(scip->stat);
}

/** calculate memory size for dynamically allocated arrays
 *
 *  @return the memory size for dynamically allocated arrays.
 */
int SCIPcalcMemGrowSize(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);

   return SCIPsetCalcMemGrowSize(scip->set, num);
}

/** extends a dynamically allocated block memory array to be able to store at least the given number of elements;
 *  use SCIPensureBlockMemoryArray() define to call this method!
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPensureBlockMemoryArray_call(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                arrayptr,           /**< pointer to dynamically sized array */
   size_t                elemsize,           /**< size in bytes of each element in array */
   int*                  arraysize,          /**< pointer to current array size */
   int                   minsize             /**< required minimal array size */
   )
{
   assert(scip != NULL);
   assert(arrayptr != NULL);
   assert(elemsize > 0);
   assert(arraysize != NULL);

   if( minsize > *arraysize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(scip->set, minsize);
      SCIP_ALLOC( BMSreallocBlockMemorySize(SCIPblkmem(scip), arrayptr, *arraysize * elemsize, newsize * elemsize) );
      *arraysize = newsize;
   }

   return SCIP_OKAY;
}

/** prints output about used memory */
void SCIPprintMemoryDiagnostic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);

   BMSdisplayMemory();

   SCIPmessagePrintInfo(scip->messagehdlr, "\nParameter Block Memory (%p):\n", (void*)scip->mem->setmem);
   BMSdisplayBlockMemory(scip->mem->setmem);

   SCIPmessagePrintInfo(scip->messagehdlr, "\nSolution Block Memory (%p):\n", (void*)scip->mem->probmem);
   BMSdisplayBlockMemory(scip->mem->probmem);

   SCIPmessagePrintInfo(scip->messagehdlr, "\nMemory Buffers:\n");
   BMSprintBufferMemory(SCIPbuffer(scip));

   SCIPmessagePrintInfo(scip->messagehdlr, "\nClean Memory Buffers:\n");
   BMSprintBufferMemory(SCIPcleanbuffer(scip));
}
