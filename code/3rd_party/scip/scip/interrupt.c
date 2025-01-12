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

/**@file   interrupt.c
 * @ingroup OTHER_CFILES
 * @brief  methods and datastructures for catching the user CTRL-C interrupt
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <sys/types.h>
#include <stdlib.h>
#include <signal.h>

#include "scip/def.h"
#include "scip/pub_message.h"
#include "blockmemshell/memory.h"
#include "scip/interrupt.h"


static volatile
int                      ninterrupts = 0;    /**< static variable counting the number of CTRL-C interrupts */
static volatile
int                      nterms = 0;         /**< static variable counting the number of times that the process received a SIGTERM signal */


#ifdef SCIP_NO_SIGACTION
typedef void (*SigHdlr)(int);

/** CTRL-C interrupt data */
struct SCIP_Interrupt
{
   SigHdlr               oldsighdlr;         /**< old CTRL-C interrupt handler */
   int                   nuses;              /**< number of times, the interrupt is captured */
};

#else

/** CTRL-C interrupt data */
struct SCIP_Interrupt
{
   struct sigaction      oldsigaction;       /**< old CTRL-C interrupt handler */
   int                   nuses;              /**< number of times, the interrupt is captured */
};
#endif

/** interrupt handler for CTRL-C interrupts */
static
void interruptHandler(
   int                   signum              /**< interrupt signal number */
   )
{
   SCIP_UNUSED(signum);

   ninterrupts++;
   if( ninterrupts >= 5 )
   {
      printf("pressed CTRL-C %d times. forcing termination.\n", ninterrupts);
      exit(1);
   }
   else
   {
      printf("pressed CTRL-C %d times (5 times for forcing termination)\n", ninterrupts);
   }
}

/** creates a CTRL-C interrupt data */
SCIP_RETCODE SCIPinterruptCreate(
   SCIP_INTERRUPT**      interrupt           /**< pointer to store the CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);

   SCIP_ALLOC( BMSallocMemory(interrupt) );
   (*interrupt)->nuses = 0;

   return SCIP_OKAY;
}

/** frees a CTRL-C interrupt data */
void SCIPinterruptFree(
   SCIP_INTERRUPT**      interrupt           /**< pointer to the CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);

   BMSfreeMemory(interrupt);
}

/** captures the CTRL-C interrupt to call the SCIP's own interrupt handler */
void SCIPinterruptCapture(
   SCIP_INTERRUPT*       interrupt           /**< CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);
   assert(interrupt->nuses >= 0);

   if( interrupt->nuses == 0 )
   {
#ifdef SCIP_NO_SIGACTION
      interrupt->oldsighdlr = signal(SIGINT, interruptHandler);
#else
      struct sigaction newaction;

      /* initialize new signal action */
      newaction.sa_handler = interruptHandler;
      newaction.sa_flags = 0;
      (void)sigemptyset(&newaction.sa_mask);

      /* set new signal action, and remember old one */
      (void)sigaction(SIGINT, &newaction, &interrupt->oldsigaction);
#endif

      ninterrupts = 0;
      nterms = 0;
   }
   interrupt->nuses++;
}

/** releases the CTRL-C interrupt and restores the old interrupt handler */
void SCIPinterruptRelease(
   SCIP_INTERRUPT*       interrupt           /**< CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);
   assert(interrupt->nuses >= 1);

   interrupt->nuses--;
   if( interrupt->nuses == 0 )
   {
#ifdef SCIP_NO_SIGACTION
      (void)signal(SIGINT, interrupt->oldsighdlr);
#else
      (void)sigaction(SIGINT, &interrupt->oldsigaction, NULL);
#endif
   }
}

/** returns whether the user interrupted by pressing CTRL-C */
SCIP_Bool SCIPinterrupted(
   void
   )
{
   return (ninterrupts > 0);
}

/** returns whether a process termination signal was received */
SCIP_Bool SCIPterminated(
   void
   )
{
   return (nterms > 0);
}

/** sends a termination signal to all SCIP processes so that they try to terminate as soon as possible
 *
 *  @note For terminating a specific SCIP process use SCIPinterruptSolve().
 */
void SCIPtryTerminate(
   void
   )
{
   nterms++;
}

/** resets the number of interrupts to 0 */
void SCIPresetInterrupted(
   void
   )
{
   ninterrupts = 0;
   nterms = 0;
}

