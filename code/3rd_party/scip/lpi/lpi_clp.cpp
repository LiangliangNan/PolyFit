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

/**@file   lpi_clp.cpp
 * @ingroup LPIS
 * @brief  LP interface for Clp
 * @author Stefan Heinz
 * @author Marc Pfetsch
 * @author John Forrest
 *
 *
 * Notes on this interface:
 *
 * - Currently, Clp (Version 1.16) supports two ways of adding rows/columns from arrays: One uses a
 *   length array that for each row/column specifies the number of nonzeros to be added. The second
 *   uses the @p beg array that gives the starting index for each row/column. We use the latter
 *   variant. Since for LPI there should be no gaps in the corresponding arrays, i.e., every entry in
 *   @p val and @a ind gives a nonzero entry, one can switch between the two formats. With the current
 *   Clp implementation both formats involve an overhead:
 *
 *    - For the @p beg variant, Clp gets the end of the array from the last position in @p beg
 *      (i.e., the entry one after the last row/column) and we have to copy and extend @p beg for this
 *      purpose. In the matrix implementation a length information is then again computed.
 *
 *    - For the @p length variant, Clp computes the number of elements from this length variant and
 *      there exists no matrix implementation that uses the length information, i.e., it is recomputed
 *      again.
 *
 *    Concluding: the implementation of Clp/CoinPackeMatrix could be improved. The functions
 *    affected by this are SCIPlpiLoadColLP(), SCIPlpiAddCols(), SCIPlpiAddRows()
 *
 * - In former versions Clp used an "auxiliary model" that allows to save time when the model is
 *   scaled. This is discarded from version higher than 1.8.2.
 *
 * - Clp needs the setting of several special flags to make sure that necessary data (e.g., rays etc.)
 *   are available. If the FASTMIP option in SCIP is true, some settings that are supposed to be faster
 *   are used. We tried to use the best settings, while still working correctly. These settings probably
 *   have to be adapted to future Clp versions. Maybe more possibilities will appear.
 *
 * - At several places this interface corrects the return value of some Clp functions, e.g.,
 *   isProvenPrimalInfeasible(). Currently (Version 1.16) no change in the Clp functions will be made,
 *   but this might change in the future.
 *
 * @todo Check about using fastDual().
 */
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <ClpSimplex.hpp>
#include <ClpPrimalColumnSteepest.hpp>
#include <ClpDualRowSteepest.hpp>
#include <CoinIndexedVector.hpp>
#include <ClpConfig.h>
#ifndef CLP_VERSION
#include <config_clp.h>
#define CLP_VERSION VERSION
#endif

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#include "lpi/lpi.h"
#include "scip/bitencode.h"
#include "scip/pub_message.h"

/* do defines for windows directly her to make the lpi more independent*/
#if defined(_WIN32) || defined(_WIN64)
#define snprintf _snprintf
#endif

/* for debugging: alternatingly write files "debug_[p|d]_[0|1].mps" after each run - use with care! */
#ifdef LPI_CLP_DEBUG_WRITE_FILES
static int fileNr = 0;
#endif

/* bound for accepting primal or dual sum of infeasibilities */
#define SUMINFEASBOUND   1.0e-3

/** LP interface for Clp */
struct SCIP_LPi
{
   ClpSimplex*           clp;                        /**< Clp simiplex solver class */
   int*                  cstat;                      /**< array for storing column basis status */
   int*                  rstat;                      /**< array for storing row basis status */
   int                   cstatsize;                  /**< size of cstat array */
   int                   rstatsize;                  /**< size of rstat array */
   bool                  startscratch;               /**< start from scratch? */
   SCIP_PRICING          pricing;                    /**< SCIP pricing setting  */
   bool                  validFactorization;         /**< whether we have a valid factorization in clp */
   SCIP_Bool             solved;                     /**< was the current LP solved? */
   bool                  setFactorizationFrequency;  /**< store whether the factorization frequency is set */
   SCIP_Bool             fastmip;                    /**< are fast mip settings turned on */
   int                   lastalgorithm;              /**< type of last algorithm call (0 = none, 1 = primal, -1 = dual, 2 = barrier) */
};






/** Definitions for storing basis status  (copied from lpi_spx.cpp) */
typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE

/** LPi state stores basis information */
struct SCIP_LPiState
{
   int                   ncols;              /**< number of LP columns */
   int                   nrows;              /**< number of LP rows */
   COLPACKET*            packcstat;          /**< column basis status in compressed form */
   ROWPACKET*            packrstat;          /**< row basis status in compressed form */
};




/*
 * dynamic memory arrays
 */

/** resizes cstat array to have at least num entries */
static
SCIP_RETCODE ensureCstatMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->cstatsize )
   {
      int newsize;

      newsize = MAX(2*lpi->cstatsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->cstat, newsize) );
      lpi->cstatsize = newsize;
   }
   assert(num <= lpi->cstatsize);

   return SCIP_OKAY;
}

/** resizes rstat array to have at least num entries */
static
SCIP_RETCODE ensureRstatMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->rstatsize )
   {
      int newsize;

      newsize = MAX(2*lpi->rstatsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rstat, newsize) );
      lpi->rstatsize = newsize;
   }
   assert(num <= lpi->rstatsize);

   return SCIP_OKAY;
}




/*
 * LPi state methods
 */

/** returns the number of packets needed to store column packet information */
static
int colpacketNum(
   int                   ncols               /**< number of columns to store */
   )
{
   return (ncols+(int)COLS_PER_PACKET-1)/(int)COLS_PER_PACKET;
}

/** returns the number of packets needed to store row packet information */
static
int rowpacketNum(
   int                   nrows               /**< number of rows to store */
   )
{
   return (nrows+(int)ROWS_PER_PACKET-1)/(int)ROWS_PER_PACKET;
}

/** store row and column basis status in a packed LPi state object */
static
void lpistatePack(
   SCIP_LPISTATE*        lpistate,           /**< pointer to LPi state data */
   const int*            cstat,              /**< basis status of columns in unpacked format */
   const int*            rstat               /**< basis status of rows in unpacked format */
   )
{
   assert(lpistate != 0);
   assert(lpistate->packcstat != 0);
   assert(lpistate->packrstat != 0);

   SCIPencodeDualBit(cstat, lpistate->packcstat, lpistate->ncols);
   SCIPencodeDualBit(rstat, lpistate->packrstat, lpistate->nrows);
}

/** unpacks row and column basis status from a packed LPi state object */
static
void lpistateUnpack(
   const SCIP_LPISTATE*  lpistate,           /**< pointer to LPi state data */
   int*                  cstat,              /**< buffer for storing basis status of columns in unpacked format */
   int*                  rstat               /**< buffer for storing basis status of rows in unpacked format */
   )
{
   assert(lpistate != 0);
   assert(lpistate->packcstat != 0);
   assert(lpistate->packrstat != 0);

   SCIPdecodeDualBit(lpistate->packcstat, cstat, lpistate->ncols);
   SCIPdecodeDualBit(lpistate->packrstat, rstat, lpistate->nrows);
}

/** creates LPi state information object */
static
SCIP_RETCODE lpistateCreate(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   ncols,              /**< number of columns to store */
   int                   nrows               /**< number of rows to store */
   )
{
   assert(lpistate != 0);
   assert(blkmem != 0);
   assert(ncols >= 0);
   assert(nrows >= 0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, lpistate) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packcstat, colpacketNum(ncols)) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packrstat, rowpacketNum(nrows)) );

   return SCIP_OKAY;
}

/** frees LPi state information */
static
void lpistateFree(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state information (like basis information) */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(blkmem != 0);
   assert(lpistate != 0);
   assert(*lpistate != 0);

   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packcstat, colpacketNum((*lpistate)->ncols));
   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packrstat, rowpacketNum((*lpistate)->nrows));
   BMSfreeBlockMemory(blkmem, lpistate);
}





/*
 * local methods
 */

/** marks the current LP to be unsolved */
static
void invalidateSolution(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   lpi->solved = FALSE;
}

/** set factorization frequency */
static
void setFactorizationFrequency(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   /* set the factorization frequency only once */
   if ( lpi->setFactorizationFrequency )
      return;

   lpi->clp->defaultFactorizationFrequency();
   lpi->setFactorizationFrequency = true;
}

/** set fastmip parameters of Clp */
static
void setFastmipClpParameters(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   lpi->fastmip = TRUE;

   /* Perturbation:
    *  50  - switch on perturbation
    *  100 - auto perturb if takes too long (1.0e-6 largest nonzero)
    *  101 - we are perturbed
    *  102 - don't try perturbing again
    *  - default is 100
    *  - others are for playing
    *
    * for Clp 1.8 stable: 50 seems to be 10% faster than 100
    */
   lpi->clp->setPerturbation(50);

   /* Special options description from ClpModell.hpp:
    *         1 - Don't keep changing infeasibility weight
    *         2 - Keep nonLinearCost round solves
    *         4 - Force outgoing variables to exact bound (primal)
    *         8 - Safe to use dense initial factorization
    *        16 - Just use basic variables for operation if column generation
    *        32 - Create ray even in BAB
    *        64 - Treat problem as feasible until last minute (i.e. minimize infeasibilities)
    *       128 - Switch off all matrix sanity checks
    *       256 - No row copy
    *       512 - If not in values pass, solution guaranteed, skip as much as possible
    *      1024 - In branch and bound
    *      2048 - Don't bother to re-factorize if < 20 iterations
    *      4096 - Skip some optimality checks
    *      8192 - Do Primal when cleaning up primal
    *     16384 - In fast dual (so we can switch off things)
    *     32768 - called from Osi
    *     65536 - keep arrays around as much as possible (also use maximumR/C)
    *    131072 - transposeTimes is -1.0 and can skip basic and fixed
    *    262144 - extra copy of scaled matrix
    *    524288 - Clp fast dual
    *   1048576 - don't need to finish dual (can return 3)
    *   2097152 - ray even if >2 pivots AND if problem is "crunched"
    *   4194304 - don't scale integer variables
    *   8388608 - Idiot when not really sure about it
    *  16777216 - zero costs!
    * 0x1000000 - is Cbc (and in branch and bound)
    * 0x2000000 - is in a different branch and bound
    *
    *  Comments:
    *         2 - Nonlinear costs are used in primal for infeasibility weight.
    *         4 - In anti-degeneracy operations can move variables just off a bound.
    *         8 - Means dense nucleus in factorization - normally not safe in first factorization as
    *             singularity handling is not useful. Is switched on if going from dual to primal or vv.
    *        16 - Used for "real" column generation.
    *        32 - Currently unclear, does not lead to produce ray
    *        64 - Good idea, since in B&B most problems are feasible.
    *       128 - Assumes user will not create tiny or duplicate elements.
    *       256 - Normally Clp keeps a scaled row copy for speed. For very large problems you might want to turn it off.
    *       512 - Means nonbasic variables should be at bounds and basis will be reasonable.
    *      1024 - In branch and bound - makes some rays available?
    *      2048 - Unclear.
    *      4096 - Skip some optimality checks.
    *      8192 - If the primal has a perturbed problem and needs to clean up, it normally uses dual - but in some cases can be better to use primal.
    *     16384 - Used internally.
    *     32768 - Just switches off some messages, e.g., empty problem.
    *     65536 - Unclear.
    *    131072 - Used internally.
    *    262144 - Normally Clp has unscaled column copy of matrix - this makes an extra scaled copy.
    *    524288 - Used internally.
    *   1048576 - Only set by fastDual (used internally).
    *   2097152 - This was fixed in 03/2018 to make sure that a ray is produced.
    *   4194304 - Not needed for us.
    *   8388608 - Unclear.
    * 0x1000000 - main point: does allow use of disaster handler, but also other decisions in code
    * 0x2000000 - main point: does allow use of disaster handler, but also other decisions in code
    *
    * Cbc seems to use the following special options:
    * lpi->clp->setSpecialOptions(64|128|1024|2048|4096|32768|262144|0x01000000);
    * Sometimes 512+8192 and 8192 or 8 are used as well.
    */

   // default settings:        32|64|128|1024|32768|262144|2097152|0x2000000
   // Additional flags: 512, 2048, 4096 in order to speed up things.
   // lpi->clp->setSpecialOptions(32|64|128|512|1024|2048|4096|32768|262144|2097152|0x2000000);

   // set default options
   lpi->clp->setSpecialOptions(32|64|128|1024|32768|262144|2097152|0x2000000);

   // Do not change moreSpecialOptions().

   // let memory grow only (do not shrink) - [needs specialOptions & 65536 != 0]
   // does not seem to work
   // lpi->clp->setPersistenceFlag(1);
}

/** unset fastmip parameters of Clp (reset to default parameters) */
static
void unsetFastmipClpParameters(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   lpi->fastmip = FALSE;

   // reset to default value:
   lpi->clp->setPerturbation(100);

   // set default special options (see SCIPlpiCreate())
   lpi->clp->setSpecialOptions(32|64|128|1024|32768|262144|2097152|0x2000000);

   // set default more special options
   lpi->clp->setMoreSpecialOptions(8192);

   // turn off memory enlargement
   lpi->clp->setPersistenceFlag(0);
}


/*
 * LP Interface Methods
 */


/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   // Currently Clp has no function to get version, so we hard code it ...
   return "Clp " CLP_VERSION;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "COIN-OR Linear Programming Solver developed by J. Forrest et.al. (projects.coin-or.org/Clp)";
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   return (void*) lpi->clp;
}

/** pass integrality information to LP solver */
SCIP_RETCODE SCIPlpiSetIntegralityInformation(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   ncols,              /**< length of integrality array */
   int*                  intInfo             /**< integrality array (0: continuous, 1: integer). May be NULL iff ncols is 0.  */
   )
{
   assert(lpi != NULL);

   SCIPerrorMessage("SCIPlpiSetIntegralityInformation() has not been implemented yet.\n");

   return SCIP_LPERROR;
}

/** informs about availability of a primal simplex solving method */
SCIP_Bool SCIPlpiHasPrimalSolve(
   void
   )
{
   return TRUE;
}

/** informs about availability of a dual simplex solving method */
SCIP_Bool SCIPlpiHasDualSolve(
   void
   )
{
   return TRUE;
}

/** informs about availability of a barrier solving method */
SCIP_Bool SCIPlpiHasBarrierSolve(
   void
   )
{
   return TRUE;
}


/**@} */




/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(
   SCIP_LPI**            lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   assert(lpi != NULL);
   assert(name != NULL);

   SCIPdebugMessage("calling SCIPlpiCreate()\n");

   // create lpi object
   SCIP_ALLOC( BMSallocMemory(lpi) );
   (*lpi)->clp = new ClpSimplex();
   (*lpi)->cstat = 0;
   (*lpi)->rstat = 0;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->startscratch = true;
   (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*lpi)->validFactorization = false;
   (*lpi)->setFactorizationFrequency = false;
   (*lpi)->fastmip = FALSE;
   (*lpi)->lastalgorithm = 0;
   invalidateSolution(*lpi);

   // if you want to use saveModel()
   // (*lpi)->clp->setLengthNames(255);

   // set pricing routines

   // for primal:
   // 0 is exact devex,
   // 1 full steepest,
   // 2 is partial exact devex
   // 3 switches between 0 and 2 depending on factorization
   // 4 starts as partial dantzig/devex but then may switch between 0 and 2.
   // - currently (Clp 1.8stable) default is 3
   ClpPrimalColumnSteepest primalSteepest;
   (*lpi)->clp->setPrimalColumnPivotAlgorithm(primalSteepest);

   // for dual:
   // 0 is uninitialized,
   // 1 full,
   // 2 is partial uninitialized,
   // 3 starts as 2 but may switch to 1.
   // - currently (Clp 1.8stable) default is 3
   ClpDualRowSteepest dualSteepest;
   (*lpi)->clp->setDualRowPivotAlgorithm(dualSteepest);

   // set problem name
   (*lpi)->clp->setStrParam(ClpProbName, std::string(name) );

   // set objective sense: SCIP values are the same as the ones for Clp
   (*lpi)->clp->setOptimizationDirection(objsen);

   // turn off output by default
   (*lpi)->clp->setLogLevel(0);

   // turn on auto scaling by default
   (*lpi)->clp->scaling(3);

   // set default special options (similar to Cbc):
   //        32 - Create ray even in BAB
   //        64 - good idea to be fast
   //       128 - Assumes user will not create tiny or duplicate elements.
   //      1024 - In branch and bound.
   //     32768 - Just switches off some messages, e.g., empty problem.
   //    262144 - extra copy of scaled matrix
   //   2097152 - ray even if >2 pivots AND if problem is "crunched"
   // 0x2000000 - is in a different branch and bound
   (*lpi)->clp->setSpecialOptions(32|64|128|1024|32768|262144|2097152|0x2000000);

   /* More special options:
    *        1 bit - if presolve says infeasible in ClpSolve return
    *        2 bit - if presolved problem infeasible return
    *        4 bit - keep arrays like upper_ around
    *        8 bit - no free or superBasic variables
    *       16 bit - if checking replaceColumn accuracy before updating
    *       32 bit - say optimal if primal feasible!
    *       64 bit - give up easily in dual (and say infeasible)
    *      128 bit - no objective, 0-1 and in B&B
    *      256 bit - in primal from dual or vice versa
    *      512 bit - alternative use of solveType_
    *     1024 bit - don't do row copy of factorization
    *     2048 bit - perturb in complete fathoming
    *     4096 bit - try more for complete fathoming
    *     8192 bit - don't even think of using primal if user asks for dual (and vv)
    *    16384 bit - in initialSolve so be more flexible
    *    32768 bit - don't swap algorithms from dual if small infeasibility
    *    65536 bit - perturb in postsolve cleanup (even if < 10000 rows)
    *   131072 bit - initial stateDualColumn
    *   524288 bit - stop when primal feasible
    *  1048576 bit - don't perturb even if long time
    *  2097152 bit - no primal in fastDual2 if feasible
    *  4194304 bit - tolerances have been changed by code
    *  8388608 bit - tolerances are dynamic (at first)
    * 16777216 bit - if factorization kept can still declare optimal at once
    */

   // set default more special options:
   (*lpi)->clp->setMoreSpecialOptions(8192);

   // set default pricing
   SCIP_CALL( SCIPlpiSetIntpar(*lpi, SCIP_LPPAR_PRICING, (int)(*lpi)->pricing) );

   return SCIP_OKAY;
}


/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(*lpi != 0);
   assert((*lpi)->clp != 0);

   SCIPdebugMessage("calling SCIPlpiFree()\n");

   /* free LP */
   delete (*lpi)->clp;

   /* free memory */
   BMSfreeMemoryArrayNull(&(*lpi)->cstat);
   BMSfreeMemoryArrayNull(&(*lpi)->rstat);
   BMSfreeMemory(lpi);

   return SCIP_OKAY;
}

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies LP data with column matrix into LP solver */
SCIP_RETCODE SCIPlpiLoadColLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   const SCIP_Real*      obj,                /**< objective function values of columns */
   const SCIP_Real*      lb,                 /**< lower bounds of columns */
   const SCIP_Real*      ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array */
   const int*            ind,                /**< row indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
#ifndef NDEBUG
   {
      int j;
      for( j = 0; j < nnonz; j++ )
         assert( val[j] != 0 );
   }
#endif

   SCIPdebugMessage("calling SCIPlpiLoadColLP()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(beg != NULL);
   assert(ind != NULL);
   assert(val != NULL);

   assert( nnonz > beg[ncols-1] );

   invalidateSolution(lpi);

   ClpSimplex* clp = lpi->clp;

   // copy beg-array
   int* mybeg = NULL;
   SCIP_ALLOC( BMSallocMemoryArray(&mybeg, ncols + 1) );
   BMScopyMemoryArray(mybeg, beg, ncols);
   mybeg[ncols] = nnonz;   // add additional entry at end

   // load problem
   clp->loadProblem(ncols, nrows, mybeg, ind, val, lb, ub, obj, lhs, rhs);
   BMSfreeMemoryArray( &mybeg );

   // set objective sense
   clp->setOptimizationDirection(objsen);

   // copy column and rownames if necessary
   if ( colnames || rownames )
   {
      std::vector<std::string> columnNames(ncols);
      std::vector<std::string> rowNames(nrows);
      if (colnames)
      {
         for (int j = 0; j < ncols; ++j)
            columnNames[j].assign(colnames[j]);
      }
      if (rownames)
      {
         for (int i = 0; i < ncols; ++i)
            rowNames[i].assign(rownames[i]);
      }
      clp->copyNames(rowNames, columnNames);
   }

   return SCIP_OKAY;
}


/** adds columns to the LP */
SCIP_RETCODE SCIPlpiAddCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const SCIP_Real*      obj,                /**< objective function values of new columns */
   const SCIP_Real*      lb,                 /**< lower bounds of new columns */
   const SCIP_Real*      ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or 0 */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or 0 if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or 0 if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or 0 if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddCols()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);
   assert(nnonz >= 0);
   assert(ncols >= 0);

   invalidateSolution(lpi);

   // store number of columns for later
   int numCols = lpi->clp->getNumCols();

   // copy beg-array (if not 0)
   int* mybeg = NULL;
   SCIP_ALLOC( BMSallocMemoryArray(&mybeg, ncols + 1) );

   // if columns are not empty
   if ( nnonz != 0 )
   {
#ifndef NDEBUG
      {
         int j;
         for( j = 0; j < nnonz; j++ )
         {
            assert( val[j] != 0.0 );
            /* perform check that no new rows are added - this is forbidden */
            assert( 0 <= ind[j] /*&& ind[j] < lpi->nrows*/ );
         }
      }
#endif
      BMScopyMemoryArray(mybeg, beg, ncols);
      mybeg[ncols] = nnonz;   // add additional entry at end

      // add columns
      lpi->clp->addColumns(ncols, lb, ub, obj, mybeg, ind, val);
   }
   else
   {
      for (int j = 0; j <= ncols; ++j)
         mybeg[j] = 0;

      // add empty columns
      lpi->clp->addColumns(ncols, lb, ub, obj, mybeg, 0, 0);
   }
   BMSfreeMemoryArray(&mybeg);

   // copy columnnames if necessary
   if ( colnames )
   {
      std::vector<std::string> columnNames(ncols);
      for (int j = 0; j < ncols; ++j)
         columnNames[j].assign(colnames[j]);
      lpi->clp->copyColumnNames(columnNames, numCols, numCols + ncols);
   }

   return SCIP_OKAY;
}


/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelCols()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());

   invalidateSolution(lpi);

   // Current Clp version (1.8) can't delete a range of columns; we have to use deleteColumns (see SCIPlpiDelColset)
   int num = lastcol-firstcol+1;
   int* which = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &which, num) );;

   // fill array with interval
   for (int j = firstcol; j <= lastcol; ++j)
      which[j - firstcol] = j;

   lpi->clp->deleteColumns(num, which);
   BMSfreeMemoryArray( &which );

   return SCIP_OKAY;
}


/** deletes columns from SCIP_LPI; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelColset()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(dstat != NULL);

   invalidateSolution(lpi);

   // transform dstat information
   int ncols = lpi->clp->getNumCols();
   int* which = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &which, ncols) );
   int cnt = 0;
   for (int j = 0; j < ncols; ++j)
   {
      if ( dstat[j] == 1 )
         which[cnt++] = j;
   }
   lpi->clp->deleteColumns(cnt, which);
   BMSfreeMemoryArray(&which);

   // update dstat
   cnt = 0;
   for (int j = 0; j < ncols; ++j)
   {
      if ( dstat[j] == 1 )
      {
         dstat[j] = -1;
         ++cnt;
      }
      else
         dstat[j] = j - cnt;
   }

   return SCIP_OKAY;
}


/** adds rows to the LP */
SCIP_RETCODE SCIPlpiAddRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddRows()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);

   invalidateSolution(lpi);

   // store number of rows for later use
   int numRows = lpi->clp->getNumRows();

   int* mybeg = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &mybeg, nrows + 1) );

   if ( nnonz > 0 )
   {
#ifndef NDEBUG
      /* perform check that no new columns are added - this is likely to be a mistake */
      int ncols = lpi->clp->getNumCols();
      for (int j = 0; j < nnonz; ++j)
      {
         assert( val[j] != 0.0 );
         assert( 0 <= ind[j] && ind[j] < ncols );
      }
#endif

      // copy beg-array
      BMScopyMemoryArray( mybeg, beg, nrows);
      mybeg[nrows] = nnonz;   // add additional entry at end

      // add rows
      lpi->clp->addRows(nrows, lhs, rhs, mybeg, ind, val);
   }
   else
   {
      // add empty rows
      for (int i = 0; i <= nrows; ++i)
         mybeg[i] = 0;
      lpi->clp->addRows(nrows, lhs, rhs, mybeg, 0, 0);
   }
   BMSfreeMemoryArray( &mybeg );

   // copy rownames if necessary
   if ( rownames )
   {
      std::vector<std::string> rowNames(nrows);
      for (int j = 0; j < nrows; ++j)
         rowNames[j].assign(rownames[j]);
      lpi->clp->copyRowNames(rowNames, numRows, numRows + nrows);
   }

   return SCIP_OKAY;
}


/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelRows() (number: %d)\n", lastrow-firstrow+1);

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->clp->numberRows());

   invalidateSolution(lpi);

   // Current Clp version (1.8) can't delete a range of rows; we have to use deleteRows (see SCIPlpiDelRowset)
   int num = lastrow-firstrow+1;
   int* which = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &which, num) );

   // fill array with interval
   for (int i = firstrow; i <= lastrow; ++i)
      which[i - firstrow] = i;

   lpi->clp->deleteRows(num, which);

   BMSfreeMemoryArray( &which );

   return SCIP_OKAY;
}


/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelRowset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelRowset()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(dstat != 0);

   invalidateSolution(lpi);

   // transform dstat information
   int nrows = lpi->clp->getNumRows();
   int* which = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &which, nrows) );
   int cnt = 0;
   for (int i = 0; i < nrows; ++i)
   {
      if ( dstat[i] == 1 )
         which[cnt++] = i;
   }
   lpi->clp->deleteRows(cnt, which);
   BMSfreeMemoryArray( &which );

   // update dstat
   cnt = 0;
   for (int i = 0; i < nrows; ++i)
   {
      if ( dstat[i] == 1 )
      {
         dstat[i] = -1;
         ++cnt;
      }
      else
         dstat[i] = i - cnt;
   }

   return SCIP_OKAY;
}


/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiClear()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   invalidateSolution(lpi);
   lpi->lastalgorithm = 0;

   // We use the resize(0,0) to get rid of the model but keep all other settings
   lpi->clp->resize(0,0);

   return SCIP_OKAY;
}


/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiChgBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices or NULL if ncols is zero */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds or NULL if ncols is zero */
   const SCIP_Real*      ub                  /**< values for the new upper bounds or NULL if ncols is zero */
   )
{
   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(ncols == 0 || (ind != NULL && lb != NULL && ub != NULL));

   SCIPdebugMessage("calling SCIPlpiChgBounds()\n");
   if( ncols <= 0 )
      return SCIP_OKAY;

   invalidateSolution(lpi);

   ClpSimplex* clp = lpi->clp;

#if SCIP_DISABLED_CODE
   /* The following bugfix was necessary some time ago to avoid an error in Clp and can currently be disabled: the
    * solution vector is modified to be set to the corresponding bounds. Remove if Clp versions have stabilized. */
   double* sol = lpi->clp->primalColumnSolution();
   const double* colLower = lpi->clp->getColLower();
   const double* colUpper = lpi->clp->getColUpper();
#endif

   for (int j = 0; j < ncols; ++j)
   {
      if ( SCIPlpiIsInfinity(lpi, lb[j]) )
      {
         SCIPerrorMessage("LP Error: fixing lower bound for variable %d to infinity.\n", ind[j]);
         return SCIP_LPERROR;
      }
      if ( SCIPlpiIsInfinity(lpi, -ub[j]) )
      {
         SCIPerrorMessage("LP Error: fixing upper bound for variable %d to -infinity.\n", ind[j]);
         return SCIP_LPERROR;
      }

      clp->setColumnBounds(ind[j], lb[j], ub[j]);

#if SCIP_DISABLED_CODE
      /* Old bugfix, not currently needed - see above */
      if ( sol != 0 )
      {
         if( clp->statusExists() )
         {
            assert( colLower != 0 );
            assert( colUpper != 0 );
            int k = ind[j];
            switch ( clp->getColumnStatus(k) )
            {
               case ClpSimplex::isFree:
               case ClpSimplex::superBasic:
                  sol[k] = 0.0;
                  break;
               case ClpSimplex::atUpperBound:
                  sol[k] = colUpper[k];
                  assert( colUpper[k] == ub[j] );
                  break;
               case ClpSimplex::isFixed:
               case ClpSimplex::atLowerBound:
                  sol[k] = colLower[k];
                  assert( colLower[k] == lb[j] );
                  break;
               default:;
            }
         }
         else
         { /* workaround: if there is no status, we assume something */
            sol[ind[j]] = 0.0;
         }
      }
#endif
   }

   return SCIP_OKAY;
}


/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiChgSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgSides()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(ind != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   if( nrows <= 0)
      return SCIP_OKAY;

   invalidateSolution(lpi);

   ClpSimplex* clp = lpi->clp;

   for (int i = 0; i < nrows; ++i)
      clp->setRowBounds(ind[i], lhs[i], rhs[i]);

   return SCIP_OKAY;
}


/** changes a single coefficient */
SCIP_RETCODE SCIPlpiChgCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgCoef()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(0 <= row && row < lpi->clp->numberRows());
   assert(0 <= col && col < lpi->clp->numberColumns());

   invalidateSolution(lpi);

   lpi->clp->matrix()->modifyCoefficient(row, col, newval);

   return SCIP_OKAY;
}


/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgObjsen()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   invalidateSolution(lpi);

   // set objective sense: SCIP values are the same as the ones for Clp
   lpi->clp->setOptimizationDirection(objsen);

   return SCIP_OKAY;
}


/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiChgObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   const int*            ind,                /**< column indices to change objective value for */
   const SCIP_Real*      obj                 /**< new objective values for columns */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgObj()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(ind != NULL);
   assert(obj != NULL);

   invalidateSolution(lpi);

   ClpSimplex* clp = lpi->clp;

   // updates whatsChanged in Clp (bound checking in Clp)
   for( int j = 0; j < ncols; ++j )
      clp->setObjCoeff(ind[j], obj[j]);  // inlined version of clp->setObjectiveCoefficient(ind[j], obj[j]);

   return SCIP_OKAY;
}


/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIPdebugMessage("calling SCIPlpiScaleRow()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(scaleval != 0.0);
   assert(0 <= row && row <= lpi->clp->numberRows() );

   invalidateSolution(lpi);

   // Note: if the scaling should be performed because of numerical stability,
   // there are other more effective methods in Clp to adjust the scaling values
   // for each row.

   ClpSimplex* clp = lpi->clp;

   // adjust the sides
   double* lhs = clp->rowLower();
   double* rhs = clp->rowUpper();

   double lhsval = lhs[row];
   if( lhsval > -COIN_DBL_MAX )
      lhsval *= scaleval;
   else if( scaleval < 0.0 )
      lhsval = COIN_DBL_MAX;
   double rhsval = rhs[row];
   if( rhsval < COIN_DBL_MAX)
      rhsval *= scaleval;
   else if( scaleval < 0.0 )
      rhsval = -COIN_DBL_MAX;
   if( scaleval < 0.0 )
   {
      SCIP_Real oldlhs = lhsval;
      lhsval = rhsval;
      rhsval = oldlhs;
   }
   lhs[row] = lhsval;    // change values directly into Clp data!
   rhs[row] = rhsval;

   // apply scaling ...

   // WARNING: the following is quite expensive:
   // We have to loop over the matrix to find the row entries.
   // For columns we can do better, see @c SCIPlpiScaleCol.
   CoinPackedMatrix* M = clp->matrix();
   assert( M->getNumCols() == clp->numberColumns() );

   const CoinBigIndex* beg = M->getVectorStarts();
   const int* length = M->getVectorLengths();
   const int* ind = M->getIndices();
   double* val = M->getMutableElements();

   for (int j = 0; j < M->getNumCols(); ++j)
   {
      for (CoinBigIndex k = beg[j]; k < beg[j] + length[j]; ++k)
      {
	 if (ind[k] == row)
	    val[k] *= scaleval;
      }
   }

   return SCIP_OKAY;
}


/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPlpiScaleCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIPdebugMessage("calling SCIPlpiScaleCol()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(scaleval != 0.0);
   assert(0 <= col && col <= lpi->clp->numberColumns() );

   invalidateSolution(lpi);

   // Note: if the scaling should be performed because of numerical stability,
   // there are other more effective methods in Clp to adjust the scaling values
   // for each column.

   ClpSimplex* clp = lpi->clp;

   // adjust the objective coefficients
   double* objvec = clp->objective();          // we have direct access to the data of Clp!
   objvec[col] *= scaleval;                    // adjust the objective function value

   // adjust the bounds
   double* lb = clp->columnLower();
   double* ub = clp->columnUpper();
   double lbval = lb[col];
   double ubval = ub[col];

   if( lbval > -COIN_DBL_MAX )
      lbval /= scaleval;
   else if( scaleval < 0.0 )
      lbval = COIN_DBL_MAX;
   if( ubval < COIN_DBL_MAX )
      ubval /= scaleval;
   else if( scaleval < 0.0 )
      ubval = -COIN_DBL_MAX;
   if( scaleval < 0.0 )
   {
      SCIP_Real oldlb = lbval;
      lbval = ubval;
      ubval = oldlb;
   }
   lb[col] = lbval;        // directly adjust values into Clp data
   ub[col] = ubval;

   // apply scaling directly to matrix (adapted from ClpPackedMatrix::reallyScale)
   // See also ClpModel::gutsOfScaling ...
   CoinPackedMatrix* M = clp->matrix();
   assert( M->getNumCols() == clp->numberColumns() );

   const CoinBigIndex* beg = M->getVectorStarts();
   const int* length = M->getVectorLengths();
   double* val = M->getMutableElements();
   for (CoinBigIndex k = beg[col]; k < beg[col] + length[col]; ++k)
      val[k] *= scaleval;

   return SCIP_OKAY;
}



/**@} */




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiGetNRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNRows()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(nrows != NULL);

   *nrows = lpi->clp->numberRows();

   return SCIP_OKAY;
}


/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNCols()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(ncols != NULL);

   *ncols = lpi->clp->numberColumns();

   return SCIP_OKAY;
}


/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNNonz()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(nnonz != NULL);

   *nnonz = lpi->clp->getNumElements();

   return SCIP_OKAY;
}


/** gets columns from LP problem object; the arrays have to be large enough to store all values
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetCols()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());
   assert((lb != NULL && ub != NULL) || (lb == NULL && ub == NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));

   ClpSimplex* clp = lpi->clp;

   // get lower and upper bounds for the variables
   if ( lb != NULL )
   {
      const double* colLower = clp->getColLower();    // Here we can use the const versions (see SCIPchgBounds)
      const double* colUpper = clp->getColUpper();

      BMScopyMemoryArray( lb, colLower + firstcol, (lastcol - firstcol + 1));
      BMScopyMemoryArray( ub, colUpper + firstcol, (lastcol - firstcol + 1));
   }

   if ( nnonz != NULL )
   {
      CoinPackedMatrix* M = clp->matrix();
      assert( M != NULL );
      assert( M->getNumCols() == clp->numberColumns() );

      const CoinBigIndex* Mbeg = M->getVectorStarts();   // can use const versions
      const int* Mlength = M->getVectorLengths();
      const int* Mind = M->getIndices();
      const double* Mval = M->getElements();

      *nnonz = 0;
      // can we use memcpy for the whole set (requires that columns are stored sequentially)
      for (int j = firstcol; j <= lastcol; ++j)
      {
         beg[j-firstcol] = *nnonz;

         BMScopyMemoryArray( (ind + (*nnonz)), Mind + Mbeg[j], Mlength[j]);
         BMScopyMemoryArray( (val + (*nnonz)), Mval + Mbeg[j], Mlength[j]);

         (*nnonz) += Mlength[j];
      }
   }

   return SCIP_OKAY;
}


/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_Real*            lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetRows()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->clp->numberRows());
   assert((lhs != NULL && rhs != NULL) || (lhs == NULL && rhs == NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));

   ClpSimplex* clp = lpi->clp;
   if ( lhs != NULL )
   {
      const double* rowLower = clp->getRowLower();    // Here we can use the const versions (see SCIPchgSides)
      const double* rowUpper = clp->getRowUpper();

      BMScopyMemoryArray( lhs, rowLower + firstrow, (lastrow - firstrow + 1) );
      BMScopyMemoryArray( rhs, rowUpper + firstrow, (lastrow - firstrow + 1) );
   }

   if ( nnonz != NULL )
   {
      ClpMatrixBase* M = clp->rowCopy();   // get row view on matrix
      if ( M == NULL ) // can happen e.g. if no LP was solved yet ...
	 M = clp->clpMatrix()->reverseOrderedCopy();
      assert( M != NULL );
      assert( M->getNumRows() == clp->numberRows() );

      const CoinBigIndex* Mbeg = M->getVectorStarts();
      const int* Mlength = M->getVectorLengths();
      const int* Mind = M->getIndices();
      const double* Mval = M->getElements();

      *nnonz = 0;
      for( int i = firstrow; i <= lastrow; ++i )
      {
         beg[i-firstrow] = *nnonz;
         for( CoinBigIndex k = Mbeg[i]; k < Mbeg[i] + Mlength[i]; ++k )
         {
            ind[*nnonz] = Mind[k];
            val[*nnonz] = Mval[k];
            (*nnonz)++;
         }
      }
   }

   return SCIP_OKAY;
}


/** gets column names */
SCIP_RETCODE SCIPlpiGetColNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for col names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   )
{
   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(colnames != NULL || namestoragesize == 0);
   assert(namestorage != NULL || namestoragesize == 0);
   assert(namestoragesize >= 0);
   assert(storageleft != NULL);

   SCIPerrorMessage("SCIPlpiGetColNames() has not been implemented yet.\n");

   return SCIP_LPERROR;
}


/** gets row names */
SCIP_RETCODE SCIPlpiGetRowNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get name from LP */
   int                   lastrow,            /**< last row to get name from LP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for row names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   )
{
   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(rownames != NULL || namestoragesize == 0);
   assert(namestorage != NULL || namestoragesize == 0);
   assert(namestoragesize >= 0);
   assert(storageleft != NULL);

   SCIPerrorMessage("SCIPlpiGetRowNames() has not been implemented yet.\n");

   return SCIP_LPERROR;
}


/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   SCIPdebugMessage("calling SCIPlpiIgnoreInstability()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(success != NULL);

   /* Unstable situations are currently not ignored. Could fix this similar to lpi_cpx by adjusting the solution status. */
   *success = FALSE;

   return SCIP_OKAY;
}


/** gets the objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   assert( lpi != NULL );
   assert( lpi->clp != NULL );
   assert( objsen != NULL );

   // Clp direction of optimization (1 - minimize, -1 - maximize, 0 - ignore)
   if ( lpi->clp->getObjSense() < 0 )
      *objsen = SCIP_OBJSEN_MAXIMIZE;
   else
      *objsen = SCIP_OBJSEN_MINIMIZE;

   return SCIP_OKAY;
}


/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiGetObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetObj()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());
   assert(vals != NULL);

   const double* obj = lpi->clp->getObjCoefficients();    // Here we can use the const versions (see SCIPchgObj)

   BMScopyMemoryArray(vals, obj + firstcol, (lastcol - firstcol + 1) );

   return SCIP_OKAY;
}


/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective value for */
   int                   lastcol,            /**< last column to get objective value for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBounds()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());

   if ( lbs != 0 )
   {
      const double* colLower = lpi->clp->getColLower();    // Here we can use the const versions (see SCIPchgBounds)
      BMScopyMemoryArray( lbs, colLower + firstcol, (lastcol - firstcol + 1) );
   }

   if ( ubs != 0 )
   {
      const double* colUpper = lpi->clp->getColUpper();
      BMScopyMemoryArray( ubs, colUpper + firstcol, (lastcol - firstcol + 1) );
   }

   return SCIP_OKAY;
}


/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiGetSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSides()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->clp->numberRows());

   if ( lhss != 0 )
   {
      const double* rowLower = lpi->clp->getRowLower();    // Here we can use the const versions (see SCIPchgSides)
      BMScopyMemoryArray( lhss, rowLower + firstrow, (lastrow - firstrow + 1) );
   }

   if ( rhss != 0 )
   {
      const double* rowUpper = lpi->clp->getRowUpper();
      BMScopyMemoryArray( rhss,  rowUpper + firstrow, (lastrow - firstrow + 1) );
   }

   return SCIP_OKAY;
}


/** gets a single coefficient */
SCIP_RETCODE SCIPlpiGetCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetCoef()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(0 <= col && col < lpi->clp->numberColumns());
   assert(0 <= row && row < lpi->clp->numberRows());
   assert(val != NULL);

   *val = lpi->clp->matrix()->getCoefficient(row, col);

   return SCIP_OKAY;
}

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */


/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   SCIPdebugMessage("calling Clp primal(): %d cols, %d rows\n", lpi->clp->numberColumns(), lpi->clp->numberRows());

#ifdef LPI_CLP_DEBUG_WRITE_FILES
   char filename[255];
   snprintf(filename, 255, "debug_p_%d.mps", fileNr);
   fileNr = fileNr % 2;
   SCIPlpiWriteLP(lpi, filename);
   SCIPdebugMessage("Wrote file <%s>\n", filename);
#endif

   invalidateSolution(lpi);

   // initialize factorization freq. depending on model size - applied only once
   setFactorizationFrequency(lpi);

   // if we want to construct a new basis
   if ( lpi->startscratch )
   {
      lpi->clp->allSlackBasis(true);   // reset basis
      lpi->validFactorization = false;
   }

   /*  startFinishOptions - bits
    *  1 - do not delete work areas and factorization at end
    *  2 - use old factorization if same number of rows
    *  4 - skip as much initialization of work areas as possible (work in progress)
    *
    *  4 does not seem to work.
    */
   int startFinishOptions = 1;
   if ( lpi->validFactorization )
      startFinishOptions = startFinishOptions | 2;

   /* Primal algorithm */
   int status = lpi->clp->primal(0, startFinishOptions);

#ifdef LPI_CLP_DEBUG_WRITE_FILES
   char basisname[255];
   snprintf(basisname, 255, "debug_p_%d.bas", fileNr);
   SCIP_CALL( SCIPlpiWriteState(lpi, basisname) );
   SCIPdebugMessage("Wrote basis file <%s>\n", basisname);
   ++fileNr; /* not increased above! */
   fileNr = fileNr % 2;
#endif

   lpi->lastalgorithm = 1;
   lpi->validFactorization = true;
   lpi->solved = TRUE;

   // Unfortunately the status of Clp is hard coded ...
   // -1 - did not run
   //  0 - optimal
   //  1 - primal infeasible
   //  2 - dual infeasible
   //  3 - stopped on iterations or time
   //  4 - stopped due to errors
   //  5 - stopped by event handler
   assert( status != -1 );      // did not run should not occur
   assert( status != 5 );       // begin stopped by event handler should not occur

   if ( status == 4 || status == 5 || status == -1 )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}


/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   SCIPdebugMessage("calling Clp dual(): %d cols, %d rows\n", lpi->clp->numberColumns(), lpi->clp->numberRows());

#ifdef LPI_CLP_DEBUG_WRITE_FILES
   char filename[255];
   snprintf(filename, 255, "debug_d_%d.mps", fileNr);
   SCIPlpiWriteLP(lpi, filename);
   SCIPdebugMessage("Wrote file <%s>\n", filename);
   snprintf(filename, 255, "debug_d_%d.sav", fileNr);
   // lpi->clp->saveModel(filename);
   SCIPdebugMessage("Wrote file <%s>\n", filename);
#endif

   invalidateSolution(lpi);

   // intialize factorization freq. depending on model size - applied only once
   setFactorizationFrequency(lpi);

   // if we want to construct a new basis
   if( lpi->startscratch )
   {
      lpi->clp->allSlackBasis(true);   // reset basis
      lpi->validFactorization = false;
   }

   /*  startFinishOptions - bits
    *  1 - do not delete work areas and factorization at end
    *  2 - use old factorization if same number of rows
    *  4 - skip as much initialization of work areas as possible (work in progress)
    *
    *  4 does not seem to work.
    */
   int startFinishOptions = 1;
   if ( lpi->validFactorization )
      startFinishOptions = startFinishOptions | 2;

   /* Dual algorithm */
   int status = lpi->clp->dual(0, startFinishOptions);

#ifdef LPI_CLP_DEBUG_WRITE_FILES
   char basisname[255];
   snprintf(basisname, 255, "debug_d_%d.bas", fileNr);
   SCIP_CALL( SCIPlpiWriteState(lpi, basisname) );
   SCIPdebugMessage("Wrote basis file <%s>\n", basisname);
   ++fileNr; /* not increased above! */
   fileNr = fileNr % 2;
#endif

   lpi->lastalgorithm = -1;
   lpi->validFactorization = true;
   lpi->solved = TRUE;

   // Unfortunately the status of Clp is hard coded ...
   // -1 - did not run
   //  0 - optimal
   //  1 - primal infeasible
   //  2 - dual infeasible
   //  3 - stopped on iterations or time
   //  4 - stopped due to errors
   //  5 - stopped by event handler
   assert( status != -1 );      // did not run should not occur
   assert( status != 5 );       // begin stopped by event handler should not occur

   if ( status == 4 || status == 5 || status == -1 )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}


/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{
   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   SCIPdebugMessage("calling Clp barrier(): %d cols, %d rows; crossover: %u\n", lpi->clp->numberColumns(), lpi->clp->numberRows(), crossover);

   invalidateSolution(lpi);

   // Check whether we have a factorization, if yes destroy it (Clp doesn't like it ...)
   /*
   if (lpi->haveFactorization)
      lpi->clp->finish();
   */

   // call barrier
#if (CLP_VERSION_MAJOR >= 1 && CLP_VERSION_MINOR > 17) || CLP_VERSION_MAJOR >= 2
   int startFinishOptions = 1;
   int status = lpi->clp->barrier(crossover, startFinishOptions);
#else
   int status = lpi->clp->barrier(crossover);
#endif

   lpi->lastalgorithm = 2;
   lpi->solved = TRUE;

   // We may need to call ClpModel::status()

   // Unfortunately the status of Clp is hard coded ...
   // -1 - did not run
   //  0 - optimal
   //  1 - primal infeasible
   //  2 - dual infeasible
   //  3 - stopped on iterations or time
   //  4 - stopped due to errors
   //  5 - stopped by event handler
   assert( status != -1 );      // did not run should not occur
   assert( status != 5 );       // begin stopped by event handler should not occur

   if ( status == 4 || status == 5 || status == -1 )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/** start strong branching - call before any strongbranching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   // currently do nothing; in the future: use code as in OSI
   return SCIP_OKAY;
}

/** end strong branching - call after any strongbranching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   // currently do nothing; in the future: use code as in OSI
   return SCIP_OKAY;
}

/** performs strong branching iterations on one arbitrary candidate */
static
SCIP_RETCODE lpiStrongbranch(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiStrongbranch() on variable %d (%d iterations)\n", col, itlim);

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   ClpSimplex* clp = lpi->clp;

   // set up output arrays
   int ncols = clp->numberColumns();
   assert( 0 <= col && col < ncols );
   double** outputSolution = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &outputSolution, 2) );
   SCIP_ALLOC( BMSallocMemoryArray( &outputSolution[0], ncols) );
   SCIP_ALLOC( BMSallocMemoryArray( &outputSolution[1], ncols) );

   int* outputStatus = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &outputStatus, 2) );

   int* outputIterations = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &outputIterations, 2) );

   // set iteration limit
   int iterlimit = clp->maximumIterations();
   clp->setMaximumIterations(itlim);

   // store objective value
   double objval = clp->objectiveValue();

   // store special options for later reset
   int specialoptions = clp->specialOptions();

   // lpi->clp->setSpecialOptions(64|128|512|1024|2048|4096|32768|262144|0x02000000);
   // use default settings:
   lpi->clp->setSpecialOptions(32|64|128|512|1024|2048|4096|32768|262144|2097152|0x2000000);

   /* 'startfinish' options for strong branching:
    *  1 - do not delete work areas and factorization at end
    *  2 - use old factorization if same number of rows
    *  4 - skip as much initialization of work areas as possible
    *      (based on whatsChanged in clpmodel.hpp) ** work in progress
    *
    *  4 does not seem to work in strong branching ...
    */
   int startFinishOptions = 1;
   if ( lpi->validFactorization )
      startFinishOptions = startFinishOptions | 2;

   // set new lower and upper bounds for variable
   *down = EPSCEIL(psol - 1.0, 1e-06);
   *up   = EPSFLOOR(psol + 1.0, 1e-06);

   /*  For strong branching.  On input lower and upper are new bounds while
    *  on output they are change in objective function values (>1.0e50
    *  infeasible).  Return code is
    *   0 if nothing interesting,
    *  -1 if infeasible both ways and
    *  +1 if infeasible one way (check values to see which one(s))
    *  -2 if bad factorization
    * Solutions are filled in as well - even down, odd up - also status and number of iterations
    *
    * The bools are:
    *   bool stopOnFirstInfeasible
    *   bool alwaysFinish
    *
    * At the moment: we need alwaysFinish to get correct bounds.
    */
   int res = clp->strongBranching(1, &col, up, down, outputSolution, outputStatus, outputIterations, false, true, startFinishOptions);

   // reset special options
   clp->setSpecialOptions(specialoptions);

   lpi->validFactorization = true;

   *down += objval;
   *up += objval;

   // The bounds returned by CLP seem to be valid using the above options
   *downvalid = TRUE;
   *upvalid = TRUE;

   // correct iteration count
   if (iter)
      *iter = outputIterations[0] + outputIterations[1];

   // reset iteration limit
   clp->setMaximumIterations(iterlimit);

   // free local memory
   BMSfreeMemoryArray( &outputStatus );
   BMSfreeMemoryArray( &outputIterations );
   BMSfreeMemoryArray( &outputSolution[1] );
   BMSfreeMemoryArray( &outputSolution[0] );
   BMSfreeMemoryArray( &outputSolution );

   if ( res == -2 )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/** performs strong branching iterations on given arbitrary candidates */
static
SCIP_RETCODE lpiStrongbranches(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< fractional current primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiStrongbranches() on %d variables (%d iterations)\n", ncols, itlim);

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   ClpSimplex* clp = lpi->clp;

   // set up output arrays
   int n = clp->numberColumns();
   assert( 0 < ncols && ncols <= n );
   double** outputSolution = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &outputSolution, 2*ncols) );
   for (int j = 0; j < 2*ncols; ++j)
   {
      SCIP_ALLOC( BMSallocMemoryArray( &(outputSolution[j]), n) );
   }

   int* outputStatus = NULL;
   SCIP_ALLOC( BMSallocMemoryArray(&outputStatus, 2*ncols) );

   int* outputIterations = NULL;
   SCIP_ALLOC( BMSallocMemoryArray(&outputIterations, 2*ncols) );

   // set iteration limit
   int iterlimit = clp->maximumIterations();
   clp->setMaximumIterations(itlim);

   // store objective value
   double objval = clp->objectiveValue();

   // store special options for later reset
   int specialoptions = clp->specialOptions();

   // lpi->clp->setSpecialOptions(64|128|512|1024|2048|4096|32768|262144|0x02000000);
   // use default settings:
   lpi->clp->setSpecialOptions(32|64|128|512|1024|2048|4096|32768|262144|2097152|0x2000000);

   /* 'startfinish' options for strong branching:
    *  1 - do not delete work areas and factorization at end
    *  2 - use old factorization if same number of rows
    *  4 - skip as much initialization of work areas as possible
    *      (based on whatsChanged in clpmodel.hpp) ** work in progress
    *
    *  4 does not seem to work in strong branching ...
    */
   int startFinishOptions = 1;
   if ( lpi->validFactorization )
      startFinishOptions = startFinishOptions | 2;

   // set new lower and upper bounds for variables
   for (int j = 0; j < ncols; ++j)
   {
      assert( 0 <= cols[j] && cols[j] < n );
      down[j] = EPSCEIL(psols[j] - 1.0, 1e-06);
      up[j]   = EPSFLOOR(psols[j] + 1.0, 1e-06);

      // The bounds returned by CLP seem to be valid using the above options
      downvalid[j] = TRUE;
      upvalid[j] = TRUE;
   }

   /*  For strong branching.  On input lower and upper are new bounds while
    *  on output they are change in objective function values (>1.0e50
    *  infeasible).  Return code is
    *   0 if nothing interesting,
    *  -1 if infeasible both ways and
    *  +1 if infeasible one way (check values to see which one(s))
    *  -2 if bad factorization
    * Solutions are filled in as well - even down, odd up - also status and number of iterations
    *
    * The bools are:
    *   bool stopOnFirstInfeasible
    *   bool alwaysFinish
    *
    * At the moment: we need alwaysFinish to get correct bounds.
    */
   int res = clp->strongBranching(ncols, cols, up, down, outputSolution, outputStatus, outputIterations, false, true, startFinishOptions);

   // reset special options
   clp->setSpecialOptions(specialoptions);

   lpi->validFactorization = true;

   for (int j = 0; j < ncols; ++j)
   {
      down[j] += objval;
      up[j] += objval;

      // correct iteration count
      if (iter)
         *iter += outputIterations[2*j] + outputIterations[2*j+1];

      BMSfreeMemoryArray(&outputSolution[2*j]);
      BMSfreeMemoryArray(&outputSolution[2*j+1]);
   }

   // reset iteration limit
   clp->setMaximumIterations(iterlimit);

   // free local memory
   BMSfreeMemoryArray( &outputStatus );
   BMSfreeMemoryArray( &outputIterations );
   BMSfreeMemoryArray( &outputSolution );

   if ( res == -2 )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPlpiStrongbranchFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   /* pass call on to lpiStrongbranch() */
   SCIP_CALL( lpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given @b fractional candidates */
SCIP_RETCODE SCIPlpiStrongbranchesFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< fractional current primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   if ( iter != NULL )
      *iter = 0;

   /* pass call on to lpiStrongbranches() */
   SCIP_CALL( lpiStrongbranches(lpi, cols, ncols, psols, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/** performs strong branching iterations on one candidate with @b integral value */
SCIP_RETCODE SCIPlpiStrongbranchInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current integral primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   /* pass call on to lpiStrongbranch() */
   SCIP_CALL( lpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given candidates with @b integral values */
SCIP_RETCODE SCIPlpiStrongbranchesInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< current integral primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   if ( iter != NULL )
      *iter = 0;

   /* pass call on to lpiStrongbranches() */
   SCIP_CALL( lpiStrongbranches(lpi, cols, ncols, psols, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/**@} */



/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiWasSolved(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   return lpi->solved;
}

/** gets information about primal and dual feasibility of the current LP solution
 *
 *  The feasibility information is with respect to the last solving call and it is only relevant if SCIPlpiWasSolved()
 *  returns true. If the LP is changed, this information might be invalidated.
 *
 *  Note that @a primalfeasible and @a dualfeasible should only return true if the solver has proved the respective LP to
 *  be feasible. Thus, the return values should be equal to the values of SCIPlpiIsPrimalFeasible() and
 *  SCIPlpiIsDualFeasible(), respectively. Note that if feasibility cannot be proved, they should return false (even if
 *  the problem might actually be feasible).
 */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< pointer to store primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< pointer to store dual feasibility status */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSolFeasibility()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);

   if ( lpi->clp->primalFeasible() )
      *primalfeasible = TRUE;
   else
      *primalfeasible = FALSE;

   if ( lpi->clp->dualFeasible() )
      *dualfeasible = TRUE;
   else
      *dualfeasible = FALSE;

   // say feasible if deviation is small
   if (lpi->clp->status()==0 && ( ! (*primalfeasible) || ! (*dualfeasible)) )
   {
      if ( !(*primalfeasible) && lpi->clp->sumPrimalInfeasibilities() < SUMINFEASBOUND )
      {
         lpi->clp->setNumberPrimalInfeasibilities(0);
         *primalfeasible = TRUE;
      }
      if ( !(*dualfeasible) && lpi->clp->sumDualInfeasibilities() < SUMINFEASBOUND)
      {
         lpi->clp->setNumberDualInfeasibilities(0);
         *dualfeasible = TRUE;
      }
   }

   return SCIP_OKAY;
}


/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExistsPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   /* Clp usually has a primal ray whenever it concludes "dual infeasible" (status == 2)
    * (but is not necessarily primal feasible), see ClpModel::unboundedRay(). */
   return ( lpi->clp->status() == 2 );
}


/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   /* Clp usually has a primal ray whenever it concludes "dual infeasible" (status == 2)
    * (but is not necessarily primal feasible), see ClpModel::unboundedRay(). */
   if ( lpi->clp->rayExists() )
   {
      return ( lpi->clp->status() == 2 );
   }
   return FALSE;
}


/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   return ( lpi->clp->isProvenDualInfeasible() && lpi->clp->primalFeasible() );
}


/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   /* Should return ClpModel::isProvenPrimalInfeasible() (which returns "status == 1"), but the
    * following is correct (Clp will not be changed). The secondaryStatus is 1 if the dual simplex
    * detects an objective limit exceedence. The primal simplex has no such detection (will never
    * stop with objective limit exceedence). Hence we are infeasible only if status == 1 and we have
    * not stopped due to the objective limit. */
   return ( lpi->clp->status() == 1 && (lpi->clp->secondaryStatus() == 0 || lpi->clp->secondaryStatus() == 6) );
}


/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   return ( lpi->clp->primalFeasible() );
}


/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExistsDualRay()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   /* Clp usually has a dual ray whenever it concludes "primal infeasible" (but is not necessarily dual feasible), see
    * ClpModel::infeasibilityRay. Additionally check whether ray exists in order to avoid situations in which Clp cannot
    * provide a ray. SCIP often decides to resolve in such a case and the problem might go away. */
   return ( lpi->clp->status() == 1 && lpi->clp->secondaryStatus() == 0 && lpi->clp->rayExists() );
}


/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasDualRay()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   /* Clp usually has a dual ray whenever it concludes "primal infeasible" (but is not necessarily dual feasible),
    * see ClpModel::infeasibilityRay. Additionally check whether ray exists. */
   if ( lpi->clp->rayExists() )
   {
      if ( lpi->clp->status() == 1 && lpi->clp->secondaryStatus() == 0 )
         return TRUE;
   }

   return FALSE;
}


/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   /* The dual seems to be unbounded if the status is 1 (primal unbounded), the secondaryStatus is
    * not 1 (i.e., the dual simplex has not stopped because of an objective limit exceedence), and
    * the dual is feasible. */
   return ( lpi->clp->status() == 1 && lpi->clp->secondaryStatus() == 0 && lpi->clp->dualFeasible() );
}


/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   return ( lpi->clp->isProvenDualInfeasible() );
}


/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   return ( lpi->clp->dualFeasible() );
}


/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsOptimal()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   if ( SCIPlpiIsObjlimExc(lpi) )
      return FALSE;

   /* secondaryStatus == 6 means that the problem is empty */
   return( lpi->clp->isProvenOptimal() && (lpi->clp->secondaryStatus() == 0 || lpi->clp->secondaryStatus() == 6));
}


/** returns TRUE iff current LP solution is stable
 *
 *  This function should return true if the solution is reliable, i.e., feasible and optimal (or proven
 *  infeasible/unbounded) with respect to the original problem. The optimality status might be with respect to a scaled
 *  version of the problem, but the solution might not be feasible to the unscaled original problem; in this case,
 *  SCIPlpiIsStable() should return false.
 */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsStable()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   /* Return false if infeasible, but dual ray is not present and the algorithm type has changed. This is one of the
    * cases in which Clp cannot produce a ray and the hope is to get a ray by rerunning. */
   if ( lpi->clp->status() == 1 && lpi->lastalgorithm != lpi->clp->algorithm() && ! lpi->clp->rayExists() )
      return FALSE;

   /*  We first check if status is ok, i.e., is one of the following:
    *   0 - optimal
    *   1 - primal infeasible
    *   2 - dual infeasible
    *   3 - stopped on iterations or time
    *   4 - stopped due to errors
    *   5 - stopped by event handler (virtual int ClpEventHandler::event())
    *
    *  Then we check the secondary status of Clp:
    *   0 - none
    *   1 - primal infeasible because dual limit reached OR (probably primal infeasible but can't prove it  - main status was 4)
    *   2 - scaled problem optimal - unscaled problem has primal infeasibilities
    *   3 - scaled problem optimal - unscaled problem has dual infeasibilities
    *   4 - scaled problem optimal - unscaled problem has primal and dual infeasibilities
    *   5 - giving up in primal with flagged variables
    *   6 - failed due to empty problem check
    *   7 - postSolve says not optimal
    *   8 - failed due to bad element check
    *   9 - status was 3 and stopped on time
    * 100 up - translation of enum from ClpEventHandler
    */
   SCIPdebugMessage("status: %d   secondary: %d\n", lpi->clp->status(), lpi->clp->secondaryStatus());
   assert( 0 <= lpi->clp->status() && lpi->clp->status() <= 5 );

   return( (lpi->clp->status() <= 3) && (lpi->clp->secondaryStatus() <= 1 || lpi->clp->secondaryStatus() == 6 || lpi->clp->secondaryStatus() == 9) );
}


/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsObjlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   /* if status == 1 (primal infeasible) and secondaryStatus == 1 then Clp hit the dual bound */
   if ( lpi->clp->status() == 1 )
   {
      if ( lpi->clp->secondaryStatus() == 1 )
	 return TRUE;
      else
	 return FALSE;
   }

   return ( lpi->clp->isObjectiveLimitTestValid() && lpi->clp->isDualObjectiveLimitReached() );

   /* The above code is equivalent to the following:
   if ( lpi->clp->status() == 0 || (lpi->clp->status() == 1 && lpi->clp->algorithm() < 0) || (lpi->clp->status() == 2 && lpi->clp->algorithm() > 0) )
   {
      return ( lpi->clp->isPrimalObjectiveLimitReached() || lpi->clp->isDualObjectiveLimitReached() );
   }
   */
}


/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsIterlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   /* status == 3 means that Clp stopped on time or iteration limit
    * secondary status == 9 means that status was 3 and Clp stopped on time */
   return ( lpi->clp->status() == 3 && lpi->clp->secondaryStatus() != 9 );
}


/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   /* status == 3 means that Clp stopped on time or iteration limit
    * secondary status == 9 means that status was 3 and Clp stopped on time */
   return ( lpi->clp->status() == 3 && lpi->clp->secondaryStatus() == 9 );
}


/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetInternalStatus()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   return lpi->clp->status();
}


/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetObjval()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(objval != NULL);

   *objval = lpi->clp->objectiveValue();

   return SCIP_OKAY;
}


/** gets primal and dual solution vectors for feasible LPs
 *
 *  Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
 *  SCIPlpiIsOptimal() returns true.
 */
SCIP_RETCODE SCIPlpiGetSol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Real*            activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSol()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   ClpSimplex* clp = lpi->clp;
   if( objval != NULL )
      *objval = clp->objectiveValue();

   if( primsol != NULL )
   {
      const double* sol = clp->getColSolution();
      BMScopyMemoryArray( primsol, sol, clp->numberColumns() );
   }
   if( dualsol != NULL )
   {
      const double* dsol = clp->getRowPrice();
      BMScopyMemoryArray( dualsol, dsol, clp->numberRows() );
   }
   if( activity != NULL )
   {
      const double* act = clp->getRowActivity();
      BMScopyMemoryArray( activity, act, clp->numberRows() );
   }
   if( redcost != NULL )
   {
      const double* red = clp->getReducedCost();
      BMScopyMemoryArray( redcost, red, clp->numberColumns() );
   }

   return SCIP_OKAY;
}


/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(ray != NULL);

   /* Unbounded ray (NULL returned if none/wrong). Up to user to use delete [] on these arrays.  */
   const double* clpray = lpi->clp->unboundedRay();

   if ( clpray == NULL )
      return SCIP_LPERROR;

   BMScopyMemoryArray(ray, clpray, lpi->clp->numberColumns());

   delete [] clpray;

   return SCIP_OKAY;
}

/** gets dual farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual farkas row multipliers */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetDualfarkas()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(dualfarkas != NULL);

   /* Infeasibility ray (NULL returned if none/wrong). Up to user to use delete [] on these arrays. */
   const double* dualray = lpi->clp->infeasibilityRay();

   if ( dualray == NULL )
      return SCIP_LPERROR;

   /* The dual ray returned by Clp sometimes contains large numbers. We try to scale the vector. First compute maximal
    * and minimal absolute values. */
   double minabsvalue = SCIPlpiInfinity(lpi);
   double maxabsvalue = 0.0;
   double feastol = lpi->clp->primalTolerance();
   for (int j = 0; j < lpi->clp->numberRows(); ++j)
   {
      double val = fabs(dualray[j]);

      /* only consider nonzero entries */
      if ( val >= feastol )
      {
         if ( val > maxabsvalue )
            maxabsvalue = val;
         if ( val < minabsvalue )
            minabsvalue = val;
      }
   }

   /* Possibly scale and also convert sign. */
   if ( maxabsvalue > 0.0 )
   {
      assert( 0.0 < minabsvalue && minabsvalue <= maxabsvalue );

      /* We try to make the maximum absolute value to be 1.0, but if the minimal absolute value would be less than the
       * feasibility tolerance, we adjust the factor such that it will be equal to the feasibility tolerance. */
      double scalingfactor = maxabsvalue;
      if ( minabsvalue / scalingfactor < feastol )
         scalingfactor = minabsvalue / feastol;

      for (int j = 0; j < lpi->clp->numberRows(); ++j)
         dualfarkas[j] = -dualray[j]/scalingfactor;
   }
   else
   {
      /* convert sign */
      for (int j = 0; j < lpi->clp->numberRows(); ++j)
         dualfarkas[j] = -dualray[j];
   }

   delete [] dualray;

   return SCIP_OKAY;
}


/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(iterations != NULL);

   *iterations = lpi->clp->numberIterations();

   return SCIP_OKAY;
}

/** gets information about the quality of an LP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The LPI should return SCIP_INVALID for @p quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPlpiGetRealSolQuality(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPSOLQUALITY     qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   assert(lpi != NULL);
   assert(quality != NULL);

   *quality = SCIP_INVALID;

   return SCIP_OKAY;
}

/**@} */




/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiGetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBase()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   ClpSimplex* clp = lpi->clp;

   if( rstat != NULL )
   {
      for( int i = 0; i < clp->numberRows(); ++i )
      {
	 switch ( clp->getRowStatus(i) )
	 {
	 case ClpSimplex::isFree:
            rstat[i] = SCIP_BASESTAT_ZERO;
            break;
	 case ClpSimplex::basic:
            rstat[i] = SCIP_BASESTAT_BASIC;
            break;
	 case ClpSimplex::atUpperBound:
            rstat[i] = SCIP_BASESTAT_UPPER;
            break;
	 case ClpSimplex::atLowerBound:
            rstat[i] = SCIP_BASESTAT_LOWER;
            break;
	 case ClpSimplex::superBasic:
            rstat[i] = SCIP_BASESTAT_ZERO;
            break;
	 case ClpSimplex::isFixed:
	    if (clp->getRowPrice()[i] > 0.0)
	       rstat[i] = SCIP_BASESTAT_LOWER;
	    else
	       rstat[i] = SCIP_BASESTAT_UPPER;
	    break;
	 default:
            SCIPerrorMessage("invalid basis status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
	 }
      }
   }

   if( cstat != NULL )
   {
#ifndef NDEBUG
      const double* lb = clp->getColLower();
      const double* ub = clp->getColUpper();
#endif

      for( int j = 0; j < clp->numberColumns(); ++j )
      {
	 switch ( clp->getColumnStatus(j) )
	 {
	 case ClpSimplex::isFree:
            cstat[j] = SCIP_BASESTAT_ZERO;
            break;
	 case ClpSimplex::basic:
            cstat[j] = SCIP_BASESTAT_BASIC;
            break;
	 case ClpSimplex::atUpperBound:
            cstat[j] = SCIP_BASESTAT_UPPER;
            assert( ub[j] < COIN_DBL_MAX );
            break;
	 case ClpSimplex::atLowerBound:
            cstat[j] = SCIP_BASESTAT_LOWER;
            assert( lb[j] > -COIN_DBL_MAX );
            break;
	 case ClpSimplex::superBasic:
            cstat[j] = SCIP_BASESTAT_ZERO;
            break;
	 case ClpSimplex::isFixed:
	    if (clp->getReducedCost()[j] > 0.0)
            {
	       cstat[j] = SCIP_BASESTAT_LOWER;
               assert( lb[j] > -COIN_DBL_MAX );
            }
	    else
            {
	       cstat[j] = SCIP_BASESTAT_UPPER;
               assert( ub[j] < COIN_DBL_MAX );
            }
	    break;
	 default: SCIPerrorMessage("invalid basis status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
	 }
      }
   }

   return SCIP_OKAY;
}


/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const int*            cstat,              /**< array with column basis status */
   const int*            rstat               /**< array with row basis status */
   )
{
   int ncols;
   int nrows;

   SCIPdebugMessage("calling SCIPlpiSetBase()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );

   assert(rstat != NULL || lpi->clp->numberRows() == 0);
   assert(cstat != NULL || lpi->clp->numberColumns() == 0);

   invalidateSolution(lpi);

   // Adapted from OsiClpSolverInterface::setBasisStatus

   ClpSimplex* clp = lpi->clp;
   clp->createStatus();

   const double* lhs = clp->getRowLower();
   const double* rhs = clp->getRowUpper();

   for( int i = 0; i < clp->numberRows(); ++i )
   {
      int status = rstat[i];
      assert( 0 <= status && status <= 3 );
      assert( lhs[i] > -COIN_DBL_MAX || status != SCIP_BASESTAT_LOWER); // can't be at lower bound
      assert( rhs[i] < COIN_DBL_MAX  || status != SCIP_BASESTAT_UPPER); // can't be at upper bound

      switch ( status )
      {
      case SCIP_BASESTAT_ZERO:
	 if ( lhs[i] <= -COIN_DBL_MAX && rhs[i] >= COIN_DBL_MAX )
	    clp->setRowStatus(i, ClpSimplex::isFree);
	 else
	    clp->setRowStatus(i, ClpSimplex::superBasic);
	 break;
      case SCIP_BASESTAT_BASIC:
         clp->setRowStatus(i, ClpSimplex::basic);
         break;
      case SCIP_BASESTAT_UPPER:
         clp->setRowStatus(i, ClpSimplex::atUpperBound);
         break;
      case SCIP_BASESTAT_LOWER:
	 if ( EPSEQ(rhs[i], lhs[i], 1e-6) )   // if bounds are equal
	    clp->setRowStatus(i, ClpSimplex::isFixed);
	 else
	    clp->setRowStatus(i, ClpSimplex::atLowerBound);
	 break;
      default:
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   const double* lb = clp->getColLower();
   const double* ub = clp->getColUpper();

   for( int j = 0; j < clp->numberColumns(); ++j )
   {
      int status = cstat[j];
      assert( 0 <= status && status <= 3 );
      assert( lb[j] > -COIN_DBL_MAX || status != SCIP_BASESTAT_LOWER); // can't be at lower bound
      assert( ub[j] < COIN_DBL_MAX  || status != SCIP_BASESTAT_UPPER); // can't be at upper bound

      switch ( status )
      {
      case SCIP_BASESTAT_ZERO:
	 if ( lb[j] <= -COIN_DBL_MAX && ub[j] >= COIN_DBL_MAX )
	    clp->setColumnStatus(j, ClpSimplex::isFree);
	 else
	    clp->setColumnStatus(j, ClpSimplex::superBasic);
	 break;
      case SCIP_BASESTAT_BASIC:
         clp->setColumnStatus(j, ClpSimplex::basic);
         break;
      case SCIP_BASESTAT_UPPER:
         clp->setColumnStatus(j, ClpSimplex::atUpperBound);
         break;
      case SCIP_BASESTAT_LOWER:
	 if ( EPSEQ(ub[j], lb[j], 1e-6) )
	    clp->setColumnStatus(j, ClpSimplex::isFixed);
	 else
	    clp->setColumnStatus(j, ClpSimplex::atLowerBound);
	 break;
      default:
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   /* Whats changed since last solve.
    *  Is only used when startFinishOptions used in dual or primal.
    * Bit 1 - number of rows/columns has not changed (so work arrays valid)
    *     2 - matrix has not changed
    *     4 - if matrix has changed only by adding rows
    *     8 - if matrix has changed only by adding columns
    *    16 - row lbs not changed
    *    32 - row ubs not changed
    *    64 - column objective not changed
    *   128 - column lbs not changed
    *   256 - column ubs not changed
    *	512 - basis not changed (up to user to set this to 0)
    *	      top bits may be used internally
    */
   clp->setWhatsChanged(clp->whatsChanged() & (~512));

   return SCIP_OKAY;
}


/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBasisInd()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(bind != 0);

   ClpSimplex* clp = lpi->clp;
   int nrows = clp->numberRows();
   int ncols = clp->numberColumns();

   int* idx = NULL;
   SCIP_ALLOC( BMSallocMemoryArray(&idx, nrows) );

   /* If secondaryStatus == 6, clp says the LP is empty. Mose likely this happened, because the
    * matrix is empty, i.e., all rows were redundant/empty. In this case, we construct a basis
    * consisting of slack variables. */
   if ( clp->secondaryStatus() == 6 )
   {
      assert( clp->getNumElements() == 0 );
      for (int i = 0; i < nrows; ++i)
	 idx[i] = ncols + i;
   }
   else
      clp->getBasics(idx);

   for (int i = 0; i < nrows; ++i)
   {
      if ( idx[i] < ncols )
         bind[i] = idx[i];
      else
         bind[i] = -1 - (idx[i] - ncols);
   }

   BMSfreeMemoryArray(&idx);

   return SCIP_OKAY;
}


/** get row of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 *
 *  @todo check that the result is in terms of the LP interface definition
 */
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvRow()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(coef != NULL);
   assert( 0 <= r && r <= lpi->clp->numberRows() );

   /* can only return dense result */
   if ( ninds != NULL )
      *ninds = -1;

   ClpSimplex* clp = lpi->clp;
   clp->getBInvRow(r, coef);

   return SCIP_OKAY;
}


/** get column of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 *
 *  @todo check that the result is in terms of the LP interface definition
 */
SCIP_RETCODE SCIPlpiGetBInvCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvCol()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(coef != NULL);
   assert( 0 <= c && c <= lpi->clp->numberRows() ); /* basis matrix is nrows * nrows */

   /* can only return dense result */
   if ( ninds != NULL )
      *ninds = -1;

   ClpSimplex* clp = lpi->clp;
   clp->getBInvCol(c, coef);

   return SCIP_OKAY;
}

/** get row of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 *
 *  @todo check that the result is in terms of the LP interface definition
 */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            coef,               /**< vector to return coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvARow()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(coef != NULL);
   assert( 0 <= r && r <= lpi->clp->numberRows() );

   /* can only return dense result */
   if ( ninds != NULL )
      *ninds = -1;

   ClpSimplex* clp = lpi->clp;
   clp->getBInvARow(r, coef, 0);

   return SCIP_OKAY;
}

/** get column of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 *
 *  @todo check that the result is in terms of the LP interface definition
 */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef,               /**< vector to return coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvACol()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert( coef != 0 );
   assert( 0 <= c && c <= lpi->clp->numberColumns() );

   /* can only return dense result */
   if ( ninds != NULL )
      *ninds = -1;

   ClpSimplex* clp = lpi->clp;
   clp->getBInvACol(c, coef);

   return SCIP_OKAY;
}


/**@} */




/*
 * LP State Methods
 */

/**@name LP State Methods */
/**@{ */

/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetState()\n");

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(lpistate != NULL);

   int ncols = lpi->clp->numberColumns();
   int nrows = lpi->clp->numberRows();
   assert(ncols >= 0);
   assert(nrows >= 0);

   /* allocate lpistate data */
   SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows) );

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* get unpacked basis information */
   SCIP_CALL( SCIPlpiGetBase(lpi, lpi->cstat, lpi->rstat) );

   /* pack LPi state data */
   (*lpistate)->ncols = ncols;
   (*lpistate)->nrows = nrows;
   lpistatePack(*lpistate, lpi->cstat, lpi->rstat);

   return SCIP_OKAY;
}


/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetState()
 */
SCIP_RETCODE SCIPlpiSetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const SCIP_LPISTATE*  lpistate            /**< LPi state information (like basis information), or NULL */
   )
{  /* lint --e{715} */
   int lpncols;
   int lpnrows;
   int i;

   SCIPdebugMessage("calling SCIPlpiSetState()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(blkmem != NULL);

   /* if there was no basis information available, the LPI state was not stored */
   if( lpistate == NULL )
      return SCIP_OKAY;

   lpncols = lpi->clp->numberColumns();
   lpnrows = lpi->clp->numberRows();
   assert(lpistate->ncols <= lpncols);
   assert(lpistate->nrows <= lpnrows);

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, lpncols) );
   SCIP_CALL( ensureRstatMem(lpi, lpnrows) );

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->cstat, lpi->rstat);

   /* extend the basis to the current LP beyond the previously existing columns */
   for( i = lpistate->ncols; i < lpncols; ++i )
   {
      SCIP_Real bnd = (lpi->clp->getColLower())[i];
      if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         bnd = (lpi->clp->getColUpper())[i];
         if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
            lpi->cstat[i] = SCIP_BASESTAT_ZERO;  /* variable is free */
         else
            lpi->cstat[i] = SCIP_BASESTAT_UPPER; /* use finite upper bound */
      }
      else
         lpi->cstat[i] = SCIP_BASESTAT_LOWER;    /* use finite lower bound */
   }
   for( i = lpistate->nrows; i < lpnrows; ++i )
      lpi->rstat[i] = SCIP_BASESTAT_BASIC;

   /* load basis information */
   SCIP_CALL( SCIPlpiSetBase(lpi, lpi->cstat, lpi->rstat) );

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiClearState()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   lpi->clp->allSlackBasis(true);
   lpi->validFactorization = false;

   return SCIP_OKAY;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   SCIPdebugMessage("calling SCIPlpiFreeState()\n");

   assert(lpi != NULL);
   assert(lpistate != NULL);
   assert(blkmem != NULL);

   if ( *lpistate != NULL )
      lpistateFree(lpistate, blkmem);

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information), or NULL*/
   )
{
   assert(lpi != NULL);
   return (lpistate != NULL);
}

/** reads LP state (like basis information) from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadState()\n");
   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(fname != NULL);

   /*  Read a basis from the given filename,
    *  returns -1 on file error, 0 if no values, 1 if values
    */
   if ( lpi->clp->readBasis(fname) < 0 )
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** writes LPi state (i.e. basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiWriteState()\n");
   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(fname != NULL);

   /*  Write the basis in MPS format to the specified file.
    *  If writeValues true, writes values of structurals
    *  (and adds VALUES to end of NAME card)
    *
    *  parameters:
    *  - filename
    *  - bool writeValues
    *  - int formatType  (0 - normal, 1 - extra accuracy, 2 - IEEE hex)
    */
   if ( lpi->clp->writeBasis(fname, false, 0) )
      return SCIP_WRITEERROR;

   return SCIP_OKAY;
}

/**@} */




/*
 * LP Pricing Norms Methods
 */

/**@name LP Pricing Norms Methods */
/**@{ */

/** stores LPi pricing norms information
 *  @todo should we store norm information?
 */
SCIP_RETCODE SCIPlpiGetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information */
   )
{
   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpinorms != NULL);

   (*lpinorms) = NULL;

   return SCIP_OKAY;
}

/** loads LPi pricing norms into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetNorms()
 */
SCIP_RETCODE SCIPlpiSetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const SCIP_LPINORMS*  lpinorms            /**< LPi pricing norms information, or NULL */
   )
{
   assert(lpinorms == NULL);

   /* no work necessary */
   return SCIP_OKAY;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPlpiFreeNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information, or NULL */
   )
{
   assert(lpinorms == NULL);

   /* no work necessary */
   return SCIP_OKAY;
}

/**@} */




/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiGetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetIntpar()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(ival != 0);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      *ival = lpi->startscratch;
      break;
   case SCIP_LPPAR_SCALING:
      if( lpi->clp->scalingFlag() != 0 )  // 0 -off, 1 equilibrium, 2 geometric, 3 auto, 4 dynamic(later)
         *ival = TRUE;
      else
         *ival = FALSE;
      break;
   case SCIP_LPPAR_PRICING:
      *ival = (int)lpi->pricing;          // store pricing method in LPI struct
      break;
   case SCIP_LPPAR_LPINFO:
      *ival = lpi->clp->logLevel() > 0 ? TRUE : FALSE;
      break;
   case SCIP_LPPAR_LPITLIM:
      *ival = lpi->clp->maximumIterations();
      break;
   case SCIP_LPPAR_FASTMIP:
      *ival = lpi->fastmip;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}


/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiSetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetIntpar()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   // Handle pricing separately ...
   if( type == SCIP_LPPAR_PRICING )
   {
      // for primal:
      // 0 is exact devex,
      // 1 full steepest,
      // 2 is partial exact devex
      // 3 switches between 0 and 2 depending on factorization
      // 4 starts as partial dantzig/devex but then may switch between 0 and 2.
      // - currently (Clp 1.8) default is 3

      // for dual:
      // 0 is uninitialized,
      // 1 full,
      // 2 is partial uninitialized,
      // 3 starts as 2 but may switch to 1.
      // - currently (Clp 1.8) default is 3
      lpi->pricing = (SCIP_PRICING)ival;
      int primalmode = 0;
      int dualmode = 0;
      switch( (SCIP_PRICING)ival )
      {
      case SCIP_PRICING_AUTO:
         primalmode = 3; dualmode = 3; break;
      case SCIP_PRICING_FULL:
         primalmode = 0; dualmode = 1; break;
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_STEEP:
         primalmode = 1; dualmode = 0; break;
      case SCIP_PRICING_STEEPQSTART:
         primalmode = 1; dualmode = 2; break;
      case SCIP_PRICING_DEVEX:
         primalmode = 2; dualmode = 3; break;
      default:
         SCIPerrorMessage("unkown pricing parameter %d!\n", ival);
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
      ClpPrimalColumnSteepest primalpivot(primalmode);
      lpi->clp->setPrimalColumnPivotAlgorithm(primalpivot);
      ClpDualRowSteepest dualpivot(dualmode);
      lpi->clp->setDualRowPivotAlgorithm(dualpivot);
      return SCIP_OKAY;
   }

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      lpi->startscratch = ival;
      break;
   case SCIP_LPPAR_SCALING:
      lpi->clp->scaling((ival > 0) ? 3 : 0);    // 0 -off, 1 equilibrium, 2 geometric, 3, auto, 4 dynamic(later));
      break;
   case SCIP_LPPAR_PRICING:
      /* should not happen - see above */
      SCIPABORT();
      return SCIP_LPERROR; /*lint !e527*/
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      /*  Amount of print out:
       *  0 - none
       *  1 - just final
       *  2 - just factorizations
       *  3 - as 2 plus a bit more
       *  4 - verbose
       *  above that 8,16,32 etc just for selective SCIPdebug
       */
      if ( ival )
	 lpi->clp->setLogLevel(2);      // lpi->clp->setLogLevel(63);
      else
         lpi->clp->setLogLevel(0);
      break;
   case SCIP_LPPAR_LPITLIM:
      /* ival >= 0, 0 stop immediately */
      assert( ival >= 0 );
      lpi->clp->setMaximumIterations(ival);
      break;
   case SCIP_LPPAR_FASTMIP:
      assert(ival == TRUE || ival == FALSE);
      if( ival )
         setFastmipClpParameters(lpi);
      else
         unsetFastmipClpParameters(lpi);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}


/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiGetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetRealpar()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(dval != 0);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      *dval = lpi->clp->primalTolerance();
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      *dval = lpi->clp->dualTolerance();
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      /* @todo add BARRIERCONVTOL parameter */
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_OBJLIM:
      *dval = lpi->clp->dualObjectiveLimit();
      break;
   case SCIP_LPPAR_LPTILIM:
      *dval = lpi->clp->maximumSeconds();
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiSetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetRealpar()\n");
   SCIPdebugMessage("setting parameter %d to value %g.\n", type, dval);

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      assert( dval > 0.0 );
      /* 0 < dval < 1e10 */
      if( dval > 1e+10 )
      {
         /* however dval is required to be strictly less than 1e+10 */
         dval = 9e+9;
      }

      lpi->clp->setPrimalTolerance(dval);
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      assert( dval > 0.0 );
      /* 0 < dval < 1e10 */
      if( dval > 1e+10 )
      {
         /* however dval is required to be strictly less than 1e+10 */
         dval = 9e+9;
      }

      lpi->clp->setDualTolerance(dval);
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      /* @todo add BARRIERCONVTOL parameter */
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_OBJLIM:
      /* no restriction on dval */

      lpi->clp->setDualObjectiveLimit(dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      assert( dval > 0.0 );
      /* clp poses no restrictions on dval
       * (it handles the case dval < 0 internally and sets param to -1 meaning no time limit.)
       *
       * However for consistency we assert the timelimit to be strictly positive.
       */

      lpi->clp->setMaximumSeconds(dval);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** interrupts the currently ongoing lp solve or disables the interrupt */
SCIP_RETCODE SCIPlpiInterrupt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             interrupt           /**< TRUE if interrupt should be set, FALSE if it should be disabled */
   )
{
   /*lint --e{715}*/
   assert(lpi != NULL);

   return SCIP_OKAY;
}

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /* lint --e{715} */
   assert(lpi != NULL);
   SCIPdebugMessage("calling SCIPlpiInfinity()\n");

   return COIN_DBL_MAX;
}


/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to check */
   )
{  /* lint --e{715} */
   assert(lpi != NULL);
   SCIPdebugMessage("calling SCIPlpiIsInfinity()\n");

   return (val >= COIN_DBL_MAX);
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** returns, whether the given file exists */
static
SCIP_Bool fileExists(
   const char*           filename            /**< file name */
   )
{
   FILE* f;

   f = fopen(filename, "r");
   if( f == 0 )
      return FALSE;

   fclose(f);

   return TRUE;
}

/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadLP()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(fname != NULL);

   // WARNING: can only read mps files

   if ( !fileExists(fname) )
      return SCIP_NOFILE;

   /* read file in MPS format
    * parameters:
    * filename
    * bool keepNames
    * bool ignoreErrors
    */
   if ( lpi->clp->readMps(fname, true, false) )
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiWriteLP() - %s\n", fname);

   assert(lpi != NULL);
   assert(lpi->clp != NULL);
   assert(fname != NULL);

   /*  write file in MPS format
    *  parameters:
    *  filename
    *  int formatType  (0 - normal, 1 - extra accuracy, 2 - IEEE hex)
    *  int numberAcross (1 or 2 values should be specified on every data line in the MPS file)
    *  double objSense
    */
   if ( lpi->clp->writeMps(fname, 0, 2, lpi->clp->optimizationDirection()) )
      return SCIP_WRITEERROR;

   return SCIP_OKAY;
}

/**@} */
