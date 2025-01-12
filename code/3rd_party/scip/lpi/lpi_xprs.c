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

/**@file   lpi_xprs.c
 * @ingroup LPIS
 * @brief  LP interface for Xpress-MP
 * @author Tobias Achterberg
 * @author Michael Perregaard
 * @author Livio Bertacco
 * @author Stefan Heinz
 *
 * This interface was revised for Xpress 26. Therefore, we removed all legacy code.
 *
 * Xpress requires that column and row names are unique. Since column and row names are not needed we ignore all column
 * and row names to avoid the uniqueness issue.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif

#include "xprs.h"
#include "scip/bitencode.h"
#include "scip/pub_misc.h"
#include "scip/pub_message.h"
#include "lpi/lpi.h"
#include "tinycthread/tinycthread.h"

#ifndef XPRS_LPQUICKPRESOLVE
#define XPRS_LPQUICKPRESOLVE 8207
#endif

/* For SCIP we need an extra LP status which is optimal with scaled infeasibilities. */
#define XPRS_LP_OPTIMAL_SCALEDINFEAS 16

#define CHECK_ZERO(messagehdlr, x) { int _restat_;                      \
      if( (_restat_ = (x)) != 0 )                                       \
      {                                                                 \
         SCIPmessagePrintWarning((messagehdlr), "%s:%d: LP Error: Xpress returned %d\n", __FILE__, __LINE__, _restat_); \
         return SCIP_LPERROR;                                           \
      }                                                                 \
   }

/* this macro is only called in functions returning SCIP_Bool; thus, we return retval if there is an error in optimized mode */
#define ABORT_ZERO(messagehdlr, retval, x) { int _restat_;              \
      if( (_restat_ = (x)) != 0 )                                       \
      {                                                                 \
         SCIPmessagePrintWarning((messagehdlr), "LP Error: Xpress returned %d\n", _restat_); \
         SCIPABORT();                                                   \
         return retval;                                                 \
      }                                                                 \
   }


typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE

/** LP interface */
struct SCIP_LPi
{
   XPRSprob              xprslp;             /**< Xpress LP pointer */
   char                  name[200];          /**< problem name */

   SCIP_PRICING          pricing;            /**< SCIP pricing setting  */
   int                   notfromscratch;     /**< do we not want to solve the lp from scratch */
   int                   solstat;            /**< solution status of last optimization call */
   char                  solmethod;          /**< method used to solve the LP */

   char*                 larray;             /**< array with 'L' entries for changing lower bounds */
   char*                 uarray;             /**< array with 'U' entries for changing upper bounds */
   char*                 senarray;           /**< array for storing row senses */
   SCIP_Real*            rhsarray;           /**< array for storing rhs values */
   SCIP_Real*            rngarray;           /**< array for storing range values */
   SCIP_Real*            valarray;           /**< array for storing coefficient values */
   int*                  cstat;              /**< array for storing column basis status */
   int*                  rstat;              /**< array for storing row basis status (row status w.r.t. slack columns) */
   int*                  indarray;           /**< array for storing coefficient indices */

   int                   boundchgsize;       /**< size of larray and uarray */
   int                   sidechgsize;        /**< size of senarray and rngarray */
   int                   valsize;            /**< size of valarray and indarray */
   int                   cstatsize;          /**< size of cstat array */
   int                   rstatsize;          /**< size of rstat array */

   int                   iterations;         /**< number of iterations used in the last solving call */
   SCIP_Bool             solisbasic;         /**< is current LP solution a basic solution? */
   SCIP_Bool             clearstate;         /**< should the current basis be ignored with the next LP solve */

   SCIP_Real             par_lobjlim;        /**< objective lower bound */
   SCIP_Real             par_uobjlim;        /**< objective upper bound */
   int                   par_fastlp;         /**< special meta parameter for making LP reoptimize go faster */
   int                   par_presolve;       /**< need to distinguish between the users setting and the optimizer setting of presolve */

   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehdlr handler to printing messages, or NULL */
};

/** LPi state stores basis information */
struct SCIP_LPiState
{
   int                   ncols;              /**< number of LP columns */
   int                   nrows;              /**< number of LP rows */
   COLPACKET*            packcstat;          /**< column basis status in compressed form */
   ROWPACKET*            packrstat;          /**< row basis status in compressed form (row status w.r.t. slack columns) */
};

/**@name Debug check methods
 *
 * @{
 */

#ifndef NDEBUG

/** check that the column range fits */
static
void debugCheckColrang(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   int ncols;

   (void)XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < ncols);
}

/** check that the row range fits */
static
void debugCheckRowrang(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   int nrows;

   (void)XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < nrows);
}

#else

/* in optimized mode the checks are replaced with an empty command */
#define debugCheckColrang(lpi, firstcol, lastcol) /* */
#define debugCheckRowrang(lpi, firstrow, lastrow) /* */
#endif

/**@} */


/**@name Dynamic memory arrays
 *
 * @{
 */

/** resizes larray and uarray to have at least num entries and fill it with 'L' and 'U' for the lower and upper bound
 *  markers
 */
static
SCIP_RETCODE ensureBoundchgMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->boundchgsize )
   {
      int newsize;
      int i;

      newsize = MAX(2*lpi->boundchgsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->larray, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->uarray, newsize) );
      for( i = lpi->boundchgsize; i < newsize; ++i )
      {
         lpi->larray[i] = 'L';
         lpi->uarray[i] = 'U';
      }
      lpi->boundchgsize = newsize;
   }
   assert(num <= lpi->boundchgsize);

   return SCIP_OKAY;
}

/** resizes senarray, rngarray, and rhsarray to have at least num entries */
static
SCIP_RETCODE ensureSidechgMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->sidechgsize )
   {
      int newsize;

      newsize = MAX(2*lpi->sidechgsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->senarray, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rhsarray, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rngarray, newsize) );
      lpi->sidechgsize = newsize;
   }
   assert(num <= lpi->sidechgsize);

   return SCIP_OKAY;
}

/** resizes valarray and indarray to have at least num entries */
static
SCIP_RETCODE ensureValMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->valsize )
   {
      int newsize;

      newsize = MAX(2*lpi->valsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->valarray, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->indarray, newsize) );
      lpi->valsize = newsize;
   }
   assert(num <= lpi->valsize);

   return SCIP_OKAY;
}

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

/**@} */


/**@name LPi state methods
 *
 * @{
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
   const int*            rstat               /**< basis status of rows in unpacked format (row status w.r.t. slack columns) */
   )
{
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

   SCIPencodeDualBit(cstat, lpistate->packcstat, lpistate->ncols);
   SCIPencodeDualBit(rstat, lpistate->packrstat, lpistate->nrows);
}

/** unpacks row and column basis status from a packed LPi state object */
static
void lpistateUnpack(
   const SCIP_LPISTATE*  lpistate,           /**< pointer to LPi state data */
   int*                  cstat,              /**< buffer for storing basis status of columns in unpacked format */
   int*                  rstat               /**< buffer for storing basis status of rows in unpacked format (row status w.r.t. slack columns) */
   )
{
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

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
   assert(lpistate != NULL);
   assert(blkmem != NULL);
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
   assert(blkmem != NULL);
   assert(lpistate != NULL);
   assert(*lpistate != NULL);

   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packcstat, colpacketNum((*lpistate)->ncols));
   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packrstat, rowpacketNum((*lpistate)->nrows));
   BMSfreeBlockMemory(blkmem, lpistate);
}

/**@} */


/**@name Conversion methods
 *
 * @{
 */

/** converts SCIP's objective sense into CPLEX's objective sense */
static
int xprsObjsen(
   SCIP_OBJSEN const     objsen              /**< objective sense */
   )
{
   switch( objsen )
   {
   case SCIP_OBJSEN_MAXIMIZE:
      return XPRS_OBJ_MAXIMIZE;
   case SCIP_OBJSEN_MINIMIZE:
      return XPRS_OBJ_MINIMIZE;
   default:
      SCIPerrorMessage("invalid objective sense\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }
}

/** converts SCIP's lhs/rhs pairs into Xpress' sen/rhs/rng */
static
void convertSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhss,               /**< left hand side vector */
   const SCIP_Real*      rhss                /**< right hand side vector */
   )
{
   int i;

   assert(lpi != NULL);
   assert(nrows >= 0);
   assert(lhss != NULL);
   assert(rhss != NULL);

   /* convert lhs/rhs into sen/rhs/rng */
   for( i = 0; i < nrows; ++i )
   {
      assert(lhss[i] <= rhss[i]);
      if( lhss[i] == rhss[i] ) /*lint !e777*/
      {
         assert(XPRS_MINUSINFINITY < rhss[i] && rhss[i] < XPRS_PLUSINFINITY);
         lpi->senarray[i] = 'E';
         lpi->rhsarray[i] = rhss[i];
         lpi->rngarray[i] = 0.0;
      }
      else if( lhss[i] <= XPRS_MINUSINFINITY )
      {
         lpi->senarray[i] = 'L';
         lpi->rhsarray[i] = rhss[i];
         lpi->rngarray[i] = XPRS_PLUSINFINITY;
      }
      else if( rhss[i] >= XPRS_PLUSINFINITY )
      {
         lpi->senarray[i] = 'G';
         lpi->rhsarray[i] = lhss[i];
         lpi->rngarray[i] = XPRS_PLUSINFINITY;
      }
      else
      {
         /* Xpress defines a ranged row to be within rhs-rng and rhs. */
         lpi->senarray[i] = 'R';
         lpi->rhsarray[i] = rhss[i];
         lpi->rngarray[i] = rhss[i] - lhss[i];
      }
   }
}

/** converts Xpress' sen/rhs/rng triplets into SCIP's lhs/rhs pairs */
static
void reconvertBothSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   SCIP_Real*            lhss,               /**< buffer to store the left hand side vector */
   SCIP_Real*            rhss                /**< buffer to store the right hand side vector */
   )
{
   int i;

   assert(lpi != NULL);
   assert(nrows >= 0);
   assert(lhss != NULL);
   assert(rhss != NULL);

   for( i = 0; i < nrows; ++i )
   {
      switch( lpi->senarray[i] )
      {
      case 'E':
         lhss[i] = lpi->rhsarray[i];
         rhss[i] = lpi->rhsarray[i];
         break;

      case 'L':
         lhss[i] = XPRS_MINUSINFINITY;
         rhss[i] = lpi->rhsarray[i];
         break;

      case 'G':
         lhss[i] = lpi->rhsarray[i];
         rhss[i] = XPRS_PLUSINFINITY;
         break;

      case 'R':
         assert(lpi->rngarray[i] >= 0.0);
         rhss[i] = lpi->rhsarray[i];
         lhss[i] = lpi->rhsarray[i] - lpi->rngarray[i];
         break;

      default:
         SCIPerrorMessage("invalid row sense\n");
         SCIPABORT();
      }
      assert(lhss[i] <= rhss[i]);
   }
}

/** converts Xpress' sen/rhs/rng triplets into SCIP's lhs/rhs pairs, only storing the left hand side */
static
void reconvertLhs(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   SCIP_Real*            lhss                /**< buffer to store the left hand side vector */
   )
{
   int i;

   assert(lpi != NULL);
   assert(nrows >= 0);
   assert(lhss != NULL);

   for( i = 0; i < nrows; ++i )
   {
      switch( lpi->senarray[i] )
      {
      case 'E':
         assert(lpi->rngarray[i] == 0.0);
         lhss[i] = lpi->rhsarray[i];
         break;

      case 'L':
         assert(lpi->rngarray[i] == 0.0);
         lhss[i] = XPRS_MINUSINFINITY;
         break;

      case 'G':
         assert(lpi->rngarray[i] == 0.0);
         lhss[i] = lpi->rhsarray[i];
         break;

      case 'R':
         assert(lpi->rngarray[i] >= 0.0);
         lhss[i] = lpi->rhsarray[i] - lpi->rngarray[i];
         break;

      default:
         SCIPerrorMessage("invalid row sense\n");
         SCIPABORT();
      }
   }
}

/** converts Xpress' sen/rhs/rng triplets into SCIP's lhs/rhs pairs, only storing the right hand side */
static
void reconvertRhs(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   SCIP_Real*            rhss                /**< buffer to store the right hand side vector */
   )
{
   int i;

   assert(lpi != NULL);
   assert(nrows >= 0);
   assert(rhss != NULL);

   for( i = 0; i < nrows; ++i )
   {
      switch( lpi->senarray[i] )
      {
      case 'E':
         assert(lpi->rngarray[i] == 0.0);
         rhss[i] = lpi->rhsarray[i];
         break;

      case 'L':
         assert(lpi->rngarray[i] == 0.0);
         rhss[i] = lpi->rhsarray[i];
         break;

      case 'G':
         assert(lpi->rngarray[i] == 0.0);
         rhss[i] = XPRS_PLUSINFINITY;
         break;

      case 'R':
         assert(lpi->rngarray[i] >= 0.0);
         rhss[i] = lpi->rhsarray[i];
         break;

      default:
         SCIPerrorMessage("invalid row sense\n");
         SCIPABORT();
      }
   }
}

/** converts Xpress' sen/rhs/rng triplets into SCIP's lhs/rhs pairs */
static
void reconvertSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   SCIP_Real*            lhs,                /**< buffer to store the left hand side vector, or NULL */
   SCIP_Real*            rhs                 /**< buffer to store the right hand side vector, or NULL */
   )
{
   if( lhs != NULL && rhs != NULL )
      reconvertBothSides(lpi, nrows, lhs, rhs);
   else if( lhs != NULL )
      reconvertLhs(lpi, nrows, lhs);
   else if( rhs != NULL )
      reconvertRhs(lpi, nrows, rhs);
}

/**@} */


/** marks the current LP to be unsolved */
static
void invalidateSolution(
   SCIP_LPI*             lpi
   )
{
   assert(lpi != NULL);
   lpi->solstat = -1;
}

/*
 * LP Interface Methods
 */

/**@name Miscellaneous Methods
 *
 * @{
 */

#ifdef _Thread_local
static _Thread_local char xprsname[100];
#else
static char xprsname[] = {'X', 'p', 'r', 'e', 's', 's', ' ', '0' + XPVERSION / 10, '0' + XPVERSION % 10};
#endif
/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
#ifdef _Thread_local
   char version[16];

   /* get version of Xpress */
   if( XPRSgetversion(version) == 0 )
      (void) sprintf(xprsname, "Xpress %s", version);
   else
      (void) sprintf(xprsname, "Xpress %d", XPVERSION);
#endif
   return xprsname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "Linear Programming Solver developed by FICO (www.fico.com/xpress)";
}

/** gets pointer for LP solver - use only with great care
 *
 *  Here we return the pointer to the LP environment.
 */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{ /*lint --e{715}*/
   return (void*) lpi->xprslp;
}

/** pass integrality information to LP solver */
SCIP_RETCODE SCIPlpiSetIntegralityInformation(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   ncols,              /**< length of integrality array */
   int*                  intInfo             /**< integrality array (0: continuous, 1: integer). May be NULL iff ncols is 0.  */
   )
{ /*lint --e{715}*/
   assert(lpi != NULL);
   assert(ncols >= 0);
   assert(ncols == 0 || intInfo != NULL);

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


/**@name LPI Creation and Destruction Methods
 *
 * @{
 */

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(
   SCIP_LPI**            lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   int zero = 0;

   assert(sizeof(SCIP_Real) == sizeof(double)); /*lint !e506*/ /* Xpress only works with doubles as floating points */
   assert(sizeof(SCIP_Bool) == sizeof(int));    /*lint !e506*/ /* Xpress only works with ints as bools */
   assert(lpi != NULL);
   assert(name != NULL);

   SCIPdebugMessage("SCIPlpiCreate()\n");

   /* the interface is revised for Xpress 26 or higher */
   if( XPVERSION < 26 ) /*lint !e506 !e774*/
   {
      SCIPmessagePrintWarning(messagehdlr, "Please use Xpress version 26 or higher, you are using %d\n", XPVERSION);
      return SCIP_LPERROR;
   }

   /* initialize the Xpress library (licensing) */
   CHECK_ZERO( messagehdlr, XPRSinit(NULL) );

   /* create LPi data structure */
   SCIP_ALLOC( BMSallocMemory(lpi) );

   /* copy the problem name */
   (void)strncpy((*lpi)->name, name, 199);
   (*lpi)->name[199] = '\0';

   (*lpi)->larray = NULL;
   (*lpi)->uarray = NULL;
   (*lpi)->senarray = NULL;
   (*lpi)->rhsarray = NULL;
   (*lpi)->rngarray = NULL;
   (*lpi)->indarray = NULL;
   (*lpi)->valarray = NULL;
   (*lpi)->cstat = NULL;
   (*lpi)->rstat = NULL;
   (*lpi)->boundchgsize = 0;
   (*lpi)->sidechgsize = 0;
   (*lpi)->valsize = 0;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->iterations = 0;
   (*lpi)->solisbasic = TRUE;
   (*lpi)->clearstate = FALSE;
   (*lpi)->solmethod = ' ';
   (*lpi)->par_lobjlim = -1e+40;
   (*lpi)->par_uobjlim = +1e+40;
   (*lpi)->par_fastlp = 0;
   (*lpi)->par_presolve = 0;
   (*lpi)->messagehdlr = messagehdlr;

   CHECK_ZERO( messagehdlr, XPRScreateprob(&(*lpi)->xprslp) );
   invalidateSolution(*lpi);

   /* turn logging off until the user explicitly turns it on; this should prevent any unwanted Xpress output from
    * appearing in the SCIP log.
    */
   CHECK_ZERO( messagehdlr, XPRSsetintcontrol((*lpi)->xprslp, XPRS_OUTPUTLOG, 0) );

   /* we need to create an empty LP in this prob since SCIP might attempt to add rows or columns to it */
   CHECK_ZERO( messagehdlr, XPRSloadlp((*lpi)->xprslp, (*lpi)->name, 0, 0, NULL, NULL, NULL, NULL, &zero, NULL, NULL, NULL, NULL, NULL) );

   /* set objective sense */
   SCIP_CALL( SCIPlpiChgObjsen(*lpi, objsen) );

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(*lpi != NULL);
   assert((*lpi)->xprslp != NULL);

   SCIPdebugMessage("SCIPlpiFree()\n");

   /* free LP */
   CHECK_ZERO( (*lpi)->messagehdlr, XPRSdestroyprob(((*lpi)->xprslp)) );

   /* free environment */
   CHECK_ZERO( (*lpi)->messagehdlr, XPRSfree() );

   /* free memory */
   BMSfreeMemoryArrayNull(&(*lpi)->larray);
   BMSfreeMemoryArrayNull(&(*lpi)->uarray);
   BMSfreeMemoryArrayNull(&(*lpi)->senarray);
   BMSfreeMemoryArrayNull(&(*lpi)->rhsarray);
   BMSfreeMemoryArrayNull(&(*lpi)->rngarray);
   BMSfreeMemoryArrayNull(&(*lpi)->indarray);
   BMSfreeMemoryArrayNull(&(*lpi)->valarray);
   BMSfreeMemoryArrayNull(&(*lpi)->cstat);
   BMSfreeMemoryArrayNull(&(*lpi)->rstat);
   BMSfreeMemory(lpi);

   return SCIP_OKAY;
}

/**@} */


/**@name Modification Methods
 *
 * @{
 */

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
{ /*lint --e{715}*/
   int c;

#ifndef NDEBUG
   {
      int j;
      for( j = 0; j < nnonz; j++ )
         assert( val[j] != 0 );
   }
#endif

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(beg != NULL);
   assert(ind != NULL);
   assert(val != NULL);
   SCIP_UNUSED(colnames);
   SCIP_UNUSED(rownames);

   SCIPdebugMessage("loading LP in column format into Xpress: %d cols, %d rows\n", ncols, nrows);

   invalidateSolution(lpi);

   /* ensure that the temporary arrays for the side conversion are long enough */
   SCIP_CALL( ensureSidechgMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples the sen/rhs/range are stored in the temporary arrays in lpi structure */
   convertSides(lpi, nrows, lhs, rhs);

   /* ensure that the temporary arrays are large enough */
   SCIP_CALL( ensureValMem(lpi, ncols) );

   /* calculate column lengths */
   for( c = 0; c < ncols-1; ++c )
   {
      lpi->indarray[c] = beg[c+1] - beg[c];
      assert(lpi->indarray[c] >= 0);
   }
   lpi->indarray[ncols-1] = nnonz - beg[ncols-1];
   assert(lpi->indarray[ncols-1] >= 0);

   /* copy data into Xpress */
   CHECK_ZERO( lpi->messagehdlr, XPRSloadlp(lpi->xprslp, lpi->name, ncols, nrows, lpi->senarray, lpi->rhsarray,
         lpi->rngarray, obj, beg, lpi->indarray, ind, val, lb, ub) );

   /* set objective sense */
   SCIP_CALL( SCIPlpiChgObjsen(lpi, objsen) );

   return SCIP_OKAY;
}

/** adds columns to the LP */
SCIP_RETCODE SCIPlpiAddCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const SCIP_Real*      obj,                /**< objective function values of new columns */
   const SCIP_Real*      lb,                 /**< lower bounds of new columns */
   const SCIP_Real*      ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{ /*lint --e{715}*/
   int c;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(ncols > 0);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(nnonz >= 0);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);
   SCIP_UNUSED(colnames);

   SCIPdebugMessage("adding %d columns with %d nonzeros to Xpress\n", ncols, nnonz);

   invalidateSolution(lpi);

   /* ensure that the temporary arrays are large enough */
   SCIP_CALL( ensureValMem(lpi, ncols+1) );

#ifndef NDEBUG
   {
      /* perform check that no new rows are added - this is forbidden */
      int nrows;
      int j;

      CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
      for (j = 0; j < nnonz; ++j)
      {
         assert( val[j] != 0.0 );
         assert( 0 <= ind[j] && ind[j] < nrows );
      }
   }
#endif

   /* only collect the start array if we have at least one non-zero */
   if( nnonz > 0 )
   {
      /* we need ncol+1 entries in the start array for Xpress */
      for( c = 0; c < ncols; c++ )
         lpi->indarray[c] = beg[c];
      lpi->indarray[ncols] = nnonz;
   }

   /* add the columns with (potential) non-zeros to the Xpress */
   CHECK_ZERO( lpi->messagehdlr, XPRSaddcols(lpi->xprslp, ncols, nnonz, obj, lpi->indarray, ind, val, lb, ub) );

   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   int c;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   debugCheckColrang(lpi, firstcol, lastcol);

   SCIPdebugMessage("deleting %d columns from Xpress\n", lastcol - firstcol + 1);

   invalidateSolution(lpi);

   /* ensure that the temporary arrays are large enough */
   SCIP_CALL( ensureValMem(lpi, lastcol-firstcol+1) );

   /* collect the columns indices to be deleted */
   for( c = firstcol; c <= lastcol; c++ )
      lpi->indarray[c-firstcol] = c;

   CHECK_ZERO( lpi->messagehdlr, XPRSdelcols(lpi->xprslp, lastcol-firstcol+1, lpi->indarray) );

   return SCIP_OKAY;
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   int nkeptcols;
   int ndelcols;
   int ncols;
   int c;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(dstat != NULL);

   SCIPdebugMessage("deleting a column set from Xpress\n");

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );

   nkeptcols = 0;
   ndelcols = 0;

   /* ensure that the temporary arrays are large enough */
   SCIP_CALL( ensureValMem(lpi, ncols) );

   /* collect the column indecies which should be deleted and create a the new column ordering */
   for( c = 0; c < ncols; c++ )
   {
      if( dstat[c] == 1 )
      {
         dstat[c] = -1;
         lpi->indarray[ndelcols] = c;
         ndelcols++;
      }
      else
      {
         dstat[c] = nkeptcols;
         nkeptcols++;
      }
   }

   CHECK_ZERO( lpi->messagehdlr, XPRSdelcols(lpi->xprslp, ndelcols, lpi->indarray) );

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
{ /*lint --e{715}*/
   int r;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(nrows >= 0);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(nnonz >= 0);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);
   SCIP_UNUSED(rownames);

   SCIPdebugMessage("adding %d rows with %d nonzeros to Xpress\n", nrows, nnonz);

   invalidateSolution(lpi);

#ifndef NDEBUG
   {
      /* perform check that no new cols are added - this is forbidden */
      int ncols;
      int j;

      CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );
      for (j = 0; j < nnonz; ++j)
      {
         assert( val[j] != 0.0 );
         assert( 0 <= ind[j] && ind[j] < ncols );
      }
   }
#endif

   /* ensure that the temporary arrays are large enough */
   SCIP_CALL( ensureSidechgMem(lpi, nrows) );
   SCIP_CALL( ensureValMem(lpi, nrows+1) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   convertSides(lpi, nrows, lhs, rhs);

   /* only collect the start array if we have at least one non-zero */
   if( nnonz > 0 )
   {
      for( r = 0; r < nrows; r++ )
         lpi->indarray[r] = beg[r];
      lpi->indarray[nrows] = nnonz;
   }

   CHECK_ZERO( lpi->messagehdlr, XPRSaddrows(lpi->xprslp, nrows, nnonz, lpi->senarray, lpi->rhsarray, lpi->rngarray, lpi->indarray, ind, val) );

     return SCIP_OKAY;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   int r;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   debugCheckRowrang(lpi, firstrow, lastrow);

   SCIPdebugMessage("deleting %d rows from Xpress\n", lastrow - firstrow + 1);

   invalidateSolution(lpi);

   /* ensure that the temporary arrays are large enough */
   SCIP_CALL( ensureValMem(lpi, lastrow-firstrow+1) );

   for( r = firstrow; r <= lastrow; r++ )
      lpi->indarray[r-firstrow] = r;

   CHECK_ZERO( lpi->messagehdlr, XPRSdelrows(lpi->xprslp, lastrow-firstrow+1, lpi->indarray) );

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
   int nkeptrows;
   int ndelrows;
   int nrows;
   int r;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("deleting a row set from Xpress\n");

   invalidateSolution(lpi);

   nkeptrows = 0;
   ndelrows = 0;

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );

   /* ensure that the temporary arrays are large enough */
   SCIP_CALL( ensureValMem(lpi, nrows) );

   /* collect the row indecies which should be deleted and create a the new row ordering */
   for( r = 0; r < nrows; r++ )
   {
      if( dstat[r] == 1 )
      {
         dstat[r] = -1;
         lpi->indarray[ndelrows] = r;
         ndelrows++;
      }
      else
      {
         dstat[r] = nkeptrows;
         nkeptrows++;
      }
   }

   CHECK_ZERO( lpi->messagehdlr, XPRSdelrows(lpi->xprslp, ndelrows, lpi->indarray) );

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int zero = 0;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("clearing Xpress LP\n");

   invalidateSolution(lpi);

   /* create an empty LP in this */
   CHECK_ZERO( lpi->messagehdlr, XPRSloadlp(lpi->xprslp, lpi->name, 0, 0, NULL, NULL, NULL, NULL, &zero, NULL, NULL, NULL, NULL, NULL) );

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
   int j;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(ncols == 0 || (ind != NULL && lb != NULL && ub != NULL));

   SCIPdebugMessage("changing %d bounds in Xpress\n", ncols);
   if( ncols <= 0 )
      return SCIP_OKAY;

   invalidateSolution(lpi);

   for (j = 0; j < ncols; ++j)
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
   }

   /* ensure that the temporary arrays are large enough */
   SCIP_CALL( ensureBoundchgMem(lpi, ncols) );

   CHECK_ZERO( lpi->messagehdlr, XPRSchgbounds(lpi->xprslp, ncols, ind, lpi->larray, (SCIP_Real*)lb) );
   CHECK_ZERO( lpi->messagehdlr, XPRSchgbounds(lpi->xprslp, ncols, ind, lpi->uarray, (SCIP_Real*)ub) );

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(ind != NULL);

   SCIPdebugMessage("changing %d sides in Xpress\n", nrows);
   if( nrows <= 0 )
      return SCIP_OKAY;

   invalidateSolution(lpi);

   /* ensure that the temporary arrays are large enough */
   SCIP_CALL( ensureSidechgMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   convertSides(lpi, nrows, lhs, rhs);

   /* change row sides */
   CHECK_ZERO( lpi->messagehdlr, XPRSchgrowtype(lpi->xprslp, nrows, ind, lpi->senarray) );
   CHECK_ZERO( lpi->messagehdlr, XPRSchgrhs(lpi->xprslp, nrows, ind, lpi->rhsarray) );
   CHECK_ZERO( lpi->messagehdlr, XPRSchgrhsrange(lpi->xprslp, nrows, ind, lpi->rngarray) );

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("changing coefficient row %d, column %d in Xpress to %g\n", row, col, newval);

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, XPRSchgcoef(lpi->xprslp, row, col, newval) );

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsense            /**< new objective sense */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("changing objective sense in Xpress to %d\n", objsense);

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, XPRSchgobjsense(lpi->xprslp, xprsObjsen(objsense)) );

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(ind != NULL);
   assert(obj != NULL);

   SCIPdebugMessage("changing %d objective values in Xpress\n", ncols);

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, XPRSchgobj(lpi->xprslp, ncols, ind, obj) );

   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nnonz;
   int ncols;
   int i;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(scaleval != 0.0);

   SCIPdebugMessage("scaling row %d with factor %g in Xpress\n", row, scaleval);

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );
   SCIP_CALL( ensureValMem(lpi, ncols) );

   /* get the row */
   SCIP_CALL( SCIPlpiGetSides(lpi, row, row, &lhs, &rhs) );
   CHECK_ZERO( lpi->messagehdlr, XPRSgetrows(lpi->xprslp, NULL, lpi->indarray, lpi->valarray, ncols, &nnonz, row, row) );
   assert(nnonz <= ncols);

   /* scale row coefficients */
   for( i = 0; i < nnonz; ++i )
   {
      SCIP_CALL( SCIPlpiChgCoef(lpi, row, lpi->indarray[i], lpi->valarray[i] * scaleval) );
   }

   /* scale row sides */
   if( lhs > XPRS_MINUSINFINITY )
      lhs *= scaleval;
   else if( scaleval < 0.0 )
      lhs = XPRS_PLUSINFINITY;
   if( rhs < XPRS_PLUSINFINITY )
      rhs *= scaleval;
   else if( scaleval < 0.0 )
      rhs = XPRS_MINUSINFINITY;

   if( scaleval > 0.0 )
   {
      SCIP_CALL( SCIPlpiChgSides(lpi, 1, &row, &lhs, &rhs) );
   }
   else
   {
      SCIP_CALL( SCIPlpiChgSides(lpi, 1, &row, &rhs, &lhs) );
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
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   int nnonz;
   int nrows;
   int i;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(scaleval != 0.0);

   SCIPdebugMessage("scaling column %d with factor %g in Xpress\n", col, scaleval);

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   SCIP_CALL( ensureValMem(lpi, nrows) );

   /* get the column */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetlb(lpi->xprslp, &lb, col, col) );
   CHECK_ZERO( lpi->messagehdlr, XPRSgetub(lpi->xprslp, &ub, col, col) );
   CHECK_ZERO( lpi->messagehdlr, XPRSgetcols(lpi->xprslp, NULL, lpi->indarray, lpi->valarray, nrows, &nnonz, col, col) );
   assert(nnonz <= nrows);

   /* get objective coefficient */
   SCIP_CALL( SCIPlpiGetObj(lpi, col, col, &obj) );

   /* scale column coefficients */
   for( i = 0; i < nnonz; ++i )
   {
      SCIP_CALL( SCIPlpiChgCoef(lpi, lpi->indarray[i], col, lpi->valarray[i] * scaleval) );
   }

   /* scale objective value */
   obj *= scaleval;
   SCIP_CALL( SCIPlpiChgObj(lpi, 1, &col, &obj) );

   /* scale column bounds */
   if( lb > XPRS_MINUSINFINITY )
      lb /= scaleval;
   else if( scaleval < 0.0 )
      lb = XPRS_PLUSINFINITY;
   if( ub < XPRS_PLUSINFINITY )
      ub /= scaleval;
   else if( scaleval < 0.0 )
      ub = XPRS_MINUSINFINITY;

   if( scaleval > 0.0 )
   {
      SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &col, &lb, &ub) );
   }
   else
   {
      SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &col, &ub, &lb) );
   }

   return SCIP_OKAY;
}

/**@} */


/**@name Data Accessing Methods
 *
 * @{
 */

/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiGetNRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(nrows != NULL);

   SCIPdebugMessage("getting number of rows\n");

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, nrows) );

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(ncols != NULL);

   SCIPdebugMessage("getting number of columns\n");

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_COLS, ncols) );

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(nnonz != NULL);

   SCIPdebugMessage("getting number of non-zeros\n");

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ELEMS, nnonz) );

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert((lb != NULL && ub != NULL) || (lb == NULL && ub == NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));

   debugCheckColrang(lpi, firstcol, lastcol);

   SCIPdebugMessage("getting columns %d to %d\n", firstcol, lastcol);

   if( lb != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, XPRSgetlb(lpi->xprslp, lb, firstcol, lastcol) );
      CHECK_ZERO( lpi->messagehdlr, XPRSgetub(lpi->xprslp, ub, firstcol, lastcol) );
   }

   if( nnonz != NULL )
   {
      int ntotalnonz;
      int c;

      /* ensure that the temporary buffer array is large enough */
      SCIP_CALL( ensureValMem(lpi, lastcol-firstcol+2) );

      /* get number of nonzero in the whole problem; needed to pass a proper size to XPRSgetcols() function call
       *
       * @note We are assuming that the arrays given by SCIP are large enough. Otherwise we are getting invalid writes
       */
      SCIP_CALL( SCIPlpiGetNNonz(lpi, &ntotalnonz) );

      /* get matrix entries */
      CHECK_ZERO( lpi->messagehdlr, XPRSgetcols(lpi->xprslp, lpi->indarray, ind, val, ntotalnonz, nnonz, firstcol, lastcol) );
      assert(*nnonz <= ntotalnonz);
      assert(lpi->indarray[lastcol-firstcol+1] == *nnonz);

      assert(beg != NULL); /* for lint */
      for( c = 0; c < lastcol-firstcol+1; c++ )
         beg[c] = lpi->indarray[c];
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
   SCIP_Real*            lhss,               /**< buffer to store left hand side vector, or NULL */
   SCIP_Real*            rhss,               /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert((lhss != NULL && rhss != NULL) || (lhss == NULL && rhss == NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));

   debugCheckRowrang(lpi, firstrow, lastrow);

   SCIPdebugMessage("getting rows %d to %d\n", firstrow, lastrow);

   if( lhss != NULL )
   {
      /* get left and right sides */
      SCIP_CALL( SCIPlpiGetSides(lpi, firstrow, lastrow, lhss, rhss) );
   }

   if( nnonz != NULL )
   {
      int ntotalnonz;
      int r;

      /* ensure that the temporary buffer array is large enough */
      SCIP_CALL( ensureValMem(lpi, lastrow-firstrow+2) );

      /* get number of nonzero in the whole problem; needed to pass a proper size to XPRSgetrows() function call
       *
       * @note We are assuming that the arrays given by SCIP are large enough. Otherwise we are getting invalid writes
       */
      SCIP_CALL( SCIPlpiGetNNonz(lpi, &ntotalnonz) );

      /* get matrix entries */
      CHECK_ZERO( lpi->messagehdlr, XPRSgetrows(lpi->xprslp, lpi->indarray, ind, val, ntotalnonz, nnonz, firstrow, lastrow) );
      assert(*nnonz <= ntotalnonz);
      assert(lpi->indarray[lastrow-firstrow+1] == *nnonz);

      assert(beg != NULL); /* for lint */
      for( r = 0; r < lastrow-firstrow+1; r++ )
         beg[r] = lpi->indarray[r];
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
{ /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(colnames != NULL || namestoragesize == 0);
   assert(namestorage != NULL || namestoragesize == 0);
   assert(namestoragesize >= 0);
   assert(storageleft != NULL);
   assert(0 <= firstcol && firstcol <= lastcol);

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
{ /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(rownames != NULL || namestoragesize == 0);
   assert(namestorage != NULL || namestoragesize == 0);
   assert(namestoragesize >= 0);
   assert(storageleft != NULL);
   assert(0 <= firstrow && firstrow <= lastrow);

   SCIPerrorMessage("SCIPlpiGetRowNames() has not been implemented yet.\n");
   return SCIP_LPERROR;
}

/** gets the objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   double xprsobjsen;
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(objsen != NULL);

   /* check the objective sense attribute for the current objective sense set in  Xpress */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetdblattrib(lpi->xprslp, XPRS_OBJSENSE, &xprsobjsen) );

   /* convert the Xpress objective sense attribute to a SCIP objective sense */
   if( xprsobjsen < 0.0 )
      (*objsen) = SCIP_OBJSEN_MAXIMIZE;
   else
      (*objsen) = SCIP_OBJSEN_MINIMIZE;

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(firstcol <= lastcol);
   assert(vals != NULL);

   SCIPdebugMessage("getting objective values %d to %d\n", firstcol, lastcol);

   CHECK_ZERO( lpi->messagehdlr, XPRSgetobj(lpi->xprslp, vals, firstcol, lastcol) );

   return SCIP_OKAY;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(firstcol <= lastcol);

   SCIPdebugMessage("getting bounds %d to %d\n", firstcol, lastcol);

   if( lbs != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, XPRSgetlb(lpi->xprslp, lbs, firstcol, lastcol) );
   }

   if( ubs != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, XPRSgetub(lpi->xprslp, ubs, firstcol, lastcol) );
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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(firstrow <= lastrow);

   SCIPdebugMessage("getting row sides %d to %d\n", firstrow, lastrow);

   /* ensure the array size of the temporary buffers */
   SCIP_CALL( ensureSidechgMem(lpi, lastrow - firstrow + 1) );

   /* get row sense, rhs, and ranges */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetrowtype(lpi->xprslp, lpi->senarray, firstrow, lastrow) );
   CHECK_ZERO( lpi->messagehdlr, XPRSgetrhs(lpi->xprslp, lpi->rhsarray, firstrow, lastrow) );
   CHECK_ZERO( lpi->messagehdlr, XPRSgetrhsrange(lpi->xprslp, lpi->rngarray, firstrow, lastrow) );

   /* convert sen/rhs/range into lhs/rhs tuples */
   reconvertSides(lpi, lastrow - firstrow + 1, lhss, rhss);

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(val != NULL);

   /* get the coefficient of the column in the corresponding row */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetcoef(lpi->xprslp, row, col, val) );

   return SCIP_OKAY;
}

/**@} */


/**@name Solving Methods
 *
 * @{
 */

/** solve LP */
static SCIP_RETCODE lpiSolve(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           method              /**< indicates the method to use ('p' - primal, 'd' - dual, 'b' - barrier) */
   )
{
   int primalinfeasible;
   int dualinfeasible;
   int state;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   invalidateSolution(lpi);

   /* check if the current basis should be ignored */
   if( lpi->clearstate )
   {
      CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_KEEPBASIS, 0) );
      lpi->clearstate = FALSE;
   }

   CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_PRESOLVE, 0) );
   CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_LPQUICKPRESOLVE, (lpi->par_presolve) ?  1 : 0) );

   if( lpi->par_fastlp )
   {
      /* Don't refactorize at the end of the solve. */
      CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_REFACTOR, 0) );
   }
   else
   {
      /* Use default settings for solving an lp (hopefully) robustly. */
      CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_REFACTOR, 1) );
   }

   /* solve the LP */
   CHECK_ZERO( lpi->messagehdlr, XPRSlpoptimize(lpi->xprslp, method) );

   /* evaluate the result */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_LPSTATUS, &lpi->solstat) );

   /* Make sure the LP is postsolved in case it was interrupted. */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_PRESOLVESTATE, &state) );

   if( state & (2|4) )
   {
      /* Problem is in a presolve state - postsolve it. */
      CHECK_ZERO( lpi->messagehdlr, XPRSpostsolve(lpi->xprslp) );
   }

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_SIMPLEXITER, &lpi->iterations) );
   lpi->solisbasic = TRUE;

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_PRIMALINFEAS, &primalinfeasible) );
   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_DUALINFEAS, &dualinfeasible) );
   SCIPdebugMessage(" -> Xpress returned solstat=%d, pinfeas=%d, dinfeas=%d (%d iterations)\n",
      lpi->solstat, primalinfeasible, dualinfeasible, lpi->iterations);

   /* Make sure that always a primal / dual ray exists */
   if( lpi->solstat == XPRS_LP_INFEAS  || lpi->solstat == XPRS_LP_UNBOUNDED )
   {
      int hasray;
      int presolving;

      /* check whether a dual ray exists, in that case we don't need to resolve the LP w/o presolving */
      CHECK_ZERO( lpi->messagehdlr, XPRSgetdualray(lpi->xprslp, NULL, &hasray) );

      if( hasray == 1 )
         goto TERMINATE;

      /* get the current presolving setting */
      CHECK_ZERO( lpi->messagehdlr, XPRSgetintcontrol(lpi->xprslp, XPRS_LPQUICKPRESOLVE, &presolving) );

      if( presolving != 0 )
      {
         int iterations;

         /* maybe the preprocessor solved the problem; but we need a solution, so solve again without preprocessing */
         SCIPdebugMessage("presolver may have solved the problem -> calling Xpress %s again without presolve\n",
               strcmp(method, "p") == 0 ? "primal simplex" : strcmp(method, "d") == 0 ? "dual simplex" : "barrier");

         /* switch off preprocessing */
         CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_LPQUICKPRESOLVE, 0) );

         /* resolve w/o presolving */
         CHECK_ZERO( lpi->messagehdlr, XPRSlpoptimize(lpi->xprslp, method) );

         /* evaluate the result */
         CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_LPSTATUS, &lpi->solstat) );

         CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_SIMPLEXITER, &iterations) );
         lpi->iterations += iterations;
         lpi->solisbasic = TRUE;

         CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_PRIMALINFEAS, &primalinfeasible) );
         CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_DUALINFEAS, &dualinfeasible) );
         SCIPdebugMessage(" -> Xpress returned solstat=%d, pinfeas=%d, dinfeas=%d (%d iterations)\n",
            lpi->solstat, primalinfeasible, dualinfeasible, lpi->iterations);

         /* reinstall the previous setting */
         CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_LPQUICKPRESOLVE, presolving) );
      }
   }

  TERMINATE:
   if( (lpi->solstat == XPRS_LP_OPTIMAL) && (primalinfeasible || dualinfeasible) )
      lpi->solstat = XPRS_LP_OPTIMAL_SCALEDINFEAS;

   return SCIP_OKAY;
}

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   lpi->solmethod = 'p';
   return lpiSolve(lpi, "p");
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   lpi->solmethod = 'd';
   return lpiSolve(lpi, "d");
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{
   SCIP_RETCODE retval;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   lpi->solmethod = 'b';

   /* enable or disable cross over */
   CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_CROSSOVER, crossover == TRUE ? -1 : 0) );

   retval = lpiSolve(lpi, "b");
   lpi->solisbasic = crossover;

   return retval;
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{ /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   /* currently do nothing */
   return SCIP_OKAY;
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   /* currently do nothing */
   return SCIP_OKAY;
}

/** performs strong branching iterations on one candidate */
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
   SCIP_OBJSEN  objsen;
   double dbndval[2];
   double dobjval[2];
   char cbndtype[2];
   int mbndind[2];
   int mstatus[2];

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   SCIPdebugMessage("calling Xpress strong branching on variable %d (%d iterations)\n", col, itlim);

   /* results of Xpress are valid in any case */
   *downvalid = TRUE;
   *upvalid = TRUE;

   SCIPdebugMessage(" -> strong branching on integral variable\n");

   if( iter != NULL )
      *iter = 0;

   /* get objective sense of the current LP */
   SCIP_CALL( SCIPlpiGetObjsen(lpi, &objsen) );

   /* Set the branching bounds (down first, up second). */
   mbndind[0]  = col;
   dbndval[0]  = EPSCEIL(psol-1.0, 1e-06);
   cbndtype[0] = 'U';
   mbndind[1]  = col;
   dbndval[1]  = EPSFLOOR(psol+1.0, 1e-06);
   cbndtype[1] = 'L';

   /* Apply strong branching to the two branches. */
   CHECK_ZERO( lpi->messagehdlr, XPRSstrongbranch(lpi->xprslp, 2, mbndind, cbndtype, dbndval, itlim, dobjval, mstatus) );

   /* Get the objective of the down branch. */
   if( (mstatus[0] == XPRS_LP_INFEAS) || (mstatus[0] == XPRS_LP_CUTOFF_IN_DUAL) )
      *down = objsen == SCIP_OBJSEN_MINIMIZE  ? 1e+40 : -1e+40;
   else if( (mstatus[0] == XPRS_LP_OPTIMAL) || (mstatus[0] == XPRS_LP_UNFINISHED) )
      *down = dobjval[0];
   else
   {
      /* Something weird happened. */
      *downvalid = FALSE;
   }

   /* Get the objective of the up branch. */
   if( (mstatus[1] == XPRS_LP_INFEAS) || (mstatus[1] == XPRS_LP_CUTOFF_IN_DUAL) )
      *up = objsen == SCIP_OBJSEN_MINIMIZE ? 1e+40 : -1e+40;
   else if( (mstatus[1] == XPRS_LP_OPTIMAL) || (mstatus[1] == XPRS_LP_UNFINISHED) )
      *up = dobjval[1];
   else
   {
      /* Something weird happened. */
      *upvalid = FALSE;
   }

   /* When using the XPRSstrongbranch function we are unable to provide an iteration count */
   if( iter != NULL )
      *iter = -1;

   return SCIP_OKAY;
}

/** performs strong branching iterations on given candidates */
static
SCIP_RETCODE lpiStrongbranches(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< current primal solution values of columns (might be integral) */
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
   double* dbndval;
   double* dobjval;
   char*   cbndtype;
   int*    mbndind;
   int*    mstatus;
   SCIP_OBJSEN objsen;
   int nbranches;
   int j;

   assert( lpi != NULL );
   assert( lpi->xprslp != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   SCIPdebugMessage("calling Xpress strong branching on %d variables (%d iterations)\n", ncols, itlim);

   if( iter != NULL )
      *iter = 0;

   /* compute the number of branches; for each column we have 2 branches */
   nbranches = 2*ncols;

   /* get objective sense of the current LP */
   SCIP_CALL( SCIPlpiGetObjsen(lpi, &objsen) );

   /* Set the branching bounds (down first, up second). */
   SCIP_ALLOC( BMSallocMemoryArray(&mbndind, nbranches) );
   SCIP_ALLOC( BMSallocMemoryArray(&dbndval, nbranches) );
   SCIP_ALLOC( BMSallocMemoryArray(&cbndtype, nbranches) );
   SCIP_ALLOC( BMSallocMemoryArray(&dobjval, nbranches) );
   SCIP_ALLOC( BMSallocMemoryArray(&mstatus, nbranches) );

   /* construct the bounds for the strong branches */
   for( j = 0; j < ncols; ++j )
   {
      mbndind[2*j] = cols[j];
      dbndval[2*j] = EPSCEIL(psols[j] - 1.0, 1e-06);
      cbndtype[2*j] = 'U';

      mbndind[2*j+1] = cols[j];
      dbndval[2*j+1] = EPSFLOOR(psols[j] + 1.0, 1e-06);
      cbndtype[2*j+1] = 'L';
   }

   /* apply strong branching to the 2*ncols branches. */
   CHECK_ZERO( lpi->messagehdlr, XPRSstrongbranch(lpi->xprslp, nbranches, mbndind, cbndtype, dbndval, itlim, dobjval, mstatus) );

   for( j = 0; j < ncols; ++j )
   {
      upvalid[j]   = TRUE;
      downvalid[j] = TRUE;

      /* Get the objective of the down branch. */
      if( (mstatus[2*j] == XPRS_LP_INFEAS) || (mstatus[2*j] == XPRS_LP_CUTOFF_IN_DUAL) )
         down[j] = objsen == SCIP_OBJSEN_MINIMIZE ? 1e+40 : -1e+40;
      else if( (mstatus[2*j] == XPRS_LP_OPTIMAL) || (mstatus[2*j] == XPRS_LP_UNFINISHED) )
         down[j] = dobjval[2*j];
      else
      {
         /* Something weird happened. */
         downvalid[j] = FALSE;
      }

      /* Get the objective of the up branch. */
      if( (mstatus[2*j+1] == XPRS_LP_INFEAS) || (mstatus[2*j+1] == XPRS_LP_CUTOFF_IN_DUAL) )
         up[j] = objsen == SCIP_OBJSEN_MINIMIZE ? 1e+40 : -1e+40;
      else if( (mstatus[2*j+1] == XPRS_LP_OPTIMAL) || (mstatus[2*j+1] == XPRS_LP_UNFINISHED) )
         up[j] = dobjval[2*j+1];
      else
      {
         /* Something weird happened. */
         upvalid[j] = FALSE;
      }
   }

   /* When using the XPRSstrongbranch function we are unable to provide
    * an iteration count.
    */
   if( iter != NULL )
      *iter = -1;

   BMSfreeMemoryArray(&mstatus);
   BMSfreeMemoryArray(&dobjval);
   BMSfreeMemoryArray(&cbndtype);
   BMSfreeMemoryArray(&dbndval);
   BMSfreeMemoryArray(&mbndind);

   return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPlpiStrongbranchFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< fractional current primal solution value of column */
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
   /* pass call on to lpiStrongbranches() */
   SCIP_CALL( lpiStrongbranches(lpi, cols, ncols, psols, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/**@} */


/**@name Solution Information Methods
 *
 * @{
 */

/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiWasSolved(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   return (lpi->solstat != -1);
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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);

   SCIPdebugMessage("getting solution feasibility\n");

   *primalfeasible = SCIPlpiIsPrimalFeasible(lpi);
   *dualfeasible = SCIPlpiIsDualFeasible(lpi);

   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == XPRS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int hasRay;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   /* check if the LP solution status is unbounded and that primal was solving the LP */
   if (lpi->solstat != XPRS_LP_UNBOUNDED || lpi->solmethod != 'p')
      return FALSE;

   /* check if we can construct a primal ray */
   ABORT_ZERO( lpi->messagehdlr, FALSE, XPRSgetprimalray(lpi->xprslp, NULL, &hasRay) );

   return (SCIP_Bool)hasRay;
}

/** returns TRUE iff LP is proven to be primal feasible and unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for primal unboundedness\n");

   /* If the solution status of Xpress is XPRS_LP_UNBOUNDED, it only means, there is an unbounded ray,
    * but not necessarily a feasible primal solution. If problem is declared LP_UNBOUNDED by dual,
    * we have no way to decide primal feasibility.
    */

   return lpi->solstat == XPRS_LP_UNBOUNDED && lpi->solmethod == 'p';
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for primal infeasibility\n");

   return (lpi->solstat == XPRS_LP_INFEAS);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int nInfeasible;
   int nIter;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for primal feasibility\n");

   /* check if problem is solved to optimality */
   if (lpi->solstat == XPRS_LP_OPTIMAL || lpi->solstat == XPRS_LP_OPTIMAL_SCALEDINFEAS)
     return TRUE;

   /* check if problem is unbounded (found by primal) */
   if (lpi->solstat == XPRS_LP_UNBOUNDED && lpi->solmethod == 'p')
     return TRUE;

   /* get number of primal infeasibilities and number of simplex iterations */
   ABORT_ZERO( lpi->messagehdlr, FALSE, XPRSgetintattrib(lpi->xprslp, XPRS_PRIMALINFEAS, &nInfeasible) );
   ABORT_ZERO( lpi->messagehdlr, FALSE, XPRSgetintattrib(lpi->xprslp, XPRS_SIMPLEXITER, &nIter) );

   /* check if the number of primal infeasibilities is zero
    * We need to make sure that the LP was indeed solved by primal, otherwise infeasibility might have been found
    * in setup (e.g. if conflicting bounds x >= 1, x <= 0 are present),
    */
   if (nInfeasible == 0  && nIter > 0 && lpi->solmethod == 'p')
     return TRUE;

   return FALSE;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == XPRS_LP_INFEAS);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int hasRay;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   ABORT_ZERO( lpi->messagehdlr, FALSE, XPRSgetdualray(lpi->xprslp, NULL, &hasRay) );

   return (SCIP_Bool) hasRay;
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for dual unboundedness\n");

   return ((lpi->solstat == XPRS_LP_INFEAS) && (lpi->solmethod == 'd'));
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for dual infeasibility\n");

   return (lpi->solstat == XPRS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int nInfeasible;
   int nIter;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for dual feasibility\n");

   /* check if problem solved to optimality */
   if (lpi->solstat == XPRS_LP_OPTIMAL || lpi->solstat == XPRS_LP_OPTIMAL_SCALEDINFEAS)
     return TRUE;

   /* check if problem infeasibility detected by dual */
   if (lpi->solstat == XPRS_LP_INFEAS && lpi->solmethod == 'd')
     return TRUE;

   /* get number of dual infeasibilities and number of simplex iterations */
   ABORT_ZERO( lpi->messagehdlr, FALSE, XPRSgetintattrib(lpi->xprslp, XPRS_DUALINFEAS, &nInfeasible) );
   ABORT_ZERO( lpi->messagehdlr, FALSE, XPRSgetintattrib(lpi->xprslp, XPRS_SIMPLEXITER, &nIter) );

   /* check if the number of dual infeasibilities is zero
    * We need to make sure that the LP was indeed solved by primal, otherwise infeasibility might have been found
    * in setup (e.g. if conflicting bounds x >= 1, x <= 0 are present),
    */
   if (nInfeasible == 0 && nIter > 0 && lpi->solmethod == 'd')
      return TRUE;

   return FALSE;
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == XPRS_LP_OPTIMAL) || (lpi->solstat == XPRS_LP_OPTIMAL_SCALEDINFEAS);
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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for stability: Xpress solstat = %d\n", lpi->solstat);

#ifdef SCIP_DISABLED_CODE
   /* The following workaround is not needed anymore for SCIP, since it tries to heuristically construct a feasible
    * solution or automatically resolves the problem if the status is "unbounded"; see SCIPlpGetUnboundedSol().
    */

   /* If the solution status of Xpress is XPRS_LP_UNBOUNDED, it only means, there is an unbounded ray,
    * but not necessarily a feasible primal solution. If primalfeasible == FALSE, we interpret this
    * result as instability, s.t. the problem is resolved from scratch
    */
   if( lpi->solstat == XPRS_LP_UNBOUNDED )
   {
      int retcode;
      int pinfeas;

      retcode = XPRSgetintattrib(lpi->xprslp, XPRS_PRIMALINFEAS, &pinfeas);

      if( retcode != 0 || pinfeas )
         return FALSE;
   }
#endif

   if( lpi->solstat == XPRS_LP_OPTIMAL_SCALEDINFEAS )
   {
      /* presolved problem was solved to optimality but infeasibilities were introduced by postsolve */
      return FALSE;
   }

   return TRUE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == XPRS_LP_CUTOFF || lpi->solstat == XPRS_LP_CUTOFF_IN_DUAL);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int lpiter;
   int lpiterlimit;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   ABORT_ZERO( lpi->messagehdlr, TRUE, XPRSgetintattrib(lpi->xprslp, XPRS_SIMPLEXITER, &lpiter) );
   ABORT_ZERO( lpi->messagehdlr, TRUE, XPRSgetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, &lpiterlimit) );

   if( (lpi->solstat == XPRS_LP_UNFINISHED) && (lpiter >= lpiterlimit) )
      return TRUE;
   else
      return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int lpiter;
   int lpiterlimit;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   ABORT_ZERO( lpi->messagehdlr, TRUE, XPRSgetintattrib(lpi->xprslp, XPRS_SIMPLEXITER, &lpiter) );
   ABORT_ZERO( lpi->messagehdlr, TRUE, XPRSgetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, &lpiterlimit) );

   if( (lpi->solstat == XPRS_LP_UNFINISHED) && (lpiter < lpiterlimit) )
      return TRUE;
   else
      return FALSE;
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   return lpi->solstat;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(success != NULL);

   /* Nothing to do here for Xpress. */
   *success = TRUE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(objval != NULL);

   SCIPdebugMessage("getting solution's objective value\n");

   CHECK_ZERO( lpi->messagehdlr, XPRSgetdblattrib(lpi->xprslp, XPRS_LPOBJVAL, objval) );

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("getting solution\n");

   CHECK_ZERO( lpi->messagehdlr, XPRSgetsol(lpi->xprslp, primsol, activity, dualsol, redcost) );

   if( objval != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, XPRSgetdblattrib(lpi->xprslp, XPRS_LPOBJVAL, objval) );
   }

   if( activity != NULL )
   {
      /* Convert the slack values into activity values. */
      int nrows;
      int r;

      CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );

      SCIP_CALL( ensureSidechgMem(lpi, nrows) );

      CHECK_ZERO( lpi->messagehdlr, XPRSgetrhs(lpi->xprslp, lpi->rhsarray, 0, nrows-1) );

      for( r = 0; r < nrows; r++ )
         activity[r] = lpi->rhsarray[r] - activity[r];
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
   int hasRay;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(ray != NULL);
   assert(lpi->solstat >= 0);

   CHECK_ZERO( lpi->messagehdlr, XPRSgetprimalray(lpi->xprslp, ray, &hasRay) );

   if( !hasRay )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
   )
{
   int hasRay;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);
   assert(dualfarkas != NULL);

   /**@note The Farkas proof might be numerically questionable which is indicated by "hasRay" use SCIPlpiHasDualRay() to
    *       check that!
    */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetdualray(lpi->xprslp, dualfarkas, &hasRay) );

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(iterations != NULL);

   *iterations = lpi->iterations;

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
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(quality != NULL);
   SCIP_UNUSED(qualityindicator);

   *quality = SCIP_INVALID;

   return SCIP_OKAY;
}

/**@} */


/**@name LP Basis Methods
 *
 * @{
 */

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiGetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL (the status is need for the row and not for the slack column) */
   )
{
   int nrows;
   int r;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   /*lint --e{506}*/
   assert((int) SCIP_BASESTAT_LOWER == 0);
   assert((int) SCIP_BASESTAT_BASIC == 1);
   assert((int) SCIP_BASESTAT_UPPER == 2);

   SCIPdebugMessage("saving Xpress basis into %p/%p\n", (void*)rstat, (void*)cstat);

   /* get the basis status */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetbasis(lpi->xprslp, rstat, cstat) );

   /* get the number of rows */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );

   /* XPRSgetbasis collects the basis status for the column and the slack column, since SCIP request the basis for
    * columns and rows we need to convert slack column status to row status
    */
   for( r = 0; r < nrows; ++r )
   {
     if (rstat[r] == (int) SCIP_BASESTAT_LOWER)
         rstat[r] = (int) SCIP_BASESTAT_UPPER;
      else if (rstat[r] == (int) SCIP_BASESTAT_UPPER)
         rstat[r] = (int) SCIP_BASESTAT_LOWER;
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
   int* slackstats;
   int ncols;
   int nrows;
   int r;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   /*  get the number of rows/columns */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );

   assert(cstat != NULL || ncols == 0);
   assert(rstat != NULL || nrows == 0);

   /*lint --e{506}*/
   assert((int) SCIP_BASESTAT_LOWER == 0);
   assert((int) SCIP_BASESTAT_BASIC == 1);
   assert((int) SCIP_BASESTAT_UPPER == 2);

   SCIPdebugMessage("loading basis %p/%p into Xpress\n", (void*)rstat, (void*)cstat);

   invalidateSolution(lpi);

   SCIP_ALLOC( BMSallocMemoryArray(&slackstats, nrows) );

   /* XPRSloadbasis expects the basis status for the column and the slack column, since SCIP has the basis status for
    * columns and rows we need to convert row status to slack column status
    */
   for( r = 0; r < nrows; ++r )
   {
      if (rstat[r] == (int) SCIP_BASESTAT_LOWER) /*lint !e613*/
         slackstats[r] = (int) SCIP_BASESTAT_UPPER;
      else if (rstat[r] == (int) SCIP_BASESTAT_UPPER) /*lint !e613*/
         slackstats[r] = (int) SCIP_BASESTAT_LOWER;
      else
         slackstats[r] = rstat[r]; /*lint !e613*/
   }

   /* load basis information into Xpress
    *
    * @note Xpress expects the row status w.r.t. slack columns!
    */
   CHECK_ZERO( lpi->messagehdlr, XPRSloadbasis(lpi->xprslp, slackstats, cstat) );

   BMSfreeMemoryArray(&slackstats);

   lpi->clearstate = FALSE;

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   )
{
   int irspace;
   int nrows;
   int r;

   /* In the basis methods we assume that xprs basis flags coincide with scip, so assert it */
   /*lint --e{506}*/
   assert((int) SCIP_BASESTAT_LOWER == 0);
   assert((int) SCIP_BASESTAT_BASIC == 1);
   assert((int) SCIP_BASESTAT_UPPER == 2);
   assert((int) SCIP_BASESTAT_ZERO == 3);

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(bind != NULL);

   SCIPdebugMessage("getting basis information\n");

   CHECK_ZERO( lpi->messagehdlr, XPRSgetpivotorder(lpi->xprslp, bind) );

   /* Reindex variables to match those of SCIP. */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_SPAREROWS, &irspace) );
   irspace += nrows;

   for( r = 0; r < nrows; r++ )
   {
      if( bind[r] < nrows )
         bind[r] = -bind[r]-1;
      else
      {
         assert(bind[r] >= irspace);
         bind[r] = bind[r] - irspace;
      }
   }

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
   int                   row,                /**< row number */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{  /*lint --e{715}*/
   int nrows;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(coef != NULL);
   SCIP_UNUSED(inds);

   SCIPdebugMessage("getting binv-row %d\n", row);

   /* can only return dense result */
   if ( ninds != NULL )
      *ninds = -1;

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   BMSclearMemoryArray(coef, nrows);
   coef[row] = 1.0;
   CHECK_ZERO( lpi->messagehdlr, XPRSbtran(lpi->xprslp, coef) );

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
{  /*lint --e{715}*/
   int nrows;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(coef != NULL);
   SCIP_UNUSED(inds);

   SCIPdebugMessage("getting binv-col %d\n", c);

   /* can only return dense result */
   if ( ninds != NULL )
      *ninds = -1;

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   BMSclearMemoryArray(coef, nrows);
   coef[c] = 1.0;
   CHECK_ZERO( lpi->messagehdlr, XPRSftran(lpi->xprslp, coef) );

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
{  /*lint --e{715}*/
   SCIP_Real* binv;
   SCIP_Real* buffer;
   int ncols;
   int nrows;
   int nnonz;
   int c;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(coef != NULL);

   SCIPdebugMessage("getting binva-row %d\n", r);

   /* can only return dense result */
   if ( ninds != NULL )
      *ninds = -1;

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );

   buffer = NULL;

   /* get (or calculate) the row in B^-1 */
   if( binvrow == NULL )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&buffer, nrows) );
      SCIP_CALL( SCIPlpiGetBInvRow(lpi, r, buffer, inds, ninds) );
      binv = buffer;
   }
   else
      binv = (double*) binvrow;

   /* We need space to extract a single column. */
   SCIP_CALL( ensureValMem(lpi, nrows) );

   for( c = 0; c < ncols; c++ )
   {
      int i;

      coef[c] = 0;

      /* Extract the column. */
      CHECK_ZERO( lpi->messagehdlr, XPRSgetcols(lpi->xprslp, NULL, lpi->indarray, lpi->valarray, nrows, &nnonz, c, c) );
      assert(nnonz <= nrows);

      /* Price out the column. */
      for( i = 0; i < nnonz; i++ )
         coef[c] += binv[lpi->indarray[i]] * lpi->valarray[i];
   }

   /* Free allocated memory. */
   BMSfreeMemoryArrayNull(&buffer);

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
{  /*lint --e{715}*/
   int nrows;
   int nnonz;
   int i;

   /* Ftran */

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(coef != NULL);
   SCIP_UNUSED(inds);

   SCIPdebugMessage("getting binv-col %d\n", c);

   /* can only return dense result */
   if ( ninds != NULL )
      *ninds = -1;

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );

   /* We need space to extract the column. */
   SCIP_CALL( ensureValMem(lpi, nrows) );

   /* Get the column to transform. */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetcols(lpi->xprslp, NULL, lpi->indarray, lpi->valarray, nrows, &nnonz, c, c) );
   assert(nnonz <= nrows);

   /* Transform the column. */
   BMSclearMemoryArray(coef, nrows);
   for( i = 0; i < nnonz; i++ )
      coef[lpi->indarray[i]] = lpi->valarray[i];

   CHECK_ZERO( lpi->messagehdlr, XPRSftran(lpi->xprslp, coef) );

   return SCIP_OKAY;
}

/**@} */


/**@name LP State Methods
 *
 * @{
 */

/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   int ncols;
   int nrows;

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpistate != NULL);

   /* if there is no basis information available (e.g. after barrier without crossover), or no state can be saved; if
    * SCIPlpiClearState() has been called, do not return the state
    */
   if( !lpi->solisbasic || lpi->clearstate )
   {
      *lpistate = NULL;
      return SCIP_OKAY;
   }

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );
   assert(ncols >= 0);
   assert(nrows >= 0);

   /* allocate lpistate data */
   SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows) );

   SCIPdebugMessage("storing Xpress LPI state in %p (%d cols, %d rows)\n", (void*)*lpistate, ncols, nrows);

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* get unpacked basis information from Xpress
    *
    * @note The row status is w.r.t. slack columns!
    */
   CHECK_ZERO( lpi->messagehdlr, XPRSgetbasis(lpi->xprslp, lpi->rstat, lpi->cstat) );

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
{
   int nrows;
   int ncols;
   int i;

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   /* if there was no basis information available, the LPI state was not stored */
   if( lpistate == NULL )
      return SCIP_OKAY;

   if( lpistate->ncols == 0 || lpistate->nrows == 0 )
      return SCIP_OKAY;

   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   CHECK_ZERO( lpi->messagehdlr, XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );

   /* the dimension of the lpi state should not be larger than the current problem; it might be that columns and rows
    * are added since the saving of the lpi state
    */
   assert(lpistate->ncols <= ncols);
   assert(lpistate->nrows <= nrows);

   SCIPdebugMessage("loading LPI state %p (%d cols, %d rows) into Xpress\n", (void*)lpistate, lpistate->ncols, lpistate->nrows);

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* unpack LPi state data
    *
    * @note The row status is w.r.t. slack column!
    */
   lpistateUnpack(lpistate, lpi->cstat, lpi->rstat);

   /* extend the basis to the current LP beyond the previously existing columns */
   for( i = lpistate->ncols; i < ncols; ++i )
   {
      SCIP_Real bnd;

      CHECK_ZERO( lpi->messagehdlr, XPRSgetlb(lpi->xprslp, &bnd, i, i) );

      if( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         CHECK_ZERO( lpi->messagehdlr, XPRSgetub(lpi->xprslp, &bnd, i, i) );

         if( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
            lpi->cstat[i] = (int) SCIP_BASESTAT_ZERO;  /* variable is free */
         else
            lpi->cstat[i] = (int) SCIP_BASESTAT_UPPER; /* use finite upper bound */
      }
      else
        lpi->cstat[i] = (int) SCIP_BASESTAT_LOWER;    /* use finite lower bound */
   }
   for( i = lpistate->nrows; i < nrows; ++i )
      lpi->rstat[i] = (int) SCIP_BASESTAT_BASIC;

   /* load basis information into Xpress
    *
    * @note Xpress expects the row status w.r.t. slack columns!
    */
   CHECK_ZERO( lpi->messagehdlr, XPRSloadbasis(lpi->xprslp, lpi->rstat, lpi->cstat) );

   lpi->clearstate = FALSE;

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   /* set KEEPBASIS to 0 for the next solve */
   lpi->clearstate = TRUE;

   return SCIP_OKAY;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   assert(lpi != NULL);
   assert(lpistate != NULL);
   assert(blkmem != NULL);

   if( *lpistate != NULL )
   {
      lpistateFree(lpistate, blkmem);
   }

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information), or NULL */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   return (lpistate != NULL);
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(fname != NULL);

   SCIPdebugMessage("reading LP state from file <%s>\n", fname);

   CHECK_ZERO( lpi->messagehdlr, XPRSreadbasis(lpi->xprslp, fname, "") );

   return SCIP_OKAY;
}

/** writes LPi state (i.e. basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(fname != NULL);

   SCIPdebugMessage("writing LP state to file <%s>\n", fname);

   CHECK_ZERO( lpi->messagehdlr, XPRSwritebasis(lpi->xprslp, fname, "") );

   return SCIP_OKAY;
}

/**@} */


/**@name LP Pricing Norms Methods
 *
 * @{
 */

/** stores LPi pricing norms information
 *  @todo should we store norm information?
 */
SCIP_RETCODE SCIPlpiGetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information */
   )
{ /*lint --e{715}*/
   assert(lpi != NULL);
   assert(blkmem != NULL);
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
{ /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpinorms == NULL);
   SCIP_UNUSED(blkmem);

   /* no work necessary */
   return SCIP_OKAY;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPlpiFreeNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information, or NULL */
   )
{ /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpinorms == NULL);
   SCIP_UNUSED(blkmem);

   /* no work necessary */
   return SCIP_OKAY;
}

/**@} */


/**@name Parameter Methods
 *
 * @{
 */

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiGetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{
   int ictrlval;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(ival != NULL);

   SCIPdebugMessage("getting int parameter %d\n", type);

   switch( type )
   {
   case SCIP_LPPAR_PRICING:
      *ival = (int)lpi->pricing;
      break;
   case SCIP_LPPAR_FROMSCRATCH:
#if  1
      *ival = (lpi->notfromscratch == 0);
#else
      CHECK_ZERO( lpi->messagehdlr, XPRSgetintcontrol(lpi->xprslp, XPRS_KEEPBASIS, &ictrlval) );
      *ival = (ictrlval == 0);
#endif
      break;
   case SCIP_LPPAR_SCALING:
      CHECK_ZERO( lpi->messagehdlr, XPRSgetintcontrol(lpi->xprslp, XPRS_SCALING, &ictrlval) );
      if( ictrlval == 0 )
         *ival = 0;
      else if( ictrlval == 16 )
         *ival = 2;
      else
         *ival = 1;
      break;
   case SCIP_LPPAR_PRESOLVING:
      *ival = lpi->par_presolve;
      break;
   case SCIP_LPPAR_LPINFO:
      CHECK_ZERO( lpi->messagehdlr, XPRSgetintcontrol(lpi->xprslp, XPRS_OUTPUTLOG, &ictrlval) );
      *ival = (ictrlval != 0);
      break;
   case SCIP_LPPAR_LPITLIM:
      CHECK_ZERO( lpi->messagehdlr, XPRSgetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, &ictrlval) );
      *ival = ictrlval;
      if( *ival >= XPRS_MAXINT )
         *ival = XPRS_MAXINT;
      break;
   case SCIP_LPPAR_THREADS:
      CHECK_ZERO( lpi->messagehdlr, XPRSgetintcontrol(lpi->xprslp, XPRS_THREADS, &ictrlval) );
      *ival = ictrlval;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiSetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("setting int parameter %d to %d\n", type, ival);

   switch( type )
   {
   case SCIP_LPPAR_PRICING:
      /* every non-supported pricing strategy is promoted to the default pricing strategy */
      lpi->pricing = (SCIP_PRICING)ival; /* store pricing method in LPI struct */
      switch( lpi->pricing )
      {
      case SCIP_PRICING_PARTIAL:
         CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_PRICINGALG, XPRS_PRICING_PARTIAL) );
         break;
      case SCIP_PRICING_DEVEX:
         CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_PRICINGALG, XPRS_PRICING_DEVEX) );
         break;
      case SCIP_PRICING_AUTO:
      case SCIP_PRICING_FULL:
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_STEEP:
      case SCIP_PRICING_STEEPQSTART:
      default:
         CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_PRICINGALG, XPRS_PRICING_DEFAULT) );
         break;
      }
      break;
   case SCIP_LPPAR_FROMSCRATCH:
      assert(ival == TRUE || ival == FALSE);
      lpi->notfromscratch = (int)(!ival);
      CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_KEEPBASIS, (ival == FALSE) ? 1 : 0) );
      break;
   case SCIP_LPPAR_SCALING:
      assert(ival >= 0 && ival <= 2);
      if( ival == 0 )
      {
         CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_SCALING, 0) );
      }
      else if( ival == 1 )
      {
         CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_SCALING, 163) );
      }
      else
      {
         CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_SCALING, 16) );
      }

      break;
   case SCIP_LPPAR_PRESOLVING:
      assert(ival == TRUE || ival == FALSE);
      lpi->par_presolve = ival;
      break;
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_OUTPUTLOG, (ival == TRUE) ? 1 : 0) );
      break;
   case SCIP_LPPAR_LPITLIM:
      assert( ival >= 0 );
      /* 0 <= ival, 0 stopping immediately */
      ival = MIN(ival, XPRS_MAXINT); /*lint !e685*//*lint !e2650*//*lint !e587*/
      CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, ival) );
      break;
   case SCIP_LPPAR_THREADS:
      CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_THREADS, ival) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiGetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   int ictrlval;
   double dctrlval;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(dval != NULL);

   SCIPdebugMessage("getting real parameter %d\n", type);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      CHECK_ZERO( lpi->messagehdlr, XPRSgetdblcontrol(lpi->xprslp, XPRS_FEASTOL, &dctrlval) );
      *dval = dctrlval;
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      CHECK_ZERO( lpi->messagehdlr, XPRSgetdblcontrol(lpi->xprslp, XPRS_OPTIMALITYTOL, &dctrlval) );
      *dval = dctrlval;
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      CHECK_ZERO( lpi->messagehdlr, XPRSgetdblcontrol(lpi->xprslp, XPRS_BARGAPSTOP, &dctrlval) );
      *dval = dctrlval;
      break;
   case SCIP_LPPAR_LPTILIM:
      CHECK_ZERO( lpi->messagehdlr, XPRSgetintcontrol(lpi->xprslp, XPRS_MAXTIME, &ictrlval) );
      /* ictrlval is the negative of the timelimit (see SCIPlpiSetRealpar) */
      *dval = (double) -ictrlval;
      break;
   case SCIP_LPPAR_MARKOWITZ:
      CHECK_ZERO( lpi->messagehdlr, XPRSgetdblcontrol(lpi->xprslp, XPRS_MARKOWITZTOL, &dctrlval) );
      *dval = dctrlval;
      break;
   case SCIP_LPPAR_OBJLIM:
      CHECK_ZERO( lpi->messagehdlr, XPRSgetdblcontrol(lpi->xprslp, XPRS_MIPABSCUTOFF, &dctrlval) );
      *dval = dctrlval;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiSetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("setting real parameter %d to %g\n", type, dval);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      /* Xpress does not pose any restriction on dval, its absolute value is used as tolerance.
       * For consistency we assert it to be strictly positive.
       */
      assert( dval > 0.0 );
      CHECK_ZERO( lpi->messagehdlr, XPRSsetdblcontrol(lpi->xprslp, XPRS_FEASTOL, dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      /* Xpress does not pose any restriction on dval,
       * however for consistency we assert it to be strictly positive.
       */
      assert( dval > 0.0 );
      CHECK_ZERO( lpi->messagehdlr, XPRSsetdblcontrol(lpi->xprslp, XPRS_OPTIMALITYTOL, dval) );
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      assert( dval >= 0.0 );
      /* Xpress poses no restriction on dval
       * However for consistency we assert it to be nonnegative.
       */
      CHECK_ZERO( lpi->messagehdlr, XPRSsetdblcontrol(lpi->xprslp, XPRS_BARGAPSTOP, dval) );
      break;
   case SCIP_LPPAR_LPTILIM:
   {
     int ival;

     /* From the Xpress documentation:
      * dval>0   If an integer solution has been found, stop MIP search after dval seconds,
      *   otherwise continue until an integer solution is finally found.
      * dval<0   Stop in LP or MIP search after dval seconds.
      * dval=0   No time limit
      */
      assert( dval > 0.0 );
      if( dval >= INT_MAX )
         ival = 0;
      else
         ival = (int) -floor(dval);

      CHECK_ZERO( lpi->messagehdlr, XPRSsetintcontrol(lpi->xprslp, XPRS_MAXTIME, ival) );
      break;
   }
   case SCIP_LPPAR_MARKOWITZ:
      /* no restriction on dval */
      CHECK_ZERO( lpi->messagehdlr, XPRSsetdblcontrol(lpi->xprslp, XPRS_MARKOWITZTOL, dval) );
      break;
   case SCIP_LPPAR_OBJLIM:
      /* no restriction on dval */
      CHECK_ZERO( lpi->messagehdlr, XPRSsetdblcontrol(lpi->xprslp, XPRS_MIPABSCUTOFF, dval) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

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


/**@name Numerical Methods
 *
 * @{
 */

/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   return XPRS_PLUSINFINITY;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   return (val >= XPRS_PLUSINFINITY);
}

/**@} */


/**@name File Interface Methods
 *
 * @{
 */

/** reads LP from a file
 *
 * The file extension defines the format. That can be lp or mps. Any given file name needs to have one of these two
 * extension. If not nothing is read and a SCIP_READERROR is returned.
 */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIP_RETCODE retcode = SCIP_OKAY;

   char* basename = NULL;
   char* compression = NULL;
   char* extension = NULL;
   char* filename = NULL;
   char* path = NULL;
   char* xpressfilename = NULL;

   int size;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(fname != NULL);

   SCIPdebugMessage("reading LP from file <%s>\n", fname);

   /* get the length of the file name */
   size = (int)strlen(fname)+1;

   /* check that the file name is longer than Xpress can handle */
   if (size > XPRS_MAXPROBNAMELENGTH)
     return SCIP_WRITEERROR;

   /* get char array for the file name we pass to Xpress */
   SCIP_ALLOC( BMSallocMemoryArray(&xpressfilename, size) );

   /* copy filename to be able to split it into its components */
   SCIP_ALLOC( BMSduplicateMemoryArray(&filename, fname, size) );

   /* get path, base file name, extension, and compression of the given file name */
   SCIPsplitFilename(filename, &path, &basename, &extension, &compression);

   /* construct file name without extension */
   if (path != NULL)
     (void) SCIPsnprintf(xpressfilename, size, "%s/%s", path, basename);
   else
     (void) SCIPsnprintf(xpressfilename, size, "%s", basename);

   /* check that the file name did not has a compression extension, has an lp or mps extension, and actually a base name */
   if (compression != NULL || extension == NULL || basename == NULL)
     retcode = SCIP_READERROR;
   if (strcasecmp(extension, "mps") == 0) {
     CHECK_ZERO( lpi->messagehdlr, XPRSreadprob(lpi->xprslp, xpressfilename, "") );
   }
   else if (strcasecmp(extension, "lp") == 0) {
     CHECK_ZERO( lpi->messagehdlr, XPRSreadprob(lpi->xprslp, xpressfilename, "l") );
   }
   else
     retcode = SCIP_READERROR;

   /* free array */
   BMSfreeMemoryArrayNull(&filename);
   BMSfreeMemoryArrayNull(&xpressfilename);

   return retcode;
}

/** writes LP to a file
 *
 * The file extension defines the format. That can be lp or mps. Any given file name needs to have one of these two
 * extension. If not nothing is written and a SCIP_WRITEERROR is returned.
 */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIP_RETCODE retcode = SCIP_OKAY;

   char* basename = NULL;
   char* compression = NULL;
   char* extension = NULL;
   char* filename = NULL;
   char* path = NULL;
   char* xpressfilename = NULL;

   int size;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(fname != NULL);

   SCIPdebugMessage("writing LP to file <%s>\n", fname);

   /* get the length of the file name */
   size = (int)strlen(fname)+1;

   /* check that the file name is longer than Xpress can handle */
   if (size > XPRS_MAXPROBNAMELENGTH)
     return SCIP_WRITEERROR;

   /* get char array for the file name we pass to Xpress */
   SCIP_ALLOC( BMSallocMemoryArray(&xpressfilename, size) );

   /* copy filename to be able to split it into its components */
   SCIP_ALLOC( BMSduplicateMemoryArray(&filename, fname, size) );

   /* get path, base file name, extension, and compression of the given file name */
   SCIPsplitFilename(filename, &path, &basename, &extension, &compression);

   /* construct file name without extension */
   if (path != NULL)
     (void) SCIPsnprintf(xpressfilename, size, "%s/%s", path, basename);
   else
     (void) SCIPsnprintf(xpressfilename, size, "%s", basename);

   /* check that the file name did not has a compression extension, has an lp or mps extension, and actually a base name */
   if (compression != NULL || extension == NULL || basename == NULL)
     retcode = SCIP_WRITEERROR;
   if (strcasecmp(extension, "mps") == 0) {
     CHECK_ZERO( lpi->messagehdlr, XPRSwriteprob(lpi->xprslp, xpressfilename, "p") );
   }
   else if (strcasecmp(extension, "lp") == 0) {
     CHECK_ZERO( lpi->messagehdlr, XPRSwriteprob(lpi->xprslp, xpressfilename, "lp") );
   }
   else
     retcode = SCIP_WRITEERROR;

   /* free array */
   BMSfreeMemoryArrayNull(&filename);
   BMSfreeMemoryArrayNull(&xpressfilename);

   return retcode;
}

/**@} */
