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

/**@file   lpi_grb.c
 * @ingroup LPIS
 * @brief  LP interface for Gurobi
 * @author Marc Pfetsch
 * @author Tobias Achterberg
 * @author Michael Winkler
 *
 * This LPI only works with Gurobi versions >= 7.0.2.
 *
 * @todo Try quad-precision and concurrent runs.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "gurobi_c.h"
#include "lpi/lpi.h"
#include "scip/pub_message.h"
#include "scip/pub_misc_sort.h"
#include "tinycthread/tinycthread.h"

#ifdef _WIN32
#define snprintf _snprintf
#endif

#if ( GRB_VERSION_MAJOR < 6 || ( GRB_VERSION_MAJOR == 7 && GRB_VERSION_MINOR == 0 && GRB_VERSION_TECHNICAL < 2 ) )
#error "The Gurobi intreface only works for Gurobi versions at least 7.0.2"
#endif

/* define infinity value of Gurobi */
#define GRB_INFBOUND 1e+20

/* macro for checking return codes of Gurobi */
#define CHECK_ZERO(messagehdlr, x) do { int _restat_;                   \
      if( (_restat_ = (x)) != 0 )                                       \
      {                                                                 \
         SCIPmessagePrintWarning((messagehdlr), "Gurobi error %d: %s\n", _restat_, GRBgeterrormsg(lpi->grbenv)); \
         return SCIP_LPERROR;                                           \
      }                                                                 \
   } while(0)

/* variant of macro for checking return codes of Gurobi */
#define CHECK_ZERO_STAR(messagehdlr, x) do { int _restat_;              \
      if( (_restat_ = (x)) != 0 )                                       \
      {                                                                 \
         SCIPmessagePrintWarning((messagehdlr), "Gurobi error %d: %s\n", _restat_, GRBgeterrormsg((*lpi)->grbenv)); \
         return SCIP_LPERROR;                                           \
      }                                                                 \
   } while(0)

#ifndef SVECTOR
#define SVECTOR GRBsvec
#endif

typedef unsigned int SCIP_DUALPACKET;        /**< storing bit pairs in packed format */
#define SCIP_DUALPACKETSIZE   (sizeof(SCIP_DUALPACKET)*4)   /**< each entry needs two bits of information */

typedef SCIP_DUALPACKET COLPACKET;           /**< each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /**< each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE


/* At several places we need to guarantee to have a factorization of an optimal basis and call the simplex to produce
 * it. In a numerical perfect world, this should need no iterations. However, due to numerical inaccuracies after
 * refactorization, it might be necessary to do a few extra pivot steps. */
#define GRB_REFACTORMAXITERS     50          /**< maximal number of iterations allowed for producing a refactorization of the basis */


/** number of Gurobi integer parameters that can be changed */
#define NUMINTPARAM 6

static const char* intparam[NUMINTPARAM] =
{
   GRB_INT_PAR_SCALEFLAG,
   GRB_INT_PAR_PRESOLVE,
   GRB_INT_PAR_SIMPLEXPRICING,
   GRB_INT_PAR_OUTPUTFLAG,
   GRB_INT_PAR_THREADS,
   GRB_INT_PAR_SEED
};

/** number of Gurobi double parameters that can be changed */
#define NUMDBLPARAM 7

static const char* dblparam[NUMDBLPARAM] =
{
   GRB_DBL_PAR_FEASIBILITYTOL,
   GRB_DBL_PAR_OPTIMALITYTOL,
   GRB_DBL_PAR_BARCONVTOL,
   GRB_DBL_PAR_CUTOFF,
   GRB_DBL_PAR_TIMELIMIT,
   GRB_DBL_PAR_ITERATIONLIMIT,
   GRB_DBL_PAR_MARKOWITZTOL
};

/** minimal values for double parameters */
static const double dblparammin[NUMDBLPARAM] =
{
   +1e-09,               /* GRB_DBL_PAR_FEASIBILITYTOL */
   +1e-09,               /* GRB_DBL_PAR_OPTIMALITYTOL */
   0.0,                  /* GRB_DBL_PAR_BARCONVTOL */
   -GRB_INFINITY,        /* GRB_DBL_PAR_CUTOFF */
   0,                    /* GRB_DBL_PAR_TIMELIMIT */
   0,                    /* GRB_DBL_PAR_ITERATIONLIMIT */
   1e-04                 /* GRB_DBL_PAR_MARKOWITZTOL */
};

/** Gurobi parameter settings */
struct GRBParam
{
   int                   intparval[NUMINTPARAM]; /**< integer parameter values */
   double                dblparval[NUMDBLPARAM]; /**< double parameter values */
};
typedef struct GRBParam GRBPARAM;


/** LP interface */
struct SCIP_LPi
{
   GRBenv*               grbenv;             /**< environment corresponding to model */
   GRBmodel*             grbmodel;           /**< Gurobi model pointer */
   int                   solstat;            /**< solution status of last optimization call */
   GRBPARAM              defparam;           /**< default parameter values */
   GRBPARAM              curparam;           /**< current parameter values stored in Gurobi LP */
   GRBPARAM              grbparam;           /**< current parameter values for this LP */
   char*                 senarray;           /**< array for storing row senses */
   SCIP_Real*            rhsarray;           /**< array for storing rhs values */
   SCIP_Real*            rngarray;           /**< array for storing range values */
   int*                  rngidxarray;        /**< array for storing the indices of ranged rows in sen/rhs/rngarray */
   SCIP_Real*            valarray;           /**< array for storing coefficient values */
   int*                  cstat;              /**< array for storing column basis status */
   int*                  rstat;              /**< array for storing row basis status */
   int*                  indarray;           /**< array for storing coefficient indices */
   int                   sidechgsize;        /**< size of senarray */
   int                   valsize;            /**< size of valarray and indarray */
   int                   cstatsize;          /**< size of cstat array */
   int                   rstatsize;          /**< size of rstat array */
   int                   iterations;         /**< number of iterations used in the last solving call */
   SCIP_Bool             solisbasic;         /**< is current LP solution a basic solution? */
   SCIP_Bool             fromscratch;        /**< should each solve be performed without previous basis state? */
   SCIP_PRICING          pricing;            /**< SCIP pricing setting  */
   SCIP_Real             conditionlimit;     /**< maximum condition number of LP basis counted as stable (-1.0: no limit) */
   SCIP_Bool             checkcondition;     /**< should condition number of LP basis be checked for stability? */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehdlr handler to printing messages, or NULL */
   int*                  rngrowmap;          /**< maps row id to rngrows array position, or -1 if not a ranged row
                                              *   (can be NULL, which means that no ranged rows exist) */
   int*                  rngrows;            /**< indices of ranged rows */
   SCIP_Real*            rngvals;            /**< range values of ranged rows */
   int                   rngrowmapsize;      /**< size of rngrowmap array */
   int                   nrngrows;           /**< number of ranged rows in the LP */
   int                   rngrowssize;        /**< size of rngrows and rngvals arrays */
   SCIP_Bool             rngvarsadded;       /**< did we add the range variables to the Gurobi model? */
};

/** LPi state stores basis information */
struct SCIP_LPiState
{
   int                   ncols;              /**< number of LP columns */
   int                   nrows;              /**< number of LP rows */
   int                   nrngrows;           /**< number of ranged rows in LP */
   COLPACKET*            packcstat;          /**< column basis status in compressed form */
   ROWPACKET*            packrstat;          /**< row basis status in compressed form */
};

/** LPi norms stores pricing norms */
struct SCIP_LPiNorms
{
   int                   ncols;              /**< number of columns for which dual norm is stored */
   int                   nrows;              /**< number of rows for which dual norm is stored */
   double*               colnorm;            /**< dual norms for columns */
   double*               rownorm;            /**< dual norms for rows */
};


#ifdef SCIP_THREADSAFE
   #if defined(_Thread_local)
      /* Use thread local environment in order to not create a new environment for each new LP. */
      static _Thread_local GRBenv*    reusegrbenv = NULL; /**< thread local Gurobi environment */
      static _Thread_local int        numlp = 0;          /**< number of open LP objects */
      #define SCIP_REUSEENV
   #endif
#else
   /* Global Gurobi environment in order to not create a new environment for each new LP. This is not thread safe. */
   static GRBenv*           reusegrbenv = NULL; /**< global Gurobi environment */
   static int               numlp = 0;          /**< number of open LP objects */
   #define SCIP_REUSEENV
#endif


/*
 * dynamic memory arrays
 */

/** resizes senarray to have at least num entries */
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
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rngidxarray, newsize) );
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

/** resizes rngrowmap array to have at least num entries */
static
SCIP_RETCODE ensureRngrowmapMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->rngrowmapsize )
   {
      int newsize;
      int r;

      newsize = MAX(2*lpi->rngrowmapsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rngrowmap, newsize) );
      for (r = lpi->rngrowmapsize; r < newsize; r++)
         lpi->rngrowmap[r] = -1;
      lpi->rngrowmapsize = newsize;
   }
   assert(num <= lpi->rngrowmapsize);

   return SCIP_OKAY;
}

/** resizes rngrows and rngvals arrays to have at least num entries */
static
SCIP_RETCODE ensureRngrowsMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->rngrowssize )
   {
      int newsize;

      newsize = MAX(2*lpi->rngrowssize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rngrows, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rngvals, newsize) );
      lpi->rngrowssize = newsize;
   }
   assert(num <= lpi->rngrowssize);

   return SCIP_OKAY;
}

/** stores current basis in internal arrays of LPI data structure */
static
SCIP_RETCODE getBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< whether basis information has successfully been obtained */
   )
{
   int ncols;
   int nrows;
   int res;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

   SCIPdebugMessage("getBase()\n");
   if ( success != NULL )
      *success = TRUE;

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* get unpacked basis information from Gurobi */
   res = GRBgetintattrarray(lpi->grbmodel, GRB_INT_ATTR_VBASIS, 0, ncols, lpi->cstat);
   if ( res == GRB_ERROR_DATA_NOT_AVAILABLE )
   {
      /* if the model is infeasible, Gurobi does not currently return basis information */
      if ( success != NULL )
         *success = FALSE;
      return SCIP_OKAY;
   }
   else if ( res != 0 )
   {
      SCIPerrorMessage("Gurobi error %d: %s\n", res, GRBgeterrormsg(lpi->grbenv));
      return SCIP_LPERROR;
   }

   res = GRBgetintattrarray(lpi->grbmodel, GRB_INT_ATTR_CBASIS, 0, nrows, lpi->rstat);
   if ( res == GRB_ERROR_DATA_NOT_AVAILABLE )
   {
      /* if the model is infeasible Gurobi does not currently return basis information */
      if ( success != NULL )
         *success = FALSE;
      return SCIP_OKAY;
   }
   else if ( res != 0 )
   {
      SCIPerrorMessage("Gurobi error %d: %s\n", res, GRBgeterrormsg(lpi->grbenv));
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** loads basis stored in internal arrays of LPI data structure into Gurobi */
static
SCIP_RETCODE setBase(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int ncols;
   int nrows;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );

   SCIPdebugMessage("setBase()\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );

   /* load basis information into Gurobi */
   CHECK_ZERO( lpi->messagehdlr, GRBsetintattrarray(lpi->grbmodel, GRB_INT_ATTR_VBASIS, 0, ncols, lpi->cstat) );
   CHECK_ZERO( lpi->messagehdlr, GRBsetintattrarray(lpi->grbmodel, GRB_INT_ATTR_CBASIS, 0, nrows, lpi->rstat) );

   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

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


/* The basis information for Gurobi is negative. So we cannot use the functions in bitencode.h/c. The functions below are a modified copy. */

/** encode a negated dual bit vector into packed format */
static
void SCIPencodeDualBitNeg(
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
   int i;
   int rest;
   int nfull;

   assert(inp != NULL || count == 0);
   assert(out != NULL || count == 0);
   assert(count >= 0);
   assert(SCIP_DUALPACKETSIZE == 16); /*lint !e506*/

   rest = count % (int)SCIP_DUALPACKETSIZE;
   nfull = count - rest;

   for( i = 0; i < nfull; i += (int)SCIP_DUALPACKETSIZE, inp += (int)SCIP_DUALPACKETSIZE ) /*lint !e679*/
   {
      assert(inp != NULL);
      assert(out != NULL);

#ifndef NDEBUG
      {
         unsigned int j;
         for( j = 0; j < SCIP_DUALPACKETSIZE; ++j )
            assert(0 <= -inp[j] && -inp[j] <= 3);
      }
#endif
      *out++ =
         mask[0][-inp[0]] | mask[1][-inp[1]] | mask[2][-inp[2]] | mask[3][-inp[3]]
         | mask[4][-inp[4]] | mask[5][-inp[5]] | mask[6][-inp[6]]
         | mask[7][-inp[7]] | mask[8][-inp[8]] | mask[9][-inp[9]]
         | mask[10][-inp[10]] | mask[11][-inp[11]] | mask[12][-inp[12]]
         | mask[13][-inp[13]] | mask[14][-inp[14]] | mask[15][-inp[15]];
   }

   if( rest > 0 )
   {
      SCIP_DUALPACKET m = (SCIP_DUALPACKET) 0u;

      assert(inp != NULL);
      assert(out != NULL);

      for( i = 0; i < rest; i++ )
         m |= mask[i][-inp[i]];  /*lint !e661*/
      *out = m;
   }
}

/** decode a packed dual bit vector into negated unpacked format */
static
void SCIPdecodeDualBitNeg(
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

      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      assert(m >> 2 == 0);
   }

   if( rest > 0 )
   {
      assert(inp != NULL);
      assert(out != NULL);

      m = *inp;
      for( i = 0; i < rest; i++ )
      {
         *out++ = -(int)(m & 3);
         m >>= 2;
      }
   }
}

/** store row and column basis status in a packed LPi state object */
static
void lpistatePack(
   SCIP_LPISTATE*        lpistate,           /**< pointer to LPi state data */
   const int*            cstat,              /**< basis status of columns in unpacked format */
   const int*            rstat               /**< basis status of rows in unpacked format */
   )
{
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

   SCIPencodeDualBitNeg(cstat, lpistate->packcstat, lpistate->ncols + lpistate->nrngrows);
   SCIPencodeDualBitNeg(rstat, lpistate->packrstat, lpistate->nrows);
}

/** unpacks row and column basis status from a packed LPi state object */
static
void lpistateUnpack(
   const SCIP_LPISTATE*  lpistate,           /**< pointer to LPi state data */
   int*                  cstat,              /**< buffer for storing basis status of columns in unpacked format */
   int*                  rstat               /**< buffer for storing basis status of rows in unpacked format */
   )
{
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

   SCIPdecodeDualBitNeg(lpistate->packcstat, cstat, lpistate->ncols + lpistate->nrngrows);
   SCIPdecodeDualBitNeg(lpistate->packrstat, rstat, lpistate->nrows);
}

/** creates LPi state information object */
static
SCIP_RETCODE lpistateCreate(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   ncols,              /**< number of columns to store */
   int                   nrows,              /**< number of rows to store */
   int                   nrngrows            /**< number of ranged rows */
   )
{
   assert(lpistate != NULL);
   assert(blkmem != NULL);
   assert(ncols >= 0);
   assert(nrows >= 0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, lpistate) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packcstat, colpacketNum(ncols + nrngrows)) );
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

   BMSfreeBlockMemoryArrayNull(blkmem, &(*lpistate)->packcstat, colpacketNum((*lpistate)->ncols + (*lpistate)->nrngrows));
   BMSfreeBlockMemoryArrayNull(blkmem, &(*lpistate)->packrstat, rowpacketNum((*lpistate)->nrows));
   BMSfreeBlockMemory(blkmem, lpistate);
}



/*
 * local methods
 */

/** gets all Gurobi parameters used in LPI */
static
SCIP_RETCODE getParameterValues(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   GRBPARAM*             grbparam            /**< Gurobi parameters */
   )
{
   int i;

   assert( lpi != NULL );
   assert( lpi->grbenv != NULL );
   assert( grbparam != NULL );

   SCIPdebugMessage("getParameterValues()\n");

   for( i = 0; i < NUMINTPARAM; ++i )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, intparam[i], &(grbparam->intparval[i])) );
   }
   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, dblparam[i], &(grbparam->dblparval[i])) );
   }

   return SCIP_OKAY;
}

/** in debug mode, checks validity of Gurobi parameters */
static
SCIP_RETCODE checkParameterValues(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
#ifndef NDEBUG
   GRBPARAM par;
   int i;

   SCIP_CALL( getParameterValues(lpi, &par) );
   for (i = 0; i < NUMINTPARAM; ++i)
      assert( lpi->curparam.intparval[i] == par.intparval[i] );
   for (i = 0; i < NUMDBLPARAM; ++i)
      assert(MAX(lpi->curparam.dblparval[i], dblparammin[i]) == par.dblparval[i]); /*lint !e777*/
#endif

   return SCIP_OKAY;
}

/** sets all Gurobi parameters used in LPI */
static
SCIP_RETCODE setParameterValues(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   GRBPARAM*             grbparam            /**< Gurobi parameters */
   )
{
   int i;

   assert( lpi != NULL );
   assert( lpi->grbenv != NULL );
   assert( grbparam != NULL );

   SCIPdebugMessage("setParameterValues()\n");

   for( i = 0; i < NUMINTPARAM; ++i )
   {
      if( lpi->curparam.intparval[i] != grbparam->intparval[i] )
      {
         SCIPdebugMessage("setting Gurobi int parameter %s from %d to %d\n",
            intparam[i], lpi->curparam.intparval[i], grbparam->intparval[i]);
         lpi->curparam.intparval[i] = grbparam->intparval[i];
         CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, intparam[i], lpi->curparam.intparval[i]) );
      }
   }
   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      if( lpi->curparam.dblparval[i] != grbparam->dblparval[i] ) /*lint !e777*/
      {
         SCIPdebugMessage("setting Gurobi dbl parameter %s from %g to %g\n",
            dblparam[i], lpi->curparam.dblparval[i], MAX(grbparam->dblparval[i], dblparammin[i]));
         lpi->curparam.dblparval[i] = MAX(grbparam->dblparval[i], dblparammin[i]);
         CHECK_ZERO( lpi->messagehdlr, GRBsetdblparam(lpi->grbenv, dblparam[i], lpi->curparam.dblparval[i]) );
      }
   }

   SCIP_CALL( checkParameterValues(lpi) );

   return SCIP_OKAY;
}

/** copies Gurobi parameters from source to dest */
static
void copyParameterValues(
   GRBPARAM*             dest,               /**< destination Gurobi parameters */
   const GRBPARAM*       source              /**< original Gurobi parameters */
   )
{
   int i;

   for( i = 0; i < NUMINTPARAM; ++i )
      dest->intparval[i] = source->intparval[i];
   for( i = 0; i < NUMDBLPARAM; ++i )
      dest->dblparval[i] = source->dblparval[i];
}

/** gets a single integer parameter value */
static
SCIP_RETCODE getIntParam(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           param,              /**< parameter name */
   int*                  p                   /**< value of parameter */
   )
{
   int i;

   assert( lpi != NULL );

   for( i = 0; i < NUMINTPARAM; ++i )
   {
      if( strcmp(intparam[i], param) == 0 )
      {
         *p = lpi->grbparam.intparval[i];
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("unknown Gurobi integer parameter <%s>.\n", param);
   return SCIP_LPERROR;
}

/** sets a single integer parameter value */
static
SCIP_RETCODE setIntParam(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           param,              /**< parameter name */
   int                   parval              /**< value of parameter */
   )
{
   int i;

   assert( lpi != NULL );

   for( i = 0; i < NUMINTPARAM; ++i )
   {
      if( strcmp(intparam[i], param) == 0 )
      {
         lpi->grbparam.intparval[i] = parval;
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("unknown Gurobi integer parameter <%s>.\n", param);
   return SCIP_LPERROR;
}

/** gets a single double parameter value */
static
SCIP_RETCODE getDblParam(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           param,              /**< parameter name */
   double*               p                   /**< value of parameter */
   )
{
   int i;

   assert(lpi != NULL);

   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      if( strcmp(dblparam[i], param) == 0 )
      {
         *p = lpi->grbparam.dblparval[i];
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("unknown Gurobi double parameter <%s>.\n", param);
   return SCIP_LPERROR;
}

/** sets a single double parameter value */
static
SCIP_RETCODE setDblParam(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           param,              /**< parameter name */
   double                parval              /**< value of parameter */
   )
{
   int i;

   assert( lpi != NULL );

   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      if( strcmp(dblparam[i], param) == 0 )
      {
         lpi->grbparam.dblparval[i] = parval;
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("unknown Gurobi double parameter <%s>.\n", param);
   return SCIP_LPERROR;
}

/** marks the current LP to be unsolved */
static
void invalidateSolution(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   lpi->solstat = -1;
}

/** converts SCIP's lhs/rhs pairs into Gurobi's sen/rhs */
static
SCIP_RETCODE convertSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand side vector */
   const SCIP_Real*      rhs,                /**< right hand side vector */
   int*                  rngcount            /**< number of ranged rows found */
   )
{
   int i;

   assert(lpi != NULL);
   assert(nrows >= 0);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(rngcount != NULL);

   /* convert lhs/rhs into sen/rhs */
   *rngcount = 0;
   for( i = 0; i < nrows; ++i )
   {
      assert(lhs[i] <= rhs[i]);

      if( lhs[i] == rhs[i] ) /*lint !e777*/
      {
         assert(-GRB_INFINITY < rhs[i] && rhs[i] < GRB_INFINITY);
         lpi->senarray[i] = GRB_EQUAL;
         lpi->rhsarray[i] = rhs[i];
         lpi->rngarray[i] = 0.0;
      }
      else if( lhs[i] <= -GRB_INFINITY )
      {
         /* this includes the case of free rows */
         assert(-GRB_INFINITY < rhs[i]);
         lpi->senarray[i] = GRB_LESS_EQUAL;
         lpi->rhsarray[i] = rhs[i];
      }
      else if( rhs[i] >= GRB_INFINITY )
      {
         assert(-GRB_INFINITY < lhs[i] && lhs[i] < GRB_INFINITY);
         lpi->senarray[i] = GRB_GREATER_EQUAL;
         lpi->rhsarray[i] = lhs[i];
      }
      else
      {
         /* we treat ranged rows as equations with an extra slack variable */
         assert(-GRB_INFINITY < lhs[i] && lhs[i] < GRB_INFINITY);
         assert(-GRB_INFINITY < rhs[i] && rhs[i] < GRB_INFINITY);
         lpi->senarray[i] = GRB_EQUAL;
         lpi->rhsarray[i] = lhs[i];
         lpi->rngarray[i] = rhs[i] - lhs[i];
         lpi->rngidxarray[(*rngcount)++] = i;
      }
   }
   return SCIP_OKAY;
}

/** converts Gurobi's sen/rhs pairs into SCIP's lhs/rhs pairs */
static
SCIP_RETCODE reconvertSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhs,                /**< buffer to store the left hand side vector, or NULL */
   SCIP_Real*            rhs                 /**< buffer to store the right hand side vector, or NULL */
   )
{
   int nrows;
   int i;

   nrows = lastrow-firstrow+1;

   assert(lpi != NULL);
   assert(nrows >= 0);

   for (i = 0; i < nrows; ++i)
   {
      switch( lpi->senarray[i] )
      {
      case GRB_EQUAL:
         if ( lhs != NULL )
            lhs[i] = lpi->rhsarray[i];
         if ( rhs != NULL )
         {
            int row;

            rhs[i] = lpi->rhsarray[i];
            row = firstrow+i;
            if ( lpi->rngrowmap != NULL && lpi->rngrowmap[row] >= 0 )
            {
               assert(lpi->rngrowmap[row] < lpi->nrngrows);
               rhs[i] += lpi->rngvals[lpi->rngrowmap[row]];
            }
         }
         break;

      case GRB_LESS_EQUAL:
         if ( lhs != NULL )
            lhs[i] = -GRB_INFINITY;
         if ( rhs != NULL )
            rhs[i] = lpi->rhsarray[i];
         break;

      case GRB_GREATER_EQUAL:
         if ( lhs != NULL )
            lhs[i] = lpi->rhsarray[i];
         if ( rhs != NULL )
            rhs[i] = GRB_INFINITY;
         break;

      default:
         SCIPerrorMessage("invalid row sense\n");
         SCIPABORT();
         return SCIP_LPERROR; /*lint !e527*/
      }
      assert(lhs == NULL || rhs == NULL || lhs[i] <= rhs[i]);
   }
   return SCIP_OKAY;
}

/** after restoring old LP data, need to resolve the LP to be able to retrieve correct information */
static
SCIP_RETCODE restoreLPData(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );

   /* set dual simplex */
   CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL) );
   CHECK_ZERO( lpi->messagehdlr, GRBoptimize(lpi->grbmodel) );

#ifndef NDEBUG
   {
      double cnt;

      /* modifying the LP, restoring the old LP, and loading the old basis is not enough for Gurobi to be able to return
       * the basis -> we have to resolve the LP;
       *
       * In a numerical perfect world, GRB_REFACTORMAXITERS below should be zero. However, due to numerical inaccuracies
       * after refactorization, it might be necessary to do a few extra pivot steps.
       */
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
      if ( cnt > (double) GRB_REFACTORMAXITERS )
         SCIPmessagePrintWarning(lpi->messagehdlr, "Gurobi needed %d iterations to restore optimal basis.\n", (int) cnt);
   }
#endif

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** verifies in debug mode that ranged row information is consistent */
static
void checkRangeInfo(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi->rngrowssize >= lpi->nrngrows);

   if ( lpi->nrngrows > 0 )
   {
      int nrngrows = 0;
      int nrows;
      int i;

      assert(lpi->rngrowmap != NULL);
      assert(lpi->rngrows != NULL);
      assert(lpi->rngvals != NULL);

      SCIP_CALL_ABORT( SCIPlpiGetNRows(lpi, &nrows) );

      assert(lpi->rngrowmapsize >= nrows);

      for (i = 0; i < nrows; i++)
      {
         int rngrow;

         rngrow = lpi->rngrowmap[i];
         assert(-1 <= rngrow && rngrow < lpi->nrngrows);
         if ( rngrow >= 0 )
         {
            assert(lpi->rngrows[rngrow] == i);
            assert(lpi->rngvals[rngrow] > 0.0);
            nrngrows++;
         }
      }
      assert(lpi->nrngrows == nrngrows);
   }
}
#else
#define checkRangeInfo(lpi) /**/
#endif

/** adds range variables to Gurobi LP */
static
SCIP_RETCODE addRangeVars(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int i;

   assert(!lpi->rngvarsadded);
   assert(lpi->nrngrows > 0);
   assert(lpi->rngrowmap != NULL);
   assert(lpi->rngrows != NULL);

   for (i = 0; i < lpi->nrngrows; i++)
   {
      double coeff = -1.0;
      int row;

      row = lpi->rngrows[i];

      CHECK_ZERO( lpi->messagehdlr, GRBaddvar(lpi->grbmodel, 1, &row, &coeff, 0.0, 0.0, lpi->rngvals[i], GRB_CONTINUOUS, NULL) );
   }

   /* flush model changes */
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   lpi->rngvarsadded = TRUE;

   return SCIP_OKAY;
}

/** deletes range variables from Gurobi LP */
static
SCIP_RETCODE delRangeVars(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int* which;
   int ncols;
   int i;

   assert(lpi->rngvarsadded);
   assert(lpi->nrngrows > 0);

   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );

   /* Gurobi can't delete a range of columns, we have to set up an index array */
   SCIP_ALLOC( BMSallocMemoryArray(&which, lpi->nrngrows) );
   for (i = 0; i < lpi->nrngrows; i++)
      which[i] = ncols+i;

   CHECK_ZERO( lpi->messagehdlr, GRBdelvars(lpi->grbmodel, lpi->nrngrows, which) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   BMSfreeMemoryArray( &which );

   lpi->rngvarsadded = FALSE;

   return SCIP_OKAY;
}

/** clear ranged row information */
static
void clearRangeInfo(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(!lpi->rngvarsadded);

   BMSfreeMemoryArrayNull(&lpi->rngrowmap);
   BMSfreeMemoryArrayNull(&lpi->rngrows);
   BMSfreeMemoryArrayNull(&lpi->rngvals);

   lpi->nrngrows = 0;
   lpi->rngrowssize = 0;
   lpi->rngrowmapsize = 0;
}

/** creates or updates maps for ranged rows after new rows have been added */
static
SCIP_RETCODE addRangeInfo(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   rngcount,           /**< number of ranged rows added */
   int                   firstrow            /**< index of first row that was added */
   )
{
   int ncols;
   int nrows;
   int r;
   int i;

   assert( lpi != NULL );

   /* get rid of range variables */
   if ( lpi->rngvarsadded )
   {
      SCIP_CALL( delRangeVars(lpi) );
   }
   assert( !lpi->rngvarsadded );

   /* query problem size in terms of SCIP's view */
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );

   /* set up and extend rngrowmap array */
   SCIP_CALL( ensureRngrowmapMem(lpi, nrows) );
   for (r = firstrow; r < nrows; r++)
      lpi->rngrowmap[r] = -1;

   /* extend rngrows and rngvals arrays */
   SCIP_CALL( ensureRngrowsMem(lpi, lpi->nrngrows + rngcount) );

   /* update maps for ranged rows */
   for (i = 0; i < rngcount; i++)
   {
      int pos;
      int row;

      pos = lpi->rngidxarray[i];
      row  = firstrow + pos;

      lpi->rngrowmap[row] = lpi->nrngrows;
      lpi->rngrows[lpi->nrngrows] = row;
      lpi->rngvals[lpi->nrngrows] = lpi->rngarray[pos];
      lpi->nrngrows++;
   }

   return SCIP_OKAY;
}



/*
 * LP Interface Methods
 */


/*
 * Miscellaneous Methods
 */

static const char grbname[] = {'G', 'u', 'r', 'o', 'b', 'i', ' ',
#if GRB_VERSION_MAJOR < 10
   GRB_VERSION_MAJOR + '0',
#else
   (GRB_VERSION_MAJOR/10) + '0', (GRB_VERSION_MAJOR%10) + '0',
#endif
   '.', GRB_VERSION_MINOR + '0', '.', GRB_VERSION_TECHNICAL + '0'};

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   return grbname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "Linear Programming Solver developed by Gurobi Optimization (www.gurobi.com)";
}

/** gets pointer for LP solver - use only with great care
 *
 *  Here we return the pointer to the model.
 */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   return (void*) lpi->grbmodel;
}

/** pass integrality information to LP solver */
SCIP_RETCODE SCIPlpiSetIntegralityInformation(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   ncols,              /**< length of integrality array */
   int*                  intInfo             /**< integrality array (0: continuous, 1: integer). May be NULL iff ncols is 0. */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( ncols >= 0 );
   assert( ncols == 0 || intInfo != NULL );

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
   assert(sizeof(SCIP_Real) == sizeof(double)); /*lint !e506*/ /* Gurobi only works with doubles as floating points */
   assert(sizeof(SCIP_Bool) == sizeof(int));    /*lint !e506*/ /* Gurobi only works with ints as bools */
   assert(lpi != NULL);
   assert(name != NULL);
#ifdef SCIP_REUSEENV
   assert(numlp >= 0);
#endif

   SCIPdebugMessage("SCIPlpiCreate()\n");

   /* create empty LPI */
   SCIP_ALLOC( BMSallocMemory(lpi) );

   /* create environment */
#ifdef SCIP_REUSEENV
   /* temporarily set environment for error messages (might be NULL) */
   (*lpi)->grbenv = reusegrbenv;

   /* Try to reuse Gurobi environment (either thread local or not being thread safe). */
   if ( reusegrbenv == NULL )
   {
      int restat;

      assert( numlp == 0 );

      /* create evironment */
      restat = GRBloadenv(&reusegrbenv, NULL);
      if ( restat != 0 )
      {
         SCIPmessagePrintWarning(messagehdlr, "Gurobi error %d: Something went wrong with creating the environment.\n", restat);
         return SCIP_LPERROR;
      }

      /* turn off output for all models */
      CHECK_ZERO_STAR( messagehdlr, GRBsetintparam(reusegrbenv, GRB_INT_PAR_OUTPUTFLAG, 0) );

      /* turn on that basis information for infeasible and unbounded models is available */
      CHECK_ZERO_STAR( messagehdlr, GRBsetintparam(reusegrbenv, GRB_INT_PAR_INFUNBDINFO, 1) );
   }

   /* create empty model */
   CHECK_ZERO_STAR( messagehdlr, GRBnewmodel(reusegrbenv, &(*lpi)->grbmodel, name, 0, NULL, NULL, NULL, NULL, NULL) );

   /* replace by local copy of environment */
   (*lpi)->grbenv = GRBgetenv((*lpi)->grbmodel);
   ++numlp;

#else

   /* Create new environment for each new instaniation; note that this involves additional work and
    * uses a new license for each new instantiation. */
   CHECK_ZERO_STAR( messagehdlr, GRBloadenv(&(*lpi)->grbenv, NULL) );

   /* turn off output for all models */
   CHECK_ZERO_STAR( messagehdlr, GRBsetintparam((*lpi)->grbenv, GRB_INT_PAR_OUTPUTFLAG, 0) );

   /* turn on that basis information for infeasible and unbounded models is available */
   CHECK_ZERO_STAR( messagehdlr, GRBsetintparam((*lpi)->grbenv, GRB_INT_PAR_INFUNBDINFO, 1) );

   /* create empty model */
   CHECK_ZERO_STAR( messagehdlr, GRBnewmodel((*lpi)->grbenv, &(*lpi)->grbmodel, name, 0, NULL, NULL, NULL, NULL, NULL) );

#endif
   assert( (*lpi)->grbenv != NULL );

   (*lpi)->senarray = NULL;
   (*lpi)->rhsarray = NULL;
   (*lpi)->rngarray = NULL;
   (*lpi)->rngidxarray = NULL;
   (*lpi)->valarray = NULL;
   (*lpi)->cstat = NULL;
   (*lpi)->rstat = NULL;
   (*lpi)->indarray = NULL;
   (*lpi)->rngrowmap = NULL;
   (*lpi)->rngrows = NULL;
   (*lpi)->rngvals = NULL;
   (*lpi)->sidechgsize = 0;
   (*lpi)->valsize = 0;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->rngrowmapsize = 0;
   (*lpi)->nrngrows = 0;
   (*lpi)->rngrowssize = 0;
   (*lpi)->rngvarsadded = FALSE;
   (*lpi)->iterations = 0;
   (*lpi)->solisbasic = FALSE;
   (*lpi)->fromscratch = FALSE;
   (*lpi)->conditionlimit = -1.0;
   (*lpi)->checkcondition = FALSE;
   (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*lpi)->messagehdlr = messagehdlr;
   invalidateSolution(*lpi);

   /* get default parameter values */
   SCIP_CALL( getParameterValues((*lpi), &((*lpi)->defparam)) );
   copyParameterValues(&((*lpi)->curparam), &((*lpi)->defparam));
   copyParameterValues(&((*lpi)->grbparam), &((*lpi)->defparam));

   /* set objective sense */
   SCIP_CALL( SCIPlpiChgObjsen(*lpi, objsen) );

   /* set default pricing */
   SCIP_CALL( SCIPlpiSetIntpar(*lpi, SCIP_LPPAR_PRICING, (int) (*lpi)->pricing) );

   checkRangeInfo(*lpi);

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(*lpi != NULL);
   assert((*lpi)->grbenv != NULL);

   SCIPdebugMessage("SCIPlpiFree()\n");

   /* free model */
   CHECK_ZERO_STAR( (*lpi)->messagehdlr, GRBfreemodel((*lpi)->grbmodel) );

   /* free memory */
   BMSfreeMemoryArrayNull(&(*lpi)->senarray);
   BMSfreeMemoryArrayNull(&(*lpi)->rhsarray);
   BMSfreeMemoryArrayNull(&(*lpi)->rngarray);
   BMSfreeMemoryArrayNull(&(*lpi)->rngidxarray);
   BMSfreeMemoryArrayNull(&(*lpi)->cstat);
   BMSfreeMemoryArrayNull(&(*lpi)->rstat);
   BMSfreeMemoryArrayNull(&(*lpi)->rngrowmap);
   BMSfreeMemoryArrayNull(&(*lpi)->rngrows);
   BMSfreeMemoryArrayNull(&(*lpi)->rngvals);
   BMSfreeMemoryArrayNull(&(*lpi)->indarray);
   BMSfreeMemoryArrayNull(&(*lpi)->valarray);

   /* free environment */
#ifdef SCIP_REUSEENV
   --numlp;
   if( numlp == 0 )
   {
      /* free reused environment */
      GRBfreeenv(reusegrbenv);
      reusegrbenv = NULL;
   }
#else
   /* free local environment */
   GRBfreeenv((*lpi)->grbenv);
#endif

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
   int grbobjsen;
   int* cnt;
   int rngcount;
   int c;

#ifndef NDEBUG
   {
      int j;
      for( j = 0; j < nnonz; j++ )
         assert( val[j] != 0 );
   }
#endif

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->grbenv != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(beg != NULL);
   assert(ind != NULL);
   assert(val != NULL);

   assert(objsen == SCIP_OBJSEN_MAXIMIZE || objsen == SCIP_OBJSEN_MINIMIZE);

   SCIPdebugMessage("loading LP in column format into Gurobi: %d cols, %d rows\n", ncols, nrows);

   invalidateSolution(lpi);

   SCIP_CALL( ensureSidechgMem(lpi, nrows) );

   /* convert objective sense */
   grbobjsen = (objsen == SCIP_OBJSEN_MINIMIZE) ? GRB_MINIMIZE : GRB_MAXIMIZE;

   /* convert lhs/rhs into sen/rhs/range tuples */
   SCIP_CALL( convertSides(lpi, nrows, lhs, rhs, &rngcount) );

   /* calculate column lengths */
   SCIP_ALLOC( BMSallocMemoryArray(&cnt, ncols) );
   for( c = 0; c < ncols-1; ++c )
   {
      cnt[c] = beg[c+1] - beg[c];
      assert(cnt[c] >= 0);
   }
   cnt[ncols-1] = nnonz - beg[ncols-1];
   assert(cnt[ncols-1] >= 0);

   /* load model - all variables are continuous */
   CHECK_ZERO( lpi->messagehdlr, GRBloadmodel(lpi->grbenv, &(lpi->grbmodel), NULL, ncols, nrows, grbobjsen, 0.0, (SCIP_Real*)obj,
         lpi->senarray, lpi->rhsarray, (int*)beg, cnt, (int*)ind, (SCIP_Real*)val, (SCIP_Real*)lb, (SCIP_Real*)ub, NULL, colnames, rownames) );

   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   /* free temporary memory */
   BMSfreeMemoryArray(&cnt);

   /* update maps for ranged rows */
   if ( rngcount > 0 )
   {
      SCIP_CALL( addRangeInfo(lpi, rngcount, 0) );
   }

#ifndef NDEBUG
   {
      int temp;

      SCIP_CALL( SCIPlpiGetNCols(lpi, &temp) );
      assert(temp == ncols);

      SCIP_CALL( SCIPlpiGetNRows(lpi, &temp) );
      assert(temp == nrows);

      SCIP_CALL( SCIPlpiGetNNonz(lpi, &temp) );
      assert(temp == nnonz);
   }
#endif

   checkRangeInfo(lpi);

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
{

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);
   assert(nnonz >= 0);
   assert(ncols >= 0);

   SCIPdebugMessage("adding %d columns with %d nonzeros to Gurobi\n", ncols, nnonz);

   invalidateSolution(lpi);

#ifndef NDEBUG
   if ( nnonz > 0 )
   {
      /* perform check that no new rows are added - this is forbidden */
      int nrows;
      int j;

      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
      for (j = 0; j < nnonz; ++j)
      {
         assert( 0 <= ind[j] && ind[j] < nrows );
         assert( val[j] != 0.0 );
      }
   }
#endif

   /* delete range variables from Gurobi LP, so that structural variables always come first */
   if ( lpi->nrngrows > 0 && lpi->rngvarsadded )
   {
      /**@todo Save and restore basis - currently, the basis is destroyed if we discard (and later re-add) range variables */
      SCIP_CALL( delRangeVars(lpi) );
   }

   /* add columns - all new variables are continuous */
   CHECK_ZERO( lpi->messagehdlr, GRBaddvars(lpi->grbmodel, ncols, nnonz, (int*)beg, (int*)ind, (SCIP_Real*)val,
      (SCIP_Real*)obj, (SCIP_Real*)lb, (SCIP_Real*)ub, NULL, colnames) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   checkRangeInfo(lpi);

   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   int ndelcols;
   int* which;
   int j;

   ndelcols = lastcol-firstcol+1;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
#ifndef NDEBUG
   {
      int temp;

      SCIP_CALL( SCIPlpiGetNCols(lpi, &temp) );
      assert(0 <= firstcol && firstcol <= lastcol && lastcol < temp);
   }
#endif

   SCIPdebugMessage("deleting %d columns from Gurobi\n", lastcol - firstcol + 1);

   invalidateSolution(lpi);

   /* Gurobi can't delete a range of columns, we have to set up an index array */
   SCIP_ALLOC( BMSallocMemoryArray(&which, ndelcols) );
   for( j = firstcol; j <= lastcol; ++j )
      which[j - firstcol] = j;

   CHECK_ZERO( lpi->messagehdlr, GRBdelvars(lpi->grbmodel, ndelcols, which) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   BMSfreeMemoryArray( &which );

   checkRangeInfo(lpi);

   return SCIP_OKAY;
}

/** deletes columns from LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   int* which;
   int ncols;
   int num;
   int j;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(dstat != NULL);

   SCIPdebugMessage("deleting a column set from Gurobi\n");

   invalidateSolution(lpi);

   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );

   /* Gurobi can't delete a set of marked columns, we have to set up an index array */
   SCIP_ALLOC( BMSallocMemoryArray(&which, ncols) );
   num = 0;
   for( j = 0; j < ncols; ++j )
   {
      if( dstat[j] )
         which[num++] = j;
   }

   CHECK_ZERO( lpi->messagehdlr, GRBdelvars(lpi->grbmodel, num, which) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   /* update dstat */
   num = 0;
   for( j = 0; j < ncols; ++j )
   {
      if( dstat[j] )
      {
         dstat[j] = -1;
         ++num;
      }
      else
         dstat[j] = j - num;
   }

   BMSfreeMemoryArray( &which );

   checkRangeInfo(lpi);

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
   int rngcount;
   int oldnrows = -1;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert((lpi->nrngrows > 0) == (lpi->rngrowmap != NULL));
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);

   SCIPdebugMessage("adding %d rows with %d nonzeros to Gurobi\n", nrows, nnonz);

   invalidateSolution(lpi);

#ifndef NDEBUG
   if ( nnonz > 0 )
   {
      /* perform check that no new cols are added - this is forbidden */
      int ncols;
      int j;

      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
      for (j = 0; j < nnonz; ++j) {
         assert( 0 <= ind[j] && ind[j] < ncols );
         assert( val[j] != 0.0 );
      }
   }
#endif

   SCIP_CALL( ensureSidechgMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   SCIP_CALL( convertSides(lpi, nrows, lhs, rhs, &rngcount) );
   if ( lpi->nrngrows > 0 || rngcount > 0 )
   {
      SCIP_CALL( SCIPlpiGetNRows(lpi, &oldnrows) );
   }

   /* add rows to LP */
   CHECK_ZERO( lpi->messagehdlr, GRBaddconstrs(lpi->grbmodel, nrows, nnonz, (int*)beg, (int*)ind, (SCIP_Real*)val, lpi->senarray, lpi->rhsarray, rownames) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   /* update maps for ranged rows */
   if ( rngcount > 0 )
   {
      SCIP_CALL( addRangeInfo(lpi, rngcount, oldnrows) );
   }
   else if ( lpi->nrngrows > 0 )
   {
      int r;

      /* extend existing rngrowmap array */
      assert(lpi->rngrowmap != NULL);
      assert(lpi->rngrows != NULL);
      SCIP_CALL( ensureRngrowmapMem(lpi, oldnrows+nrows) );
      for (r = oldnrows; r < oldnrows+nrows; r++)
         lpi->rngrowmap[r] = -1;
   }

   checkRangeInfo(lpi);

   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   int ndelrows;
   int* which;
   int i;

   ndelrows = lastrow-firstrow+1;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
#ifndef NDEBUG
   {
      int nrows;
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
      assert(0 <= firstrow && firstrow <= lastrow && lastrow < nrows);
   }
#endif

   SCIPdebugMessage("deleting %d rows from Gurobi\n", ndelrows);

   invalidateSolution(lpi);

   /* Gurobi can't delete a range of rows, we have to set up an index array */
   SCIP_ALLOC( BMSallocMemoryArray(&which, ndelrows) );
   for( i = firstrow; i <= lastrow; ++i )
      which[i - firstrow] = i;

   CHECK_ZERO( lpi->messagehdlr, GRBdelconstrs(lpi->grbmodel, ndelrows, which) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   BMSfreeMemoryArray( &which );

   /* update ranged row info */
   if ( lpi->nrngrows > 0 )
   {
      int nrngrows;
      int nrows;

      assert(lpi->rngrowmap != NULL);
      assert(lpi->rngrows != NULL);

      /* find first ranged row that has been deleted */
      for (i = 0; i < lpi->nrngrows; i++)
      {
         if ( lpi->rngrows[i] >= firstrow )
            break;
      }
      nrngrows = i;

      /* skip all deleted ranged rows */
      for (; i < lpi->nrngrows; i++)
      {
         if ( lpi->rngrows[i] > lastrow )
            break;
      }

      /* move remaining ranged rows to the front */
      for (; i < lpi->nrngrows; i++)
      {
         int oldrow = lpi->rngrows[i];
         lpi->rngrowmap[oldrow] = nrngrows; /* store at old place for now */
         lpi->rngrows[nrngrows] = oldrow - ndelrows;
         lpi->rngvals[nrngrows] = lpi->rngvals[i];
         nrngrows++;
      }

      if ( nrngrows < lpi->nrngrows && lpi->rngvarsadded )
      {
         /* For simplicity, just delete all range variables from Gurobi LP - it would suffice to only delete those
          * corresponding to deleted ranged rows, but this should not matter much. */
         SCIP_CALL( delRangeVars(lpi) );
      }

      lpi->nrngrows = nrngrows;

      if ( nrngrows == 0 )
         clearRangeInfo(lpi);
      else
      {
         /* move rngrowmap entries */
         SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
         for (i = firstrow; i < nrows; i++)
         {
            lpi->rngrowmap[i] = lpi->rngrowmap[i+ndelrows];
            assert(-1 <= lpi->rngrowmap[i] && lpi->rngrowmap[i] < lpi->nrngrows);
         }
      }
   }

   checkRangeInfo(lpi);

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
   int i;
   int num = 0;
   int nrows;
   int* which;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("deleting a row set from Gurobi\n");

   invalidateSolution(lpi);

   /* Gurobi can't delete a range of rows, we have to set up an index array */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&which, nrows) );
   for( i = 0; i < nrows; ++i )
   {
      if( dstat[i] )
         which[num++] = i;
   }
   CHECK_ZERO( lpi->messagehdlr, GRBdelconstrs(lpi->grbmodel, num, which) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   /* update dstat */
   num = 0;
   for( i = 0; i < nrows; ++i )
   {
      if( dstat[i] )
      {
         dstat[i] = -1;
         ++num;
      }
      else
         dstat[i] = i - num;
   }

   /* update ranged row info */
   if ( lpi->nrngrows > 0 )
   {
      int nrngrows = 0;

      assert(lpi->rngrowmap != NULL);
      assert(lpi->rngrows != NULL);

      for (i = 0; i < lpi->nrngrows; i++)
      {
         int oldrow = lpi->rngrows[i];
         int newrow = dstat[oldrow];
         if ( newrow >= 0 )
         {
            lpi->rngrowmap[oldrow] = nrngrows; /* store at old place for now */
            lpi->rngrows[nrngrows] = newrow;
            lpi->rngvals[nrngrows] = lpi->rngvals[i];
            nrngrows++;
         }
      }

      if ( nrngrows < lpi->nrngrows && lpi->rngvarsadded )
      {
         /* for simplicity, just delete all range variables from
          * Gurobi LP - it would suffice to only delete those
          * corresponding to deleted ranged rows, but this should
          * not matter much
          */
         SCIP_CALL( delRangeVars(lpi) );
      }

      lpi->nrngrows = nrngrows;

      if ( nrngrows == 0 )
         clearRangeInfo(lpi);
      else
      {
         /* move rngrowmap entries */
         for (i = 0; i < nrows; i++)
         {
            int newrow = dstat[i];
            assert(newrow <= i);
            if ( newrow >= 0 )
            {
               lpi->rngrowmap[newrow] = lpi->rngrowmap[i];
               assert(-1 <= lpi->rngrowmap[newrow] && lpi->rngrowmap[newrow] < lpi->nrngrows);
            }
         }
      }
   }

   BMSfreeMemoryArray( &which );

   checkRangeInfo(lpi);

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int nrows;
   int ncols;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

   SCIPdebugMessage("clearing Gurobi LP\n");

   invalidateSolution(lpi);

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );

   if ( nrows >= 1 )
   {
      SCIP_CALL( SCIPlpiDelRows(lpi, 0, nrows-1) );
   }
   if ( ncols >= 1 )
   {
      SCIP_CALL( SCIPlpiDelCols(lpi, 0, ncols-1) );
   }

#ifdef SCIP_DISABLED_CODE
   /* the following seems to be slower */
   CHECK_ZERO( lpi->messagehdlr, GRBfreemodel(lpi->grbmodel) );
   CHECK_ZERO( lpi->messagehdlr, GRBnewmodel(lpi->grbenv, &(lpi->grbmodel), "", 0, NULL, NULL, NULL, NULL, NULL) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );
#endif

   /* clear ranged row info */
   clearRangeInfo(lpi);

   checkRangeInfo(lpi);

   return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiChgBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices or NULL if ncols is zero */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds or NULL if ncols is zero */
   const SCIP_Real*      ub                  /**< values for the new upper bounds or NULL if ncols is zero*/
   )
{
   int i;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(ncols == 0 || (ind != NULL && lb != NULL && ub != NULL));

   SCIPdebugMessage("changing %d bounds in Gurobi\n", ncols);
   if( ncols <= 0 )
      return SCIP_OKAY;

   for (i = 0; i < ncols; ++i)
   {
      SCIPdebugPrintf("  col %d: [%g,%g]\n", ind[i], lb[i], ub[i]);

      if ( SCIPlpiIsInfinity(lpi, lb[i]) )
      {
         SCIPerrorMessage("LP Error: fixing lower bound for variable %d to infinity.\n", ind[i]);
         return SCIP_LPERROR;
      }
      if ( SCIPlpiIsInfinity(lpi, -ub[i]) )
      {
         SCIPerrorMessage("LP Error: fixing upper bound for variable %d to -infinity.\n", ind[i]);
         return SCIP_LPERROR;
      }
   }

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrlist(lpi->grbmodel, GRB_DBL_ATTR_LB, ncols, (int*)ind, (SCIP_Real*)lb) );
   CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrlist(lpi->grbmodel, GRB_DBL_ATTR_UB, ncols, (int*)ind, (SCIP_Real*)ub) );

   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   checkRangeInfo(lpi);

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
   int rngcount;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(ind != NULL);

   SCIPdebugMessage("changing %d sides in Gurobi\n", nrows);
   if( nrows <= 0)
      return SCIP_OKAY;

   invalidateSolution(lpi);

   /* convert lhs/rhs into sen/rhs/range tuples */
   SCIP_CALL( ensureSidechgMem(lpi, nrows) );
   SCIP_CALL( convertSides(lpi, nrows, lhs, rhs, &rngcount) );

   /* change row sides */
   CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrlist(lpi->grbmodel, GRB_DBL_ATTR_RHS, nrows, (int*)ind, lpi->rhsarray) );
   CHECK_ZERO( lpi->messagehdlr, GRBsetcharattrlist(lpi->grbmodel, GRB_CHAR_ATTR_SENSE, nrows, (int*)ind, lpi->senarray) );

   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   /* update ranged row info */
   if ( rngcount > 0 || lpi->nrngrows > 0 )
   {
      int modified = 0;
      int nnewrngrows = 0;
      int ntotrows;
      int ncols;
      int i;

      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
      SCIP_CALL( SCIPlpiGetNRows(lpi, &ntotrows) );

      SCIP_CALL( ensureRngrowmapMem(lpi, ntotrows) );

      for (i = 0; i < nrows; i++)
      {
         int rngrowidx;
         int row;

         row = ind[i];
         rngrowidx = lpi->rngrowmap[row];

         assert(-1 <= rngrowidx && rngrowidx < lpi->nrngrows);

         if ( lpi->senarray[i] == GRB_EQUAL && lpi->rngarray[i] > 0.0 )
         {
            /* row is (now) a ranged row */
            if ( rngrowidx >= 0 )
            {
               /* row was already a ranged row: just update rngval and ub of associated column */
               lpi->rngvals[rngrowidx] = lpi->rngarray[i];
               if ( !modified && lpi->rngvarsadded )
               {
                  CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_UB, ncols+rngrowidx, lpi->rngvals[rngrowidx]) );
               }
            }
            else
            {
               /* row was not ranged before: we need to reset range variables */
               modified = 1;

               /* for now, add row to end of rngrows/rngvals arrays */
               SCIP_CALL( ensureRngrowsMem(lpi, lpi->nrngrows + nnewrngrows + 1) );
               lpi->rngrowmap[row] = lpi->nrngrows + nnewrngrows;
               lpi->rngrows[lpi->nrngrows + nnewrngrows] = row;
               lpi->rngvals[lpi->nrngrows + nnewrngrows] = lpi->rngarray[i];
               nnewrngrows++;
            }
         }
         else
         {
            /* row is not (no longer) a ranged row */
            if ( rngrowidx >= 0 )
            {
               /* row was a ranged row before: we need to reset range variables */
               modified = 1;
               lpi->rngrowmap[row] = -1;
            }
         }
      }

      if ( modified )
      {
         int nrngrows = 0;

         /* the range status of at least one row changed: discard range variables */
         if ( lpi->rngvarsadded )
         {
            /**@todo Save and restore basis - currently, the basis is destroyed if we discard (and later re-add) range variables */
            SCIP_CALL( delRangeVars(lpi) );
         }
         assert(!lpi->rngvarsadded);

         if ( nnewrngrows > 0 )
         {
            /* integrate new ranged rows into arrays */
            lpi->nrngrows += nnewrngrows;
            SCIPsortIntReal(lpi->rngrows, lpi->rngvals, lpi->nrngrows);
         }

         /* update rngrowmap and discard rows that are no longer ranged */
         for (i = 0; i < lpi->nrngrows; i++)
         {
            int row = lpi->rngrows[i];
            if ( lpi->rngrowmap[row] >= 0 )
            {
               lpi->rngrowmap[row] = nrngrows;
               lpi->rngrows[nrngrows] = row;
               lpi->rngvals[nrngrows] = lpi->rngvals[i];
               nrngrows++;
            }
         }
         lpi->nrngrows = nrngrows;

         /* discard ranged row info if no ranged rows remain */
         if ( nrngrows == 0 )
            clearRangeInfo(lpi);
      }
   }

   checkRangeInfo(lpi);

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
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("changing coefficient row %d, column %d in Gurobi to %g\n", row, col, newval);

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, GRBchgcoeffs(lpi->grbmodel, 1, &row, &col, &newval) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   int grbobjsen;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(objsen == SCIP_OBJSEN_MAXIMIZE || objsen == SCIP_OBJSEN_MINIMIZE);

   /* convert objective sense */
   grbobjsen = (objsen == SCIP_OBJSEN_MINIMIZE) ? GRB_MINIMIZE : GRB_MAXIMIZE;

   SCIPdebugMessage("changing objective sense in Gurobi to %d\n", grbobjsen);

   invalidateSolution(lpi);

   /* The objective sense of Gurobi and SCIP are equal */
   CHECK_ZERO( lpi->messagehdlr, GRBsetintattr(lpi->grbmodel, GRB_INT_ATTR_MODELSENSE, grbobjsen) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

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
   assert(lpi->grbmodel != NULL);
   assert(ind != NULL);
   assert(obj != NULL);

   SCIPdebugMessage("changing %d objective values in Gurobi\n", ncols);

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrlist(lpi->grbmodel, GRB_DBL_ATTR_OBJ, ncols, (int*)ind, (SCIP_Real*)obj) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

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
   int beg;
   int i;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(scaleval != 0.0);

   SCIPdebugMessage("scaling row %d with factor %g in Gurobi\n", row, scaleval);

   invalidateSolution(lpi);

   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   SCIP_CALL( ensureValMem(lpi, ncols+1) ); /* +1 for range variable */

   /* get the row */
   SCIP_CALL( SCIPlpiGetRows(lpi, row, row, &lhs, &rhs, &nnonz, &beg, lpi->indarray, lpi->valarray) );

   /* scale row coefficients */
   for ( i = 0; i < nnonz; ++i )
   {
      SCIP_CALL( SCIPlpiChgCoef(lpi, row, lpi->indarray[i], lpi->valarray[i] * scaleval) );
   }

   /* scale row sides */
   if( lhs > -GRB_INFBOUND )
      lhs *= scaleval;
   else if( scaleval < 0.0 )
      lhs = GRB_INFBOUND;
   if( rhs < GRB_INFBOUND )
      rhs *= scaleval;
   else if( scaleval < 0.0 )
      rhs = -GRB_INFBOUND;
   if( scaleval > 0.0 )
   {
      SCIP_CALL( SCIPlpiChgSides(lpi, 1, &row, &lhs, &rhs) );
   }
   else
   {
      SCIP_CALL( SCIPlpiChgSides(lpi, 1, &row, &rhs, &lhs) );
   }

   checkRangeInfo(lpi);

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
   int beg;
   int i;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(scaleval != 0.0);

   SCIPdebugMessage("scaling column %d with factor %g in Gurobi\n", col, scaleval);

   invalidateSolution(lpi);

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( ensureValMem(lpi, nrows) );

   /* get the column */
   SCIP_CALL( SCIPlpiGetCols(lpi, col, col, &lb, &ub, &nnonz, &beg, lpi->indarray, lpi->valarray) );

   /* get objective coefficient */
   SCIP_CALL( SCIPlpiGetObj(lpi, col, col, &obj) );

   /* scale column coefficients */
   for(  i = 0; i < nnonz; ++i )
   {
      SCIP_CALL( SCIPlpiChgCoef(lpi, lpi->indarray[i], col, lpi->valarray[i] * scaleval) );
   }

   /* scale objective value */
   obj *= scaleval;
   SCIP_CALL( SCIPlpiChgObj(lpi, 1, &col, &obj) );

   /* scale column bounds */
   if( lb > -GRB_INFINITY )
      lb /= scaleval;
   else if( scaleval < 0.0 )
      lb = GRB_INFINITY;
   if( ub < GRB_INFINITY )
      ub /= scaleval;
   else if( scaleval < 0.0 )
      ub = -GRB_INFINITY;
   if( scaleval > 0.0 )
   {
      SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &col, &lb, &ub) );
   }
   else
   {
      SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &col, &ub, &lb) );
   }

   checkRangeInfo(lpi);

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
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(nrows != NULL);

   SCIPdebugMessage("getting number of rows\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, nrows) );

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(ncols != NULL);

   SCIPdebugMessage("getting number of columns\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, ncols) );

   /* subtract number of ranged rows, as these are the LPI internal columns */
   if ( lpi->rngvarsadded )
     (*ncols) -= lpi->nrngrows;

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(nnonz != NULL);

   SCIPdebugMessage("getting number of non-zeros\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMNZS, nnonz) );

   /* subtract number of ranged rows, as these are non-zeros for the LPI internal columns */
   if ( lpi->rngvarsadded )
      (*nnonz) -= lpi->nrngrows;

   return SCIP_OKAY;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values;
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
   assert(lpi->grbmodel != NULL);
   assert((lb != NULL && ub != NULL) || (lb == NULL && ub == NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));
#ifndef NDEBUG
   {
      int ncols;
      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
      assert(0 <= firstcol && firstcol <= lastcol && lastcol < ncols);
   }
#endif

   SCIPdebugMessage("getting columns %d to %d\n", firstcol, lastcol);

   if( lb != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_LB, firstcol, lastcol-firstcol+1, lb) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_UB, firstcol, lastcol-firstcol+1, ub) );
   }

   if( nnonz != NULL )
   {
      /* get matrix entries */
      CHECK_ZERO( lpi->messagehdlr, GRBgetvars(lpi->grbmodel, nnonz, beg, ind, val, firstcol, lastcol-firstcol+1) );
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
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert((lhs == NULL && rhs == NULL) || (rhs != NULL && lhs != NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));

#ifndef NDEBUG
   {
      int nrows;
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
      assert(0 <= firstrow && firstrow <= lastrow && lastrow < nrows);
   }
#endif

   SCIPdebugMessage("getting rows %d to %d\n", firstrow, lastrow);

   if( lhs != NULL )
   {
      /* get row sense and rhs */
      SCIP_CALL( ensureSidechgMem(lpi, lastrow - firstrow + 1) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_RHS, firstrow, lastrow-firstrow+1, lpi->rhsarray) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetcharattrarray(lpi->grbmodel, GRB_CHAR_ATTR_SENSE, firstrow, lastrow-firstrow+1, lpi->senarray) );

      /* convert sen and rhs into lhs/rhs tuples */
      SCIP_CALL( reconvertSides(lpi, firstrow, lastrow, lhs, rhs) );
   }

   if( nnonz != NULL )
   {
      assert(beg != NULL && ind != NULL && val != NULL); /* for lint */

      /* get matrix entries */
      CHECK_ZERO( lpi->messagehdlr, GRBgetconstrs(lpi->grbmodel, nnonz, beg, ind, val, firstrow, lastrow-firstrow+1) );

      if ( lpi->rngvarsadded )
      {
         int i;

         assert(lpi->rngrowmap != NULL);
         assert(lpi->rngrows != NULL);

         /* remove non-zeros for range variables from rows */
         for (i = firstrow; i <= lastrow; i++)
         {
            assert(-1 <= lpi->rngrowmap[i] && lpi->rngrowmap[i] < lpi->nrngrows);
            if ( lpi->rngrowmap[i] >= 0 )
               break;
         }
         if ( i <= lastrow )
         {
            /* skip last non-zero of this first ranged row */
            int newnz = (i < lastrow ? beg[i - firstrow +1]-1 : (*nnonz)-1); /*lint !e661*/

            /* process remaining rows, moving non-zeros to the front */
            for (; i <= lastrow; i++)
            {
               int thebeg;
               int theend;

               thebeg = beg[i - firstrow]; /*lint !e661*/
               theend = (i < lastrow ? beg[i - firstrow +1] : *nnonz);

               assert(-1 <= lpi->rngrowmap[i] && lpi->rngrowmap[i] < lpi->nrngrows);
               if ( lpi->rngrowmap[i] >= 0 )
                  theend--;

               assert( theend >= thebeg );
               memmove(&ind[newnz], &ind[thebeg], ((size_t) (theend - thebeg)) * sizeof(*ind)); /*lint !e776 !e571*/
               memmove(&val[newnz], &val[thebeg], ((size_t) (theend - thebeg)) * sizeof(*val)); /*lint !e776 !e571*/
               beg[i - firstrow] = newnz; /*lint !e661*/
               newnz += theend - thebeg;
            }
            assert(newnz < *nnonz);
            *nnonz = newnz;
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
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
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
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
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
   int grbobjsen;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( objsen != NULL );

   SCIPdebugMessage("getting objective sense\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_MODELSENSE, &grbobjsen) );
   assert(grbobjsen == GRB_MINIMIZE || grbobjsen == GRB_MAXIMIZE);

   *objsen = (grbobjsen == GRB_MINIMIZE) ? SCIP_OBJSEN_MINIMIZE : SCIP_OBJSEN_MAXIMIZE;

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
   assert(lpi->grbmodel != NULL);
   assert(firstcol <= lastcol);
   assert(vals != NULL);

   SCIPdebugMessage("getting objective values %d to %d\n", firstcol, lastcol);

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_OBJ, firstcol, lastcol-firstcol+1, vals) );

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
   assert(lpi->grbmodel != NULL);
#ifndef NDEBUG
   {
      int ncols;
      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
      assert(0 <= firstcol && firstcol <= lastcol && lastcol < ncols);
   }
#endif

   SCIPdebugMessage("getting bounds %d to %d\n", firstcol, lastcol);

   if( lbs != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_LB, firstcol, lastcol-firstcol+1, lbs) );
   }

   if( ubs != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_UB, firstcol, lastcol-firstcol+1, ubs) );
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
   assert(lpi->grbmodel != NULL);
   assert(firstrow <= lastrow);

   SCIPdebugMessage("getting row sides %d to %d\n", firstrow, lastrow);

   /* get row sense, rhs, and ranges */
   SCIP_CALL( ensureSidechgMem(lpi, lastrow - firstrow + 1) );

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_RHS, firstrow, lastrow-firstrow+1, lpi->rhsarray) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetcharattrarray(lpi->grbmodel, GRB_CHAR_ATTR_SENSE, firstrow, lastrow-firstrow+1, lpi->senarray) );

   /* convert sen and rhs into lhs/rhs tuples */
   SCIP_CALL( reconvertSides(lpi, firstrow, lastrow, lhss, rhss) );

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
   assert(lpi->grbmodel != NULL);
   assert(val != NULL);

   SCIPdebugMessage("getting coefficient of row %d col %d\n", row, col);

   CHECK_ZERO( lpi->messagehdlr, GRBgetcoeff(lpi->grbmodel, row, col, val) );

   return SCIP_OKAY;
}

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP
 *
 *  @todo Check concurrent (GRB_METHOD_CONCURRENT or GRB_METHOD_DETERMINISTIC_CONCURRENT)
 */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   double cnt;
   int retval;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

#ifdef SCIP_DEBUG
   {
      int ncols, nrows;
      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
      SCIPdebugMessage("calling Gurobi primal simplex: %d cols, %d rows\n", ncols, nrows);
   }
#endif

   invalidateSolution(lpi);

   if ( lpi->fromscratch )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBresetmodel(lpi->grbmodel) );
   }

   SCIPdebugMessage("calling GRBoptimize() - primal\n");

   /* set primal simplex */
   SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );
   CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, GRB_METHOD_PRIMAL) );

   /* add range variables */
   if ( lpi->nrngrows > 0 && !lpi->rngvarsadded )
   {
      SCIP_CALL( addRangeVars(lpi) );
   }

   retval = GRBoptimize(lpi->grbmodel);
   switch( retval  )
   {
   case 0:
      break;
   case GRB_ERROR_OUT_OF_MEMORY:
      return SCIP_NOMEMORY;
   default:
      return SCIP_LPERROR;
   }

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
   lpi->iterations = (int) cnt;

   lpi->solisbasic = TRUE;
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );

   SCIPdebugMessage("Gurobi primal simplex needed %d iterations to gain LP status %d\n", (int) cnt, lpi->solstat);

   /* maybe the preprocessor solved the problem; but we need a solution, so solve again without preprocessing */
   assert( lpi->solstat != GRB_INF_OR_UNBD );
   if( lpi->solstat == GRB_INFEASIBLE )
   {
      int presolve;

      CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, GRB_INT_PAR_PRESOLVE, &presolve) );

      if( presolve != GRB_PRESOLVE_OFF )
      {
         SCIPdebugMessage("presolver may have solved the problem -> calling Gurobi primal simplex again without presolve\n");

         /* switch off preprocessing */
         CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF) );

         retval = GRBoptimize(lpi->grbmodel);
         switch( retval  )
         {
         case 0:
            break;
         case GRB_ERROR_OUT_OF_MEMORY:
            return SCIP_NOMEMORY;
         default:
            return SCIP_LPERROR;
         }

         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
         lpi->iterations += (int) cnt;
         CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );
         SCIPdebugMessage(" -> Gurobi returned solstat=%d (%d iterations)\n", lpi->solstat, lpi->iterations);

         /* reset parameters */
         CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_PRESOLVE, presolve) );
      }

      if( lpi->solstat == GRB_INF_OR_UNBD )
      {
         /* preprocessing was not the problem; issue a warning message and treat LP as infeasible */
         SCIPerrorMessage("Gurobi primal simplex returned GRB_INF_OR_UNBD after presolving was turned off\n");
         return SCIP_LPERROR;
      }
   }
   else if ( lpi->solstat == GRB_UNBOUNDED )
   {
      /* Unbounded means that there exists an unbounded primal ray. However, this does not state whether the problem is
       * feasible. Thus, we temporarily set the objective to 0 and solve again. */
      SCIP_Real* zeroobjcoefs;
      SCIP_Real* objcoefs;
      SCIP_Real oldobjcutoff;
      int grbobjsen;
      int status;
      int ncols;

      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
      SCIP_ALLOC( BMSallocMemoryArray(&objcoefs, ncols) );
      SCIP_ALLOC( BMSallocClearMemoryArray(&zeroobjcoefs, ncols) );

      /* preserve objective coefficients */
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_OBJ, 0, ncols, objcoefs) );

      /* set objective to 0 */
      CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_OBJ, 0, ncols, zeroobjcoefs) );

      /* disable cutoff */
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, &oldobjcutoff) );

      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_MODELSENSE, &grbobjsen) );
      if ( grbobjsen == GRB_MINIMIZE )
      {
         CHECK_ZERO( lpi->messagehdlr, GRBsetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, GRB_INFINITY) );
      }
      else
      {
         CHECK_ZERO( lpi->messagehdlr, GRBsetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, -GRB_INFINITY) );
         assert( grbobjsen == GRB_MAXIMIZE );
      }

      /* solve problem again */
      CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );
      CHECK_ZERO( lpi->messagehdlr, GRBoptimize(lpi->grbmodel) );

      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
      lpi->iterations += (int) cnt;

      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &status) );

      /* restore objective */
      CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_OBJ, 0, ncols, objcoefs) );
      CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

      /* restore objective limit */
      CHECK_ZERO( lpi->messagehdlr, GRBsetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, oldobjcutoff) );

      BMSfreeMemoryArray(&zeroobjcoefs);
      BMSfreeMemoryArray(&objcoefs);

      /* possibly correct status */
      switch ( status )
      {
      case GRB_INF_OR_UNBD:
      case GRB_INFEASIBLE:
         lpi->solstat = GRB_INFEASIBLE;
         break;

      case GRB_OPTIMAL:
         /* We again have to solve the problem to restore possible unbounded rays. */
         CHECK_ZERO( lpi->messagehdlr, GRBoptimize(lpi->grbmodel) );
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
         lpi->iterations += (int) cnt;
         break;

      case GRB_ITERATION_LIMIT:
      case GRB_TIME_LIMIT:
         /* do nothing */
         break;

         /* GRB_LOADED, GRB_NODE_LIMIT, GRB_CUTOFF, GRB_SOLUTION_LIMIT, GRB_INTERRUPTED, GRB_NUMERIC, GRB_SUBOPTIMAL, GRB_INPROGRESS, GRB_USER_OBJ_LIMIT */
      default:
         SCIPerrorMessage("Gurobi returned wrong status %d.\n", status);
         return SCIP_LPERROR;
      }
   }

   checkRangeInfo(lpi);

   return SCIP_OKAY;
}

/** calls dual simplex to solve the LP
 *
 *  @todo Check concurrent (GRB_METHOD_CONCURRENT or GRB_METHOD_DETERMINISTIC_CONCURRENT)
 */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int oldprimdual = 0;
   int oldpresolve = GRB_PRESOLVE_OFF;
   int retval;
   double cnt;
   double itlim;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

#ifdef SCIP_DEBUG
   {
      int ncols, nrows;
      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
      SCIPdebugMessage("calling Gurobi dual simplex: %d cols, %d rows\n", ncols, nrows);
   }
#endif

   invalidateSolution(lpi);

   if ( lpi->fromscratch )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBresetmodel(lpi->grbmodel) );
   }

   SCIPdebugMessage("calling GRBoptimize() - dual\n");

   SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );

   /* set dual simplex */
   CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL) );

   /* add range variables */
   if ( lpi->nrngrows > 0 && !lpi->rngvarsadded )
   {
      SCIP_CALL( addRangeVars(lpi) );
   }

   SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_ITERATIONLIMIT, &itlim) );
   if ( itlim < GRB_INFINITY )
   {
      /* turn off primal-dual switching for an LP solve that might be a strong branching LP solve */
      CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, "GURO_PAR_PRIMDUALSWITCH", &oldprimdual) );
      if ( oldprimdual != 0 )
      {
         CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, "GURO_PAR_PRIMDUALSWITCH", 0) );
      }

      /* turn off presolve to avoid the case where the iteration limit is reached
       * and we do not get a valid dual bound installed for the original model */
      CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, GRB_INT_PAR_PRESOLVE, &oldpresolve) );
      if ( oldpresolve != GRB_PRESOLVE_OFF )
      {
         CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF) );
      }
   }

   retval = GRBoptimize(lpi->grbmodel);
   switch( retval  )
   {
   case 0:
      break;
   case GRB_ERROR_OUT_OF_MEMORY:
      return SCIP_NOMEMORY;
   default:
      return SCIP_LPERROR;
   }

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
   lpi->iterations = (int) cnt;

   lpi->solisbasic = TRUE;
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );

   SCIPdebugMessage("Gurobi dual simplex needed %d iterations to gain LP status %d\n", (int) cnt, lpi->solstat);

   if( lpi->solstat == GRB_INF_OR_UNBD )
   {
      int presolve;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_PAR_PRESOLVE, &presolve) );

      if( presolve != GRB_PRESOLVE_OFF )
      {
         /* maybe the preprocessor solved the problem; but we need a solution, so solve again without preprocessing */
         SCIPdebugMessage("presolver may have solved the problem -> calling Gurobi dual simplex again without presolve\n");

         /* switch off preprocessing */
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattr(lpi->grbmodel, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF) );
         SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );

         retval = GRBoptimize(lpi->grbmodel);
         switch( retval  )
         {
         case 0:
            break;
         case GRB_ERROR_OUT_OF_MEMORY:
            return SCIP_NOMEMORY;
         default:
            return SCIP_LPERROR;
         }

         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
         lpi->iterations += (int) cnt;
         CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );
         SCIPdebugMessage(" -> Gurobi returned solstat=%d (%d iterations)\n", lpi->solstat, lpi->iterations);

         /* switch on preprocessing again */
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattr(lpi->grbmodel, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_AUTO) );
      }

      if( lpi->solstat == GRB_INF_OR_UNBD )
      {
         /* preprocessing was not the problem; issue a warning message and treat LP as infeasible */
         SCIPerrorMessage("Gurobi dual simplex returned GRB_INF_OR_UNBD after presolving was turned off.\n");
         return SCIP_LPERROR;
      }
   }

   checkRangeInfo(lpi);

   /* reset parameters to their original values */
   if ( oldprimdual != 0 )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, "GURO_PAR_PRIMDUALSWITCH", oldprimdual) );
   }
   if ( oldpresolve != GRB_PRESOLVE_OFF )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_PRESOLVE, oldpresolve) );
   }

   return SCIP_OKAY;
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{
   int retval;
   double cnt;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

#ifdef SCIP_DEBUG
   {
      int ncols, nrows;
      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
      SCIPdebugMessage("calling Gurobi barrier: %d cols, %d rows\n", ncols, nrows);
   }
#endif

   invalidateSolution(lpi);

   if ( lpi->fromscratch )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBresetmodel(lpi->grbmodel) );
   }

   SCIPdebugMessage("calling GRBoptimize() - barrier\n");

   /* set barrier */
   SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );

   if( crossover )
   {
      /* turn on crossover to automatic setting (-1) */
      CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_CROSSOVER, -1) );
   }
   else
   {
      /* turn off crossover */
      CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_CROSSOVER, 0) );
   }

   CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER) );

   /* add range variables */
   if ( lpi->nrngrows > 0 && !lpi->rngvarsadded )
   {
      SCIP_CALL( addRangeVars(lpi) );
   }

   retval = GRBoptimize(lpi->grbmodel);
   switch( retval  )
   {
   case 0:
      break;
   case GRB_ERROR_OUT_OF_MEMORY:
      return SCIP_NOMEMORY;
   default:
      return SCIP_LPERROR;
   }

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
   lpi->iterations = (int) cnt;

   lpi->solisbasic = crossover;
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );

   SCIPdebugMessage("Gurobi barrier needed %d iterations to gain LP status %d\n", (int) cnt, lpi->solstat);

   if( lpi->solstat == GRB_INF_OR_UNBD )
   {
      int presolve;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_PAR_PRESOLVE, &presolve) );

      if( presolve != GRB_PRESOLVE_OFF )
      {
         /* maybe the preprocessor solved the problem; but we need a solution, so solve again without preprocessing */
         SCIPdebugMessage("presolver may have solved the problem -> calling Gurobi barrier again without presolve\n");

         /* switch off preprocessing */
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattr(lpi->grbmodel, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF) );
         SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );

         retval = GRBoptimize(lpi->grbmodel);
         switch( retval  )
         {
         case 0:
            break;
         case GRB_ERROR_OUT_OF_MEMORY:
            return SCIP_NOMEMORY;
         default:
            return SCIP_LPERROR;
         }

         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
         lpi->iterations += (int) cnt;
         CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );
         SCIPdebugMessage(" -> Gurobi returned solstat=%d (%d iterations)\n", lpi->solstat, lpi->iterations);

         /* switch on preprocessing again */
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattr(lpi->grbmodel, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_AUTO) );
      }

      if( lpi->solstat == GRB_INF_OR_UNBD )
      {
         /* preprocessing was not the problem; issue a warning message and treat LP as infeasible */
         SCIPerrorMessage("Gurobi barrier returned GRB_INF_OR_UNBD after presolving was turned off\n");
         return SCIP_LPERROR;
      }
   }

   checkRangeInfo(lpi);

   return SCIP_OKAY;
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

   /* currently do nothing */
   return SCIP_OKAY;
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

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
   SCIP_Real oldlb;
   SCIP_Real oldub;
   SCIP_Real newlb;
   SCIP_Real newub;
   SCIP_Real olditlim;
   SCIP_Bool error = FALSE;
   SCIP_Bool success;
   int it;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   SCIPdebugMessage("performing strong branching on variable %d (%d iterations)\n", col, itlim);

   SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );

   *downvalid = FALSE;
   *upvalid = FALSE;
   if( iter != NULL )
      *iter = 0;

   /* save current LP basis and bounds*/
   SCIP_CALL( getBase(lpi, &success) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_LB, col, &oldlb) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_UB, col, &oldub) );

   if ( lpi->fromscratch )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBresetmodel(lpi->grbmodel) );
   }

   /* save old iteration limit and set iteration limit to strong branching limit */
   if( itlim < 0 )
      itlim = INT_MAX;

   SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_ITERATIONLIMIT, &olditlim) );
   SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_ITERATIONLIMIT, (double) itlim) );

   /* add range variables */
   if ( lpi->nrngrows > 0 && !lpi->rngvarsadded )
   {
      SCIP_CALL( addRangeVars(lpi) );
   }

   /* down branch */
   newub = EPSCEIL(psol-1.0, 1e-06);
   if( newub >= oldlb - 0.5 )
   {
      SCIPdebugMessage("strong branching down (%g) on x%d (%g) with %d iterations\n", newub, col, psol, itlim);

      CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_UB, col, newub) );

      SCIP_CALL( SCIPlpiSolveDual(lpi) );
      /* when iteration limit was reached the objective value is not computed */
      if( SCIPlpiIsOptimal(lpi) ) /*|| SCIPlpiIsIterlimExc(lpi) ) */
      {
         SCIP_CALL( SCIPlpiGetObjval(lpi, down) );
         *downvalid = TRUE;
      }
      else if( SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) )
      {
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, down) );
      }
      else if( !SCIPlpiIsIterlimExc(lpi) )
         error = TRUE;

      if( iter != NULL )
      {
         SCIP_CALL( SCIPlpiGetIterations(lpi, &it) );
         *iter += it;
      }
      SCIPdebugMessage(" -> down (x%d <= %g): %g\n", col, newub, *down);

      CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_UB, col, oldub) );
      CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );
#ifdef SCIP_DEBUG
      {
         double b;
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_UB, col, &b) );
         assert( b == oldub );
      }
#endif

      if ( success )
      {
         SCIP_CALL( setBase(lpi) );
      }
   }
   else
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, down) );
      *downvalid = TRUE;
   }

   /* up branch */
   if( !error )
   {
      newlb = EPSFLOOR(psol+1.0, 1e-06);
      if( newlb <= oldub + 0.5 )
      {
         SCIPdebugMessage("strong branching  up (%g) on x%d (%g) with %d iterations\n", newlb, col, psol, itlim);

         CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_LB, col, newlb) );

         SCIP_CALL( SCIPlpiSolveDual(lpi) );
         /* when iteration limit was reached the objective value is not computed */
         if( SCIPlpiIsOptimal(lpi) ) /*|| SCIPlpiIsIterlimExc(lpi) ) */
         {
            SCIP_CALL( SCIPlpiGetObjval(lpi, up) );
            *upvalid = TRUE;
         }
         else if( SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) )
         {
            CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, up) );
         }
         else if( !SCIPlpiIsIterlimExc(lpi) )
            error = TRUE;

         if( iter != NULL )
         {
            SCIP_CALL( SCIPlpiGetIterations(lpi, &it) );
            *iter += it;
         }
         SCIPdebugMessage(" -> up  (x%d >= %g): %g\n", col, newlb, *up);

         CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_LB, col, oldlb) );
         CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );
#ifdef SCIP_DEBUG
         {
            double b;
            CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_LB, col, &b) );
            assert( b == oldlb );
         }
#endif

         if ( success )
         {
            SCIP_CALL( setBase(lpi) );
         }
      }
      else
      {
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, up) );
         *upvalid = TRUE;
      }
   }

   /* reset iteration limit */
   SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_ITERATIONLIMIT, olditlim) );
   /* CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) ); */

   if( error )
   {
      SCIPerrorMessage("LP error in strong branching.\n");
      return SCIP_LPERROR;
   }

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

   checkRangeInfo(lpi);

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
   int j;

   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if( iter != NULL )
      *iter = 0;

   for( j = 0; j < ncols; ++j )
   {
      /* pass call on to lpiStrongbranch() */
      SCIP_CALL( lpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter) );
   }

   checkRangeInfo(lpi);

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

   checkRangeInfo(lpi);

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
   int j;

   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if( iter != NULL )
      *iter = 0;

   for( j = 0; j < ncols; ++j )
   {
      /* pass call on to lpiStrongbranch() */
      SCIP_CALL( lpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter) );
   }

   checkRangeInfo(lpi);

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
   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 1 );
   assert( primalfeasible != NULL );
   assert( dualfeasible != NULL );

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
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_UNBOUNDED);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int algo;
   int res;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 0 );

   res = GRBgetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, &algo);
   if ( res != 0 )
   {
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }

   return (lpi->solstat == GRB_UNBOUNDED && algo == GRB_METHOD_PRIMAL);
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int algo;
   int res;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 0 );

   res = GRBgetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, &algo);
   if ( res != 0 )
   {
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }

   /* GRB_UNBOUNDED means that there exists a primal ray. SCIPlpiSolvePrimal() will determine whether the problem is
    * actually infeasible or (feasible and) unbounded. In the latter case, the status will be GRB_UNBOUNDED.
    */
   return (lpi->solstat == GRB_UNBOUNDED && algo == GRB_METHOD_PRIMAL);
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for primal infeasibility\n");

   assert( lpi->solstat != GRB_INF_OR_UNBD );
   return (lpi->solstat == GRB_INFEASIBLE);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int algo;
   int res;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 0 );

   SCIPdebugMessage("checking for primal feasibility\n");

   if ( lpi->solstat == GRB_OPTIMAL )
      return TRUE;

   res = GRBgetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, &algo);
   if ( res != 0 )
   {
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }
   if ( algo != GRB_METHOD_PRIMAL )
      return FALSE;

   if( lpi->solstat == GRB_ITERATION_LIMIT )
   {
      double consviol;
      double boundviol;
      double eps;

      /* get feasibility tolerance */
      res = GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_FEASIBILITYTOL, &eps);
      if ( res != 0 )
      {
         SCIPABORT();
         return FALSE; /*lint !e527*/
      }
      res = GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_CONSTR_VIO, &consviol);
      if ( res != 0 )
      {
         /* If Gurobi cannot return the constraint violation, there is no feasible solution available. */
         return FALSE;
      }
      res = GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_BOUND_VIO, &boundviol);
      if ( res != 0 )
      {
         /* If Gurobi cannot return the bound violation, there is no feasible solution available. */
         return FALSE;
      }

      if ( consviol <= eps && boundviol <= eps )
         return TRUE;
   }

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
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_INFEASIBLE);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int algo;
   int res;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 0 );

   res = GRBgetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, &algo);
   if ( res != 0 )
   {
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }

   return (lpi->solstat == GRB_INFEASIBLE && algo == GRB_METHOD_DUAL);
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int algo;
   int res;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 0 );

   SCIPdebugMessage("checking for dual unboundedness\n");

   res = GRBgetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, &algo);
   if ( res != 0 )
   {
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }

   return (lpi->solstat == GRB_INFEASIBLE && algo == GRB_METHOD_DUAL);
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->solstat >= 0 );

   SCIPdebugMessage("checking for dual infeasibility\n");

   return (lpi->solstat == GRB_UNBOUNDED);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int algo;
   int res;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 0 );

   SCIPdebugMessage("checking for dual feasibility\n");

   res = GRBgetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, &algo);
   if ( res != 0 )
   {
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }

   return (lpi->solstat == GRB_OPTIMAL                                      ||
           (lpi->solstat == GRB_INFEASIBLE      && algo == GRB_METHOD_DUAL) ||
           (lpi->solstat == GRB_ITERATION_LIMIT && algo == GRB_METHOD_DUAL)   );
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_OPTIMAL);
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
   double consviol;
   double boundviol;
   double dualviol;
   double feastol;
   double optimalitytol;
   int res;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for stability: Gurobi solstat = %d\n", lpi->solstat);

   /* If the condition number of the basis should be checked, everything above the specified threshold is counted as
    * instable. */
   if ( lpi->checkcondition && (SCIPlpiIsOptimal(lpi) || SCIPlpiIsObjlimExc(lpi)) )
   {
      SCIP_Real kappa;
      SCIP_RETCODE retcode;

      retcode = SCIPlpiGetRealSolQuality(lpi, SCIP_LPSOLQUALITY_ESTIMCONDITION, &kappa);
      if ( retcode != SCIP_OKAY ) /*lint !e774*/
      {
         SCIPABORT();
         return FALSE; /*lint !e527*/
      }

      /* if the kappa could not be computed (e.g., because we do not have a basis), we cannot check the condition */
      if ( kappa != SCIP_INVALID || kappa > lpi->conditionlimit ) /*lint !e777*/
         return FALSE;
   }

   /* test whether we have unscaled infeasibilities */
   if ( SCIPlpiIsOptimal(lpi) )
   {
      /* first get tolerances */
      res = GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_FEASIBILITYTOL, &feastol);
      if ( res != 0 )
      {
         SCIPABORT();
         return FALSE; /*lint !e527*/
      }
      res = GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_OPTIMALITYTOL, &optimalitytol);
      if ( res != 0 )
      {
         SCIPABORT();
         return FALSE; /*lint !e527*/
      }

      /* next get constraint, bound, and reduced cost violations */
      res = GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_CONSTR_VIO, &consviol);
      if ( res != 0 )
      {
         SCIPABORT();
         return FALSE; /*lint !e527*/
      }
      res = GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_BOUND_VIO, &boundviol);
      if ( res != 0 )
      {
         SCIPABORT();
         return FALSE; /*lint !e527*/
      }
      res = GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_DUAL_VIO, &dualviol);
      if ( res != 0 )
      {
         SCIPABORT();
         return FALSE; /*lint !e527*/
      }

      return ( consviol <= feastol && boundviol <= feastol && dualviol <= optimalitytol );
   }

   return (lpi->solstat != GRB_NUMERIC);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_CUTOFF);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_ITERATION_LIMIT);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_TIME_LIMIT);
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   return lpi->solstat;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(success != NULL);

   *success = FALSE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   int ret;

#ifndef NDEBUG
   double oval = GRB_INFINITY;
   double obnd = -GRB_INFINITY;
#endif

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(objval != NULL);

   SCIPdebugMessage("getting solution's objective value\n");

#ifndef NDEBUG
   (void)GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_OBJVAL, &oval);
   (void)GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_OBJBOUND, &obnd);

   assert(lpi->solstat != GRB_OPTIMAL || oval == obnd); /*lint !e777*/
#endif

   /* obtain objective value */
   ret = GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_OBJVAL, objval);
   assert( ret == 0 || ret == GRB_ERROR_DATA_NOT_AVAILABLE );
   SCIP_UNUSED(ret);

   /* return minus infinity if value not available and we reached the iteration limit (see lpi_cpx) */
   if( lpi->solstat == GRB_ITERATION_LIMIT )
   {
      (void)GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_OBJBOUND, objval);
   }
   else if( lpi->solstat == GRB_CUTOFF )
   {
      SCIP_Real cutoff;

      /* if we reached the cutoff, then OBJBOUND seems to be -infinity; we set the value to the cutoff in this case */
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_CUTOFF, &cutoff) );
      *objval = cutoff;

#ifdef SCIP_DISABLED_CODE
      /**@todo The following is some kind of hack which works with the current SCIP implementation and should be fixed.  In
       * the case that the LP status is GRB_CUTOFF it might be that certain attributes cannot be queried (e.g., objval,
       * primal and dual solution), in this case we just return the installed cutoff value minus some epsilon. This is some
       * kind of hack for the code in conflict.c:7595 were some extra code handles CPLEX' FASTMIP case that is similar to
       * this case.
       */
      SCIP_Real dval;
      SCIP_OBJSEN objsense;

      SCIP_CALL( SCIPlpiGetObjsen(lpi, &objsense) );
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_CUTOFF, &dval) );

      if( objsense == SCIP_OBJSEN_MINIMIZE )
         *objval = dval - 1e-06;
      else
         *objval = dval + 1e-06;
#endif
   }

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
   int ncols;
   int nrows;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("getting solution\n");

   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   assert( ncols >= 0 && nrows >= 0 );

   if( objval != NULL )
   {
      SCIP_CALL( SCIPlpiGetObjval(lpi, objval) );
   }

   if( primsol != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_X, 0, ncols, primsol) );
   }

   if( dualsol != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_PI, 0, nrows, dualsol) );
   }

   if( activity != NULL )
   {
      int i;

      /* first get the values of the slack variables */
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_SLACK, 0, nrows, activity) );

      SCIP_CALL( ensureSidechgMem(lpi, nrows) );

      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_RHS, 0, nrows, lpi->rhsarray) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetcharattrarray(lpi->grbmodel, GRB_CHAR_ATTR_SENSE, 0, nrows, lpi->senarray) );

      for( i = 0; i < nrows; ++i )
      {
         switch(lpi->senarray[i])
         {
         case GRB_EQUAL:
            if ( lpi->rngrowmap != NULL && lpi->rngrowmap[i] >= 0 )
            {
               /* get solution value of range variable */
               SCIP_Real solval;
               assert(lpi->rngrowmap[i] < lpi->nrngrows);
               CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_X, ncols + lpi->rngrowmap[i], &solval) );
               activity[i] = lpi->rhsarray[i] + solval;
            }
            else
            {
               activity[i] = lpi->rhsarray[i] - activity[i];
            }
            break;
         case GRB_LESS_EQUAL:
            activity[i] = lpi->rhsarray[i] - activity[i];
            break;
         case GRB_GREATER_EQUAL:
            activity[i] = lpi->rhsarray[i] - activity[i];
            break;
         default:
            SCIPerrorMessage("Unkown sense %c.\n", lpi->senarray[i]);
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }

   if( redcost != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_RC, 0, ncols, redcost) );
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
   int ncols;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);
   assert(ray != NULL);

   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   assert( ncols >= 0 );

   SCIPdebugMessage("calling Gurobi get primal ray: %d cols\n", ncols);

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_UNBDRAY, 0, ncols, ray) );

   return SCIP_OKAY;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
   )
{
   int nrows;
   int i;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);
   assert(dualfarkas != NULL);

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   assert( nrows >= 0 );

   SCIPdebugMessage("calling Gurobi dual Farkas: %d rows\n", nrows);

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_FARKASDUAL, 0, nrows, dualfarkas) );

   /* correct sign of ray */
   for (i = 0; i < nrows; ++i)
      dualfarkas[i] *= -1.0;

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
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
   const char* what;
   int ret;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(quality != NULL);

   SCIPdebugMessage("requesting solution quality from Gurobi: quality %d\n", qualityindicator);

   switch( qualityindicator )
   {
   case SCIP_LPSOLQUALITY_ESTIMCONDITION:
      what = GRB_DBL_ATTR_KAPPA;
      break;

   case SCIP_LPSOLQUALITY_EXACTCONDITION:
      what = GRB_DBL_ATTR_KAPPA_EXACT;
      break;

   default:
      SCIPerrorMessage("Solution quality %d unknown.\n", qualityindicator);
      return SCIP_INVALIDDATA;
   }

   ret = GRBgetdblattr(lpi->grbmodel, what, quality);
   if( ret != 0 )
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
   int nrows;
   int ncols;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("saving Gurobi basis into %p/%p\n", (void*) cstat, (void*) rstat);

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );

   if( rstat != NULL )
   {
      int i;

      SCIP_CALL( ensureSidechgMem(lpi, nrows) );

      CHECK_ZERO( lpi->messagehdlr, GRBgetintattrarray(lpi->grbmodel, GRB_INT_ATTR_CBASIS, 0, nrows, rstat) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetcharattrarray(lpi->grbmodel, GRB_CHAR_ATTR_SENSE, 0, nrows, lpi->senarray) );

      for( i = 0; i < nrows; ++i )
      {
         if ( lpi->rngrowmap != NULL && lpi->rngrowmap[i] >= 0 && rstat[i] != GRB_BASIC )
         {
            int idx;

            /* get range row basis status from corresponding range variable */
            idx = ncols + lpi->rngrowmap[i];
            assert(lpi->rngrowmap[i] < lpi->nrngrows);
            CHECK_ZERO( lpi->messagehdlr, GRBgetintattrelement(lpi->grbmodel, GRB_INT_ATTR_VBASIS, idx, &rstat[i]) );

            switch( rstat[i] )
            {
            case GRB_BASIC:
               rstat[i] = (int) SCIP_BASESTAT_BASIC;
               break;

            case GRB_NONBASIC_LOWER:
               rstat[i] = (int) SCIP_BASESTAT_LOWER;
               break;

            case GRB_NONBASIC_UPPER:
               rstat[i] = (int) SCIP_BASESTAT_UPPER;
               break;

               /*lint -fallthrough*/
            case GRB_SUPERBASIC:
            default:
               SCIPerrorMessage("invalid basis status %d for ranged row.\n", rstat[i]);
               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }
         }
         else
         {
            /* Slack variables can only be basic or at their lower bounds in Gurobi. */
            switch( rstat[i] )
            {
            case GRB_BASIC:
               rstat[i] = (int) SCIP_BASESTAT_BASIC;
               break;

            case GRB_NONBASIC_LOWER:
               if ( lpi->senarray[i] == '>' || lpi->senarray[i] == '=' )
                  rstat[i] = (int) SCIP_BASESTAT_LOWER;
               else
               {
                  assert( lpi->senarray[i] == '<' );
                  rstat[i] = (int) SCIP_BASESTAT_UPPER;
               }
               break;

               /*lint -fallthrough*/
            case GRB_NONBASIC_UPPER:
            case GRB_SUPERBASIC:
            default:
               SCIPerrorMessage("invalid basis status %d for row.\n", rstat[i]);
               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }
         }
      }
   }

   if( cstat != 0 )
   {
      int j;

      CHECK_ZERO( lpi->messagehdlr, GRBgetintattrarray(lpi->grbmodel, GRB_INT_ATTR_VBASIS, 0, ncols, cstat) );

      for( j = 0; j < ncols; ++j )
      {
         switch( cstat[j] )
         {
         case GRB_BASIC:
            cstat[j] = (int) SCIP_BASESTAT_BASIC;
            break;

         case GRB_NONBASIC_LOWER:
            cstat[j] = (int) SCIP_BASESTAT_LOWER;
            break;

         case GRB_NONBASIC_UPPER:
            cstat[j] = (int) SCIP_BASESTAT_UPPER;
            break;

         case GRB_SUPERBASIC:
            cstat[j] = (int) SCIP_BASESTAT_ZERO;
            break;

         default:
            SCIPerrorMessage("invalid basis status %d for column.\n", cstat[j]);
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
   int i, j;
   int nrows, ncols;
#ifndef NDEBUG
   int nrngsfound = 0;
#endif

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );

   assert(cstat != NULL || ncols == 0);
   assert(rstat != NULL || nrows == 0);

   SCIPdebugMessage("loading basis %p/%p into Gurobi\n", (void*) cstat, (void*) rstat);

   invalidateSolution(lpi);

   SCIP_CALL( ensureCstatMem(lpi, ncols+lpi->nrngrows) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   for( i = 0; i < nrows; ++i )
   {
      if ( lpi->rngrowmap != NULL && lpi->rngrowmap[i] >= 0 )
      {
         int idx;

         /* set basis status of corresponding range variable; ranged row is always non-basic */
         idx = ncols + lpi->rngrowmap[i];
         assert(lpi->rngrowmap[i] < lpi->nrngrows);
         lpi->cstat[idx] = lpi->rstat[i];
         lpi->rstat[i] = GRB_NONBASIC_LOWER;
#ifndef NDEBUG
         nrngsfound++;
#endif
      }
      else
      {
         switch( rstat[i] ) /*lint !e613*/
         {
         case SCIP_BASESTAT_BASIC:
            lpi->rstat[i] = GRB_BASIC;
            break;

         case SCIP_BASESTAT_UPPER:
         {
#ifndef NDEBUG
            char sense;
            CHECK_ZERO( lpi->messagehdlr, GRBgetcharattrarray(lpi->grbmodel, GRB_CHAR_ATTR_SENSE, i, 1, &sense) );
            assert( sense == '<' );
#endif
            /* Slack variables can only be basic or at their lower bounds in Gurobi. */
            lpi->rstat[i] = GRB_NONBASIC_LOWER;
            break;
         }

         case SCIP_BASESTAT_LOWER:
         {
#ifndef NDEBUG
            char sense;
            CHECK_ZERO( lpi->messagehdlr, GRBgetcharattrarray(lpi->grbmodel, GRB_CHAR_ATTR_SENSE, i, 1, &sense) );
            assert( sense == '>' || sense == '=' );
#endif
            lpi->rstat[i] = GRB_NONBASIC_LOWER;
            break;
         }

         case SCIP_BASESTAT_ZERO:
         default:
            SCIPerrorMessage("invalid basis status %d for row.\n", rstat[i]); /*lint !e613*/
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }

   for( j = 0; j < ncols; ++j )
   {
      switch( cstat[j] ) /*lint !e613*/
      {
      case SCIP_BASESTAT_BASIC:
         lpi->cstat[j] = GRB_BASIC;
         break;

      case SCIP_BASESTAT_LOWER:
         lpi->cstat[j] = GRB_NONBASIC_LOWER;
         break;

      case SCIP_BASESTAT_UPPER:
         lpi->cstat[j] = GRB_NONBASIC_UPPER;
         break;

      case SCIP_BASESTAT_ZERO:
         lpi->cstat[j] = GRB_SUPERBASIC;
         break;

      default:
         SCIPerrorMessage("invalid basis status %d\n", cstat[j]); /*lint !e613*/
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

#ifndef NDEBUG
   assert(nrngsfound == lpi->nrngrows);
#endif

   CHECK_ZERO( lpi->messagehdlr, GRBsetintattrarray(lpi->grbmodel, GRB_INT_ATTR_CBASIS, 0, nrows, lpi->rstat) );
   CHECK_ZERO( lpi->messagehdlr, GRBsetintattrarray(lpi->grbmodel, GRB_INT_ATTR_VBASIS, 0, ncols+lpi->nrngrows, lpi->cstat) );

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   )
{
   int i;
   int nrows;
   int ncols;
   int ngrbcols;
   int* bhead;
   int status;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(bind != NULL);

   SCIPdebugMessage("getting basis information\n");

   /* check whether we have to reoptimize */
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &status) );
   if ( status == GRB_LOADED || status == GRB_INTERRUPTED || status == GRB_INPROGRESS )
   {
      SCIP_CALL_QUIET( restoreLPData(lpi) );
   }

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ngrbcols) );

   /**@todo avoid memory allocation by using bind directly */
   /* get space for bhead */
   SCIP_ALLOC( BMSallocMemoryArray(&bhead, nrows) );

   /* get basis indices */
   CHECK_ZERO( lpi->messagehdlr, GRBgetBasisHead(lpi->grbmodel, bhead) );

   for (i = 0; i < nrows; ++i)
   {
      /* entries >= ncols refer to slack variables */
      if ( bhead[i] < ncols )
         bind[i] = bhead[i];
      else if ( bhead[i] < ngrbcols )
      {
         /* a range variable: use corresponding ranged row */
         int rngrow = bhead[i]-ncols;
         assert(rngrow < lpi->nrngrows);
         assert(lpi->rngrowmap != NULL);
         assert(lpi->rngrows != NULL);
         assert(lpi->rngrowmap[lpi->rngrows[rngrow]] == rngrow);
         bind[i] = -1 - lpi->rngrows[rngrow];
      }
      else
      {
         /* a regular slack variable */
         bind[i] = -1 - (bhead[i] - ngrbcols);
      }
   }
   BMSfreeMemoryArray(&bhead);

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
   SVECTOR x;
   SVECTOR b;
   int nrows;
   double val;
   int ind;
   int status;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(coef != NULL);

   SCIPdebugMessage("getting binv-row %d\n", r);

   /* check whether we have to reoptimize */
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &status) );
   if ( status == GRB_LOADED || status == GRB_INTERRUPTED || status == GRB_INPROGRESS )
   {
      SCIP_CALL_QUIET( restoreLPData(lpi) );
   }

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );

   /* set up solution vector */
   x.len = 0;
   SCIP_ALLOC( BMSallocMemoryArray(&(x.ind), nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&(x.val), nrows) );

   /* get basis indices, temporarily using memory of x.ind */
   SCIP_CALL( SCIPlpiGetBasisInd(lpi, x.ind) );

   /* set up rhs */
   b.len = 1;
   ind = r;
   val = (x.ind)[r] >= 0 ? 1.0 : -1.0;
   b.ind = &ind;
   b.val = &val;

   /* solve B^T x = e_r, which results in the r-th row of the basis inverse */
   CHECK_ZERO( lpi->messagehdlr, GRBBSolve(lpi->grbmodel, &b, &x) );

   /* size should be at most the number of rows */
   assert( x.len <= nrows );

   /* check whether we require a dense or sparse result vector */
   if ( ninds != NULL && inds != NULL )
   {
      int idx;
      int i;

      /* copy sparse solution */
      for (i = 0; i < x.len; ++i)
      {
         idx = (x.ind)[i];
         assert( idx >= 0 && idx < nrows );
         inds[i] = idx;
         coef[idx] = (x.val)[i];
      }
      *ninds = x.len;
   }
   else
   {
      int idx;
      int i;

      /* copy solution to dense vector */
      BMSclearMemoryArray(coef, nrows);
      for (i = 0; i < x.len; ++i)
      {
         idx = (x.ind)[i];
         assert( idx >= 0 && idx < nrows );
         coef[idx] = (x.val)[i];
      }
   }

   /* free solution space */
   BMSfreeMemoryArray(&(x.val));
   BMSfreeMemoryArray(&(x.ind));

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
   SVECTOR x;
   SVECTOR b;
   int* bind;
   int nrows;
   double val;
   int ind;
   int status;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(coef != NULL);

   SCIPdebugMessage("getting binv-col %d\n", c);

   /* check whether we have to reoptimize */
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &status) );
   if ( status == GRB_LOADED || status == GRB_INTERRUPTED || status == GRB_INPROGRESS )
   {
      SCIP_CALL_QUIET( restoreLPData(lpi) );
   }

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );

   /* set up solution vector */
   x.len = 0;
   SCIP_ALLOC( BMSallocMemoryArray(&(x.ind), nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&(x.val), nrows) );

   /* set up rhs */
   b.len = 1;
   ind = c;
   val = 1.0;
   b.ind = &ind;
   b.val = &val;

   /* solve B x = e_c, which results in the c-th columns of the basis inverse */
   CHECK_ZERO( lpi->messagehdlr, GRBFSolve(lpi->grbmodel, &b, &x) );

   /* size should be at most the number of rows */
   assert( x.len <= nrows );

   /* get basis indices: entries that correspond to slack variables with coefficient -1 must be negated */
   SCIP_ALLOC( BMSallocMemoryArray(&bind, nrows) );
   SCIP_CALL( SCIPlpiGetBasisInd(lpi, bind) );

   /* check whether we require a dense or sparse result vector */
   if ( ninds != NULL && inds != NULL )
   {
      int idx;
      int i;

      /* copy sparse solution */
      for (i = 0; i < x.len; ++i)
      {
         idx = (x.ind)[i];
         assert( idx >= 0 && idx < nrows );
         inds[i] = idx;
         coef[idx] = (x.val)[i];
         if( bind[idx] < 0 )
            coef[idx] *= -1.0;
      }
      *ninds = x.len;
   }
   else
   {
      int idx;
      int i;

      /* copy solution to dense vector */
      BMSclearMemoryArray(coef, nrows);
      for (i = 0; i < x.len; ++i)
      {
         idx = (x.ind)[i];
         assert( idx >= 0 && idx < nrows );
         coef[idx] = (x.val)[i];
         if( bind[idx] < 0 )
            coef[idx] *= -1.0;
      }
   }

   /* free solution space and basis index array */
   BMSfreeMemoryArray(&bind);
   BMSfreeMemoryArray(&(x.val));
   BMSfreeMemoryArray(&(x.ind));

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
   SVECTOR x;
   int nrows;
   int ncols;
   int ngrbcols;
   int status;
   SCIP_Bool isslackvar;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(coef != NULL);
   SCIP_UNUSED( binvrow );

   SCIPdebugMessage("getting binv-row %d\n", r);

   /* check whether we have to reoptimize */
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &status) );
   if ( status == GRB_LOADED || status == GRB_INTERRUPTED || status == GRB_INPROGRESS )
   {
      SCIP_CALL_QUIET( restoreLPData(lpi) );
   }

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ngrbcols) );
   assert( r >= 0 && r < nrows );

   x.len = 0;
   SCIP_ALLOC( BMSallocMemoryArray(&(x.ind), ngrbcols + nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&(x.val), ngrbcols + nrows) );

   /* get basis indices, temporarily using memory of x.ind: if r corresponds to a slack variable with coefficient -1 we
    * have to negate all values
    */
   SCIP_CALL( SCIPlpiGetBasisInd(lpi, x.ind) );
   isslackvar = ((x.ind)[r] < 0);

   /* retrieve row */
   CHECK_ZERO( lpi->messagehdlr, GRBBinvRowi(lpi->grbmodel, r, &x) );

   /* size should be at most the number of columns plus rows for slack variables */
   assert( x.len <= ngrbcols + nrows );

   /* check whether we require a dense or sparse result vector */
   if ( ninds != NULL && inds != NULL )
   {
      int idx;
      int k;
      int j;

      /* Copy sparse solution: Column indices ngrbcols and larger correspond to slack variables artificially introduced
       * by Gurobi; column indices ncols, ncols+1, ..., ngrbcols-1 correspond to slack variables introduced by the LPI
       * implementation. Both must simply be ignored.
       */
      k = 0;
      for (j = 0; j < x.len; ++j)
      {
         idx = (x.ind)[j];
         assert( idx >= 0 && idx < ngrbcols+nrows );
         if ( idx < ncols )
         {
            inds[k++] = idx;
            coef[idx] = (x.val)[j];
            if( isslackvar )
               coef[idx] *= -1.0;
         }
      }
      *ninds = k;
   }
   else
   {
      int idx;
      int j;

      /* Copy dense solution (see comment above). */
      BMSclearMemoryArray(coef, ncols);
      for (j = 0; j < x.len; ++j)
      {
         idx = (x.ind)[j];
         assert( idx >= 0 && idx < ngrbcols+nrows );
         if ( idx < ncols )
         {
            coef[idx] = (x.val)[j];
            if( isslackvar )
               coef[idx] *= -1.0;
         }
      }
   }

   /* free solution space */
   BMSfreeMemoryArray(&(x.val));
   BMSfreeMemoryArray(&(x.ind));

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
   SVECTOR x;
   int* bind;
   int nrows;
   int status;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(coef != NULL);

   SCIPdebugMessage("getting binv-col %d\n", c);

   /* check whether we have to reoptimize */
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &status) );
   if ( status == GRB_LOADED || status == GRB_INTERRUPTED || status == GRB_INPROGRESS )
   {
      SCIP_CALL_QUIET( restoreLPData(lpi) );
   }

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );

   x.len = 0;
   SCIP_ALLOC( BMSallocMemoryArray(&(x.ind), nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&(x.val), nrows) );

   CHECK_ZERO( lpi->messagehdlr, GRBBinvColj(lpi->grbmodel, c, &x) );

   /* size should be at most the number of rows */
   assert( x.len <= nrows );

   /* get basis indices: entries that correspond to slack variables with coefficient -1 must be negated */
   SCIP_ALLOC( BMSallocMemoryArray(&bind, nrows) );
   SCIP_CALL( SCIPlpiGetBasisInd(lpi, bind) );

   /* check whether we require a dense or sparse result vector */
   if ( ninds != NULL && inds != NULL )
   {
      int idx;
      int j;

      /* copy sparse solution */
      for (j = 0; j < x.len; ++j)
      {
         idx = (x.ind)[j];
         assert( idx >= 0 && idx < nrows );
         inds[j] = idx;
         coef[idx] = (x.val)[j];
         if( bind[idx] < 0 )
            coef[idx] *= -1.0;
      }
      *ninds = x.len;
   }
   else
   {
      int idx;
      int j;

      /* copy dense solution */
      BMSclearMemoryArray(coef, nrows);
      for (j = 0; j < x.len; ++j)
      {
         idx = (x.ind)[j];
         assert( idx >= 0 && idx < nrows );
         coef[idx] = (x.val)[j];
         if( bind[idx] < 0 )
            coef[idx] *= -1.0;
      }
   }

   /* free solution space and basis index array */
   BMSfreeMemoryArray(&bind);
   BMSfreeMemoryArray(&(x.val));
   BMSfreeMemoryArray(&(x.ind));

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
   SCIP_Bool success;
   int ncols;
   int nrows;

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpistate != NULL);

   /* if there is no basis information available, no state can be saved */
   if( !lpi->solisbasic )
   {
      *lpistate = NULL;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   assert(ncols >= 0);
   assert(nrows >= 0);

   /* get unpacked basis information from Gurobi */
   SCIP_CALL( getBase(lpi, &success) );

   if ( success )
   {
      /* allocate lpistate data */
      SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows, lpi->nrngrows) );
      (*lpistate)->ncols = ncols;
      (*lpistate)->nrows = nrows;
      (*lpistate)->nrngrows = lpi->nrngrows;

      SCIPdebugMessage("stored Gurobi LPI state in %p (%d cols, %d rows, %d ranged rows)\n",
         (void*) *lpistate, ncols, nrows, lpi->nrngrows);

      /* pack LPi state data */
      lpistatePack(*lpistate, lpi->cstat, lpi->rstat);
   }
   else
   {
      /* In this case no basis information is available. Since SCIP expects the information to work in any case, we
       * allocate the lpistate, but do not use the packed information. This might happen if the model is infeasible,
       * since Gurobi currently does not return basis information in this case. */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, lpistate) );
      (*lpistate)->ncols = ncols;
      (*lpistate)->nrows = nrows;
      (*lpistate)->nrngrows = lpi->nrngrows;
      (*lpistate)->packrstat = NULL;
      (*lpistate)->packcstat = NULL;
   }

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
   int ncols;
   int nrows;
   int i;

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   /* if there was no basis information available, the LPI state was not stored */
   if( lpistate == NULL || lpistate->packrstat == NULL || lpistate->packcstat == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   assert(lpistate->ncols <= ncols);
   assert(lpistate->nrows <= nrows);
   assert(lpistate->nrngrows <= lpi->nrngrows);

   SCIPdebugMessage("loading LPI state %p (%d cols, %d rows, %d ranged rows) into Gurobi LP with %d cols, %d rows, and %d ranged rows\n",
      (void*) lpistate, lpistate->ncols, lpistate->nrows, lpistate->nrngrows, ncols, nrows, lpi->nrngrows);

   if( lpistate->ncols == 0 || lpistate->nrows == 0 )
      return SCIP_OKAY;

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols + lpi->nrngrows) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->cstat, lpi->rstat);

   if ( lpistate->nrngrows > 0 && lpistate->ncols < ncols )
   {
      /* New columns have been added: need to move range variable information */
      memmove(&lpi->cstat[ncols], &lpi->cstat[lpistate->ncols], (size_t) lpistate->nrngrows * sizeof(*lpi->cstat)); /*lint !e571*/
   }

   /* extend the basis to the current LP beyond the previously existing columns */
   for( i = lpistate->ncols; i < ncols; ++i )
   {
      SCIP_Real bnd;
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_LB, i, &bnd) );
      if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_UB, i, &bnd) );
         if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
            lpi->cstat[i] = (int) SCIP_BASESTAT_ZERO;  /* variable is free */
         else
            lpi->cstat[i] = (int) SCIP_BASESTAT_UPPER; /* use finite upper bound */
      }
      else
         lpi->cstat[i] = (int) SCIP_BASESTAT_LOWER;    /* use finite lower bound */
   }
   for( i = lpistate->nrngrows; i < lpi->nrngrows; ++i )
      lpi->cstat[ncols + i] = (int) SCIP_BASESTAT_LOWER;
   for( i = lpistate->nrows; i < nrows; ++i )
      lpi->rstat[i] = (int) SCIP_BASESTAT_BASIC;

   /* load basis information into Gurobi */
   SCIP_CALL( setBase(lpi) );

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   CHECK_ZERO( lpi->messagehdlr, GRBresetmodel(lpi->grbmodel) );

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
      lpistateFree(lpistate, blkmem);

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information), or NULL */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   return (lpistate != NULL && lpistate->packcstat != NULL);
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   size_t l;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(fname != NULL);

   SCIPdebugMessage("reading LP state from file <%s>\n", fname);

   /* gurobi reads a basis if the extension is ".bas" */
   l = strlen(fname);
   if ( l > 4 && fname[l-4] == '.' && fname[l-3] == 'b' && fname[l-2] == 'a' && fname[l-1] == 's' )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBread(lpi->grbmodel, fname) );
   }
   else
   {
      SCIPerrorMessage("To read a basis with gurobi, the extension has to be '.bas'.\n");
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** writes LPi state (i.e. basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   size_t l;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(fname != NULL);

   SCIPdebugMessage("writing basis state to file <%s>\n", fname);

   /* gurobi writes the basis if the extension is ".bas" */
   l = strlen(fname);
   if ( l > 4 && fname[l-4] == '.' && fname[l-3] == 'b' && fname[l-2] == 'a' && fname[l-1] == 's' )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBwrite(lpi->grbmodel, fname) );
   }
   else
   {
      char name[SCIP_MAXSTRLEN];

      /* force extension to be ".bas" */
      if ( strlen(fname) > SCIP_MAXSTRLEN-4)
      {
         SCIPerrorMessage("Basis file name too long.\n");
         return SCIP_LPERROR;
      }
      (void) snprintf(name, SCIP_MAXSTRLEN, "%s.bas", fname);
      CHECK_ZERO( lpi->messagehdlr, GRBwrite(lpi->grbmodel, fname) );
   }

   return SCIP_OKAY;
}

/**@} */




/*
 * LP Pricing Norms Methods
 */

/**@name LP Pricing Norms Methods */
/**@{ */

/** stores LPi pricing norms information */
SCIP_RETCODE SCIPlpiGetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information */
   )
{  /*lint --e{715}*/
   int hasnorm;
   int ncols;
   int nrows;

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpinorms != NULL);

   *lpinorms = NULL;

   /* if there is no basis information available (e.g. after barrier without crossover), norms cannot be saved */
   if( !lpi->solisbasic )
      return SCIP_OKAY;

   /* check if dual norms are available:
    *  value 0: no basis, so no norms available
    *  value 1: basis exists, so norms can be computed
    *  value 2: norms are available
    */
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_HASDUALNORM, &hasnorm) );
   if( hasnorm <= 1 )
      return SCIP_OKAY;

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );

   /* allocate lpinorms data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, lpinorms) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpinorms)->colnorm, ncols) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpinorms)->rownorm, nrows) );
   (*lpinorms)->ncols = ncols;
   (*lpinorms)->nrows = nrows;

   /* query dual norms from Gurobi */
   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_VDUALNORM, 0, ncols, (*lpinorms)->colnorm) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_CDUALNORM, 0, nrows, (*lpinorms)->rownorm) );

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
{  /*lint --e{715}*/
   int error;

   assert(blkmem != NULL);
   assert(lpi != NULL);

   /* if there was no pricing norms information available, the LPI norms were not stored */
   if( lpinorms == NULL )
      return SCIP_OKAY;

   /* store dual norms in Gurobi */
   error = GRBsetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_VDUALNORM, 0, lpinorms->ncols, lpinorms->colnorm);
   /* it can fail to set the norms if no basis was previously set, e.g.,
    * this can happen if flushing an LP did not change anything and
    * therefore no basis was set, as a result Gurobi has no extra user
    * warmstart information and cannot set norms */
#ifdef SCIP_DEBUG
   if( error )
      SCIPmessagePrintWarning(lpi->messagehdlr, "Warning: setting dual variable norms failed with Gurobi error %d\n", error);
#else
   (void)error;
#endif

   error = GRBsetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_CDUALNORM, 0, lpinorms->nrows, lpinorms->rownorm);
   /* it can fail to set the norms if no basis was previously set, e.g.,
    * this can happen if flushing an LP did not change anything and
    * therefore no basis was set, as a result Gurobi has no extra user
    * warmstart information and cannot set norms */
#ifdef SCIP_DEBUG
   if( error )
      SCIPmessagePrintWarning(lpi->messagehdlr, "Warning: setting dual constraint norms failed with Gurobi error %d\n", error);
#else
   (void)error;
#endif

   return SCIP_OKAY;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPlpiFreeNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information, or NULL */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpinorms != NULL);

   if ( *lpinorms != NULL )
   {
      BMSfreeBlockMemoryArray(blkmem, &(*lpinorms)->colnorm, (*lpinorms)->ncols);
      BMSfreeBlockMemoryArray(blkmem, &(*lpinorms)->rownorm, (*lpinorms)->nrows);
      BMSfreeBlockMemory(blkmem, lpinorms);
   }

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
   int temp;
   SCIP_Real dtemp;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(ival != NULL);

   SCIPdebugMessage("getting int parameter %d\n", type);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      *ival = (int) lpi->fromscratch;
      break;
   case SCIP_LPPAR_FASTMIP:
      /* maybe set perturbation */
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_SCALING:
      SCIP_CALL( getIntParam(lpi, GRB_INT_PAR_SCALEFLAG, &temp) );
      assert(temp >= -1 && temp <= 3);
      if( temp == -1 )
         *ival = 1;
      else
         *ival = temp;
      break;
   case SCIP_LPPAR_PRESOLVING:
      SCIP_CALL( getIntParam(lpi, GRB_INT_PAR_PRESOLVE, &temp) );
      assert( temp == GRB_PRESOLVE_AUTO || temp == GRB_PRESOLVE_OFF || temp == GRB_PRESOLVE_CONSERVATIVE || temp == GRB_PRESOLVE_AGGRESSIVE );
      *ival = (temp == GRB_PRESOLVE_OFF) ? FALSE : TRUE;
      break;
   case SCIP_LPPAR_PRICING:
      *ival = (int) lpi->pricing;
      break;
   case SCIP_LPPAR_LPINFO:
      SCIP_CALL( getIntParam(lpi, GRB_INT_PAR_OUTPUTFLAG, &temp) );
      assert( temp == 0 || temp == 1 );
      *ival = (temp == 1) ? TRUE : FALSE;
      break;
   case SCIP_LPPAR_LPITLIM:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_ITERATIONLIMIT, &dtemp) );
      assert( dtemp >= 0.0 );
      if( dtemp >= INT_MAX )
         *ival = INT_MAX;
      else
         *ival = (int) dtemp;
      break;
   case SCIP_LPPAR_THREADS:
      SCIP_CALL( getIntParam(lpi, GRB_INT_PAR_THREADS, ival) );
      break;
   case SCIP_LPPAR_RANDOMSEED:
      SCIP_CALL( getIntParam(lpi, GRB_INT_PAR_SEED, ival) );
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
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("setting int parameter %d to %d\n", type, ival);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      assert(ival == TRUE || ival == FALSE);
      lpi->fromscratch = (SCIP_Bool) ival;
      break;
   case SCIP_LPPAR_FASTMIP:
      assert(ival == TRUE || ival == FALSE);
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_SCALING:
      if( ival == 1 )
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SCALEFLAG, -1) );
      else
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SCALEFLAG, ival) );
      break;
   case SCIP_LPPAR_PRESOLVING:
      assert(ival == TRUE || ival == FALSE);
      if( ival )
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_AUTO) );
      else
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF) );
      break;
   case SCIP_LPPAR_PRICING:
      lpi->pricing = (SCIP_PRICING)ival;
      switch( (SCIP_PRICING)ival )
      {
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_AUTO:
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_AUTO) );
         break;
      case SCIP_PRICING_FULL:
         /* full does not seem to exist -> use auto */
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_AUTO) );
         break;
      case SCIP_PRICING_PARTIAL:
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_PARTIAL) );
         break;
      case SCIP_PRICING_STEEP:
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_STEEPEST_EDGE) );
         break;
      case SCIP_PRICING_STEEPQSTART:
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_STEEPEST_QUICK) );
         break;
      case SCIP_PRICING_DEVEX:
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_DEVEX) );
         break;
      default:
         return SCIP_PARAMETERUNKNOWN;
      }
      break;
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      if( ival )
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_OUTPUTFLAG, 1) );
      else
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_OUTPUTFLAG, 0) );
      break;
   case SCIP_LPPAR_LPITLIM:
      assert( ival >= 0 );
      /* 0 <= ival, 0 stopping immediately */
      {
         double itlim;
         itlim = (ival >= INT_MAX ? GRB_INFINITY : ival);
         SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_ITERATIONLIMIT, itlim) );
      }
      break;
   case SCIP_LPPAR_THREADS:
      assert( ival >= 0 );
      SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_THREADS, ival) );
      break;
   case SCIP_LPPAR_RANDOMSEED:
      assert( ival >= 0 );
      SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SEED, ival) );
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
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(dval != NULL);

   SCIPdebugMessage("getting real parameter %d\n", type);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_FEASIBILITYTOL, dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_OPTIMALITYTOL, dval) );
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_BARCONVTOL, dval) );
      break;
   case SCIP_LPPAR_OBJLIM:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_CUTOFF, dval) );
      break;
   case SCIP_LPPAR_LPTILIM:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_TIMELIMIT, dval) );
      break;
   case SCIP_LPPAR_MARKOWITZ:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_MARKOWITZTOL, dval) );
      break;
   case SCIP_LPPAR_CONDITIONLIMIT:
      *dval = lpi->conditionlimit;
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
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("setting real parameter %d to %g\n", type, dval);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      assert( dval > 0.0 );
      /* 1e-9 <= dval <= 1e-2 */
      if( dval < 1e-9 )
         dval = 1e-9;
      else if( dval > 1e-2 )
         dval = 1e-2;

      SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_FEASIBILITYTOL, dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      assert( dval > 0.0 );
      /* 1e-9 <= dval <= 1e-2 */
      if (dval < 1e-9)
         dval = 1e-9;
      else if( dval > 1e-2 )
         dval = 1e-2;

      SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_OPTIMALITYTOL, dval) );
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      /* 0 <= dval <= 1 */
      assert( dval >= 0.0 );
      if( dval > 1.0 )
         dval = 1.0;

      SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_BARCONVTOL, dval) );
      break;
   case SCIP_LPPAR_OBJLIM:
      /* no restriction on dval */

      SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_CUTOFF, dval) );
      break;
   case SCIP_LPPAR_LPTILIM:
      assert( dval > 0.0 );
      /* gurobi requires 0 <= dval
       *
       * However for consistency we assert the timelimit to be strictly positive.
       */

      SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_TIMELIMIT, dval) );
      break;
   case SCIP_LPPAR_MARKOWITZ:
      /* 1e-4 <= dval <= 0.999 */
      if( dval < 1e-4 )
         dval = 1e-4;
      else if( dval > 0.999 )
         dval = 0.999;

      SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_MARKOWITZTOL, dval) );
      break;
   case SCIP_LPPAR_CONDITIONLIMIT:
      lpi->conditionlimit = dval;
      lpi->checkcondition = (dval >= 0.0) ? TRUE : FALSE;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** interrupts the currently ongoing lp solve or disables the interrupt */ /*lint -e{715}*/
SCIP_RETCODE SCIPlpiInterrupt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             interrupt           /**< TRUE if interrupt should be set, FALSE if it should be disabled */
   )
{  /*lint --e{715}*/
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
{  /*lint --e{715}*/
   assert(lpi != NULL);
   return GRB_INFINITY;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   return (val >= GRB_INFINITY);
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(fname != NULL);

   SCIPdebugMessage("reading LP from file <%s>\n", fname);

   CHECK_ZERO( lpi->messagehdlr, GRBreadmodel(lpi->grbenv, fname, &lpi->grbmodel) );

   /* the model name seems to be empty, use filename */
   CHECK_ZERO( lpi->messagehdlr, GRBsetstrattr(lpi->grbmodel, GRB_STR_ATTR_MODELNAME, fname) );

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(fname != NULL);

   SCIPdebugMessage("writing LP to file <%s>\n", fname);

   /* if range rows were not added, add, print and remove them; otherwise, just print */
   if ( lpi->nrngrows > 0 && !lpi->rngvarsadded )
   {
      SCIP_CALL( addRangeVars(lpi) );
      CHECK_ZERO( lpi->messagehdlr, GRBwrite(lpi->grbmodel, fname) );
      SCIP_CALL( delRangeVars(lpi) );
   }
   else
   {
      CHECK_ZERO( lpi->messagehdlr, GRBwrite(lpi->grbmodel, fname) );
   }

   return SCIP_OKAY;
}

/**@} */
