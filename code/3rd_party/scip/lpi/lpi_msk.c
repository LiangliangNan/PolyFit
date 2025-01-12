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

/**@file   lpi_msk.c
 * @ingroup LPIS
 * @brief  LP interface for MOSEK
 * @author Bo Jensen
 * @author Tristan Gally
 * @author Marc Pfetsch
 *
 * @todo Check whether MSK_IPAR_{SIM_DUAL|PRIMAL}_RESTRICT_SELECTION should be used if problem is solved from scratch or
 *       if no basis is available.
 * @todo Revise handling of the MSK_RES_TRM_MAX_NUM_SETBACKS return value: Remove it form the check of MOSEK_CALL and
 *       include it in filterTRMrescode().
 * @todo Check whether SCIPlpiGetSolFeasibility() should also return primal/dual feasible if the status is
 *       MSK_SOL_STA_NEAR_PRIM_FEAS, MSK_SOL_STA_NEAR_DUAL_FEAS.
 * @todo Check why it can happen that the termination code is MSK_RES_OK, but the solution status is MSK_SOL_STA_UNKNOWN.
 */

/*lint -e750*/
/*lint -e830*/

#include <assert.h>

#define MSKCONST const    /* this define is needed for older MOSEK versions */
#include "mosek.h"

#include "lpi/lpi.h"
#include "scip/bitencode.h"
#include "scip/pub_message.h"
#include <string.h>
#include "tinycthread/tinycthread.h"

/* do defines for windows directly here to make the lpi more independent */
#if defined(_WIN32) || defined(_WIN64)
#define snprintf _snprintf
#endif

#if ( MSK_VERSION_MAJOR < 7 )
#error "The MOSEK interface only works for MOSEK versions 7.0.0.0 and newer"
#endif

#define scipmskobjsen MSKobjsensee
#define SENSE2MOSEK(objsen) (((objsen)==SCIP_OBJSEN_MINIMIZE)?(MSK_OBJECTIVE_SENSE_MINIMIZE):(MSK_OBJECTIVE_SENSE_MAXIMIZE))

typedef enum MSKoptimizertype_enum MSKoptimizertype;

#define MOSEK_CALL(x)  do                                               \
                       {  /*lint --e{641}*/                                                                   \
                          MSKrescodee _restat_;                                                               \
                          _restat_ = (x);                                                                     \
                          if( (_restat_) != MSK_RES_OK && (_restat_ ) != MSK_RES_TRM_MAX_NUM_SETBACKS )       \
                          {                                                                                   \
                             SCIPerrorMessage("LP Error: MOSEK returned %d.\n", (int)_restat_);               \
                             return SCIP_LPERROR;                                                             \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

/* this macro is only called in functions returning SCIP_Bool; thus, we return FALSE if there is an error in optimized mode */
#define SCIP_ABORT_FALSE(x) do                                          \
   {                                                                    \
      SCIP_RETCODE _restat_;                                            \
      if( (_restat_ = (x)) != SCIP_OKAY )                               \
      {                                                                 \
         SCIPerrorMessage("LP Error: MOSEK returned %d.\n", (int)_restat_); \
         SCIPABORT();                                                   \
         return FALSE;                                                  \
      }                                                                 \
   }                                                                    \
   while( FALSE )

#define IS_POSINF(x) ((x) >= MSK_INFINITY)
#define IS_NEGINF(x) ((x) <= -MSK_INFINITY)

#ifdef SCIP_THREADSAFE
   #if defined(_Thread_local)
      /* Use thread local environment in order to not create a new environment for each new LP. */
      static _Thread_local MSKenv_t reusemosekenv =     NULL;
      static _Thread_local int numlp         =           0;
      #define SCIP_REUSEENV
   #endif
#else
   /* Global Mosek environment in order to not create a new environment for each new LP. This is not thread safe. */
   static MSKenv_t reusemosekenv =     NULL;
   static int numlp         =           0;
   #define SCIP_REUSEENV
#endif

#if MSK_VERSION_MAJOR >= 9
#define NEAR_REL_TOLERANCE           1.0     /* MOSEK will multiply all tolerances with this factor after stalling */
#endif
#define DEBUG_PRINT_STAT             0
#define DEBUG_PARAM_SETTING          0
#define DEBUG_CHECK_DATA             0
#define DEBUG_EASY_REPRODUCE         0
#define DEBUG_DO_INTPNT_FEAS_CHECK   0
#define DEBUG_CHECK_STATE_TOL        1e-5
#define SHOW_ERRORS                  0
#define SHOW_RELATIVE_OPTIMAL_GAP    0
#define ASSERT_ON_NUMERICAL_TROUBLES 0
#define ASSERT_ON_WARNING            0
#define FORCE_MOSEK_LOG              0       /* note that changing this AND setting lpinfo will lead to asserts in lpCheckIntpar */
#define FORCE_MOSEK_SUMMARY          0
#define FORCE_NO_MAXITER             0
#define SETBACK_LIMIT                250
#define STRONGBRANCH_PRICING         MSK_SIM_SELECTION_SE
#define SUPRESS_NAME_ERROR           1
#define WRITE_DUAL                   0
#define WRITE_PRIMAL                 0
#define WRITE_INTPNT                 0
#if WRITE_DUAL > 0 || WRITE_PRIMAL > 0 || WRITE_INTPNT > 0 || FORCE_MOSEK_LOG > 0 || FORCE_MOSEK_SUMMARY > 0
#define WRITE_ABOVE                  0
#endif
#define DEGEN_LEVEL                  MSK_SIM_DEGEN_FREE
#define ALWAYS_SOLVE_PRIMAL_FORM     1
#if DEBUG_PRINT_STAT > 0
static int numstrongbranchmaxiterup =  0;
static int numstrongbranchmaxiterdo =  0;
static int numprimalmaxiter         =  0;
static int numdualmaxiter           =  0;
static int numstrongbranchobjup     =  0;
static int numstrongbranchobjdo     =  0;
static int numprimalobj             =  0;
static int numdualobj               =  0;
#endif

#if DEBUG_PRINT_STAT > 0
static int numstrongbranchmaxiterup =  0;
static int numstrongbranchmaxiterdo =  0;
static int numprimalmaxiter         =  0;
static int numdualmaxiter           =  0;
static int numstrongbranchobjup     =  0;
static int numstrongbranchobjdo     =  0;
static int numprimalobj             =  0;
static int numdualobj               =  0;
#endif


/** internal data for Mosek LPI */
struct SCIP_LPi
{
   MSKenv_t              mosekenv;           /**< Mosek environment */
   MSKtask_t             task;               /**< Mosek task */
   int                   optimizecount;      /**< optimization counter (mainly for debugging) */
   MSKrescodee           termcode;           /**< termination code of last optimization run */
   int                   itercount;          /**< iteration count of last optimization run */
   SCIP_PRICING          pricing;            /**< SCIP pricing setting */
   int                   lpid;               /**< id for LP within same task */
   MSKoptimizertype      lastalgo;           /**< algorithm type of last solving call */
   MSKstakeye*           skx;                /**< basis status for columns */
   MSKstakeye*           skc;                /**< basis status for rows */
   MSKboundkeye*         bkx;                /**< bound keys for columns */
   MSKboundkeye*         bkc;                /**< bound keys for rows */
   MSKint32t*            aptre;              /**< row or column end pointers */
   int                   skxsize;            /**< size of skx array */
   int                   skcsize;            /**< size of skc array */
   int                   bkxsize;            /**< size of bkx */
   int                   bkcsize;            /**< size of bkx */
   int                   aptresize;          /**< size of aptre */
   MSKsoltypee           lastsolvetype;      /**< Which solver was called last and which solution should be returned? */
   SCIP_Bool             solved;             /**< was the current LP solved? */
   SCIP_Bool             fromscratch;        /**< Shall solves be performed with MSK_IPAR_SIM_HOTSTART turned off? */
   SCIP_Bool             clearstate;         /**< Shall next solve be performed with MSK_IPAR_SIM_HOTSTART turned off? */
   SCIP_Bool             lpinfo;             /**< Should LP solver output information to the screen? */
   int                   restrictselectdef;  /**< default value for MSK_IPAR_SIM_DUAL_RESTRICT_SELECTION */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehdlr handler to printing messages, or NULL */
};

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE

/** basis status */
struct SCIP_LPiState
{
   int                   num;
   MSKsolstae            solsta;             /**< solution status */
   int                   ncols;              /**< number of columns */
   int                   nrows;              /**< number of rows */
   COLPACKET*            skx;                /**< basis status for columns */
   ROWPACKET*            skc;                /**< basis status for rows */
};


/*
 * Local functions
 */

/** gives problem and solution status for a Mosek Task
 *
 *  With Mosek 7.0, the routine MSK_getsolutionstatus was replaced by MSK_getprosta and MSK_getsolsta.
 */
static
MSKrescodee MSK_getsolutionstatus(
   MSKtask_t             task,               /**< Mosek Task */
   MSKsoltypee           whichsol,           /**< for which type of solution a status is requested */
   MSKprostae*           prosta,             /**< buffer to store problem status, or NULL if not needed */
   MSKsolstae*           solsta              /**< buffer to store solution status, or NULL if not needed */
   )
{
   MSKrescodee res;

   if( prosta != NULL )
   {
      res = MSK_getprosta(task, whichsol, prosta);
      if ( res != MSK_RES_OK )
         return res;
   }
   if( solsta != NULL )
   {
      res = MSK_getsolsta(task, whichsol, solsta);
      if ( res != MSK_RES_OK )
         return res;
   }

   return MSK_RES_OK;
}

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

/** print string using message handler of SCIP */
static
void MSKAPI printstr(
   MSKuserhandle_t       handle,             /**< error handle */
   const char*           str                 /**< string that contains string on output */
   )
{  /*lint --e{715}*/
#if SUPRESS_NAME_ERROR
   char errstr[32];
   (void) snprintf(errstr, 32, "MOSEK Error %d", MSK_RES_ERR_DUP_NAME);
   if (0 == strncmp(errstr, str, strlen(errstr)))
      return;
#endif

   SCIPmessagePrintInfo((SCIP_MESSAGEHDLR *) handle, "MOSEK: %s", str);
}

#if DEBUG_CHECK_DATA > 0
/** check data */
static
SCIP_RETCODE scip_checkdata(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   const char*           functionname        /**< function name */
   )
{
   int i;
   int numcon;
   int numvar;
   int gotbasicsol;
   MSKboundkeye* tbkc;
   MSKboundkeye* tbkx;
   MSKstakeye *tskc;
   MSKstakeye* tskx;
   double* tblc;
   double* tbuc;
   double* tblx;
   double* tbux;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   MOSEK_CALL( MSK_solutiondef(lpi->task, MSK_SOL_BAS, &gotbasicsol) );

   MOSEK_CALL( MSK_getnumvar(lpi->task, &numvar) );
   MOSEK_CALL( MSK_getnumcon(lpi->task, &numcon) );

   /* allocate memory */
   SCIP_ALLOC( BMSallocMemoryArray(&tbkc, numcon) );
   SCIP_ALLOC( BMSallocMemoryArray(&tskc, numcon) );
   SCIP_ALLOC( BMSallocMemoryArray(&tblc, numcon) );
   SCIP_ALLOC( BMSallocMemoryArray(&tbuc, numcon) );

   SCIP_ALLOC( BMSallocMemoryArray(&tbkx, numvar) );
   SCIP_ALLOC( BMSallocMemoryArray(&tskx, numvar) );
   SCIP_ALLOC( BMSallocMemoryArray(&tblx, numvar) );
   SCIP_ALLOC( BMSallocMemoryArray(&tbux, numvar) );

   /* Check bounds */
   if( gotbasicsol )
   {
      MOSEK_CALL( MSK_getsolution(lpi->task, MSK_SOL_BAS, NULL, NULL, tskc, tskx,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
   }

   for( i = 0; i < numvar; i++ )
   {
      MOSEK_CALL( MSK_getbound(lpi->task, MSK_ACC_VAR, i, &tbkx[i], &tblx[i], &tbux[i]) );
   }

   for( i = 0; i < numcon; i++ )
   {
      MOSEK_CALL( MSK_getbound(lpi->task, MSK_ACC_CON, i, &tbkc[i], &tblc[i], &tbuc[i]) );
   }

   for( i = 0; i < numcon; ++i )
   {
      if( gotbasicsol )
      {
         if( ( tskc[i] == MSK_SK_FIX && tbkc[i] != MSK_BK_FX ) ||
            ( tskc[i] == MSK_SK_LOW && !(tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_LO || tbkc[i] == MSK_BK_RA ) ) ||
            ( tskc[i] == MSK_SK_UPR && !(tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_UP || tbkc[i] == MSK_BK_RA ) ) )
         {
            SCIPerrorMessage("STATUS KEY ERROR i %d bkc %d skc %d %s\n", i, tbkc[i], tskc[i], functionname);
         }
      }

      if( tbkc[i] == MSK_BK_LO || tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_RA )
      {
         if( isnan(tblc[i]) )
         {
            SCIPdebugMessage("nan in blc : %s\n", functionname);
         }
      }

      if( tbkc[i] == MSK_BK_UP || tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_RA )
      {
         if( isnan(tbuc[i]) )
         {
            SCIPdebugMessage("nan in bux : %s\n", functionname);
         }
      }
   }

   for( i = 0; i < numvar; ++i )
   {
      if( tbkx[i] == MSK_BK_LO || tbkx[i] == MSK_BK_FX || tbkx[i] == MSK_BK_RA )
      {
         if( isnan(tblx[i]) )
         {
            SCIPdebugMessage("nan in blx : %s\n", functionname);
         }
      }

      if( tbkx[i] == MSK_BK_UP || tbkx[i] == MSK_BK_FX || tbkx[i] == MSK_BK_RA )
      {
         if( isnan(tbux[i]) )
         {
            SCIPdebugMessage("nan in bux : %s\n", functionname);
            getchar();
         }
      }
   }

   BMSfreeMemoryArray(&tbkc);
   BMSfreeMemoryArray(&tskc);
   BMSfreeMemoryArray(&tblc);
   BMSfreeMemoryArray(&tbuc);
   BMSfreeMemoryArray(&tbkx);
   BMSfreeMemoryArray(&tskx);
   BMSfreeMemoryArray(&tblx);
   BMSfreeMemoryArray(&tbux);

   return SCIP_OKAY;
}
#endif

/** resizes bound keys array bkx to have at least ncols entries */
static
SCIP_RETCODE ensureBkxMem(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   ncols               /**< number of columns */
   )
{
   if ( lpi->bkxsize < ncols )
   {
      int newsize;
      newsize = MAX(2*lpi->bkxsize, ncols);

      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->bkx), newsize) );
      lpi->bkxsize = newsize;
   }

   return SCIP_OKAY;
}

/** resizes bound keys array bkc to have at least nrows entries */
static
SCIP_RETCODE ensureBkcMem(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   nrows               /**< number of rows */
   )
{
   if ( lpi->bkcsize < nrows )
   {
      int newsize;
      newsize = MAX(2*lpi->bkcsize, nrows);

      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->bkc), newsize) );
      lpi->bkcsize = newsize;
   }

   return SCIP_OKAY;
}

/** resizes aptre array to have at least n entries */
static
SCIP_RETCODE ensureAptreMem(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   n                   /**< number of entries */
   )
{
   if ( lpi->aptresize < n )
   {
      int newsize;
      newsize = MAX(2*lpi->aptresize, n);

      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->aptre), newsize) );
      lpi->aptresize = newsize;
   }

   return SCIP_OKAY;
}

/** marks the current LP to be unsolved */
static
void invalidateSolution(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);

   lpi->solved = FALSE;
}

/** compute boundkeys to inform MOSEK about fixed/free/ranged/lower bounded/upper bounded variables or constraints */
static
void generateMskBoundkeys(
   int                   n,                  /**< array size */
   const double*         lb,                 /**< lower bounds of variables or left-hand sides of ranged rows */
   const double*         ub,                 /**< upper bounds of variables or right-hand sides of ranged rows */
   MSKboundkeye*         bk                  /**< pointer to store boundkeys to inform MOSEK about status of var/row */
   )
{
   int i;

   assert(lb != NULL);
   assert(ub != NULL);
   assert(bk != NULL);

   for( i = 0; i < n; i++ )
   {
      if (IS_NEGINF(lb[i]))
      {
         if (IS_POSINF(ub[i]))
         {
            bk[i] = MSK_BK_FR;
         }
         else
         {
            assert(!IS_NEGINF(ub[i]));
            bk[i] = MSK_BK_UP;
         }
      }
      else
      {
         assert(!IS_POSINF(lb[i]));
         if (IS_POSINF(ub[i]))
         {
            bk[i] = MSK_BK_LO;
         }
         else if (lb[i] == ub[i])/*lint !e777*/  /* No epsilon-test since MOSEK will also test for exact equality */
         {
            assert(lb[i] - ub[i] == 0);
            assert(ub[i] - lb[i] == 0);
            bk[i] = MSK_BK_FX;
         }
         else
         {
            assert(lb[i] < ub[i]);
            bk[i] = MSK_BK_RA;
         }
      }
   }
}

/** get end pointers of arrays */
static
SCIP_RETCODE getEndptrs(
   int                   n,                  /**< array size */
   const int*            beg,                /**< array of beginning indices */
   int                   nnonz,              /**< number of nonzeros */
   MSKint32t*            aptre               /**< array to store the result */
   )
{
   int i;

   assert(beg != NULL || nnonz == 0);

   if (nnonz > 0)
   {
      assert(beg != NULL);
      for(i = 0; i < n-1; i++)
      {
         aptre[i] = beg[i+1];
         assert(aptre[i] >= beg[i]);
      }

      aptre[n-1] = nnonz;
      assert(aptre[n-1] >= beg[n-1]);
   }
   else
   {
      for( i = 0; i < n; i++ )
         aptre[i] = 0;
   }

   return SCIP_OKAY;
}

/** compute indices from range */
static
SCIP_RETCODE getIndicesRange(
   int                   first,              /**< first index */
   int                   last,               /**< last index */
   int**                 sub                 /**< pointer to store the indices ranges */
   )
{
   int i;

   assert(first <= last);

   SCIP_ALLOC( BMSallocMemoryArray(sub, (last - first + 1)) );

   for( i = first; i <= last; i++ )
   {
      (*sub)[i-first] = i;
   }

   return SCIP_OKAY;
}

/** compute indices from dense array */
static
SCIP_RETCODE getIndicesFromDense(
   int*                  dstat,              /**< array */
   int                   n,                  /**< size of array */
   int*                  count,              /**< array of counts (sizes) */
   int**                 sub                 /**< pointer to store array of indices */
   )
{
   int i;
   int j;

   assert(dstat != NULL);
   assert(count != NULL);

   *count = 0;
   for( i = 0; i < n; i++ )
   {
      if (dstat[i] == 1)
      {
         (*count)++;
      }
   }

   if( (*count) > 0 )
   {
      SCIP_ALLOC( BMSallocMemoryArray(sub, (*count)) );
   }
   else
      return SCIP_OKAY;

   j = 0;
   for( i = 0; i < n; i++ )
   {
      if (dstat[i] == 1)
      {
         (*sub)[j++] = i;
      }
   }

   return SCIP_OKAY;
}

/** scale a vector */
static
void scale_vec(
   int                   len,                /**< length of vector */
   double*               vec,                /**< vector to be scaled */
   double                s                   /**< scaling factor */
   )
{
   int i;
   for( i = 0; i < len; i++ )
   {
      vec[i] *= s;
   }
}

/** scale lower and upper bound */
static
void scale_bound(
   MSKboundkeye*         bk,                 /**< pointer to store boundkeys to inform MOSEK about status of var/row */
   double*               bl,                 /**< lower bound */
   double*               bu,                 /**< upper bound */
   double                s                   /**< scaling factor */
   )
{
   switch( *bk )
   {
   case MSK_BK_LO:
      *bl *= s;
      if ( s < 0.0 )
         *bk = MSK_BK_UP;
      break;
   case MSK_BK_UP:
      *bu *= s;
      if ( s < 0.0 )
         *bk = MSK_BK_LO;
      break;
   case MSK_BK_FX:
   case MSK_BK_RA:
      *bl *= s;
      *bu *= s;
      break;
   case MSK_BK_FR:
      break;
   default:
      SCIPABORT();
      break;
   }  /*lint !e788*/

   /* switch bounds if scaling is negative */
   if ( s < 0.0 )
   {
      double tmp;
      tmp = *bl;
      *bl = *bu;
      *bu = tmp;
   }
}

/** resizes state arrays to have at least ncols/nrows entries */
static
SCIP_RETCODE ensureStateMem(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   ncols,              /**< number of columns */
   int                   nrows               /**< number of rows */
   )
{
   if ( lpi->skxsize < ncols )
   {
      int newsize;
      newsize = MAX(2*lpi->skxsize, ncols);

      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->skx), newsize) );
      lpi->skxsize = newsize;
   }

   if ( lpi->skcsize < nrows )
   {
      int newsize;
      newsize = MAX(2*lpi->skcsize, nrows);

      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->skc), newsize) );
      lpi->skcsize = newsize;
   }

   return SCIP_OKAY;
}

/** get base and store in skc/skx arrays */
static
SCIP_RETCODE getbase(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   ncols,              /**< number of columns */
   int                   nrows               /**< number of rows */
   )
{
   assert(lpi->lastsolvetype == MSK_SOL_BAS);

   SCIPdebugMessage("Calling getbase (%d)\n", lpi->lpid);

   SCIP_CALL( ensureStateMem(lpi, ncols, nrows) );
   MOSEK_CALL( MSK_getsolution(lpi->task, MSK_SOL_BAS, NULL, NULL, lpi->skc, lpi->skx,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}

/** set base to the values given in skc/skx arrays */
static
SCIP_RETCODE setbase(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   SCIPdebugMessage("Calling setbase (%d)\n", lpi->lpid);

   lpi->lastsolvetype = MSK_SOL_BAS;
   lpi->solved = FALSE;

   MOSEK_CALL( MSK_putsolution(lpi->task, MSK_SOL_BAS, lpi->skc, lpi->skx, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}



/*
 * Miscellaneous Methods
 */

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#if MSK_VERSION_MAJOR < 9
   #define mskname "MOSEK " STR(MSK_VERSION_MAJOR) "." STR(MSK_VERSION_MINOR) "." STR(MSK_VERSION_BUILD) "." STR(MSK_VERSION_REVISION)
#else
   #define mskname "MOSEK " STR(MSK_VERSION_MAJOR) "." STR(MSK_VERSION_MINOR) "." STR(MSK_VERSION_REVISION)
#endif

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   return mskname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "Linear Programming Solver developed by MOSEK Optimization Software (www.mosek.com)";
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   return (void*) lpi->task;
}

/** pass integrality information to LP solver */
SCIP_RETCODE SCIPlpiSetIntegralityInformation(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   ncols,              /**< length of integrality array */
   int*                  intInfo             /**< integrality array (0: continuous, 1: integer). May be NULL iff ncols is 0.  */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( ncols >= 0 );
   assert( intInfo != NULL );

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

   SCIPdebugMessage("Calling SCIPlpiCreate\n");

   SCIP_ALLOC( BMSallocMemory(lpi) );

#ifdef SCIP_REUSEENV
   if ( reusemosekenv == NULL )
   {
      assert(numlp == 0);
      MOSEK_CALL( MSK_makeenv(&reusemosekenv, NULL) );
      MOSEK_CALL( MSK_linkfunctoenvstream(reusemosekenv, MSK_STREAM_LOG, (MSKuserhandle_t) messagehdlr, printstr) );
#if MSK_VERSION_MAJOR < 8
      MOSEK_CALL( MSK_initenv(reusemosekenv) );
#endif
   }
   (*lpi)->mosekenv = reusemosekenv;
   (*lpi)->lpid = numlp++;

#else

   MOSEK_CALL( MSK_makeenv(&(*lpi)->mosekenv, NULL) );
   MOSEK_CALL( MSK_linkfunctoenvstream((*lpi)->mosekenv, MSK_STREAM_LOG, (MSKuserhandle_t) messagehdlr, printstr) );
#if MSK_VERSION_MAJOR < 8
   MOSEK_CALL( MSK_initenv((*lpi)->mosekenv) );
#endif
#endif

   MOSEK_CALL( MSK_makeemptytask((*lpi)->mosekenv, &((*lpi)->task)) );

   MOSEK_CALL( MSK_linkfunctotaskstream((*lpi)->task, MSK_STREAM_LOG, (MSKuserhandle_t) messagehdlr, printstr) );

   MOSEK_CALL( MSK_putobjsense((*lpi)->task, SENSE2MOSEK(objsen)) );
   MOSEK_CALL( MSK_putintparam((*lpi)->task, MSK_IPAR_SIM_MAX_NUM_SETBACKS, SETBACK_LIMIT) );
   MOSEK_CALL( MSK_putintparam((*lpi)->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_FREE_SIMPLEX) );
   MOSEK_CALL( MSK_putintparam((*lpi)->task, MSK_IPAR_SIM_DEGEN, DEGEN_LEVEL) );
   MOSEK_CALL( MSK_putintparam((*lpi)->task, MSK_IPAR_SIM_SWITCH_OPTIMIZER, MSK_ON) );
   MOSEK_CALL( MSK_puttaskname((*lpi)->task, (char*) name) );

   /* disable errors for huge values */
   MOSEK_CALL( MSK_putdouparam((*lpi)->task, MSK_DPAR_DATA_TOL_AIJ_HUGE, MSK_INFINITY * 2)); /* not clear why the *2 is needed */
   MOSEK_CALL( MSK_putdouparam((*lpi)->task, MSK_DPAR_DATA_TOL_C_HUGE, MSK_INFINITY));

   /* disable warnings for large values */
   MOSEK_CALL( MSK_putdouparam((*lpi)->task, MSK_DPAR_DATA_TOL_AIJ_LARGE, MSK_INFINITY * 2));
   MOSEK_CALL( MSK_putdouparam((*lpi)->task, MSK_DPAR_DATA_TOL_CJ_LARGE, MSK_INFINITY));

   /* disable warnings for large bounds */
   MOSEK_CALL( MSK_putdouparam((*lpi)->task, MSK_DPAR_DATA_TOL_BOUND_WRN, MSK_INFINITY));

   (*lpi)->termcode = MSK_RES_OK;
   (*lpi)->itercount = 0;
   (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*lpi)->lastalgo = MSK_OPTIMIZER_FREE;
   (*lpi)->skx = NULL;
   (*lpi)->skc = NULL;
   (*lpi)->bkx = NULL;
   (*lpi)->bkc = NULL;
   (*lpi)->aptre = NULL;
   (*lpi)->skxsize = 0;
   (*lpi)->skcsize = 0;
   (*lpi)->bkxsize = 0;
   (*lpi)->bkcsize = 0;
   (*lpi)->aptresize = 0;
   (*lpi)->lastsolvetype = (MSKsoltypee) -1;
   (*lpi)->lpinfo = FALSE;
   (*lpi)->restrictselectdef = 50;
   (*lpi)->fromscratch = FALSE;
   (*lpi)->clearstate = FALSE;
   (*lpi)->messagehdlr = messagehdlr;

   invalidateSolution(*lpi);

   MOSEK_CALL( MSK_putintparam((*lpi)->task, MSK_IPAR_LOG, MSK_OFF) );
   MOSEK_CALL( MSK_putintparam((*lpi)->task, MSK_IPAR_LOG_SIM, MSK_OFF) );

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(*lpi != NULL);

   SCIPdebugMessage("Calling SCIPlpiFree (%d)\n", (*lpi)->lpid);

   MOSEK_CALL( MSK_deletetask(&(*lpi)->task) );

   BMSfreeMemoryArrayNull(&(*lpi)->aptre);
   BMSfreeMemoryArrayNull(&(*lpi)->bkx);
   BMSfreeMemoryArrayNull(&(*lpi)->bkc);
   BMSfreeMemoryArrayNull(&(*lpi)->skx);
   BMSfreeMemoryArrayNull(&(*lpi)->skc);

#ifdef SCIP_REUSEENV
   assert(numlp > 0);
   numlp--;
   if ( numlp == 0 )
   {
      MOSEK_CALL( MSK_deleteenv(&reusemosekenv) );
      reusemosekenv = NULL;
   }
#else
   MOSEK_CALL( MSK_deleteenv(&(*lpi)->mosekenv) );
#endif

   BMSfreeMemory(lpi);

   return SCIP_OKAY;
}

/*
 * Modification Methods
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
{
#ifndef NDEBUG
   {
      int j;
      for( j = 0; j < nnonz; j++ )
         assert( val[j] != 0.0 );
   }
#endif

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(beg != NULL);
   assert(ind != NULL);
   assert(val != NULL);

   SCIPdebugMessage("Calling SCIPlpiLoadColLP (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiLoadColLP") );
#endif

   if (nrows > 0)
   {
      SCIP_CALL( ensureBkcMem(lpi, nrows) );
      generateMskBoundkeys(nrows, lhs, rhs, lpi->bkc);
   }

   if (ncols > 0)
   {
      SCIP_CALL( ensureBkxMem(lpi, ncols) );
      generateMskBoundkeys(ncols, lb, ub, lpi->bkx);

      SCIP_CALL( ensureAptreMem(lpi, ncols) );
      SCIP_CALL( getEndptrs(ncols, beg, nnonz, lpi->aptre) );
   }

   MOSEK_CALL( MSK_inputdata(lpi->task, nrows, ncols, nrows, ncols, obj, 0.0, beg, lpi->aptre, ind, val,
         lpi->bkc, lhs, rhs, lpi->bkx, lb, ub) );

   MOSEK_CALL( MSK_putobjsense(lpi->task, SENSE2MOSEK(objsen)) );

   if( colnames != NULL )
   {
      int c;

      for( c = 0; c < ncols; c++ )
      {
         MOSEK_CALL( MSK_putvarname(lpi->task, c, colnames[c]) );
      }
   }

   if( rownames != NULL )
   {
      int r;

      for( r = 0; r < nrows; r++ )
      {
         MOSEK_CALL( MSK_putconname(lpi->task, r, rownames[r]) );
      }
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiLoadColLP") );
#endif

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
#if MSK_VERSION_MAJOR < 7
   const int* aptrb;
#endif

   int oldcols;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);
   assert(nnonz >= 0);
   assert(ncols >= 0);

   SCIPdebugMessage("Calling SCIPlpiAddCols (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiAddCols") );
#endif

   if (ncols == 0)
      return SCIP_OKAY;

   SCIP_CALL( ensureBkxMem(lpi, ncols) );
   generateMskBoundkeys(ncols, lb, ub, lpi->bkx);

   MOSEK_CALL( MSK_getnumvar(lpi->task, &oldcols) );

   MOSEK_CALL( MSK_appendvars(lpi->task, ncols) );
   MOSEK_CALL( MSK_putcslice(lpi->task, oldcols, oldcols+ncols, obj) );
   MOSEK_CALL( MSK_putvarboundslice(lpi->task, oldcols, oldcols+ncols, lpi->bkx, lb, ub) );

   if( nnonz > 0 )
   {
#ifndef NDEBUG
      /* perform check that no new rows are added - this is forbidden */
      int nrows;
      int j;

      MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );
      for (j = 0; j < nnonz; ++j)
      {
         assert( 0 <= ind[j] && ind[j] < nrows );
         assert( val[j] != 0.0 );
      }
#endif

      SCIP_CALL( ensureAptreMem(lpi, ncols) );
      SCIP_CALL( getEndptrs(ncols, beg, nnonz, lpi->aptre) );
      MOSEK_CALL( MSK_putacolslice(lpi->task, oldcols, oldcols+ncols, beg, lpi->aptre, ind, val) );
   }

   if( colnames != NULL )
   {
      int c;

      for( c = 0; c < ncols; c++ )
      {
         MOSEK_CALL( MSK_putvarname(lpi->task, c, colnames[c]) );
      }
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiAddCols") );
#endif

   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   int* sub;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiDelCols (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelCols") );
#endif

   SCIP_CALL( getIndicesRange(firstcol, lastcol, &sub) );

   MOSEK_CALL( MSK_removevars(lpi->task, lastcol-firstcol+1, sub) );

   BMSfreeMemoryArray(&sub);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelCols") );
#endif

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
   int* sub = NULL;
   int count;
   int ncols;
   int col;
   int i;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(dstat != NULL);

   SCIPdebugMessage("Calling SCIPlpiDelColset (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelColset") );
#endif

   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );

   SCIP_CALL( getIndicesFromDense(dstat, ncols, &count, &sub) );

   col = 0;
   for( i = 0; i < ncols; i++)
   {
      if (dstat[i] == 1)
      {
         dstat[i] = -1;
      }
      else
      {
         dstat[i] = col;
         col++;
      }
   }

   if (count > 0)
   {
      SCIPdebugMessage("Deleting %d vars %d, ...\n", count, sub[0]);
      MOSEK_CALL( MSK_removevars(lpi->task, count, sub) );
      BMSfreeMemoryArray(&sub);
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelColset") );
#endif

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
   int oldrows;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);

   SCIPdebugMessage("Calling SCIPlpiAddRows (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiAddRows") );
#endif

   if (nrows == 0)
      return SCIP_OKAY;

   SCIP_CALL( ensureBkcMem(lpi, nrows) );

   generateMskBoundkeys(nrows, lhs, rhs, lpi->bkc);

   MOSEK_CALL( MSK_getnumcon(lpi->task, &oldrows) );

   MOSEK_CALL( MSK_appendcons(lpi->task, nrows) );
   MOSEK_CALL( MSK_putconboundslice(lpi->task, oldrows, oldrows+nrows, lpi->bkc, lhs, rhs) );

   if( nnonz > 0 )
   {
#ifndef NDEBUG
      /* perform check that no new cols are added - this is forbidden */
      int ncols;
      int j;

      MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
      for (j = 0; j < nnonz; ++j)
      {
         assert( val[j] != 0.0 );
         assert( 0 <= ind[j] && ind[j] < ncols );
      }
#endif

      SCIP_CALL( ensureAptreMem(lpi, nrows) );
      SCIP_CALL( getEndptrs(nrows, beg, nnonz, lpi->aptre) );
      MOSEK_CALL( MSK_putarowslice(lpi->task, oldrows, oldrows+nrows, beg, lpi->aptre, ind, val) );
   }

   if( rownames != NULL )
   {
      int r;

      for( r = 0; r < nrows; r++ )
      {
         MOSEK_CALL( MSK_putconname(lpi->task, r, rownames[r]) );
      }
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiAddRows") );
#endif

   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   int* sub;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiDelRows (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelRows") );
#endif

   SCIP_CALL( getIndicesRange(firstrow, lastrow, &sub) );

   SCIPdebugMessage("Deleting cons %d to %d\n", firstrow, lastrow);

   MOSEK_CALL( MSK_removecons(lpi->task, lastrow-firstrow+1, sub) );

   BMSfreeMemoryArray(&sub);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelRows") );
#endif

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
   int* sub;
   int count;
   int nrows;
   int row;
   int i;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiDelRowset (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelRowset") );
#endif

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );

   sub = NULL;
   SCIP_CALL( getIndicesFromDense(dstat, nrows, &count, &sub) );

   row = 0;
   for( i = 0; i < nrows; i++ )
   {
      if (dstat[i] == 1)
      {
         dstat[i] = -1;
      }
      else
      {
         dstat[i] = row;
         row++;
      }
   }

   if (count > 0)
   {
      SCIPdebugMessage("Deleting %d cons %d, ...\n", count, sub[0]);
      MOSEK_CALL( MSK_removecons(lpi->task, count, sub) );
      BMSfreeMemoryArray(&sub);
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelRowset end") );
#endif

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int nrows;
   int ncols;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiClear (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );

   SCIP_CALL( SCIPlpiDelRows(lpi, 0, nrows - 1) );
   SCIP_CALL( SCIPlpiDelCols(lpi, 0, ncols - 1) );

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
   int i;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(ncols == 0 || (ind != NULL && lb != NULL && ub != NULL));

   SCIPdebugMessage("Calling SCIPlpiChgBounds (%d)\n", lpi->lpid);
   if( ncols <= 0 )
      return SCIP_OKAY;

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgBounds") );
#endif

   /* @todo This test could be integrated into generateMskBoundkeys, but then this function needs to be able to return an
    * error, which requires some rewriting. */
   for (i = 0; i < ncols; ++i)
   {
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

   SCIP_CALL( ensureBkxMem(lpi, ncols) );

   generateMskBoundkeys(ncols, lb, ub, lpi->bkx);
#if MSK_VERSION_MAJOR < 9
   MOSEK_CALL( MSK_putboundlist(lpi->task, MSK_ACC_VAR, ncols, ind, lpi->bkx, lb, ub) );
#else
   MOSEK_CALL( MSK_putvarboundlist(lpi->task, ncols, ind, lpi->bkx, lb, ub) );
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgBounds") );
#endif

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(ind != NULL);

   if( nrows <= 0 )
      return SCIP_OKAY;

   invalidateSolution(lpi);

   SCIPdebugMessage("Calling SCIPlpiChgSides (%d)\n", lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgSides") );
#endif

   SCIP_CALL( ensureBkcMem(lpi, nrows) );

   generateMskBoundkeys(nrows, lhs, rhs, lpi->bkc);
#if MSK_VERSION_MAJOR < 9
   MOSEK_CALL( MSK_putboundlist(lpi->task, MSK_ACC_CON, nrows, ind, lpi->bkc, lhs, rhs) );
#else
   MOSEK_CALL( MSK_putconboundlist(lpi->task, nrows, ind, lpi->bkc, lhs, rhs) );
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgSides") );
#endif

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiChgCoef (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgCoef") );
#endif

   MOSEK_CALL( MSK_putaij(lpi->task, row, col, newval) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgCoef") );
#endif

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiChgObjsen (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

   MOSEK_CALL( MSK_putobjsense(lpi->task, SENSE2MOSEK(objsen)) );

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(ind != NULL);
   assert(obj != NULL);

   SCIPdebugMessage("Calling SCIPlpiChgObj (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgObj") );
#endif

   MOSEK_CALL( MSK_putclist(lpi->task, ncols, ind, obj) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi,"SCIPlpiChgObj") );
#endif

   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   int nnonz;
   int* sub;
   double* val;
   MSKboundkeye bkc;
   double blc;
   double buc;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiScaleRow (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiScaleRow") );
#endif

   assert(scaleval != 0);

   MOSEK_CALL( MSK_getarownumnz(lpi->task, row, &nnonz) );

   if (nnonz != 0)
   {
      SCIP_ALLOC( BMSallocMemoryArray(&sub, nnonz) );
      SCIP_ALLOC( BMSallocMemoryArray(&val, nnonz) );

      MOSEK_CALL( MSK_getarow(lpi->task, row, &nnonz, sub, val) );
      scale_vec(nnonz, val, scaleval);
      MOSEK_CALL( MSK_putarow(lpi->task, row, nnonz, sub, val) );

      BMSfreeMemoryArray(&val);
      BMSfreeMemoryArray(&sub);
   }

#if MSK_VERSION_MAJOR < 9
   MOSEK_CALL( MSK_getbound(lpi->task, MSK_ACC_CON, row, &bkc, &blc, &buc) );
   scale_bound(&bkc, &blc, &buc, scaleval);
   MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_CON, row, bkc, blc, buc) );
#else
   MOSEK_CALL( MSK_getconbound(lpi->task, row, &bkc, &blc, &buc) );
   scale_bound(&bkc, &blc, &buc, scaleval);
   MOSEK_CALL( MSK_putconbound(lpi->task, row, bkc, blc, buc) );
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiScaleRow") );
#endif

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
   int nnonz;
   int *sub    = NULL;
   double *val = NULL;
   MSKboundkeye bkx;
   double blx, bux, c;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiScaleCol (%d)\n", lpi->lpid);

   invalidateSolution(lpi);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiScaleCol") );
#endif

   assert(scaleval != 0);
   MOSEK_CALL( MSK_getacolnumnz(lpi->task, col, &nnonz) );

   if (nnonz != 0)
   {
      SCIP_ALLOC( BMSallocMemoryArray(&sub, nnonz) );
      SCIP_ALLOC( BMSallocMemoryArray(&val, nnonz) );

      MOSEK_CALL( MSK_getacol(lpi->task, col, &nnonz, sub, val) );
      scale_vec(nnonz, val, scaleval);
      MOSEK_CALL( MSK_putacol(lpi->task, col, nnonz, sub, val) );

      BMSfreeMemoryArray(&val);
      BMSfreeMemoryArray(&sub);
   }

#if MSK_VERSION_MAJOR < 9
   MOSEK_CALL( MSK_getbound(lpi->task, MSK_ACC_VAR, col, &bkx, &blx, &bux) );
   scale_bound(&bkx, &blx, &bux, 1.0/scaleval);
   MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_VAR, col, bkx, blx, bux) );
#else
   MOSEK_CALL( MSK_getvarbound(lpi->task, col, &bkx, &blx, &bux) );
   scale_bound(&bkx, &blx, &bux, 1.0/scaleval);
   MOSEK_CALL( MSK_putvarbound(lpi->task, col, bkx, blx, bux) );
#endif

   MOSEK_CALL( MSK_getcslice(lpi->task, col, col+1, &c) );
   MOSEK_CALL( MSK_putcj(lpi->task, col, c*scaleval) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiScaleCol") );
#endif

   return SCIP_OKAY;
}


/*
 * Data Accessing Methods
 */


/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiGetNRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(nrows != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetNRows (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getnumcon(lpi->task, nrows) );

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(ncols != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetNCols (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getnumvar(lpi->task, ncols) );

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(nnonz != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetNNonz (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getnumanz(lpi->task, nnonz) );

   return SCIP_OKAY;
}

/** get a slice of a row or column */
static
SCIP_RETCODE getASlice(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             iscon,              /**< whether we are requesting a slice of a constraint or column */
   int                   first,              /**< first index */
   int                   last,               /**< last index */
   int*                  nnonz,              /**< pointer to store the number of nonzeros */
   int*                  beg,                /**< array for begins of indices/values */
   int*                  ind,                /**< array of row/column indices */
   double*               val                 /**< array of values */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(first <= last);

   SCIPdebugMessage("Calling SCIPlpiGetNNonz (%d)\n", lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "getASlice") );
#endif

   if( nnonz != 0 )
   {
      int surplus;

      assert(beg != NULL);
      assert(ind != NULL);
      assert(val != NULL);

      SCIP_CALL( ensureAptreMem(lpi, last - first + 1) );

#if MSK_VERSION_MAJOR < 9
      MOSEK_CALL( MSK_getaslicenumnz(lpi->task, iscon ? MSK_ACC_CON : MSK_ACC_VAR, first, last+1, nnonz) );
      surplus = *nnonz;
      MOSEK_CALL( MSK_getaslice(lpi->task, iscon ? MSK_ACC_CON : MSK_ACC_VAR, first, last+1, *nnonz, &surplus, beg, lpi->aptre, ind, val) );
#else
      if( iscon )
      {
         MOSEK_CALL( MSK_getarowslicenumnz(lpi->task, first, last+1, nnonz) );
         surplus = *nnonz;
#if MSK_VERSION_MAJOR == 9
         MOSEK_CALL( MSK_getarowslice(lpi->task, first, last+1, *nnonz, &surplus, beg, lpi->aptre, ind, val) );
#else
         MOSEK_CALL( MSK_getarowslice(lpi->task, first, last+1, *nnonz, beg, lpi->aptre, ind, val) );
         (void)surplus;
#endif
      }
      else
      {
         MOSEK_CALL( MSK_getacolslicenumnz(lpi->task, first, last+1, nnonz) );
         surplus = *nnonz;
#if MSK_VERSION_MAJOR == 9
         MOSEK_CALL( MSK_getacolslice(lpi->task, first, last+1, *nnonz, &surplus, beg, lpi->aptre, ind, val) );
#else
         MOSEK_CALL( MSK_getacolslice(lpi->task, first, last+1, *nnonz, beg, lpi->aptre, ind, val) );
         (void)surplus;
#endif
      }
#endif

      assert(surplus == 0);
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "getASlice") );
#endif

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert((lb != NULL && ub != NULL) || (lb == NULL && ub == NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));

#ifndef NDEBUG
   {
      int ncols;
      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
      assert(0 <= firstcol && firstcol <= lastcol && lastcol < ncols);
   }
#endif

   SCIPdebugMessage("Calling SCIPlpiGetCols (%d)\n", lpi->lpid);

   SCIP_CALL( SCIPlpiGetBounds(lpi, firstcol, lastcol, lb, ub) );
   SCIP_CALL( getASlice(lpi, FALSE, firstcol, lastcol, nnonz, beg, ind, val) );

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert((lhs != NULL && rhs != NULL) || (lhs == NULL && rhs == NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));

#ifndef NDEBUG
   {
      int nrows;
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
      assert(0 <= firstrow && firstrow <= lastrow && lastrow < nrows);
   }
#endif

   SCIPdebugMessage("Calling SCIPlpiGetRows (%d)\n", lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetRows") );
#endif

   SCIP_CALL( SCIPlpiGetSides(lpi, firstrow, lastrow, lhs, rhs) );
   SCIP_CALL( getASlice(lpi, TRUE, firstrow, lastrow, nnonz, beg, ind, val) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetRows") );
#endif

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(0 <= firstcol && firstcol <= lastcol);
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
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(0 <= firstrow && firstrow <= lastrow);
   assert(rownames != NULL || namestoragesize == 0);
   assert(namestorage != NULL || namestoragesize == 0);
   assert(namestoragesize >= 0);
   assert(storageleft != NULL);

   SCIPerrorMessage("SCIPlpiGetRowNames() has not been implemented yet.\n");

   return SCIP_LPERROR;
}

/** gets the objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   MSKobjsensee mskobjsen;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(objsen != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetObjsen (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getobjsense(lpi->task, &mskobjsen) );
   *objsen = (mskobjsen == MSK_OBJECTIVE_SENSE_MINIMIZE ? SCIP_OBJSEN_MINIMIZE : SCIP_OBJSEN_MAXIMIZE);

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(firstcol <= lastcol);
   assert(vals != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetObj (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getcslice(lpi->task, firstcol, lastcol+1, vals) );

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(firstcol <= lastcol);

   SCIPdebugMessage("Calling SCIPlpiGetBounds (%d)\n", lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetBounds") );
#endif

#if MSK_VERSION_MAJOR < 9
   MOSEK_CALL( MSK_getboundslice(lpi->task, MSK_ACC_VAR, firstcol, lastcol+1, NULL, lbs, ubs) );
#else
   MOSEK_CALL( MSK_getvarboundslice(lpi->task, firstcol, lastcol+1, NULL, lbs, ubs) );
#endif

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(firstrow <= lastrow);

   SCIPdebugMessage("Calling SCIPlpiGetSides (%d)\n", lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetSides") );
#endif

#if MSK_VERSION_MAJOR < 9
   MOSEK_CALL( MSK_getboundslice(lpi->task, MSK_ACC_CON, firstrow, lastrow+1, NULL, lhss, rhss) );
#else
   MOSEK_CALL( MSK_getconboundslice(lpi->task, firstrow, lastrow+1, NULL, lhss, rhss) );
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetSides") );
#endif

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(val != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetCoef (%d)\n", lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetCoef") );
#endif

   MOSEK_CALL( MSK_getaij(lpi->task, row, col, val) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetCoef") );
#endif

   return SCIP_OKAY;
}

/*
 * Solving Methods
 */


/** gets the internal solution status of the solver */
static
SCIP_RETCODE getSolutionStatus(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   MSKprostae*           prosta,             /**< pointer to store the problem status */
   MSKsolstae*           solsta              /**< pointer to store the solution status */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   MOSEK_CALL( MSK_getsolutionstatus(lpi->task, lpi->lastsolvetype, prosta, solsta) );

   return SCIP_OKAY;
}

/** helper method to filter out numerical problems */
static
MSKrescodee filterTRMrescode(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   MSKrescodee*          termcode,           /**< pointer to store output termination code */
   MSKrescodee           res                 /**< input result of call to Mosek function */
   )
{  /*lint --e{715}*/
   assert( termcode != NULL );

#if ASSERT_ON_NUMERICAL_TROUBLES > 0
   if ( res == MSK_RES_TRM_MAX_NUM_SETBACKS || res == MSK_RES_TRM_NUMERICAL_PROBLEM )
   {
      SCIPmessagePrintWarning(messagehdlr, "Return code %d.\n", res);
      assert(0);
      *termcode = res;
      return MSK_RES_OK;
   }
#else
   SCIP_UNUSED(messagehdlr);
#endif

   if (  res == MSK_RES_TRM_MAX_ITERATIONS || res == MSK_RES_TRM_MAX_TIME
      || res == MSK_RES_TRM_OBJECTIVE_RANGE || res == MSK_RES_TRM_STALL )
   {
      *termcode = res;
      res = MSK_RES_OK;
   }
   else
      *termcode = MSK_RES_OK;

   return res;
}

/** solve problem with the simplex algorithm */
static
SCIP_RETCODE SolveWSimplex(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int itercount_primal;
   int itercount_dual;
   int gotbasicsol;
   int presolve;
   int maxiter;
   MSKprostae prosta;
   MSKsolstae solsta;
   double pobj, dobj;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   invalidateSolution(lpi);
   lpi->lastsolvetype = MSK_SOL_BAS;

   /* store original settings */
   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_PRESOLVE_USE, &presolve) );
   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, &maxiter) );

   /* set some paramters */
#if DEBUG_EASY_REPRODUCE
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_AUTO_SORT_A_BEFORE_OPT, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_OFF) );
#else
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );
#endif

   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_AUTO_UPDATE_SOL_INFO, MSK_OFF) );

#if FORCE_MOSEK_LOG
   if( lpi->optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG_SIM, 4) );
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG_SIM_FREQ, 1) );
   }
   else
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG_SIM, 4) );
   }
#endif

   MOSEK_CALL( MSK_solutiondef(lpi->task, MSK_SOL_BAS, &gotbasicsol) );

   if( gotbasicsol )
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_OFF) );
   }
   else
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_ON) );
   }

#if ALWAYS_SOLVE_PRIMAL_FORM > 0
   /* always solve the primal formulation */
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_SOLVE_FORM, MSK_SOLVE_PRIMAL) );
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SolveWSimplex") );
#endif

   if( gotbasicsol && maxiter < 20000 )
   {
      /* Since max iter often is set, we switch off restricted pricing */
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_DUAL_RESTRICT_SELECTION, 0) );
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_PRIMAL_RESTRICT_SELECTION, 0) );
   }
   else
   {
      /* otherwise use default value */
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_DUAL_RESTRICT_SELECTION, lpi->restrictselectdef) );
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_PRIMAL_RESTRICT_SELECTION, lpi->restrictselectdef) );
   }

#if FORCE_NO_MAXITER > 0
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, 2000000000) );
#endif


#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "Begin optimize with simplex") );
#endif

#if FORCE_MOSEK_SUMMARY > 1
   if( lpi->optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_solutionsummary(lpi->task, MSK_STREAM_LOG) );
   }
#endif

   /* perform actual optimization */
   MOSEK_CALL( filterTRMrescode(lpi->messagehdlr, &lpi->termcode, MSK_optimize(lpi->task)) );

#if MSK_VERSION_MAJOR < 10
   /* resolve with aggressive scaling if the maximal number of setbacks has been reached */
   if( lpi->termcode == MSK_RES_TRM_MAX_NUM_SETBACKS )
   {
      int scaling;

      MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_SCALING, &scaling) );
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_SCALING, MSK_SCALING_AGGRESSIVE) );
      MOSEK_CALL( filterTRMrescode(lpi->messagehdlr, &lpi->termcode, MSK_optimize(lpi->task)) );
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_SCALING, scaling) );
   }
#endif

#if FORCE_MOSEK_SUMMARY
   if( lpi->optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_solutionsummary(lpi->task, MSK_STREAM_LOG) );
   }
#else
   if( lpi->lpinfo )
   {
      MOSEK_CALL( MSK_solutionsummary(lpi->task, MSK_STREAM_LOG) );
   }
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "End optimize with simplex") );
#endif

   /* set parameters to their original values */
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_PRESOLVE_USE, presolve) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, maxiter) );

   /* obtain iteration count */
   MOSEK_CALL( MSK_getintinf(lpi->task, MSK_IINF_SIM_PRIMAL_ITER, &itercount_primal) );
   MOSEK_CALL( MSK_getintinf(lpi->task, MSK_IINF_SIM_DUAL_ITER, &itercount_dual) );

   lpi->itercount = itercount_primal + itercount_dual;

   /* get solution information */
   MOSEK_CALL( MSK_getprimalobj(lpi->task, MSK_SOL_BAS, &pobj) );
   MOSEK_CALL( MSK_getdualobj(lpi->task, MSK_SOL_BAS, &dobj) );

   MOSEK_CALL( MSK_getsolutionstatus(lpi->task, MSK_SOL_BAS, &prosta, &solsta) );

   SCIPdebugMessage("maxiter = %d, termcode = %d, prosta = %d, solsta = %d, objval = %g : %g, iter = %d+%d\n",
      maxiter, lpi->termcode, prosta, solsta, pobj, dobj, itercount_primal, itercount_dual);

   switch (solsta)
   {
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_FEAS:
   case MSK_SOL_STA_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      if (lpi->termcode == MSK_RES_OK)
         lpi->solved = TRUE;
      break;

   case MSK_SOL_STA_UNKNOWN:
      /* Mosek seems to have status unknown on the following termination codes */
      assert( lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS || lpi->termcode == MSK_RES_TRM_MAX_TIME ||
         lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE || lpi->termcode == MSK_RES_TRM_STALL ||
         lpi->termcode == MSK_RES_OK );

      if ( lpi->termcode != MSK_RES_TRM_MAX_ITERATIONS && lpi->termcode != MSK_RES_TRM_MAX_TIME &&
           lpi->termcode != MSK_RES_TRM_OBJECTIVE_RANGE )
      {
         SCIPmessagePrintWarning(lpi->messagehdlr, "Numerical problem: simplex[%d] returned solsta = %d.\n", lpi->optimizecount, solsta);
         lpi->termcode = MSK_RES_TRM_NUMERICAL_PROBLEM;
#if ASSERT_ON_WARNING
         assert(0);
#endif
      }
      break;

#if MSK_VERSION_MAJOR < 9
   case MSK_SOL_STA_NEAR_OPTIMAL:
   case MSK_SOL_STA_NEAR_PRIM_FEAS:
   case MSK_SOL_STA_NEAR_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
   case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:

      assert(lpi->termcode == MSK_RES_OK);

      SCIPmessagePrintWarning(lpi->messagehdlr, "Simplex[%d] returned solsta = %d (numerical problem).\n", lpi->optimizecount, solsta);
      lpi->termcode = MSK_RES_TRM_NUMERICAL_PROBLEM;
#if ASSERT_ON_WARNING
      assert(0);
#endif
      break;
#endif

   case MSK_SOL_STA_INTEGER_OPTIMAL:
#if MSK_VERSION_MAJOR < 9
   case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
#endif
   default:
#if SHOW_ERRORS
      SCIPerrorMessage("Simplex[%d] returned solsta = %d\n", lpi->optimizecount, solsta);
#endif

#if ASSERT_ON_WARNING
      assert(0);
#endif

      return SCIP_LPERROR;
   }  /*lint !e788*/

   switch (prosta)
   {
   /* already handled above */
   case MSK_PRO_STA_PRIM_AND_DUAL_FEAS:
   case MSK_PRO_STA_PRIM_FEAS:
   case MSK_PRO_STA_DUAL_FEAS:
   case MSK_PRO_STA_PRIM_AND_DUAL_INFEAS:
   case MSK_PRO_STA_PRIM_INFEAS:
   case MSK_PRO_STA_DUAL_INFEAS:
   case MSK_PRO_STA_UNKNOWN:
      break;

#if MSK_VERSION_MAJOR < 9
   case MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS:
   case MSK_PRO_STA_NEAR_PRIM_FEAS:
   case MSK_PRO_STA_NEAR_DUAL_FEAS:
#endif
   case MSK_PRO_STA_ILL_POSED:
   case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
      assert(lpi->termcode == MSK_RES_OK);

      SCIPmessagePrintWarning(lpi->messagehdlr, "Simplex[%d] returned prosta = %d\n", lpi->optimizecount, prosta);
      lpi->termcode = MSK_RES_TRM_NUMERICAL_PROBLEM;
      invalidateSolution(lpi);
#if ASSERT_ON_WARNING
      assert(0);
#endif
      break;

   default:
#if SHOW_ERRORS
      SCIPerrorMessage("Simplex[%d] returned prosta = %d\n", lpi->optimizecount, prosta);
#endif

#if ASSERT_ON_WARNING
      assert(0);
#endif

      return SCIP_LPERROR;
   }  /*lint !e788*/

   /* todo: replace numbers by constants, e.g., tolerances */
#if SHOW_RELATIVE_OPTIMAL_GAP
   if ( solsta == MSK_SOL_STA_OPTIMAL && fabs(pobj) + fabs(dobj) > 1.0e-6 && fabs(pobj-dobj) > 0.0001*(fabs(pobj) + fabs(dobj)))
   {
      SCIPerrorMessage("Simplex[%d] returned optimal solution with different objvals %g != %g reldiff %.2g%%\n",
         lpi->optimizecount, pobj, dobj, 100.0 * fabs(pobj-dobj)/ MAX(fabs(pobj), fabs(dobj))); /*lint !e666*/
   }
#endif

   /* The optimizer terminated with an objective value outside the objective range. */
   if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
   {
      if (solsta != MSK_SOL_STA_DUAL_FEAS && solsta != MSK_SOL_STA_OPTIMAL && solsta != MSK_SOL_STA_PRIM_AND_DUAL_FEAS)
      {
         SCIPerrorMessage("[%d] Terminated on objective range without dual feasible solsta.\n", lpi->optimizecount);

         /* solve again with barrier */
         SCIP_CALL( SCIPlpiSolveBarrier(lpi, TRUE) );
      }
   }

   /* if the simplex took too many iterations, solve again with barrier */
   if (maxiter >= 2000000000)
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, maxiter) );

      if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      {
         SCIPmessagePrintWarning(lpi->messagehdlr, "Simplex[%d] failed to terminate in 10000 iterations, switching to interior point\n",
            lpi->optimizecount);

         SCIP_CALL( SCIPlpiSolveBarrier(lpi, TRUE) );
      }
   }

#if DEBUG_DO_INTPNT_FEAS_CHECK
   if (solsta == MSK_SOL_STA_PRIM_INFEAS_CER || solsta == MSK_SOL_STA_DUAL_INFEAS_CER)
   {
      SCIPdebugMessage("Checking infeasibility[%d]... ", lpi->optimizecount);

      SCIP_CALL( SCIPlpiSolveBarrier(lpi, true) );

      MOSEK_CALL( MSK_getsolutionstatus(lpi->task, MSK_SOL_BAS, &prosta, &solsta) );

      if (solsta == MSK_SOL_STA_PRIM_INFEAS_CER || solsta == MSK_SOL_STA_DUAL_INFEAS_CER)
      {
         SCIPdebugPrintf("ok\n");
      }
      else
      {
         SCIPdebugPrintf("wrong [%d] prosta = %d, solsta = %d\n", lpi->optimizecount, prosta, solsta);
      }
   }
#endif


#if DEBUG_PRINT_STAT > 0
   SCIPdebugMessage("Max iter stat    : Count %d branchup = %d branchlo = %d primal %d dual %d\n",
      lpi->optimizecount, numstrongbranchmaxiterup, numstrongbranchmaxiterdo, numprimalmaxiter, numdualmaxiter);
   SCIPdebugMessage("Objcut iter stat : Count %d branchup = %d branchlo = %d primal %d dual %d\n",
      lpi->optimizecount, numstrongbranchobjup, numstrongbranchobjdo, numprimalobj, numdualobj);
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SolveWSimplex") );
#endif

   return SCIP_OKAY;
}

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   lpi->optimizecount++;

   SCIPdebugMessage("Calling SCIPlpiSolvePrimal[%d] (%d) ", lpi->optimizecount, lpi->lpid);

   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );

   /* Set warmstarting information in MOSEK. We only have status keys (recalculate dual solution without dual superbasics) */
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART, lpi->fromscratch || lpi->clearstate ?
         MSK_SIM_HOTSTART_NONE : MSK_SIM_HOTSTART_STATUS_KEYS) );
   lpi->clearstate = FALSE;

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiSolvePrimal") );
#endif

   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_PRIMAL_SIMPLEX) );
   lpi->lastalgo = MSK_OPTIMIZER_PRIMAL_SIMPLEX;

#if WRITE_PRIMAL > 0
   if( lpi->optimizecount > WRITE_ABOVE )
   {
      char fname[40];
      snprintf(fname, 40, "primal_%d.lp", lpi->optimizecount);
      SCIPdebugMessage("\nWriting lp %s\n", fname);
      /*MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_GENERIC_NAMES, MSK_ON) );*/
      MSK_writedata(lpi->task, fname);
   }
#endif

   SCIP_CALL( SolveWSimplex(lpi) );

#ifdef SCIP_DISABLED_CODE
   /* the following code is unclear: Why should the resolve change anything ?????? */
   if ( lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE )
   {
      MSKsolstae solsta;

      MOSEK_CALL( MSK_getsolutionstatus(lpi->task, MSK_SOL_BAS, NULL, &solsta) );

      if( solsta != MSK_SOL_STA_PRIM_FEAS )
      {
         SCIP_CALL( SolveWSimplex(lpi) );
      }
   }
#endif

#if DEBUG_PRINT_STAT > 0
   if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
      ++numprimalobj;
#endif

#if DEBUG_PRINT_STAT > 0
   if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      ++numprimalmaxiter;
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiSolvePrimal") );
#endif

   return SCIP_OKAY;
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   lpi->optimizecount++;

   SCIPdebugMessage("Calling SCIPlpiSolveDual[%d] (%d)\n", lpi->optimizecount, lpi->lpid);

/* MSK_IPAR_SIM_INTEGER is removed in Mosek 8.1 */
#if (MSK_VERSION_MAJOR < 8) || (MSK_VERSION_MAJOR == 8 && MSK_VERSION_MINOR == 0)
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_INTEGER, MSK_ON) );
#endif
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );

   /* Set warmstarting information in MOSEK. We only have status keys (recalculate dual solution without dual superbasics) */
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART, (lpi->fromscratch || lpi->clearstate) ?
         MSK_SIM_HOTSTART_NONE : MSK_SIM_HOTSTART_STATUS_KEYS) );
   lpi->clearstate = FALSE;

   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_DUAL_SIMPLEX) );
   lpi->lastalgo = MSK_OPTIMIZER_DUAL_SIMPLEX;

#if WRITE_DUAL > 0
   if( lpi->optimizecount > WRITE_ABOVE )
   {
      char fname[40];
      snprintf(fname,40,"dual_%d.lp", lpi->optimizecount);
      SCIPdebugMessage("\nWriting lp %s\n", fname);
      MSK_writedata(lpi->task, fname);
   }
#endif

   SCIP_CALL( SolveWSimplex(lpi) );

#ifdef SCIP_DISABLED_CODE
   /* the following code is unclear: Why should the resolve change anything ?????? */
   if ( lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE )
   {
      MSKsolstae solsta;

      MOSEK_CALL( MSK_getsolutionstatus(lpi->task, MSK_SOL_BAS, NULL, &solsta) );

      if( solsta != MSK_SOL_STA_DUAL_FEAS )
      {
         SCIP_CALL( SolveWSimplex(lpi) );
      }
   }
#endif

#if DEBUG_PRINT_STAT > 0
   if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
      ++numdualobj;
#endif

#if DEBUG_PRINT_STAT > 0
   if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      ++numdualmaxiter;
#endif

   return SCIP_OKAY;
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{
   MSKprostae prosta;
   MSKsolstae solsta;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   lpi->optimizecount++;

   invalidateSolution(lpi);
   lpi->lastsolvetype = crossover ? MSK_SOL_BAS : MSK_SOL_ITR;

#if FORCE_MOSEK_LOG
   if( lpi->optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, 4) );
   }
   else
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, 0) );
   }
#endif

   SCIPdebugMessage("Calling SCIPlpiSolveBarrier[%d] (%d) ", lpi->optimizecount, lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiSolveBarrier") );
#endif

#ifdef SCIP_DISABLED_CODE
   /* The parameter exists in MOSEK, but as of version 8, it is not in use and the interior-point solver is never warmstarted */
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_INTPNT_HOTSTART, (lpi->fromscratch || lpi->clearstate) ?
         MSK_SIM_HOTSTART_NONE : MSK_INTPNT_HOTSTART_PRIMAL_DUAL) );
#endif
   lpi->clearstate = FALSE;

   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_INTPNT_BASIS, crossover ? MSK_BI_ALWAYS : MSK_BI_NEVER) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT) );
   lpi->lastalgo = MSK_OPTIMIZER_INTPNT;

#if MSK_VERSION_MAJOR >= 9
   MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_INTPNT_CO_TOL_NEAR_REL, NEAR_REL_TOLERANCE) );
#endif

#if WRITE_INTPNT > 0
   if( lpi->optimizecount > WRITE_ABOVE )
   {
      char fname[40];
      snprintf(fname,40,"intpnt_%d.lp", lpi->optimizecount);
      SCIPdebugMessage("\nWriting lp %s\n", fname);
      /*MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_GENERIC_NAMES, MSK_ON) );*/
      MSK_writedata(lpi->task, fname);
   }
#endif

   MOSEK_CALL( filterTRMrescode(lpi->messagehdlr, &lpi->termcode, MSK_optimize(lpi->task)) );

#if DEBUG_PRINT_STAT > 0
   if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      ++numdualmaxiter;
#endif

   MOSEK_CALL( MSK_getintinf(lpi->task, MSK_IINF_INTPNT_ITER, &lpi->itercount) );

   MOSEK_CALL( MSK_getsolutionstatus(lpi->task, lpi->lastsolvetype, &prosta, &solsta) );
   SCIPdebugMessage("termcode = %d, prosta = %d, solsta = %d, iter = %d\n",
      lpi->termcode, prosta, solsta, lpi->itercount);

   switch (solsta)
   {
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_FEAS:
   case MSK_SOL_STA_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      if (lpi->termcode == MSK_RES_OK)
         lpi->solved = TRUE;
      break;
   case MSK_SOL_STA_UNKNOWN:
#if MSK_VERSION_MAJOR < 9
   case MSK_SOL_STA_NEAR_OPTIMAL:
   case MSK_SOL_STA_NEAR_PRIM_FEAS:
   case MSK_SOL_STA_NEAR_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
   case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
#endif
      SCIPmessagePrintWarning(lpi->messagehdlr, "Barrier[%d] returned solsta = %d\n", lpi->optimizecount, solsta);

      if (lpi->termcode == MSK_RES_OK)
         lpi->termcode = MSK_RES_TRM_NUMERICAL_PROBLEM;

#if ASSERT_ON_WARNING
      assert(0);
#endif
      break;
   case MSK_SOL_STA_INTEGER_OPTIMAL:
#if MSK_VERSION_MAJOR < 9
   case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
#endif
   default:
#if SHOW_ERRORS
      SCIPerrorMessage("Barrier[%d] returned solsta = %d\n", lpi->optimizecount, solsta);
#endif

#if ASSERT_ON_WARNING
      assert(0);
#endif

      return SCIP_LPERROR;
   }  /*lint !e788*/

   switch (prosta)
   {
   case MSK_PRO_STA_PRIM_AND_DUAL_FEAS:
   case MSK_PRO_STA_PRIM_FEAS:
   case MSK_PRO_STA_DUAL_FEAS:
   case MSK_PRO_STA_PRIM_AND_DUAL_INFEAS:
   case MSK_PRO_STA_PRIM_INFEAS:
   case MSK_PRO_STA_DUAL_INFEAS:
      break;
   case MSK_PRO_STA_UNKNOWN:
#if MSK_VERSION_MAJOR < 9
   case MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS:
   case MSK_PRO_STA_NEAR_PRIM_FEAS:
   case MSK_PRO_STA_NEAR_DUAL_FEAS:
#endif
   case MSK_PRO_STA_ILL_POSED:
   case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
      SCIPmessagePrintWarning(lpi->messagehdlr, "Barrier[%d] returned prosta = %d\n", lpi->optimizecount, prosta);

      if (lpi->termcode == MSK_RES_OK)
         lpi->termcode = MSK_RES_TRM_NUMERICAL_PROBLEM;

      invalidateSolution(lpi);

#if ASSERT_ON_WARNING
      assert(0);
#endif
      break;
   default:
#if SHOW_ERRORS
      SCIPerrorMessage("Barrier[%d] returned prosta = %d\n", lpi->optimizecount, prosta);
#endif

#if ASSERT_ON_WARNING
      assert(0);
#endif

      return SCIP_LPERROR;
   }  /*lint !e788*/

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiSolveBarrier") );
#endif

   return SCIP_OKAY;
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   /* currently do nothing */
   return SCIP_OKAY;
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{ /* lint --e{715}*/
   assert(lpi != NULL);
   /* assert(MosekEnv != NULL);
   assert(lpi->task != NULL); */

   /* currently do nothing */
   return SCIP_OKAY;
}

/** performs strong branching iterations on all candidates
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
static
SCIP_RETCODE SCIPlpiStrongbranch(
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
   MSKobjsensee objsen;
   int olditerlim;
   int oldselection;
   int oldhotstart;

   double bound;
   int ncols;
   int nrows;
   MSKboundkeye bkx;
   double blx;
   double bux;
   double newub;
   double newlb;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiStrongbranch (%d)\n", lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiStrongbranch") );
#endif

   if (lpi->termcode != MSK_RES_OK)
   {
      SCIPmessagePrintWarning(lpi->messagehdlr, "SB Warning: Previous termcode is %d\n", lpi->termcode);
   }

   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );

   SCIP_CALL( getbase(lpi, ncols, nrows) );

   MOSEK_CALL( MSK_getobjsense(lpi->task, &objsen) );
   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, &olditerlim) );
   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_DUAL_SELECTION, &oldselection) );
   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_HOTSTART, &oldhotstart) );

   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, itlim) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_DUAL_SELECTION, STRONGBRANCH_PRICING) );

   if (objsen == MSK_OBJECTIVE_SENSE_MINIMIZE)
   {
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_UPPER_OBJ_CUT, &bound) );
   }
   else /* objsen == MSK_OBJECTIVE_SENSE_MAX */
   {
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_LOWER_OBJ_CUT, &bound) );
   }

#if MSK_VERSION_MAJOR < 9
   MOSEK_CALL( MSK_getbound(lpi->task, MSK_ACC_VAR, col, &bkx, &blx, &bux) );
#else
   MOSEK_CALL( MSK_getvarbound(lpi->task, col, &bkx, &blx, &bux) );
#endif

   *iter = 0;

   newub = EPSCEIL(psol-1.0, 1e-06);

   if (newub < blx - 0.5) /* infeasible */
   {
      *down = bound;
      *downvalid = TRUE;
   }
   else
   {
      MSKboundkeye newbk;

      if (IS_NEGINF(blx))
         newbk = MSK_BK_UP;
      else if (EPSEQ(blx, newub,1.0e-6))
      {
         newbk = MSK_BK_FX;
         newub = blx;
      }
      else
         newbk = MSK_BK_RA;

#if MSK_VERSION_MAJOR < 9
      MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_VAR, col, newbk, blx, newub) );
#else
      MOSEK_CALL( MSK_putvarbound(lpi->task, col, newbk, blx, newub) );
#endif

      SCIP_CALL( SCIPlpiSolveDual(lpi) );

      *iter += lpi->itercount;

      if (SCIPlpiIsStable(lpi))
         *downvalid = TRUE;
      else
         *downvalid = FALSE;

      if (SCIPlpiExistsPrimalRay(lpi))
      {
         SCIPmessagePrintWarning(lpi->messagehdlr, "SB ERROR: Lp [%d] is dual infeasible\n", lpi->optimizecount);

         *down = -1e20;
         *downvalid = FALSE;
      }
      else if (SCIPlpiExistsDualRay(lpi))
      {
         *down = bound;
      }
      else
      {
         SCIP_Bool pfeas;
         SCIP_Bool dfeas;

         SCIP_CALL( SCIPlpiGetSolFeasibility(lpi, &pfeas, &dfeas) );

         if (!dfeas)
         {
            SCIPmessagePrintWarning(lpi->messagehdlr, "SB ERROR: Lp [%d] is not dual feasible\n", lpi->optimizecount);

            *down = -1e20;
            *downvalid = FALSE;
         }
         else
         {
            MOSEK_CALL( MSK_getdualobj(lpi->task, lpi->lastsolvetype, down) );
         }
      }

#if DEBUG_PRINT_STAT > 0
      if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
         ++numstrongbranchobjup;

      if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
         ++numstrongbranchmaxiterup;
#endif
   }

   /* Reset basis solution before doing the up branch */
#if MSK_VERSION_MAJOR < 9
   MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_VAR, col, bkx, blx, bux) );
#else
   MOSEK_CALL( MSK_putvarbound(lpi->task, col, bkx, blx, bux) );
#endif
   SCIP_CALL( setbase(lpi) );

   newlb = EPSFLOOR(psol+1.0, 1e-06);
   if (newlb > bux + 0.5) /* infeasible */
   {
      *up = bound;
      *upvalid = TRUE;
   }
   else
   {
      MSKboundkeye newbk;

      if (IS_POSINF(bux))
         newbk = MSK_BK_LO;
      else if (EPSEQ(bux, newlb,1.0e-6))
      {
         newbk = MSK_BK_FX;
         newlb = bux;
      }
      else
         newbk = MSK_BK_RA;

#if MSK_VERSION_MAJOR < 9
      MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_VAR, col, newbk, newlb, bux) );
#else
      MOSEK_CALL( MSK_putvarbound(lpi->task, col, newbk, newlb, bux) );
#endif
      SCIP_CALL( SCIPlpiSolveDual(lpi) );

      *iter += lpi->itercount;

      if (SCIPlpiIsStable(lpi))
         *upvalid = TRUE;
      else
         *upvalid = FALSE;

      if (SCIPlpiExistsPrimalRay(lpi))
      {
         *up = -1e20;
         *upvalid = FALSE;
      }
      else if (SCIPlpiExistsDualRay(lpi))
      {
         *up = bound;
      }
      else
      {
         SCIP_Bool pfeas;
         SCIP_Bool dfeas;

         SCIP_CALL( SCIPlpiGetSolFeasibility(lpi, &pfeas, &dfeas) );

         if (!dfeas)
         {
            SCIPmessagePrintWarning(lpi->messagehdlr, "SB ERROR: Lp [%d] is not dual feasible\n", lpi->optimizecount);

            *up = -1e20;
            *upvalid = FALSE;
         }
         else
         {
            MOSEK_CALL( MSK_getdualobj(lpi->task, lpi->lastsolvetype, up) );
         }
      }

#if DEBUG_PRINT_STAT > 0
      if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
         ++numstrongbranchobjdo;

      if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
         ++numstrongbranchmaxiterdo;
#endif
   }

#if MSK_VERSION_MAJOR < 9
   MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_VAR, col, bkx, blx, bux) );
#else
   MOSEK_CALL( MSK_putvarbound(lpi->task, col, bkx, blx, bux) );
#endif
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, olditerlim) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_DUAL_SELECTION, oldselection) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART, oldhotstart) );

   SCIP_CALL( setbase(lpi) );

   invalidateSolution(lpi);

   lpi->termcode = MSK_RES_OK;
   lpi->itercount = 0;

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiStrongbranch") );
#endif

   SCIPdebugMessage("End SCIPlpiStrongbranch (%d)\n", lpi->lpid);

   return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
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
   SCIP_CALL( SCIPlpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given @b fractional candidates
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
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

   if ( iter != NULL )
      *iter = 0;

   for (j = 0; j < ncols; ++j)
   {
      /* pass call on to lpiStrongbranch() */
      SCIP_CALL( SCIPlpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter) );
   }
   return SCIP_OKAY;
}

/** performs strong branching iterations on one candidate with @b integral value
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
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
   SCIP_CALL( SCIPlpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given candidates with @b integral values
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
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

   if ( iter != NULL )
      *iter = 0;

   for (j = 0; j < ncols; ++j)
   {
      /* pass call on to lpiStrongbranch() */
      SCIP_CALL( SCIPlpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter) );
   }
   return SCIP_OKAY;
}


/*
 * Solution Information Methods
 */


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
   MSKprostae prosta;

   assert( lpi != NULL );
   assert( lpi->mosekenv != NULL );
   assert( lpi->task != NULL );
   assert( primalfeasible != NULL );
   assert( dualfeasible != NULL );

   SCIPdebugMessage("Calling SCIPlpiGetSolFeasibility (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getsolutionstatus(lpi->task, lpi->lastsolvetype, &prosta, NULL) );

   switch (prosta)
   {
   case MSK_PRO_STA_PRIM_AND_DUAL_FEAS:
      *primalfeasible = TRUE;
      *dualfeasible = TRUE;
      break;
   case MSK_PRO_STA_PRIM_FEAS:
      *primalfeasible = TRUE;
      *dualfeasible = FALSE;
      break;
   case MSK_PRO_STA_DUAL_FEAS:
      *primalfeasible = FALSE;
      *dualfeasible = TRUE;
      break;
   case MSK_PRO_STA_DUAL_INFEAS:
      /* assume that we have a primal solution if we used the primal simplex */
      *primalfeasible = (lpi->lastalgo == MSK_OPTIMIZER_PRIMAL_SIMPLEX);
      *dualfeasible = FALSE;
      break;
   case MSK_PRO_STA_UNKNOWN:
   case MSK_PRO_STA_PRIM_INFEAS:
   case MSK_PRO_STA_PRIM_AND_DUAL_INFEAS:
   case MSK_PRO_STA_ILL_POSED:
#if MSK_VERSION_MAJOR < 9
   case MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS:
   case MSK_PRO_STA_NEAR_PRIM_FEAS:
   case MSK_PRO_STA_NEAR_DUAL_FEAS:
#endif
   case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
      *primalfeasible = FALSE;
      *dualfeasible = FALSE;
      break;
   default:
      return SCIP_LPERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKprostae prosta;
   MSKsolstae solsta;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiExistsPrimalRay (%d)\n", lpi->lpid);

   SCIP_ABORT_FALSE( getSolutionStatus(lpi, &prosta, &solsta) );

   return ( solsta == MSK_SOL_STA_DUAL_INFEAS_CER
      || prosta == MSK_PRO_STA_DUAL_INFEAS
      || prosta == MSK_PRO_STA_PRIM_AND_DUAL_INFEAS );
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKsolstae solsta;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiHasPrimalRay (%d)\n", lpi->lpid);

   SCIP_ABORT_FALSE( getSolutionStatus(lpi, NULL, &solsta) );

   return (solsta == MSK_SOL_STA_DUAL_INFEAS_CER);
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   MSKsolstae solsta;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIP_ABORT_FALSE( getSolutionStatus(lpi, NULL, &solsta) );

   /* assume primal solution and ray is available if we used the primal simplex and the dual is proven to be infeasible */
   return (solsta == MSK_SOL_STA_DUAL_INFEAS_CER && lpi->lastalgo == MSK_OPTIMIZER_PRIMAL_SIMPLEX);
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   return SCIPlpiExistsDualRay(lpi);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKprostae prosta;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiIsPrimalFeasible (%d)\n", lpi->lpid);

   SCIP_ABORT_FALSE( getSolutionStatus(lpi, &prosta, NULL) );

   return (prosta == MSK_PRO_STA_PRIM_FEAS || prosta == MSK_PRO_STA_PRIM_AND_DUAL_FEAS || (prosta == MSK_PRO_STA_DUAL_INFEAS && lpi->lastalgo == MSK_OPTIMIZER_PRIMAL_SIMPLEX));
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKprostae prosta;
   MSKsolstae solsta;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiExistsDualRay (%d)\n", lpi->lpid);

   SCIP_ABORT_FALSE( getSolutionStatus(lpi, &prosta, &solsta) );

   return ( solsta == MSK_SOL_STA_PRIM_INFEAS_CER
      || prosta == MSK_PRO_STA_PRIM_INFEAS
      || prosta == MSK_PRO_STA_PRIM_AND_DUAL_INFEAS );
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKsolstae solsta;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiHasDualRay (%d)\n", lpi->lpid);

   SCIP_ABORT_FALSE( getSolutionStatus(lpi, NULL, &solsta) );

   return (solsta == MSK_SOL_STA_PRIM_INFEAS_CER);
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   return FALSE;
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   return SCIPlpiExistsPrimalRay(lpi);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKprostae prosta;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiIsDualFeasible (%d)\n", lpi->lpid);

   SCIP_ABORT_FALSE( getSolutionStatus(lpi, &prosta, NULL) );

   return (prosta == MSK_PRO_STA_DUAL_FEAS || prosta == MSK_PRO_STA_PRIM_AND_DUAL_FEAS);
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKsolstae solsta;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiIsOptimal (%d)\n", lpi->lpid);

   SCIP_ABORT_FALSE( getSolutionStatus(lpi, NULL, &solsta) );

   return (solsta == MSK_SOL_STA_OPTIMAL);
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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   return ( lpi->termcode == MSK_RES_OK
      || lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS
      || lpi->termcode == MSK_RES_TRM_MAX_TIME
      || lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE );
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   return ( lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE );
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   return ( lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS );
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   return ( lpi->termcode == MSK_RES_TRM_MAX_TIME );
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKsolstae solsta;
   SCIP_RETCODE retcode;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetInternalStatus (%d)\n", lpi->lpid);

   retcode = getSolutionStatus(lpi, NULL, &solsta);
   if ( retcode != SCIP_OKAY )
      return 0;

   return (int) solsta;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(success != NULL);

   SCIPdebugMessage("Calling SCIPlpiIgnoreInstability (%d)\n", lpi->lpid);

   *success = FALSE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(objval != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetObjval (%d)\n", lpi->lpid);

   /* TODO: tjek lighed med dual objektiv i de fleste tilfaelde. */

   if ( lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE )
   {
      /* if we reached the objective limit, return this value */
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_UPPER_OBJ_CUT, objval) );
   }
   else
   {
      /* otherwise get the value from Mosek */
      MOSEK_CALL( MSK_getprimalobj(lpi->task, lpi->lastsolvetype, objval) );
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
   MSKsolstae solsta;
   double* sux = NULL;
   int ncols = 0;
   int i;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetSol (%d)\n", lpi->lpid);

   if ( objval != NULL )
   {
      if ( lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE )
      {
         /* if we reached the objective limit, return this value */
         MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_UPPER_OBJ_CUT, objval) );
      }
      else
      {
         /* otherwise get the value from Mosek */
         MOSEK_CALL( MSK_getprimalobj(lpi->task, lpi->lastsolvetype, objval) );
      }
   }

   if( redcost )
   {
      MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
      SCIP_ALLOC( BMSallocMemoryArray(&sux, ncols) );
   }

   if ( primsol != NULL && lpi->lastalgo == MSK_OPTIMIZER_PRIMAL_SIMPLEX )
   {
      /* If the status shows that the dual is infeasible this is due to the primal being unbounded. In this case, we need
       * to compute a feasible solution by setting the objective to 0.
       */
      MOSEK_CALL( MSK_getsolutionstatus(lpi->task, MSK_SOL_BAS, NULL, &solsta) );
      if ( solsta == MSK_SOL_STA_DUAL_INFEAS_CER )
      {
         SCIP_Real* objcoefs;
         int j;

         MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
         SCIP_ALLOC( BMSallocMemoryArray(&objcoefs, ncols) );

         /* store old objective coefficients and set them to 0 */
         for (j = 0; j < ncols; ++j)
         {
            MOSEK_CALL( MSK_getcj(lpi->task, j, &objcoefs[j]) );
            MOSEK_CALL( MSK_putcj(lpi->task, j, 0.0) );
         }

         /* solve problem again */
         SCIP_CALL( SolveWSimplex(lpi) );

         /* At this point we assume that the problem is feasible, since we previously ran the primal simplex and it
          * produced a ray.
          */

         /* get primal solution */
         MOSEK_CALL( MSK_getsolution(lpi->task, lpi->lastsolvetype, NULL, NULL, NULL, NULL, NULL, activity,
               primsol, NULL, NULL, NULL, NULL, NULL, NULL) );

         /* restore objective */
         MOSEK_CALL( MSK_putcslice(lpi->task, 0, ncols, objcoefs) );

         /* resolve to restore original status */
         SCIP_CALL( SolveWSimplex(lpi) );

         BMSfreeMemoryArray(&objcoefs);
      }
      else
      {
         MOSEK_CALL( MSK_getsolution(lpi->task, lpi->lastsolvetype, NULL, NULL, NULL, NULL, NULL, activity,
               primsol, dualsol, NULL, NULL, redcost, sux, NULL) );
      }
   }
   else
   {
      MOSEK_CALL( MSK_getsolution(lpi->task, lpi->lastsolvetype, NULL, NULL, NULL, NULL, NULL, activity,
            primsol, dualsol, NULL, NULL, redcost, sux, NULL) );
   }

   /* the reduced costs are given by the difference of the slx and sux variables (third and second to last parameters) */
   if( redcost )
   {
      for( i = 0; i < ncols; i++ )
      {
         assert(sux != NULL);
         redcost[i] -= sux[i];
      }
      BMSfreeMemoryArray(&sux);
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(ray != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetPrimalRay (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getsolution(lpi->task, lpi->lastsolvetype, NULL, NULL, NULL, NULL, NULL, NULL, ray,
         NULL, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(dualfarkas != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetDualfarkas (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getsolution(lpi->task, lpi->lastsolvetype, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dualfarkas,
         NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(iterations != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetIterations (%d)\n", lpi->lpid);

   *iterations = lpi->itercount;

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(quality != NULL);
   SCIP_UNUSED(qualityindicator);

   *quality = SCIP_INVALID;

   return SCIP_OKAY;
}

/** handle singular basis */
static
SCIP_RETCODE handle_singular(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  basis,              /**< array of basis indices */
   MSKrescodee           res                 /**< result */
   )
{
   if (res == MSK_RES_ERR_BASIS_SINGULAR)
   {
      SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

      MOSEK_CALL( MSK_initbasissolve(lpi->task, basis) );
   }
   else
   {
      MOSEK_CALL( res );
   }

   return SCIP_OKAY;
}


/*
 * LP Basis Methods
 */

/** convert Mosek basis status to SCIP basis status
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
static
SCIP_RETCODE convertstat_mosek2scip(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             iscon,              /**< whether constraints/variables are considered */
   MSKstakeye*           sk,                 /**< status array of Mosek */
   int                   n,                  /**< size */
   int*                  stat                /**< status array of SCIP */
   )
{
   int i;

   assert(lpi->lastsolvetype == MSK_SOL_BAS);

   for( i = 0; i < n; i++ )
   {
      double sl;
      double su;

      switch (sk[i])
      {
      case MSK_SK_BAS:
         stat[i] = (int)SCIP_BASESTAT_BASIC;
         break;
      case MSK_SK_SUPBAS:
         stat[i] = (int)SCIP_BASESTAT_ZERO;
         break;
      case MSK_SK_FIX:
#if MSK_VERSION_MAJOR < 9
         MOSEK_CALL( MSK_getsolutioni(lpi->task, iscon ? MSK_ACC_CON : MSK_ACC_VAR, i, MSK_SOL_BAS, NULL, NULL, &sl, &su, NULL) );
#else
         if( iscon )
         {
            MOSEK_CALL( MSK_getslcslice(lpi->task, MSK_SOL_BAS, i, i+1, &sl ) );
            MOSEK_CALL( MSK_getsucslice(lpi->task, MSK_SOL_BAS, i, i+1, &su ) );
         }
         else
         {
            MOSEK_CALL( MSK_getslxslice(lpi->task, MSK_SOL_BAS, i, i+1, &sl ) );
            MOSEK_CALL( MSK_getsuxslice(lpi->task, MSK_SOL_BAS, i, i+1, &su ) );
         }
#endif

         if (sl < su) /* Negative reduced cost */
            stat[i] = (int)SCIP_BASESTAT_UPPER;
         else
            stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_UNK:
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_INF:
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_LOW:
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_UPR:
         stat[i] = (int)SCIP_BASESTAT_UPPER;
         break;
#if MSK_VERSION_MAJOR < 10
      case MSK_SK_END:
         break;
#endif
      default:
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** convert Mosek to SCIP basis status - slack variables
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
static
SCIP_RETCODE convertstat_mosek2scip_slack(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             iscon,              /**< whether constraints or variables are accessed */
   MSKstakeye*           sk,                 /**< Mosek basis status */
   int                   m,                  /**< size */
   int*                  stat                /**< status array */
   )
{
   int i;

   assert(lpi->lastsolvetype == MSK_SOL_BAS);

   for( i = 0; i < m; i++ )
   {
      double sl;
      double su;
      switch (sk[i])
      {
      case MSK_SK_BAS:
         stat[i] = (int)SCIP_BASESTAT_BASIC;
         break;
      case MSK_SK_SUPBAS:
         stat[i] = (int)SCIP_BASESTAT_ZERO;
         break;
      case MSK_SK_FIX:
#if MSK_VERSION_MAJOR < 9
         MOSEK_CALL( MSK_getsolutioni(lpi->task, iscon ? MSK_ACC_CON : MSK_ACC_VAR, i, MSK_SOL_BAS, NULL, NULL, &sl, &su, NULL) );
#else
         if( iscon )
         {
            MOSEK_CALL( MSK_getslcslice(lpi->task, MSK_SOL_BAS, i, i+1, &sl ) );
            MOSEK_CALL( MSK_getsucslice(lpi->task, MSK_SOL_BAS, i, i+1, &su ) );
         }
         else
         {
            MOSEK_CALL( MSK_getslxslice(lpi->task, MSK_SOL_BAS, i, i+1, &sl ) );
            MOSEK_CALL( MSK_getsuxslice(lpi->task, MSK_SOL_BAS, i, i+1, &su ) );
         }
#endif

         if (sl < su) /* Negative reduced cost */
            stat[i] = (int)SCIP_BASESTAT_UPPER;
         else
            stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_UNK:
      case MSK_SK_INF:
      case MSK_SK_UPR:
         stat[i] = (int)SCIP_BASESTAT_UPPER;
         break;
      case MSK_SK_LOW:
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
#if MSK_VERSION_MAJOR < 10
      case MSK_SK_END:
#endif
      default:
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** convert SCIP to Mosek basis status */
static
void convertstat_scip2mosek(
   const int*            stat,               /**< SCIP status array */
   int                   n,                  /**< size of array */
   MSKstakeye*           resstat             /**< resulting Mosek status array */
   )
{
   int i;
   for( i = 0; i < n; i++ )
   {
      switch (stat[i])
      {
      case SCIP_BASESTAT_LOWER:
         resstat[i] = MSK_SK_LOW;
         break;
      case SCIP_BASESTAT_BASIC:
         resstat[i] = MSK_SK_BAS;
         break;
      case SCIP_BASESTAT_UPPER:
         resstat[i] = MSK_SK_UPR;
         break;
      case SCIP_BASESTAT_ZERO:
         resstat[i] = MSK_SK_SUPBAS;
         break;
      default:
         SCIPABORT();
      }
   }
}

/** convert SCIP to Mosek basis status - slack variables */
static
void convertstat_scip2mosek_slack(
   const int*            stat,               /**< SCIP status array */
   int                   n,                  /**< size of array */
   MSKstakeye*           resstat             /**< resulting Mosek status array */
   )
{
   /* slacks are stored as -1 in Mosek, i.e., bounds are reversed compared to SCIP  */
   int i;

   for( i = 0; i < n; i++ )
   {
      switch (stat[i])
      {
      case SCIP_BASESTAT_LOWER:
         resstat[i] = MSK_SK_UPR; /* Reversed */
         break;
      case SCIP_BASESTAT_BASIC:
         resstat[i] = MSK_SK_BAS;
         break;
      case SCIP_BASESTAT_UPPER:
         resstat[i] = MSK_SK_LOW; /* Reversed */
         break;
      case SCIP_BASESTAT_ZERO:
         resstat[i] = MSK_SK_SUPBAS;
         break;
      default:
         SCIPABORT();
      }
   }
}

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
SCIP_RETCODE SCIPlpiGetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   int nrows;
   int ncols;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(lpi->lastsolvetype == MSK_SOL_BAS);

   SCIPdebugMessage("Calling SCIPlpiGetBase (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );

   SCIP_CALL( getbase(lpi, ncols, nrows) );

   if (cstat)
   {
      SCIP_CALL( convertstat_mosek2scip(lpi, FALSE, lpi->skx, ncols, cstat) );
   }

   if (rstat)
   {
      SCIP_CALL( convertstat_mosek2scip_slack(lpi, TRUE, lpi->skc, nrows, rstat) );
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
   int nrows;
   int ncols;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );

   assert(cstat != NULL || ncols == 0);
   assert(rstat != NULL || nrows == 0);

   SCIPdebugMessage("Calling SCIPlpiSetBase (%d)\n", lpi->lpid);

   SCIP_CALL( ensureStateMem(lpi, ncols, nrows) );

   convertstat_scip2mosek(cstat, ncols, lpi->skx);
   convertstat_scip2mosek_slack(rstat, nrows, lpi->skc);

   SCIP_CALL( setbase(lpi) );

   invalidateSolution(lpi);

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   )
{
   int nrows;
   int i;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(bind != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetBasisInd (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );

   SCIP_CALL( handle_singular(lpi, bind, MSK_initbasissolve(lpi->task, bind)) );

   for (i = 0; i < nrows; i++ )
   {
      if (bind[i] < nrows) /* row bind[i] is basic */
         bind[i] = -1 - bind[i];
      else                 /* column bind[i]-nrows is basic */
         bind[i] = bind[i] - nrows;
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
 *
 *  @todo check if this should invalidate the solution
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
   int nrows;
   int i;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(coef != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetBInvRow (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );

   /* set coefficient for slack variables to be 1 instead of -1 */
   MOSEK_CALL( MSK_putnaintparam(lpi->task, MSK_IPAR_BASIS_SOLVE_USE_PLUS_ONE_, MSK_ON) );

   /* prepare basis in Mosek, since we do not need the basis ourselves, we set the return parameter to NULL */
   SCIP_CALL( handle_singular(lpi, NULL, MSK_initbasissolve(lpi->task, NULL)) );

   /* initialize rhs of system to be a dense +/- unit vector (needed for MSK_solvewithbasis()) */
   for (i = 0; i < nrows; ++i)
      coef[i] = 0.0;
   coef[r] = 1.0; /* unit vector e_r */

   /* check whether we require a dense or sparse result vector */
   if ( ninds != NULL && inds != NULL )
   {
      *ninds = 1;
      inds[0]= r;

      /* solve transposed system */
#if MSK_VERSION_MAJOR < 10
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 1, ninds, inds, coef) );
#else
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 1, *ninds, inds, coef, ninds) );
#endif
      assert( *ninds <= nrows );
   }
   else
   {
      int* sub;
      int numnz;

      SCIP_ALLOC( BMSallocMemoryArray(&sub, nrows) );

      numnz = 1;
      sub[0] = r;

      /* solve transposed system */
#if MSK_VERSION_MAJOR < 10
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 1, &numnz, sub, coef) );
#else
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 1, numnz, sub, coef, &numnz) );
#endif
      assert( numnz <= nrows );

      BMSfreeMemoryArray(&sub);
   }
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );

   SCIPdebugMessage("End SCIPlpiGetBInvRow (%d)\n", lpi->lpid);

   return SCIP_OKAY;
}

/** get column of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 *
 *  @todo check that the result is in terms of the LP interface definition
 *
 *  @todo check if this should invalidate the solution
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
   int nrows;
   int i;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(coef != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetBInvCol (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );

   /* set coefficient for slack variables to be 1 instead of -1 */
   MOSEK_CALL( MSK_putnaintparam(lpi->task, MSK_IPAR_BASIS_SOLVE_USE_PLUS_ONE_, MSK_ON) );

   /* prepare basis in Mosek, since we do not need the basis ourselves, we set the return parameter to NULL */
   SCIP_CALL( handle_singular(lpi, NULL, MSK_initbasissolve(lpi->task, NULL)) );

   /* initialize rhs of system to be a dense +/- unit vector (needed for MSK_solvewithbasis()) */
   for (i = 0; i < nrows; ++i)
      coef[i] = 0.0;
   coef[c] = 1.0; /* unit vector e_c */

   /* check whether we require a dense or sparse result vector */
   if ( ninds != NULL && inds != NULL )
   {
      *ninds = 1;
      inds[0]= c;

#if MSK_VERSION_MAJOR < 10
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 0, ninds, inds, coef) );
#else
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 0, *ninds, inds, coef, ninds) );
#endif
      assert( *ninds <= nrows );
   }
   else
   {
      int* sub;
      int numnz;

      SCIP_ALLOC( BMSallocMemoryArray(&sub, nrows) );

      numnz = 1;
      sub[0]= c;

#if MSK_VERSION_MAJOR < 10
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 0, &numnz, sub, coef) );
#else
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 0, numnz, sub, coef, &numnz) );
#endif
      assert( numnz <= nrows );

      BMSfreeMemoryArray(&sub);
   }
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );

   return SCIP_OKAY;
}

/** get row of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 *
 *  @todo check that the result is in terms of the LP interface definition
 *
 *  @todo check if this should invalidate the solution
 */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            coef,               /**< vector to return coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{  /*lint --e{715}*/
   int nrows;
   int ncols;
   int numnz;
   int* csub;
   int didalloc = 0;
   double* cval;
   double* binv;
   int i;
   int k;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(coef != NULL);
   SCIP_UNUSED(inds);

   SCIPdebugMessage("Calling SCIPlpiGetBInvARow (%d)\n", lpi->lpid);

   /* currently only return dense result */
   if ( ninds != NULL )
      *ninds = -1;

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );

   SCIP_ALLOC( BMSallocMemoryArray(&csub, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&cval, nrows) );

   if( binvrow == NULL )
   {
      didalloc = 1;

      /* get dense vector */
      SCIP_ALLOC( BMSallocMemoryArray(&binv, nrows) );
      SCIP_CALL( SCIPlpiGetBInvRow(lpi, row, binv, NULL, NULL) );
   }
   else
      binv = (SCIP_Real*)binvrow;

   /* binvrow*A */
   for (i = 0; i < ncols; ++i)
   {
      MOSEK_CALL( MSK_getacol(lpi->task, i, &numnz, csub, cval) );

      /* compute dense vector */
      coef[i] = 0.0;
      for (k = 0; k < numnz; ++k)
      {
         assert( 0 <= csub[k] && csub[k] < nrows );
         coef[i] += binv[csub[k]] * cval[k];
      }
   }

   /* free memory arrays */
   BMSfreeMemoryArray(&cval);
   BMSfreeMemoryArray(&csub);

   if ( didalloc > 0 )
   {
      BMSfreeMemoryArray(&binv);
   }

   return SCIP_OKAY;
}

/** get column of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @todo check if this should invalidate the solution
 */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef,               /**< vector to return coefficients of the columnn */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{  /*lint --e{715}*/
   SCIP_Real* val;
   int nrows;
   int numnz;
   int i;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(coef != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetBInvACol (%d)\n", lpi->lpid);

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );
   MOSEK_CALL( MSK_getacolnumnz(lpi->task, c, &numnz) );
   SCIP_ALLOC( BMSallocMemoryArray(&val, numnz+1) );

   /* set coefficient for slack variables to be 1 instead of -1 */
   MOSEK_CALL( MSK_putnaintparam(lpi->task, MSK_IPAR_BASIS_SOLVE_USE_PLUS_ONE_, MSK_ON) );

   /* prepare basis in Mosek, since we do not need the basis ourselves, we set the return parameter to NULL */
   SCIP_CALL( handle_singular(lpi, NULL, MSK_initbasissolve(lpi->task, NULL)) );

   /* init coefficients */
   for (i = 0; i < nrows; ++i)
      coef[i] = 0.0;

   /* check whether we require a dense or sparse result vector */
   if ( ninds != NULL && inds != NULL )
   {
      /* fill column into dense vector */
      MOSEK_CALL( MSK_getacol(lpi->task, c, &numnz, inds, val) );
      for (i = 0; i < numnz; ++i)
      {
         assert( 0 <= inds[i] && inds[i] < nrows );
         coef[inds[i]] = val[i];
      }

      *ninds = numnz;

#if MSK_VERSION_MAJOR < 10
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 0, ninds, inds, coef) );
#else
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 0, *ninds, inds, coef, ninds) );
#endif
      assert( *ninds <= nrows );
   }
   else
   {
      int* sub;
      SCIP_ALLOC( BMSallocMemoryArray(&sub, nrows) );

      /* fill column into dense vector */
      MOSEK_CALL( MSK_getacol(lpi->task, c, &numnz, sub, val) );
      for (i = 0; i < numnz; ++i)
      {
         assert( 0 <= sub[i] && sub[i] < nrows );
         coef[sub[i]] = val[i];
      }

#if MSK_VERSION_MAJOR < 10
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 0, &numnz, sub, coef) );
#else
      MOSEK_CALL( MSK_solvewithbasis(lpi->task, 0, numnz, sub, coef, &numnz) );
#endif

      if ( ninds != NULL )
         *ninds = numnz;

      BMSfreeMemoryArray(&sub);
   }

   BMSfreeMemoryArray(&val);
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );

   return SCIP_OKAY;
}


/*
 * LP State Methods
 */

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
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->skx, colpacketNum(ncols)) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->skc, rowpacketNum(nrows)) );

   (*lpistate)->solsta = MSK_SOL_STA_UNKNOWN;
   (*lpistate)->num    = -1;
   (*lpistate)->ncols  = ncols;
   (*lpistate)->nrows  = nrows;

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

   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->skx, colpacketNum((*lpistate)->ncols));
   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->skc, rowpacketNum((*lpistate)->nrows));
   BMSfreeBlockMemory(blkmem, lpistate);
}

#ifndef NDEBUG
/** check state
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
static
SCIP_RETCODE checkState1(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   n,                  /**< number of rows or columns */
   MSKstakeye*           sk,                 /**< basis status */
   SCIP_Bool             isrow               /**< whether rows/columns are considered */
   )
{
   char xc;
   int i;

   assert(lpi != NULL);
   assert(lpi->lastsolvetype == MSK_SOL_BAS);

#ifdef SCIP_DEBUG
   if( !isrow )
      xc = 'x';
   else
      xc = 'c';
#endif

   /* printout for all except LOW, UPR, FIX and BAS with sl[xc]==su[xc] */
   for( i = 0; i < n; i++ )
   {
      double sl;
      double su;
      switch (sk[i])
      {
      case MSK_SK_UNK:
         SCIPdebugMessage("STATE[%d]: %c[%d] = unk\n", lpi->optimizecount, xc, i);
         break;
      case MSK_SK_BAS:
         /* the following function is deprecated */
#if MSK_VERSION_MAJOR < 9
         MOSEK_CALL( MSK_getsolutioni(lpi->task, isrow ? MSK_ACC_CON : MSK_ACC_VAR, i, MSK_SOL_BAS, NULL, NULL, &sl, &su, NULL) );
#else
         if( isrow )
         {
            MOSEK_CALL( MSK_getslcslice(lpi->task, MSK_SOL_BAS, i, i+1, &sl ) );
            MOSEK_CALL( MSK_getsucslice(lpi->task, MSK_SOL_BAS, i, i+1, &su ) );
         }
         else
         {
            MOSEK_CALL( MSK_getslxslice(lpi->task, MSK_SOL_BAS, i, i+1, &sl ) );
            MOSEK_CALL( MSK_getsuxslice(lpi->task, MSK_SOL_BAS, i, i+1, &su ) );
         }
#endif
         if (fabs(sl-su) > DEBUG_CHECK_STATE_TOL)
         {
            SCIPdebugMessage("STATE[%d]: %c[%d] = bas, sl%c = %g, su%c = %g\n", lpi->optimizecount, xc, i, xc, sl, xc, su);
         }
         break;
      case MSK_SK_SUPBAS:
         SCIPdebugMessage("STATE[%d]: %c[%d] = supbas\n", lpi->optimizecount, xc, i);
         break;
      case MSK_SK_LOW:
      case MSK_SK_UPR:
      case MSK_SK_FIX:
         break;
      case MSK_SK_INF:
         SCIPdebugMessage("STATE[%d]: %c[%d] = inf\n", lpi->optimizecount, xc, i);
         break;
      default:
         SCIPdebugMessage("STATE[%d]: %c[%d] = unknown status <%d>\n", lpi->optimizecount, xc, i, sk[i]);
         break;
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** check state
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
static
SCIP_RETCODE checkState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns */
   int                   nrows               /**< number of rows */
   )
{
   assert(lpi != NULL);
   assert(lpi->lastsolvetype == MSK_SOL_BAS);

   SCIP_CALL( checkState1(lpi, ncols, lpi->skx, FALSE) );
   SCIP_CALL( checkState1(lpi, nrows, lpi->skc, TRUE) );

   return SCIP_OKAY;
 }
#endif

/** store row and column basis status in a packed LPi state object
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
static
SCIP_RETCODE lpistatePack(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< pointer to LPi state data */
   )
{
   int *skxi = (int *) lpi->skx; /* Used as temp. buffer */
   int *skci = (int *) lpi->skc; /* Used as temp. buffer */

   assert(sizeof(int) == sizeof(MSKstakeye)); /*lint !e506*/
   assert(lpi != NULL);
   assert(lpistate != NULL);
   assert(lpi->lastsolvetype == MSK_SOL_BAS);

   SCIP_CALL( convertstat_mosek2scip(lpi, FALSE, lpi->skx, lpistate->ncols, skxi) );
   SCIP_CALL( convertstat_mosek2scip_slack(lpi, TRUE, lpi->skc, lpistate->nrows, skci) );

   SCIPencodeDualBit(skxi, lpistate->skx, lpistate->ncols);
   SCIPencodeDualBit(skci, lpistate->skc, lpistate->nrows);

   return SCIP_OKAY;
}

/** unpacks row and column basis status from a packed LPi state object */
static
void lpistateUnpack(
   const SCIP_LPISTATE*  lpistate,           /**< pointer to LPi state data */
   MSKstakeye*           skx,                /**< basis status for columns */
   MSKstakeye*           skc                 /**< basis status for rows */
   )
{
   assert(sizeof(int) == sizeof(MSKstakeye)); /*lint !e506*/

   SCIPdecodeDualBit(lpistate->skx, (int*) skx, lpistate->ncols);
   SCIPdecodeDualBit(lpistate->skc, (int*) skc, lpistate->nrows);

   convertstat_scip2mosek((int*) skx, lpistate->ncols, skx);
   convertstat_scip2mosek_slack((int*) skc, lpistate->nrows, skc);
}

/** stores LP state (like basis information) into lpistate object
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
SCIP_RETCODE SCIPlpiGetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
   int gotbasicsol;
   int nrows;
   int ncols;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(lpistate != NULL);
   assert(blkmem != NULL);
   assert(lpi->lastsolvetype == MSK_SOL_BAS);

   SCIPdebugMessage("Calling SCIPlpiGetState (%d)\n", lpi->lpid);

   *lpistate = NULL;

   MOSEK_CALL( MSK_solutiondef(lpi->task, MSK_SOL_BAS, &gotbasicsol) );

   if ( gotbasicsol == 0 || SCIPlpiExistsDualRay(lpi) || lpi->clearstate )
      return SCIP_OKAY;

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );

   /* allocate lpistate data */
   SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows) );

   lpistate[0]->num = lpi->optimizecount;

   MOSEK_CALL( MSK_getsolutionstatus(lpi->task, MSK_SOL_BAS, NULL, &lpistate[0]->solsta) );

   SCIP_CALL( getbase(lpi, ncols, nrows) );

#ifndef NDEBUG
   SCIP_CALL( checkState(lpi, ncols, nrows) );
#endif

   SCIP_CALL( lpistatePack(lpi, lpistate[0]) );

   SCIPdebugMessage("Store into state from iter : %d\n", lpi->optimizecount);

   /*    if (r != SCIP_OKAY)
    *    lpistateFree(lpistate, blkmem );
    */

   return SCIP_OKAY;
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetState()
 */
SCIP_RETCODE SCIPlpiSetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const SCIP_LPISTATE*  lpistate            /**< LP state information (like basis information), or NULL */
   )
{  /*lint --e{715}*/
   int nrows;
   int ncols;
   int i;

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(blkmem != NULL);
#ifdef SCIP_DISABLED_CODE
   assert(lpi->lastsolvetype == MSK_SOL_BAS);
#endif

   if (lpistate == NULL)
   {
      SCIPdebugMessage("Setting NULL state\n");
      return SCIP_OKAY;
   }

   if (lpistate->nrows == 0 || lpistate->ncols == 0)
      return SCIP_OKAY;

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
   assert(lpistate->nrows <= nrows);
   assert(lpistate->ncols <= ncols);

   SCIP_CALL( ensureStateMem(lpi, ncols, nrows) );

#ifdef SCIP_DISABLED_CODE
   SCIP_CALL( getbase(lpi, ncols, nrows) ); /* Why do we need to get the basis ????? */
#endif

   lpistateUnpack(lpistate, lpi->skx, lpi->skc);

   /* extend the basis to the current LP beyond the previously existing columns */
   for (i = lpistate->ncols; i < ncols; ++i)
   {
      SCIP_Real lb;
      SCIP_Real ub;
#if MSK_VERSION_MAJOR < 9
      MOSEK_CALL( MSK_getboundslice(lpi->task, MSK_ACC_VAR, i, i, NULL, &lb, &ub) );
#else
      MOSEK_CALL( MSK_getvarboundslice(lpi->task, i, i, NULL, &lb, &ub) );
#endif
      if ( SCIPlpiIsInfinity(lpi, REALABS(lb)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         if ( SCIPlpiIsInfinity(lpi, REALABS(ub)) )
            lpi->skx[i] = MSK_SK_SUPBAS;  /* variable is free (super basic) */
         else
            lpi->skx[i] = MSK_SK_UPR;     /* use finite upper bound */
      }
      else
         lpi->skx[i] = MSK_SK_LOW;        /* use finite lower bound */
   }
   for (i = lpistate->nrows; i < nrows; ++i)
      lpi->skc[i] = MSK_SK_BAS;

   /* load basis information into MOSEK */
   SCIP_CALL( setbase(lpi) );

   invalidateSolution(lpi);

   SCIPdebugMessage("Store from state into task iter : %d with solsta : %d\n", lpistate->num, lpistate->solsta);

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   lpi->clearstate = TRUE;

   return SCIP_OKAY;
}

/** frees LP state information */
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(lpistate != NULL);
   assert(blkmem != NULL);

   SCIPdebugMessage("Calling SCIPlpiFreeState (%d)\n", lpi->lpid);

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiHasStateBasis (%d)\n", lpi->lpid);

   return ( lpistate != NULL && lpistate->num >= 0);
}

/** reads LP state (like basis information) from a file
 *
 * @note last solve call must have been either simplex or barrier with crossover or base must have been set manually
 */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(fname != NULL);

   SCIPdebugMessage("reading LP state from file <%s>\n", fname);

   lpi->clearstate = FALSE;

   MOSEK_CALL( MSK_readsolution(lpi->task, MSK_SOL_BAS, fname) );

   return SCIP_OKAY;
}

/** writes LPi state (i.e. basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   int v;
   int nvars;
   int c = 0;
   int nconss;
   SCIP_Bool emptyname = FALSE;
   char name[SCIP_MAXSTRLEN];

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(fname != NULL);
   assert(lpi->lastsolvetype == MSK_SOL_BAS);

   SCIPdebugMessage("writing LP state to file <%s>\n", fname);

   if( lpi->clearstate )
   {
      SCIPdebugMessage("No LP state written, since it was cleared after the last solve \n");
      return SCIP_OKAY;
   }

   /* If any rows or columns have empty names, MOSEK will make up names like C1 and X1, but will no
    * longer recognize them when reading the same state file back in, therefore we return an error in
    * this case
    */
   MOSEK_CALL( MSK_getnumvar(lpi->task, &nvars) );
   for( v = 0; v < nvars; v++ )
   {
      MOSEK_CALL( MSK_getvarname(lpi->task, v, SCIP_MAXSTRLEN, name) );
      if( strcmp(name, "") == 0 )
      {
         emptyname = TRUE;
         break;
      }
   }
   if( !emptyname )
   {
      MOSEK_CALL( MSK_getnumcon(lpi->task, &nconss) );
      for( c = 0; c < nconss; c++ )
      {
         MOSEK_CALL( MSK_getconname(lpi->task, c, SCIP_MAXSTRLEN, name) );
         if( strcmp(name, "") == 0 )
         {
            emptyname = TRUE;
            break;
         }
      }
   }

   if( emptyname )
   {
      SCIPmessagePrintWarning(lpi->messagehdlr, "Writing LP state with unnamed %s %d, using default"
            " names instead. Note that this state cannot be read back in later!\n",
            v < nvars ? "variable" : "constraint", v < nvars ? v : c);
   }

   /* set parameter to be able to write */
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_SOL_HEAD, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_SOL_VARIABLES, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_SOL_CONSTRAINTS, MSK_ON) );

   MOSEK_CALL( MSK_writesolution(lpi->task, MSK_SOL_BAS, fname) );

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
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(lpinorms != NULL);
   assert(blkmem != NULL);

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
{  /*lint --e{715}*/
   assert(lpi != NULL);
   SCIP_UNUSED(blkmem);
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
{  /*lint --e{715}*/
   assert(lpi != NULL);
   SCIP_UNUSED(blkmem);
   assert(lpinorms == NULL);

   /* no work necessary */
   return SCIP_OKAY;
}

/**@} */

/*
 * Parameter Methods
 */

/** constant array containing the parameter names */
static const char* paramname[] = {
   "SCIP_LPPAR_FROMSCRATCH",                 /* solver should start from scratch at next call? */
   "SCIP_LPPAR_FASTMIP",                     /* fast mip setting of LP solver */
   "SCIP_LPPAR_SCALING",                     /* which scaling should LP solver use? */
   "SCIP_LPPAR_PRESOLVING",                  /* should LP solver use presolving? */
   "SCIP_LPPAR_PRICING",                     /* pricing strategy */
   "SCIP_LPPAR_LPINFO",                      /* should LP solver output information to the screen? */
   "SCIP_LPPAR_FEASTOL",                     /* feasibility tolerance for primal variables and slacks */
   "SCIP_LPPAR_DUALFEASTOL",                 /* feasibility tolerance for dual variables and reduced costs */
   "SCIP_LPPAR_BARRIERCONVTOL",              /* convergence tolerance used in barrier algorithm */
   "SCIP_LPPAR_OBJLIM",                      /* objective limit (stop if objective is known be larger/smaller than limit for min/max-imization) */
   "SCIP_LPPAR_LPITLIM",                     /* LP iteration limit, greater than or equal 0 */
   "SCIP_LPPAR_LPTILIM",                     /* LP time limit, positive */
   "SCIP_LPPAR_MARKOWITZ",                   /* Markowitz tolerance */
   "SCIP_LPPAR_ROWREPSWITCH",                /* simplex algorithm shall use row representation of the basis
                                              * if number of rows divided by number of columns exceeds this value */
   "SCIP_LPPAR_THREADS",                     /* number of threads used to solve the LP */
   "SCIP_LPPAR_CONDITIONLIMIT",              /* maximum condition number of LP basis counted as stable */
   "SCIP_LPPAR_TIMING",                      /* type of timer (1 - cpu, 2 - wallclock, 0 - off) */
   "SCIP_LPPAR_RANDOMSEED",                  /* inital random seed, e.g. for perturbations in the simplex (0: LP default) */
   "SCIP_LPPAR_POLISHING",                   /* set solution polishing (0 - disable, 1 - enable) */
   "SCIP_LPPAR_REFACTOR"                     /* set refactorization interval (0 - automatic) */
};

/** method mapping parameter index to parameter name */
static
const char* paramty2str(
   SCIP_LPPARAM          type
   )
{  /*lint --e{641}*/
   /* check whether the parameters are in the right order */
   /*lint --e{506}*/
   assert(SCIP_LPPAR_FROMSCRATCH == 0);      /* solver should start from scratch at next call? */
   assert(SCIP_LPPAR_FASTMIP == 1);          /* fast mip setting of LP solver */
   assert(SCIP_LPPAR_SCALING == 2);          /* which scaling should LP solver use? */
   assert(SCIP_LPPAR_PRESOLVING == 3);       /* should LP solver use presolving? */
   assert(SCIP_LPPAR_PRICING == 4);          /* pricing strategy */
   assert(SCIP_LPPAR_LPINFO == 5);           /* should LP solver output information to the screen? */
   assert(SCIP_LPPAR_FEASTOL == 6);          /* feasibility tolerance for primal variables and slacks */
   assert(SCIP_LPPAR_DUALFEASTOL == 7);      /* feasibility tolerance for dual variables and reduced costs */
   assert(SCIP_LPPAR_BARRIERCONVTOL == 8);   /* convergence tolerance used in barrier algorithm */
   assert(SCIP_LPPAR_OBJLIM == 9);           /* objective limit (stop if objective is known be larger/smaller than limit for min/max-imization) */
   assert(SCIP_LPPAR_LPITLIM == 10);         /* LP iteration limit, greater than or equal 0 */
   assert(SCIP_LPPAR_LPTILIM == 11);         /* LP time limit, positive */
   assert(SCIP_LPPAR_MARKOWITZ == 12);       /* Markowitz tolerance */
   assert(SCIP_LPPAR_ROWREPSWITCH == 13);    /* row representation switch */
   assert(SCIP_LPPAR_THREADS == 14);         /* number of threads used to solve the LP */
   assert(SCIP_LPPAR_CONDITIONLIMIT == 15);  /* maximum condition number of LP basis counted as stable */
   assert(SCIP_LPPAR_TIMING == 16);          /* type of timer (1 - cpu, 2 - wallclock, 0 - off) */
   assert(SCIP_LPPAR_RANDOMSEED == 17);      /* inital random seed, e.g. for perturbations in the simplex (0: LP default) */
   assert(SCIP_LPPAR_POLISHING == 18);       /* set solution polishing (0 - disable, 1 - enable) */
   assert(SCIP_LPPAR_REFACTOR == 19);        /* set refactorization interval (0 - automatic) */

   return paramname[type];
}

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiGetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
                              )
{  /*lint --e{641}*/
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(ival != NULL);

   SCIPdebugMessage("getting int parameter %s\n", paramty2str(type));

   switch (type)
   {
   case SCIP_LPPAR_FROMSCRATCH:               /* solver should start from scratch at next call? */
      *ival = (int) lpi->fromscratch;
      break;
   case SCIP_LPPAR_FASTMIP:                   /* fast mip setting of LP solver */
      return  SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_SCALING:                   /* should LP solver use scaling? */
      MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_SCALING, ival) );
      if( *ival == MSK_SCALING_NONE )
         *ival = 0;
      else if( *ival == MSK_SCALING_FREE )
         *ival = 1;
#if MSK_VERSION_MAJOR < 10
      else if( *ival == MSK_SCALING_AGGRESSIVE )
         *ival = 2;
#endif
      else /* MSK_SCALING_MODERATE should not be used by the interface */
         return SCIP_PARAMETERWRONGVAL;
      break;
   case SCIP_LPPAR_PRESOLVING:                /* should LP solver use presolving? */
      MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_PRESOLVE_USE, ival) );
      *ival = (*ival != MSK_PRESOLVE_MODE_OFF);
      break;
   case SCIP_LPPAR_PRICING:                   /* pricing strategy */
      *ival = lpi->pricing;
      break;
   case SCIP_LPPAR_LPINFO:                    /* should LP solver output information to the screen? */
      *ival = (int) lpi->lpinfo;
      break;
   case SCIP_LPPAR_LPITLIM:                   /* LP iteration limit */
      MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, ival) );
      break;
   case SCIP_LPPAR_THREADS:                   /* number of threads */
      MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_NUM_THREADS, ival) );
      break;
   case SCIP_LPPAR_REFACTOR:                  /* refactorization interval */
      MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_REFACTOR_FREQ, ival) );
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
   static int pricing[7] =
   {
      (int)MSK_SIM_SELECTION_SE,             /**< mosek pricing for SCIP_PRICING_LPIDEFAULT */
      (int)MSK_SIM_SELECTION_FREE,           /**< mosek pricing for SCIP_PRICING_AUTO */
      (int)MSK_SIM_SELECTION_FULL,           /**< mosek pricing for SCIP_PRICING_FULL */
      (int)MSK_SIM_SELECTION_PARTIAL,        /**< mosek pricing for SCIP_PRICING_PARTIAL */
      (int)MSK_SIM_SELECTION_SE,             /**< mosek pricing for SCIP_PRICING_STEEP */
      (int)MSK_SIM_SELECTION_ASE,            /**< mosek pricing for SCIP_PRICING_STEEPQSTART */
      (int)MSK_SIM_SELECTION_DEVEX,          /**< mosek pricing for SCIP_PRICING_DEVEX */
   };

   /*lint --e{506}*/
   assert((int)SCIP_PRICING_LPIDEFAULT == 0);
   assert((int)SCIP_PRICING_AUTO == 1);
   assert((int)SCIP_PRICING_FULL == 2);
   assert((int)SCIP_PRICING_PARTIAL == 3);
   assert((int)SCIP_PRICING_STEEP == 4);
   assert((int)SCIP_PRICING_STEEPQSTART == 5);
   assert((int)SCIP_PRICING_DEVEX == 6);

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiSetIntpar (%d) Parameter=<%s>  Value=<%d>\n", lpi->lpid, paramty2str(type), ival);

   switch (type)
   {
   case SCIP_LPPAR_FROMSCRATCH:               /* solver should start from scratch at next call? */
      lpi->fromscratch = (SCIP_Bool) ival;
      break;
   case SCIP_LPPAR_FASTMIP:                   /* fast mip setting of LP solver */
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_SCALING:                   /* should LP solver use scaling? */
#if MSK_VERSION_MAJOR < 10
      assert( ival >= 0 && ival <= 2 );
#else
      assert( ival >= 0 && ival <= 1 );
#endif
      if( ival == 0 )
      {
         MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_SCALING, MSK_SCALING_NONE) );
         MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_NONE) );
      }
#if MSK_VERSION_MAJOR < 10
      else if( ival == 1 )
      {
         MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_SCALING, MSK_SCALING_FREE) );
         MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_FREE) );
      }
      else
      {
         MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_SCALING, MSK_SCALING_AGGRESSIVE) );
         MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_AGGRESSIVE) );
      }
#else
      else
      {
         MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_SCALING, MSK_SCALING_FREE) );
         MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_FREE) );
      }
#endif

      break;
   case SCIP_LPPAR_PRESOLVING:                /* should LP solver use presolving? */
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_PRESOLVE_USE,
            ival ? MSK_PRESOLVE_MODE_FREE : MSK_PRESOLVE_MODE_OFF) );
      break;
   case SCIP_LPPAR_PRICING:                   /* pricing strategy */
      assert(ival >= 0 && ival <= (int)SCIP_PRICING_DEVEX);
      lpi->pricing = (SCIP_PRICING)ival;

      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_PRIMAL_SELECTION, pricing[ival]) );
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_DUAL_SELECTION, pricing[ival]) );

      /* for certain pricing values, do not use restricted pricing */
      if( lpi->pricing == SCIP_PRICING_PARTIAL || lpi->pricing == SCIP_PRICING_AUTO )
         lpi->restrictselectdef = 50;
      else
         lpi->restrictselectdef = 0;

      break;
   case SCIP_LPPAR_LPINFO:
      /* should LP solver output information to the screen? */
#if FORCE_MOSEK_LOG
      SCIPdebugMessage("Ignoring log setting!\n");
#else
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, ival ? 4 : MSK_OFF) );
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG_SIM, ival ? 4 : MSK_OFF) );
#endif
      lpi->lpinfo = (SCIP_Bool) ival;
      break;
   case SCIP_LPPAR_LPITLIM:                   /* LP iteration limit */
#if DEBUG_PARAM_SETTING
      if( ival )
      {
         SCIPdebugMessage("Setting max iter to : %d\n", ival);
      }
#endif
      /* 0 <= ival, 0 stopping immediately */
      assert( ival >= 0 );
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, ival) );
      break;
   case SCIP_LPPAR_THREADS:                   /* number of threads (0 => MOSEK chooses) */
      assert(ival >= 0);
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_NUM_THREADS, ival) );
      break;
   case SCIP_LPPAR_REFACTOR:                  /* refactorization interval */
      assert(ival >= 0);
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_REFACTOR_FREQ, ival) );
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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(dval != NULL);

   SCIPdebugMessage("getting real parameter %s\n", paramty2str(type));

   switch (type)
   {
   case SCIP_LPPAR_FEASTOL:                   /* feasibility tolerance for primal variables and slacks */
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_BASIS_TOL_X, dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:               /* feasibility tolerance for dual variables and reduced costs */
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_BASIS_TOL_S, dval) );
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:            /* convergence tolerance used in barrier algorithm */
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_INTPNT_TOL_REL_GAP, dval) );
      break;
   case SCIP_LPPAR_OBJLIM:                    /* objective limit */
   {
      MSKobjsensee objsen;
      MOSEK_CALL( MSK_getobjsense(lpi->task, &objsen) );
      if (objsen == MSK_OBJECTIVE_SENSE_MINIMIZE)
      {
         MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_UPPER_OBJ_CUT, dval) );
      }
      else /* objsen == MSK_OBJECTIVE_SENSE_MAX */
      {
         MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_LOWER_OBJ_CUT, dval) );
      }
      break;
   }
   case SCIP_LPPAR_LPTILIM:                   /* LP time limit */
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_OPTIMIZER_MAX_TIME, dval) );
      break;
   case SCIP_LPPAR_MARKOWITZ:                 /* Markowitz tolerance */
   default:
      return SCIP_PARAMETERUNKNOWN;
   } /*lint !e788*/

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
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("setting real parameter %s to %g\n", paramty2str(type), dval);

   /**@todo Limits shouldn't be hardcoded */

   switch (type)
   {
   case SCIP_LPPAR_FEASTOL:                   /* feasibility tolerance for primal variables and slacks */
      assert( dval > 0.0 );
      /* 1e-9 <= dval <= inf */
      if( dval < 1e-9 )
         dval = 1e-9;

      MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_BASIS_TOL_X, dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:               /* feasibility tolerance for dual variables and reduced costs */
      assert( dval > 0.0 );
      /* 1e-9 <= dval <= inf */
      if( dval < 1e-9 )
         dval = 1e-9;

      MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_BASIS_TOL_S, dval) );
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:            /* convergence tolerance used in barrier algorithm */
      /* 1e-14 <= dval <= inf */
      assert( dval >= 0.0 );
      if( dval < 1e-14 )
         dval = 1e-14;

      MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_INTPNT_TOL_REL_GAP, dval) );
      break;
   case SCIP_LPPAR_OBJLIM:                    /* objective limit */
   {
      /* no restriction on dval */
      MSKobjsensee objsen;
      MOSEK_CALL( MSK_getobjsense(lpi->task, &objsen) );
      if (objsen == MSK_OBJECTIVE_SENSE_MINIMIZE)
      {
         MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_UPPER_OBJ_CUT, dval) );
      }
      else /* objsen == MSK_OBJECTIVE_SENSE_MAX */
      {
         MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_LOWER_OBJ_CUT, dval) );
      }
      break;
   }
   case SCIP_LPPAR_LPTILIM:                   /* LP time limit */
      assert( dval > 0.0 );
      /* mosek requires 0 <= dval
       *
       * However for consistency we assert the timelimit to be strictly positive.
       */
      MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_OPTIMIZER_MAX_TIME, dval) );
      break;
   case SCIP_LPPAR_MARKOWITZ:                 /* Markowitz tolerance */
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


/*
 * Numerical Methods
 */


/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   return MSK_INFINITY;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);

   return IS_POSINF(val);
}


/*
 * File Interface Methods
 */


/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
#if MSK_VERSION_MAJOR < 9
   int olddataformat;
#endif

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(fname != NULL);

   SCIPdebugMessage("Calling SCIPlpiReadLP (%d), filename <%s>\n", lpi->lpid, fname);

#if MSK_VERSION_MAJOR < 9
   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_READ_DATA_FORMAT, &olddataformat) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_READ_DATA_FORMAT, MSK_DATA_FORMAT_LP) );
   MOSEK_CALL( MSK_readdata(lpi->task, fname) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_READ_DATA_FORMAT, olddataformat) );
#else
   MOSEK_CALL( MSK_readdataformat(lpi->task, fname, MSK_DATA_FORMAT_LP, MSK_COMPRESS_FREE) );
#endif

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
#if MSK_VERSION_MAJOR < 9
   int olddataformat;
#endif

   assert(lpi != NULL);
   assert(lpi->mosekenv != NULL);
   assert(lpi->task != NULL);
   assert(fname != NULL);
#if MSK_VERSION_MAJOR >= 9
   /* Mosek 9 derives file format from given filename */
   assert(strstr(fname, ".lp") != NULL);
#endif

   SCIPdebugMessage("Calling SCIPlpiWriteLP (%d), filename <%s>\n", lpi->lpid, fname);

#if MSK_VERSION_MAJOR < 9
   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_WRITE_DATA_FORMAT, &olddataformat) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_DATA_FORMAT, MSK_DATA_FORMAT_LP) );
   MOSEK_CALL( MSK_writedata(lpi->task, fname) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_DATA_FORMAT, olddataformat) );
#else
   MOSEK_CALL( MSK_writedata(lpi->task, fname) );
#endif

   return SCIP_OKAY;
}
