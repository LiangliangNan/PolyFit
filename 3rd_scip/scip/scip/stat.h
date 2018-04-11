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

/**@file   stat.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for problem statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STAT_H__
#define __SCIP_STAT_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_prob.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_mem.h"
#include "scip/pub_message.h"
#include "scip/concurrent.h"

#include "scip/struct_stat.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates problem statistics data */
extern
SCIP_RETCODE SCIPstatCreate(
   SCIP_STAT**           stat,               /**< pointer to problem statistics data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem, or NULL */
   SCIP_PROB*            origprob,           /**< original problem, or NULL */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** frees problem statistics data */
extern
SCIP_RETCODE SCIPstatFree(
   SCIP_STAT**           stat,               /**< pointer to problem statistics data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** diables the collection of any statistic for a variable */
extern
void SCIPstatDisableVarHistory(
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** enables the collection of statistics for a variable */
extern
void SCIPstatEnableVarHistory(
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** marks statistics to be able to reset them when solving process is freed */
extern
void SCIPstatMark(
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** reset statistics to the data before solving started */
extern
void SCIPstatReset(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem, or NULL */
   SCIP_PROB*            origprob            /**< original problem, or NULL */
   );

/** reset implication counter */
extern
void SCIPstatResetImplications(
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** reset presolving and current run specific statistics */
extern
void SCIPstatResetPresolving(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem, or NULL */
   SCIP_PROB*            origprob            /**< original problem, or NULL */
   );

/* reset primal-dual integral */
extern
void SCIPstatResetPrimalDualIntegral(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             partialreset        /**< should time and integral value be kept? (in combination with no statistical
                                              *  reset, integrals are added for each problem to be solved) */
   );

/** update the primal-dual integral statistic. method accepts + and - SCIPsetInfinity() as values for
 *  upper and lower bound, respectively
 */
extern
void SCIPstatUpdatePrimalDualIntegral(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_Real             primalbound,        /**< current primal bound in transformed problem, or infinity */
   SCIP_Real             dualbound           /**< current lower bound in transformed space, or -infinity */
   );

/** update and return the primal-dual integral statistic */
extern
SCIP_Real SCIPstatGetPrimalDualIntegral(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob            /**< original problem */
   );

/** reset current branch and bound run specific statistics */
extern
void SCIPstatResetCurrentRun(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem, or NULL */
   SCIP_PROB*            origprob,           /**< original problem, or NULL */
   SCIP_Bool             solved              /**< is problem already solved? */
   );

/** resets display statistics, such that a new header line is displayed before the next display line */
extern
void SCIPstatResetDisplay(
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** increases LP count, such that all lazy updates depending on the LP are enforced again */
extern
void SCIPstatEnforceLPUpdates(
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** depending on the current memory usage, switches mode flag to standard or memory saving mode */
extern
void SCIPstatUpdateMemsaveMode(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_MEM*             mem                 /**< block memory pools */
   );

/** returns the estimated number of bytes used by extern software, e.g., the LP solver */
extern
SCIP_Longint SCIPstatGetMemExternEstim(
   SCIP_STAT*            stat                /**< dynamic SCIP statistics */
   );

/** enables or disables all statistic clocks of \p stat concerning LP execution time, strong branching time, etc.
 *
 *  @note: The (pre-)solving time clocks which are relevant for the output during (pre-)solving
 *         are not affected by this method
 *
 *  @see: For completely disabling all timing of SCIP, consider setting the parameter timing/enabled to FALSE
 */
extern
void SCIPstatEnableOrDisableStatClocks(
   SCIP_STAT*            stat,               /**< SCIP statistics */
   SCIP_Bool             enable              /**< should the LP clocks be enabled? */
   );

/** recompute root LP best-estimate from scratch */
extern
void SCIPstatComputeRootLPBestEstimate(
   SCIP_STAT*            stat,               /**< SCIP statistics */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             rootlpobjval,       /**< root LP objective value */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars               /**< number of variables */
   );

/** update root LP best-estimate with changed variable pseudo-costs */
extern
SCIP_RETCODE SCIPstatUpdateVarRootLPBestEstimate(
   SCIP_STAT*            stat,               /**< SCIP statistics */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable with changed pseudo costs */
   SCIP_Real             oldrootpscostscore  /**< old minimum pseudo cost score of variable */
   );

#ifdef TPI_NONE
/* no TPI included so just update the stats */

#define SCIPstatUpdate(stat, set, field, val) do { \
  (stat)->field = (val); \
  } while(0)

#define SCIPstatIncrement(stat, set, field) do { \
   ++(stat)->field; \
   } while(0)

#define SCIPstatAdd(stat, set, field, val) do { \
   (stat)->field += (val); \
   } while(0)

#else
/* TPI not none, so increment deterministic time for relevant stats */

#define SCIPupdateDeterministicTimeCount(stat, set, val) do { \
        (stat)->detertimecnt += (val); \
        if( (stat)->detertimecnt > 10000.0 ) { \
           SCIP_CALL_ABORT( SCIPincrementConcurrentTime( (set)->scip, (stat)->detertimecnt ) ); \
           (stat)->detertimecnt = 0.0;                                  \
        }\
    } while(0) \

#define SCIPstatUpdate(stat, set, field, val) do { \
   switch( offsetof(SCIP_STAT, field) ) \
   { \
      default: \
         break; \
      case offsetof(SCIP_STAT, nprimalresolvelpiterations): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.00328285264101 * ((val) - (stat)->field) * (stat)->nnz ); \
         break; \
      case offsetof(SCIP_STAT, ndualresolvelpiterations): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.00531625104146 * ((val) - (stat)->field) * (stat)->nnz ); \
         break; \
      case offsetof(SCIP_STAT, nprobboundchgs): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.000738719124051 * ((val) - (stat)->field) * (stat)->nnz ); \
         break; \
      case offsetof(SCIP_STAT, nisstoppedcalls): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.0011123144764 * ((val) - (stat)->field) * (stat)->nnz ); \
   } \
   (stat)->field = (val); \
   } while(0)


#define SCIPstatIncrement(stat, set, field) do { \
   switch( offsetof(SCIP_STAT, field) ) \
   { \
      default: \
         break; \
      case offsetof(SCIP_STAT, nprimalresolvelpiterations): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.00328285264101 * (stat)->nnz ); \
         break; \
      case offsetof(SCIP_STAT, ndualresolvelpiterations): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.00531625104146 * (stat)->nnz ); \
         break; \
      case offsetof(SCIP_STAT, nprobboundchgs): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.000738719124051 * (stat)->nnz ); \
         break; \
      case offsetof(SCIP_STAT, nisstoppedcalls): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.0011123144764 * (stat)->nnz ); \
   } \
   ++(stat)->field; \
   } while(0)

#define SCIPstatAdd(stat, set, field, val) do { \
   switch( offsetof(SCIP_STAT, field) ) \
   { \
      default: \
         break; \
      case offsetof(SCIP_STAT, nprimalresolvelpiterations): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.00328285264101 * (val) * (stat)->nnz); \
         break; \
      case offsetof(SCIP_STAT, ndualresolvelpiterations): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.00531625104146 * (val) * (stat)->nnz); \
         break; \
      case offsetof(SCIP_STAT, nprobboundchgs): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.000738719124051 * (val) * (stat)->nnz ); \
         break; \
      case offsetof(SCIP_STAT, nisstoppedcalls): \
         SCIPupdateDeterministicTimeCount(stat, set, 0.0011123144764 * (val) * (stat)->nnz ); \
   } \
   (stat)->field += (val); \
   } while(0)
#endif


/* if we have a C99 compiler */
#ifdef SCIP_HAVE_VARIADIC_MACROS

/** prints a debugging message if SCIP_DEBUG flag is set */
#ifdef SCIP_DEBUG
#define SCIPstatDebugMsg(set, ...)      SCIPstatPrintDebugMessage(stat, __FILE__, __LINE__, __VA_ARGS__)
#define SCIPstatDebugMsgPrint(set, ...) SCIPstatPrintDebugMessagePrint(stat, __VA_ARGS__)
#else
#define SCIPstatDebugMsg(set, ...)      while ( FALSE ) SCIPstatPrintDebugMessage(stat, __FILE__, __LINE__, __VA_ARGS__)
#define SCIPstatDebugMsgPrint(set, ...) while ( FALSE ) SCIPstatPrintDebugMessagePrint(stat, __VA_ARGS__)
#endif

#else
/* if we do not have a C99 compiler, use a workaround that prints a message, but not the file and linenumber */

/** prints a debugging message if SCIP_DEBUG flag is set */
#ifdef SCIP_DEBUG
#define SCIPstatDebugMsg                printf("debug: "), SCIPstatDebugMessagePrint
#define SCIPstatDebugMsgPrint           SCIPstatDebugMessagePrint
#else
#define SCIPstatDebugMsg                while ( FALSE ) SCIPstatDebugMessagePrint
#define SCIPstatDebugMsgPrint           while ( FALSE ) SCIPstatDebugMessagePrint
#endif

#endif


/** prints a debug message */
EXTERN
void SCIPstatPrintDebugMessage(
   SCIP_STAT*            stat,               /**< SCIP statistics */
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline,         /**< line in the source file where the function was called */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a debug message without precode */
EXTERN
void SCIPstatDebugMessagePrint(
   SCIP_STAT*            stat,               /**< SCIP statistics */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

#ifdef __cplusplus
}
#endif

#endif
