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

/**@file   type_compr.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for tree compression
 * @author Jakob Witzig
 *
 *  This file defines the interface for tree compression implemented in C.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_COMPR_H__
#define __SCIP_TYPE_COMPR_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"
#include "scip/type_timing.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Compr SCIP_COMPR;              /**< tree compression */
typedef struct SCIP_ComprData SCIP_COMPRDATA;      /**< locally defined tree compression data */


/** copy method for compression plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - compr           : the compression technique itself
 */
#define SCIP_DECL_COMPRCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_COMPR* compr)

/** destructor of tree compression to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - compr           : the compression technique itself
 */
#define SCIP_DECL_COMPRFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_COMPR* compr)

/** initialization method of tree compression (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - compr           : the compression technique itself
 */
#define SCIP_DECL_COMPRINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_COMPR* compr)

/** deinitialization method of tree compression (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - compr           : the compression technique itself
 */
#define SCIP_DECL_COMPREXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_COMPR* compr)

/** solving process initialization method of tree compressionc (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The tree compression may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - compr           : the compression technique itself
 */
#define SCIP_DECL_COMPRINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_COMPR* compr)

/** solving process deinitialization method of tree compression (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The tree compression should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - compr           : the compression technique itself
 */
#define SCIP_DECL_COMPREXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_COMPR* compr)

/** execution method of tree compression technique
 *
 *  Try to compress the current search tree. The method is called in the node processing loop.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - compr           : the compression technique itself
 *  - result          : pointer to store the result of the heuristic call
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the tree could be compressed
 *  - SCIP_DIDNITFIND : the method could not compress the tree
 *  - SCIP_DIDNOTRUN  : the compression was skipped
 */
#define SCIP_DECL_COMPREXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_COMPR* compr, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
