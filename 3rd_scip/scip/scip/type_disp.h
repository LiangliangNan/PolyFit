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

/**@file   type_disp.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for displaying runtime statistics
 * @author Tobias Achterberg
 *
 *  This file defines the interface for display columns implemented in C.
 *
 * - \ref DISP "Instructions for implementing a display column"
 * - \ref DISPLAYS "List of available display columns"
 * - \ref scip::ObjDisp "C++ wrapper class
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_DISP_H__
#define __SCIP_TYPE_DISP_H__

#include <stdio.h>

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** display activation status of display column */
enum SCIP_DispStatus
{
   SCIP_DISPSTATUS_OFF  = 0,            /**< display column is not displayed */
   SCIP_DISPSTATUS_AUTO = 1,            /**< display column is switched on and off automatically */
   SCIP_DISPSTATUS_ON   = 2             /**< display column is displayed */
};
typedef enum SCIP_DispStatus SCIP_DISPSTATUS;

/** display activation status of display column */
enum SCIP_DispMode
{
   SCIP_DISPMODE_DEFAULT    = 0x00000001u,        /**< display column is displayed only in sequential mode */
   SCIP_DISPMODE_CONCURRENT = 0x00000002u,        /**< display column is displayed only in concurrent mode */
   SCIP_DISPMODE_ALL        = 0x00000003u         /**< display column is displayed in concurrent and sequential mode*/
};
typedef enum SCIP_DispMode SCIP_DISPMODE;

typedef struct SCIP_Disp SCIP_DISP;               /**< display column data structure */
typedef struct SCIP_DispData SCIP_DISPDATA;       /**< display column specific data */


/**  copy method for display plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 */
#define SCIP_DECL_DISPCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_DISP* disp)

/** destructor of display column to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 */
#define SCIP_DECL_DISPFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_DISP* disp)

/** initialization method of display column (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 */
#define SCIP_DECL_DISPINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_DISP* disp)

/** deinitialization method of display column (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 */
#define SCIP_DECL_DISPEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_DISP* disp)

/** solving process initialization method of display column (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The display column may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 */
#define SCIP_DECL_DISPINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_DISP* disp)

/** solving process deinitialization method of display column (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The display column should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 */
#define SCIP_DECL_DISPEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_DISP* disp)

/** output method of display column to output file stream 'file'
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 *  - file            : file stream for output
 */
#define SCIP_DECL_DISPOUTPUT(x) SCIP_RETCODE x (SCIP* scip, SCIP_DISP* disp, FILE* file)

#ifdef __cplusplus
}
#endif

#endif
