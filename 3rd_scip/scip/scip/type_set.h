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

/**@file   type_set.h
 * @brief  type definitions for global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SET_H__
#define __SCIP_TYPE_SET_H__

/**! [SnippetCodeStyleExample] */

#ifdef __cplusplus
extern "C" {
#endif

/** SCIP operation stage */
enum SCIP_Stage
{
   SCIP_STAGE_INIT         =  0,        /**< SCIP data structures are initialized, no problem exists */
   SCIP_STAGE_PROBLEM      =  1,        /**< the problem is being created and modified */
   SCIP_STAGE_TRANSFORMING =  2,        /**< the problem is being transformed into solving data space */
   SCIP_STAGE_TRANSFORMED  =  3,        /**< the problem was transformed into solving data space */
   SCIP_STAGE_INITPRESOLVE =  4,        /**< presolving is initialized */
   SCIP_STAGE_PRESOLVING   =  5,        /**< the problem is being presolved */
   SCIP_STAGE_EXITPRESOLVE =  6,        /**< presolving is exited */
   SCIP_STAGE_PRESOLVED    =  7,        /**< the problem was presolved */
   SCIP_STAGE_INITSOLVE    =  8,        /**< the solving process data is being initialized */
   SCIP_STAGE_SOLVING      =  9,        /**< the problem is being solved */
   SCIP_STAGE_SOLVED       = 10,        /**< the problem was solved */
   SCIP_STAGE_EXITSOLVE    = 11,        /**< the solving process data is being freed */
   SCIP_STAGE_FREETRANS    = 12,        /**< the transformed problem is being freed */
   SCIP_STAGE_FREE         = 13         /**< SCIP data structures are being freed */
};
typedef enum SCIP_Stage SCIP_STAGE;

/** possible settings for enabling/disabling algorithms and other features */
enum SCIP_Setting
{
   SCIP_UNDEFINED = 0,                  /**< undefined setting */
   SCIP_DISABLED  = 1,                  /**< feature is disabled */
   SCIP_AUTO      = 2,                  /**< feature is set to automatic mode */
   SCIP_ENABLED   = 3                   /**< feature is enabled */
};
typedef enum SCIP_Setting SCIP_SETTING;

typedef struct SCIP_Set SCIP_SET;                 /**< global SCIP settings */

#ifdef __cplusplus
}
#endif

/**! [SnippetCodeStyleExample] */

#endif
