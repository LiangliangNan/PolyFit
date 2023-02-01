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

/**@file   type_sepastore.h
 * @brief  type definitions for storing separated cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SEPASTORE_H__
#define __SCIP_TYPE_SEPASTORE_H__

#ifdef __cplusplus
extern "C" {
#endif

/** possible settings for specifying the solution for which cuts are selected */
enum SCIP_Efficiacychoice
{
   SCIP_EFFICIACYCHOICE_LP    = 0,           /**< use LP solution to base efficacy on */
   SCIP_EFFICIACYCHOICE_RELAX = 1,           /**< use relaxation solution to base efficacy on */
   SCIP_EFFICIACYCHOICE_NLP   = 2            /**< use NLP solution to base efficacy on */
};
typedef enum SCIP_Efficiacychoice SCIP_EFFICIACYCHOICE;

typedef struct SCIP_SepaStore SCIP_SEPASTORE;     /**< storage for separated variables */

#ifdef __cplusplus
}
#endif

#endif
