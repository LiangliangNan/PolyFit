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

/**@file   struct_compr.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for tree compression techniques
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_COMPR_H__
#define __SCIP_STRUCT_COMPR_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_compr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** tree compression data */
struct SCIP_Compr
{
   SCIP_Longint          ncalls;             /**< number of times, this compression was called */
   SCIP_Longint          nfound;              /**< number of compressions found so far by this method */
   SCIP_Real             rate;               /**< rate of the last compression */
   SCIP_Real             loi;                /**< loss of information of the last compression */
   char*                 name;               /**< name of tree compression */
   char*                 desc;               /**< description of tree compression */
   SCIP_DECL_COMPRCOPY   ((*comprcopy));     /**< copy method of tree compression or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_COMPRFREE   ((*comprfree));     /**< destructor of tree compression */
   SCIP_DECL_COMPRINIT   ((*comprinit));     /**< initialize tree compression */
   SCIP_DECL_COMPREXIT   ((*comprexit));     /**< deinitialize tree compression */
   SCIP_DECL_COMPRINITSOL ((*comprinitsol)); /**< solving process initialization method of tree compression */
   SCIP_DECL_COMPREXITSOL ((*comprexitsol)); /**< solving process deinitialization method of tree compression */
   SCIP_DECL_COMPREXEC   ((*comprexec));     /**< execution method of tree compression */
   SCIP_COMPRDATA*       comprdata;          /**< tree compression local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this compression for the next stages */
   SCIP_CLOCK*           comprclock;         /**< compression execution time */
   int                   priority;           /**< priority of the tree compression */
   int                   minnnodes;          /**< minimal number of nodes for calling compression, -1 if no threshold exists */
   int                   nnodes;             /**< number of nodes of the last compression */
   SCIP_Bool             initialized;        /**< is tree compression initialized? */
};

#ifdef __cplusplus
}
#endif

#endif
