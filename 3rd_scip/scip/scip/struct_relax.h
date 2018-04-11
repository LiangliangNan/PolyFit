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

/**@file   struct_relax.h
 * @ingroup INTERNALAPI
 * @brief  data structures for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_RELAX_H__
#define __SCIP_STRUCT_RELAX_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_relax.h"

#ifdef __cplusplus
extern "C" {
#endif

/** relaxators data */
struct SCIP_Relax
{
   SCIP_Longint          ncalls;             /**< number of times, this relaxator was called */
   SCIP_Longint          lastsolvednode;     /**< last total nodes counter, where the current relaxation was solved */
   char*                 name;               /**< name of relaxator */
   char*                 desc;               /**< description of relaxator */
   SCIP_DECL_RELAXCOPY   ((*relaxcopy));     /**< copy method of relaxator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_RELAXFREE   ((*relaxfree));     /**< destructor of relaxator */
   SCIP_DECL_RELAXINIT   ((*relaxinit));     /**< initialize relaxator */
   SCIP_DECL_RELAXEXIT   ((*relaxexit));     /**< deinitialize relaxator */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol));  /**< solving process initialization method of relaxator */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol));  /**< solving process deinitialization method of relaxator */
   SCIP_DECL_RELAXEXEC   ((*relaxexec));     /**< execution method of relaxator */
   SCIP_RELAXDATA*       relaxdata;          /**< relaxators local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this relaxator for the next stages */
   SCIP_CLOCK*           relaxclock;         /**< relaxation time */
   int                   priority;           /**< priority of the relaxator */
   int                   freq;               /**< frequency for calling relaxator */
   SCIP_Bool             initialized;        /**< is relaxator initialized? */
};

/** relaxation information data */
struct SCIP_Relaxation
{
   SCIP_Real             relaxsolobjval;
   SCIP_Bool             relaxsolvalid;
   SCIP_Bool             relaxsolincludeslp;
   SCIP_Bool             relaxsolzero;
};


#ifdef __cplusplus
}
#endif

#endif
