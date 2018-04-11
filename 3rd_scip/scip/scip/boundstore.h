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

/**@file   boundstore.h
 * @ingroup PARALLEL
 * @brief  the interface of the boundstore structure
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/def.h"
#include "scip/type_syncstore.h"
#include "scip/type_scip.h"
#include "scip/type_lp.h"
#include "scip/type_retcode.h"

#ifndef __SCIP_BOUNDSTORE_H__
#define __SCIP_BOUNDSTORE_H__

/** create bound store data structure */
EXTERN
SCIP_RETCODE SCIPboundstoreCreate(
   SCIP*                 scip,               /**< scip main datastructure */
   SCIP_BOUNDSTORE**     boundstore,         /**< pointer to store the bound store datastructure */
   int                   nvars               /**< number of variables for which bounds may be stored */
   );

/** free bound store data structure */
EXTERN
void SCIPboundstoreFree(
   SCIP*                 scip,               /**< scip main datastructure */
   SCIP_BOUNDSTORE**     boundstore          /**< pointer to the bound store datastructure */
   );

/** add bound change to bound store data structure */
EXTERN
SCIP_RETCODE SCIPboundstoreAdd(
   SCIP*                 scip,               /**< scip main datastructure */
   SCIP_BOUNDSTORE*      boundstore,         /**< the bound store datastructure */
   int                   varidx,             /**< variable index of bound change, must be smaller than the
                                              *   number of variables given during creation of bound store */
   SCIP_Real             newbound,           /**< bound value of variable */
   SCIP_BOUNDTYPE        boundtype           /**< type of new bound */
   );

/** add all bound changes of source to target */
EXTERN
SCIP_RETCODE SCIPboundstoreMerge(
   SCIP*                 scip,               /**< scip main datastructure for target boundstore   */
   SCIP_BOUNDSTORE*      target,             /**< the bound store datastructure where the bounds get merged in */
   SCIP_BOUNDSTORE*      source              /**< the bound store datastructure from which the bounds get merged in */
   );

/** remove all boundchanges from bound store */
EXTERN
void SCIPboundstoreClear(
   SCIP_BOUNDSTORE*      boundstore          /**< the bound store datastructure */
   );

/** gets variable index of the i'th stored boundchange */
EXTERN
int SCIPboundstoreGetChgVaridx(
   SCIP_BOUNDSTORE*      boundstore,         /**< the bound store datastructure */
   int                   i                   /**< the index of the bound change */
   );

/** gets the type of the i'th stored boundchange */
EXTERN
SCIP_BOUNDTYPE SCIPboundstoreGetChgType(
   SCIP_BOUNDSTORE*      boundstore,         /**< the bound store datastructure */
   int                   i                   /**< the index of the bound change */
   );

/** gets the bound value of the i'th stored boundchange */
EXTERN
SCIP_Real SCIPboundstoreGetChgVal(
   SCIP_BOUNDSTORE*      boundstore,         /**< the bound store datastructure */
   int                   i                   /**< the index of the bound change */
   );

/** gets the number of stored bound changes */
EXTERN
int SCIPboundstoreGetNChgs(
   SCIP_BOUNDSTORE*      boundstore          /**< the bound store datastructure */
   );

#endif
