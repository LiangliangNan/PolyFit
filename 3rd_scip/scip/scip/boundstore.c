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

/**@file   boundstore.c
 * @ingroup PARALLEL
 * @brief  the implementation of the bound store datastructure
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/boundstore.h"
#include "scip/struct_syncstore.h"
#include "blockmemshell/memory.h"
#include "scip/scip.h"

/** create bound store data structure */
SCIP_RETCODE SCIPboundstoreCreate(
   SCIP*                 scip,               /**< scip main datastructure */
   SCIP_BOUNDSTORE**     boundstore,         /**< pointer to store the bound store datastructure */
   int                   nvars               /**< number of variables for which bounds may be stored */
   )
{
   assert(scip != NULL);
   assert(boundstore != NULL);

   SCIP_CALL( SCIPallocMemory(scip, boundstore) );

   (*boundstore)->bndchg = NULL;
   (*boundstore)->bndchgsize = 0;
   (*boundstore)->nbndchg = 0;
   (*boundstore)->nvars = nvars;
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(*boundstore)->bndpos, nvars) );

   return SCIP_OKAY;
}

/** free bound store data structure */
void SCIPboundstoreFree(
   SCIP*                 scip,               /**< scip main datastructure */
   SCIP_BOUNDSTORE**     boundstore          /**< pointer to the bound store datastructure */
   )
{
   assert(scip != NULL);
   assert(boundstore != NULL);
   assert((*boundstore) != NULL);

   SCIPfreeBlockMemoryArray(scip, &(*boundstore)->bndpos, (*boundstore)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*boundstore)->bndchg, (*boundstore)->bndchgsize);
   SCIPfreeMemory(scip, boundstore);
}

/** add bound change to bound store data structure */
SCIP_RETCODE SCIPboundstoreAdd(
   SCIP*                 scip,               /**< scip main datastructure */
   SCIP_BOUNDSTORE*      boundstore,         /**< the bound store datastructure */
   int                   varidx,             /**< variable index of bound change, must be smaller than the
                                              *   number of variables given during creation of bound store */
   SCIP_Real             newbound,           /**< bound value of variable */
   SCIP_BOUNDTYPE        boundtype           /**< type of new bound */
   )
{
   /* check if already stored a bound of same type for this variable */
   int pos;

   assert(scip != NULL);
   assert(boundstore != NULL);

   pos = boundstore->bndpos[varidx].pos[boundtype];

   if( pos == 0 )
   {
      /* variable has no bound stored yet so store this bound */
      int i;
      i = boundstore->nbndchg++;
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, &boundstore->bndchg, &boundstore->bndchgsize, boundstore->nbndchg) );
      boundstore->bndchg[i].varidx = varidx;
      boundstore->bndchg[i].newbound = newbound;
      boundstore->bndchg[i].boundtype = boundtype;
      boundstore->bndpos[varidx].pos[boundtype] = boundstore->nbndchg;
   }
   else
   {
      /* since pos == 0 is reserved if no bound is stored
       * the index is pos-1
       */
      --pos;
      assert(boundstore->bndchg[pos].boundtype == boundtype);
      assert(boundstore->bndchg[pos].varidx == varidx);

      /* if the bound is better overwrite the old one */
      if( (boundtype == SCIP_BOUNDTYPE_LOWER && newbound > boundstore->bndchg[pos].newbound) ||
          (boundtype == SCIP_BOUNDTYPE_UPPER && newbound < boundstore->bndchg[pos].newbound) )
      {
         boundstore->bndchg[pos].newbound = newbound;
      }
   }

   return SCIP_OKAY;
}

/** add all bound changes of source to target */
SCIP_RETCODE SCIPboundstoreMerge(
   SCIP*                 scip,               /**< scip main datastructure for target boundstore */
   SCIP_BOUNDSTORE*      target,             /**< the bound store datastructure where the bounds get merged in */
   SCIP_BOUNDSTORE*      source              /**< the bound store datastructure from which the bounds get merged in */
   )
{
   int i;

   assert(scip != NULL);
   assert(source != NULL);
   assert(target != NULL);

   /* just iterate over the boundchanges in the source and add them to the target */
   for( i = 0; i < source->nbndchg; ++i )
   {
      SCIP_CALL( SCIPboundstoreAdd(scip, target, source->bndchg[i].varidx, source->bndchg[i].newbound, source->bndchg[i].boundtype) );
   }

   return SCIP_OKAY;
}

/** remove all boundchanges from bound store */
void SCIPboundstoreClear(
   SCIP_BOUNDSTORE*      boundstore          /**< the bound store datastructure */
   )
{
   assert(boundstore != NULL);

   /* clearing the positions is enough */
   if( boundstore->nbndchg > 0 )
   {
      BMSclearMemoryArray(boundstore->bndpos, boundstore->nvars);
      boundstore->nbndchg = 0;
   }
}

/** gets variable index of the i'th stored boundchange */
int SCIPboundstoreGetChgVaridx(
   SCIP_BOUNDSTORE*      boundstore,         /**< the bound store datastructure */
   int                   i                   /**< the index of the bound change */
   )
{
   assert(boundstore != NULL);
   assert(i < boundstore->nbndchg);

   return boundstore->bndchg[i].varidx;
}

/** gets the type of the i'th stored boundchange */
SCIP_BOUNDTYPE SCIPboundstoreGetChgType(
   SCIP_BOUNDSTORE*      boundstore,         /**< the bound store datastructure */
   int                   i                   /**< the index of the bound change */
   )
{
   assert(boundstore != NULL);
   assert(i < boundstore->nbndchg);

   return boundstore->bndchg[i].boundtype;
}

/** gets the bound value of the i'th stored boundchange */
SCIP_Real SCIPboundstoreGetChgVal(
   SCIP_BOUNDSTORE*      boundstore,         /**< the bound store datastructure */
   int                   i                   /**< the index of the bound change */
   )
{
   assert(boundstore != NULL);
   assert(i < boundstore->nbndchg);

   return boundstore->bndchg[i].newbound;
}

/** gets the number of stored bound changes */
int SCIPboundstoreGetNChgs(
   SCIP_BOUNDSTORE*      boundstore          /**< the bound store datastructure */
   )
{
   assert(boundstore != NULL);

   return boundstore->nbndchg;
}
