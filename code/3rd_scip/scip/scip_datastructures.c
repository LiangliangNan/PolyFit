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

/**@file   scip_datastructures.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for data structures
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/misc.h"
#include "scip/pub_message.h"
#include "scip/scip_datastructures.h"
#include "scip/scip_mem.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"

/** creates a dynamic array of real values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPcreateRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY**      realarray           /**< pointer to store the real array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarrayCreate(realarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of real values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPfreeRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY**      realarray           /**< pointer to the real array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarrayFree(realarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPextendRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarrayExtend(realarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic real array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPclearRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarrayClear(realarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array
 *
 *  @return  value of entry in dynamic array
 */
SCIP_Real SCIPgetRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(scip != NULL);

   return SCIPrealarrayGetVal(realarray, idx);
}

/** sets value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx,                /**< array index to set value for */
   SCIP_Real             val                 /**< value to set array index to */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarraySetVal(realarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, val) );

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPincRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx,                /**< array index to increase value for */
   SCIP_Real             incval              /**< value to increase array index */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarrayIncVal(realarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, incval) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements
 *
 *  @return the minimal index of all stored non-zero elements
 */
int SCIPgetRealarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(scip != NULL);

   return SCIPrealarrayGetMinIdx(realarray);
}

/** returns the maximal index of all stored non-zero elements
 *
 *  @return the maximal index of all stored non-zero elements
 */
int SCIPgetRealarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(scip != NULL);

   return SCIPrealarrayGetMaxIdx(realarray);
}

/** creates a dynamic array of int values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPcreateIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY**       intarray            /**< pointer to store the int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarrayCreate(intarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of int values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPfreeIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY**       intarray            /**< pointer to the int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarrayFree(intarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPextendIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarrayExtend(intarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic int array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPclearIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarrayClear(intarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array
 *
 *  @return value of entry in dynamic array
 */
int SCIPgetIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(scip != NULL);

   return SCIPintarrayGetVal(intarray, idx);
}

/** sets value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx,                /**< array index to set value for */
   int                   val                 /**< value to set array index to */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarraySetVal(intarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, val) );

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPincIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx,                /**< array index to increase value for */
   int                   incval              /**< value to increase array index */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarrayIncVal(intarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, incval) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements
 *
 *  @return the minimal index of all stored non-zero elements
 */
int SCIPgetIntarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(scip != NULL);

   return SCIPintarrayGetMinIdx(intarray);
}

/** returns the maximal index of all stored non-zero elements
 *
 *  @return the maximal index of all stored non-zero elements
 */
int SCIPgetIntarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(scip != NULL);

   return SCIPintarrayGetMaxIdx(intarray);
}

/** creates a dynamic array of bool values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPcreateBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY**      boolarray           /**< pointer to store the bool array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPboolarrayCreate(boolarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of bool values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPfreeBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY**      boolarray           /**< pointer to the bool array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPboolarrayFree(boolarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPextendBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPboolarrayExtend(boolarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic bool array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPclearBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPboolarrayClear(boolarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array
 *
 *  @return value of entry in dynamic array at position idx
 */
SCIP_Bool SCIPgetBoolarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(scip != NULL);

   return SCIPboolarrayGetVal(boolarray, idx);
}

/** sets value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetBoolarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx,                /**< array index to set value for */
   SCIP_Bool             val                 /**< value to set array index to */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPboolarraySetVal(boolarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, val) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements
 *
 *  @return the minimal index of all stored non-zero elements
 */
int SCIPgetBoolarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(scip != NULL);

   return SCIPboolarrayGetMinIdx(boolarray);
}

/** returns the maximal index of all stored non-zero elements
 *
 *  @return the maximal index of all stored non-zero elements
 */
int SCIPgetBoolarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(scip != NULL);

   return SCIPboolarrayGetMaxIdx(boolarray);
}

/** creates a dynamic array of pointers
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPcreatePtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY**       ptrarray            /**< pointer to store the int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPptrarrayCreate(ptrarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of pointers
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPfreePtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY**       ptrarray            /**< pointer to the int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPptrarrayFree(ptrarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPextendPtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPptrarrayExtend(ptrarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic pointer array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPclearPtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPptrarrayClear(ptrarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
void* SCIPgetPtrarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(scip != NULL);

   return SCIPptrarrayGetVal(ptrarray, idx);
}

/** sets value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetPtrarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   idx,                /**< array index to set value for */
   void*                 val                 /**< value to set array index to */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPptrarraySetVal(ptrarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, val) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements
 *
 *  @return the minimal index of all stored non-zero elements
 */
int SCIPgetPtrarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(scip != NULL);

   return SCIPptrarrayGetMinIdx(ptrarray);
}

/** returns the maximal index of all stored non-zero elements
 *
 *  @return the maximal index of all stored non-zero elements
 */
int SCIPgetPtrarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(scip != NULL);

   return SCIPptrarrayGetMaxIdx(ptrarray);
}

/** creates directed graph structure */
SCIP_RETCODE SCIPcreateDigraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH**        digraph,            /**< pointer to store the created directed graph */
   int                   nnodes              /**< number of nodes */
   )
{
   assert(scip != NULL);
   assert(digraph != NULL);

   SCIP_CALL( SCIPdigraphCreate(digraph, scip->mem->probmem, nnodes) );

   return SCIP_OKAY;
}

/** copies directed graph structure
 *
 *  The copying procedure uses the memory of the passed SCIP instance. The user must ensure that the digraph lives
 *  as most as long as the SCIP instance.
 *
 *  @note The data in nodedata is copied verbatim. This possibly has to be adapted by the user.
 */
SCIP_RETCODE SCIPcopyDigraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH**        targetdigraph,      /**< pointer to store the copied directed graph */
   SCIP_DIGRAPH*         sourcedigraph       /**< source directed graph */
   )
{
   assert(scip != NULL);
   assert(sourcedigraph != NULL);
   assert(targetdigraph != NULL);

   SCIP_CALL( SCIPdigraphCopy(targetdigraph, sourcedigraph, scip->mem->probmem) );

   return SCIP_OKAY;
}

/** creates a disjoint set (union find) structure \p uf for \p ncomponents many components (of size one) */
SCIP_RETCODE SCIPcreateDisjointset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DISJOINTSET**    djset,              /**< disjoint set (union find) data structure */
   int                   ncomponents         /**< number of components */
   )
{
   assert(scip != NULL);
   assert(djset != NULL);

   SCIP_CALL( SCIPdisjointsetCreate(djset, scip->mem->probmem, ncomponents) );

   return SCIP_OKAY;
}

/** frees the disjoint set (union find) data structure */
void SCIPfreeDisjointset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DISJOINTSET**    djset               /**< disjoint set (union find) data structure */
   )
{
   assert(scip != NULL);

   SCIPdisjointsetFree(djset, scip->mem->probmem);
}
