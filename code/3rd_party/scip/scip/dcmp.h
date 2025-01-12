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

/**@file   dcmp.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for decompositions and the decomposition store
 * @author Gregor Hendel
 *
 * @todo get a decomposition score, and compute other stuff that may be important
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef SRC_SCIP_DECOMP_H_
#define SRC_SCIP_DECOMP_H_

#include "scip/type_dcmp.h"
#include "scip/type_var.h"
#include "scip/type_cons.h"
#include "blockmemshell/memory.h"


#ifdef __cplusplus
extern "C" {
#endif

#define SCIP_DECOMPSTORE_CAPA 10             /**< hardcoded maximum capacity of decomposition store */

/** creates a decomposition storage */
SCIP_RETCODE SCIPdecompstoreCreate(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to store decomposition storage */
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   int                   nslots              /**< maximum number of decomposition slots in storage */
   );

/** frees a decomposition storage */
void SCIPdecompstoreFree(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to free decomposition storage */
   BMS_BLKMEM*           blkmem              /**< block memory data structure */
   );

/** adds decomposition to storage */
SCIP_RETCODE SCIPdecompstoreAdd(
   SCIP_DECOMPSTORE*     decompstore,        /**< decomposition storage */
   SCIP_DECOMP*          decomp              /**< decomposition to add */
   );

/** transforms all available original decompositions into transformed space */
SCIP_RETCODE SCIPtransformDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** frees all decompositions in transformed space */
void SCIPexitSolveDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets decompositions from storage */
SCIP_DECOMP** SCIPdecompstoreGetDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** gets number of decompositions in storage */
int SCIPdecompstoreGetNDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** gets decompositions in original space from storage */
SCIP_DECOMP** SCIPdecompstoreGetOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** gets number of decompositions in original space in storage */
int SCIPdecompstoreGetNOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

#ifdef __cplusplus
}
#endif

#endif
