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

/**@file   objcloneable.h
 * @brief  definition of base class for all clonable classes
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJCLONEABLE_H__
#define __SCIP_OBJCLONEABLE_H__

#include "scip/def.h"
#include "scip/scip.h"
#include "objscip/type_objcloneable.h"


namespace scip
{
   /** @brief Definition of base class for all clonable classes
    *
    * All C++ wrapper object plugins should extend this class, except constraint handlers and variable pricers. This is
    * needed to be able to copy (clone) a SCIP instance.
    */
   struct ObjCloneable
   {
      virtual ~ObjCloneable() {}

      /** assignment of polymorphic classes causes slicing and is therefore disabled. */
      ObjCloneable& operator=(const ObjCloneable& o) = delete;

      /** assignment of polymorphic classes causes slicing and is therefore disabled. */
      ObjCloneable& operator=(ObjCloneable&& o) = delete;

      /** clone method, used to copy plugins which are not constraint handlers or variable pricer plugins */
      virtual SCIP_DECL_OBJCLONEABLECLONE(ObjCloneable* clone)
      {
         return 0;
      }

      /** returns whether the plugin object is copyable */
      virtual SCIP_DECL_OBJCLONEABLEISCLONEABLE(iscloneable)
      {
         return false;
      }
   };
}

#endif
