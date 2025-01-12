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

/**@file   objprobcloneable.h
 * @brief Definition of base class for all clonable classes which define problem data
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Ambros Gleixner
 * @author Stefan Heinz
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJPROBCLONEABLE_H__
#define __SCIP_OBJPROBCLONEABLE_H__

#include "scip/def.h"
#include "scip/scip.h"
#include "objscip/type_objprobcloneable.h"



namespace scip
{
   /** @brief Definition of base class for all clonable classes which define problem data
    *
    *  Constraint handler and variable pricer C++ wrapper object plugins should extend this class
    */
   struct ObjProbCloneable
   {
      virtual ~ObjProbCloneable() {}

      /** assignment of polymorphic classes causes slicing and is therefore disabled. */
      ObjProbCloneable& operator=(const ObjProbCloneable& o) = delete;

      /** assignment of polymorphic classes causes slicing and is therefore disabled. */
      ObjProbCloneable& operator=(ObjProbCloneable&& o) = delete;

      /** clone method which will be used to copy constraint handler and variable pricer objects */
      virtual SCIP_DECL_OBJPROBCLONE(ObjProbCloneable* clone)
      {
         return 0;
      }

      /** returns whether the plugin object is copyable */
      virtual SCIP_DECL_OBJPROBISCLONEABLE(iscloneable)
      {
         return FALSE;
      }
   };
}

#endif
