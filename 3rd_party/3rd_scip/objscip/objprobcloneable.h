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
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
