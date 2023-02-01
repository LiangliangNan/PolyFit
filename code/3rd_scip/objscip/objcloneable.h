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
