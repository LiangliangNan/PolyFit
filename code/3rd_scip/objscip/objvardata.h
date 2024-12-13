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

/**@file   objvardata.h
 * @brief  C++ wrapper for user variable data
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJVARDATA_H__
#define __SCIP_OBJVARDATA_H__


#include <cassert>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for user variable data
 *
 *  This class defines the interface for user variable data implemented in C++. Each variable can be equipped with a
 *  variable data class. This data can be accessed via the function SCIPgetObjVardata() at any time after it is created
 *  and before it is deleted.
 *
 *  - \ref type_var.h "Corresponding C interface"
 */
class ObjVardata
{
public:
   /** default constructor */
   ObjVardata()
   {
   }

   /** destructor */
   virtual ~ObjVardata()
   {
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjVardata& operator=(const ObjVardata& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjVardata& operator=(ObjVardata&& o) = delete;

   /** destructor of user variable data to free original user data (called when original variable is freed)
    *
    *  If the "deleteobject" flag in the SCIPcreateObjVar() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user variable data can be done in the destructor of the user variable
    *  data object. If the "deleteobject" flag was set to FALSE, and the user variable data object stays alive
    *  after the SCIP variable is freed, this method should delete all the variable specific data that is no
    *  longer needed.
    */ /*lint -e715*/
   virtual SCIP_RETCODE scip_delorig(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_VAR*          var                 /**< original variable, the data to free is belonging to */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** creates user data of transformed variable by transforming the original user variable data
    *  (called after variable was transformed)
    *
    *  The user has two possibilities to implement this method:
    *   1. Return the pointer to the original variable data object (this) as pointer to the transformed variable data
    *      object. The user may modify some internal attributes, but he has to make sure, that these modifications are
    *      reversed in the scip_deltrans() method, such that the original variable data is restored. In this case,
    *      he should set *deleteobject to FALSE, because the variable data must not be destructed by SCIP after the
    *      solving process is terminated.
    *   2. Call the copy constructor of the variable data object and return the created copy as transformed variable
    *      data object. In this case, he probably wants to set *deleteobject to TRUE, thus letting SCIP call the
    *      destructor of the object if the transformed variable data is no longer needed.
    */ /*lint -e715*/
   virtual SCIP_RETCODE scip_trans(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_VAR*          var,                /**< transformed variable, the data to create is belonging to */
      ObjVardata**       objvardata,         /**< pointer to store the transformed variable data object */
      SCIP_Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
      )
   {  /*lint --e{715}*/
      assert(objvardata != NULL);
      assert(deleteobject != NULL);

      /* the default implementation just copies the pointer to the variable data object;
       * SCIP must not destruct the transformed variable data object, because the original variable data must stay alive
       */
      *objvardata = this;
      *deleteobject = FALSE;

      return SCIP_OKAY;
   }

   /** destructor of user variable data to free transformed user data (called when transformed variable is freed)
    *
    *  If the "*deleteobject" flag in the scip_trans() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user variable data can be done in the destructor of the user variable
    *  data object. If the "*deleteobject" flag was set to FALSE, and the user variable data object stays alive
    *  after the SCIP variable is freed, this method should delete all the variable specific data that is no
    *  longer needed.
    */ /*lint -e715*/
   virtual SCIP_RETCODE scip_deltrans(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_VAR*          var                 /**< transformed variable, the data to free is belonging to */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** copies variable data of source SCIP variable for the target SCIP variable
    *
    *  This method should copy the variable data of the source SCIP and create a target variable data for target
    *  variable. This callback is optional. If the copying process was successful, the target variable gets this variable
    *  data assigned. In case the result pointer is set to SCIP_DIDNOTRUN, the target variable will have no variable data at
    *  all.
    *
    *  The variable map and the constraint map can be used via the function SCIPgetVarCopy() and SCIPgetConsCopy(),
    *  respectively, to get for certain variables and constraints of the source SCIP the counter parts in the target
    *  SCIP. You should be very carefully in using these two methods since they could lead to infinite loop.
    *
    *  possible return values for *result:
    *  - SCIP_DIDNOTRUN  : the copying process was not performed
    *  - SCIP_SUCCESS    : the copying process was successfully performed
    */ /*lint -e715*/
   virtual SCIP_RETCODE scip_copy(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP*              sourcescip,         /**< source SCIP main data structure */
      SCIP_VAR*          sourcevar,          /**< variable of the source SCIP */
      SCIP_HASHMAP*      varmap,             /**< a hashmap which stores the mapping of source variables to corresponding
                                              *   target variables */
      SCIP_HASHMAP*      consmap,            /**< a hashmap which stores the mapping of source contraints to corresponding 
                                              *   target constraints */
      SCIP_VAR*          targetvar,          /**< variable of the (targert) SCIP (targetvar is the copy of sourcevar) */
      ObjVardata**       objvardata,         /**< pointer to store the copied variable data object */
      SCIP_RESULT*       result              /**< pointer to store the result of the call */
      )
   {  /*lint --e{715}*/
      (*objvardata) = 0;
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
};

} /* namespace scip */



/** create and capture problem variable and associates the given variable data with the variable;
 *  if variable is of integral type, fractional bounds are automatically rounded
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateObjVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   scip::ObjVardata*     objvardata,         /**< user variable data object */
   SCIP_Bool             deleteobject        /**< should the user variable data object be deleted when variable is freed? */
   );

/** gets user variable data object for given problem variable
 *  Warning! This method should only be called after a variable was created with SCIPcreateObjVar().
 *  Otherwise, a segmentation fault may arise, or an undefined pointer is returned.
 */
SCIP_EXPORT
scip::ObjVardata* SCIPgetObjVardata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

#endif
