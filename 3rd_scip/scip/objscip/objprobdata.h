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

/**@file   objprobdata.h
 * @brief  C++ wrapper for user problem data
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJPROBDATA_H__
#define __SCIP_OBJPROBDATA_H__


#include <cassert>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for user problem data
 *
 *  This class defines the interface for user problem data implemented in C++. This class can be accessed at any time
 *  using the methods SCIPgetObjProbData(). Therefore, it can be used to store data which has to be accessible within
 *  several plugins.
 *
 *  - \ref type_prob.h "Corresponding C interface"
 */
class ObjProbData
{
public:
   /** default constructor */
   ObjProbData()
   {
   }

   /** destructor */
   virtual ~ObjProbData()
   {
   }

   /** destructor of user problem data to free original user data (called when original problem is freed)
    *
    *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual SCIP_RETCODE scip_delorig(
      SCIP*              scip                /**< SCIP data structure */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** creates user data of transformed problem by transforming the original user problem data
    *  (called after problem was transformed)
    *
    *  The user has two possibilities to implement this method:
    *   1. Return the pointer to the original problem data object (this) as pointer to the transformed problem data
    *      object. The user may modify some internal attributes, but he has to make sure, that these modifications are
    *      reversed in the scip_deltrans() method, such that the original problem data is restored. In this case,
    *      he should set *deleteobject to FALSE, because the problem data must not be destructed by SCIP after the
    *      solving process is terminated.
    *   2. Call the copy constructor of the problem data object and return the created copy as transformed problem
    *      data object. In this case, he probably wants to set *deleteobject to TRUE, thus letting SCIP call the
    *      destructor of the object if the transformed problem data is no longer needed.
    */
   virtual SCIP_RETCODE scip_trans(
      SCIP*              scip,               /**< SCIP data structure */
      ObjProbData**      objprobdata,        /**< pointer to store the transformed problem data object */
      SCIP_Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
      )
   {  /*lint --e{715}*/
      assert(objprobdata != NULL);
      assert(deleteobject != NULL);

      /* the default implementation just copies the pointer to the problem data object;
       * SCIP must not destruct the transformed problem data object, because the original problem data must stay alive
       */
      *objprobdata = this;
      *deleteobject = FALSE;

      return SCIP_OKAY;
   }      

   /** destructor of user problem data to free transformed user data (called when transformed problem is freed)
    *
    *  If the "*deleteobject" flag in the scip_trans() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "*deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual SCIP_RETCODE scip_deltrans(
      SCIP*              scip                /**< SCIP data structure */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of transformed data (called before the branch and bound process begins)
    *
    *  This method is called before the branch and bound process begins and can be used to initialize user problem
    *  data that depends for example on the number of active problem variables, because these are now fixed.
    */
   virtual SCIP_RETCODE scip_initsol(
      SCIP*              scip                /**< SCIP data structure */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of transformed data (called before the branch and bound data is freed)
    *
    *  This method is called before the branch and bound data is freed and should be used to free all data that
    *  was allocated in the solving process initialization method. The user has to make sure, that all LP rows associated
    *  to the transformed user problem data are released.
    */
   virtual SCIP_RETCODE scip_exitsol(
      SCIP*              scip,                /**< SCIP data structure */
      SCIP_Bool          restart              /**< was this exit solve call triggered by a restart? */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** copies user data of source SCIP for the target SCIP
    *
    *  This method should copy the problem data of the source SCIP and create a target problem data for (target)
    *  SCIP. Implementing this callback is optional. If the copying process was successful the target SCIP gets this
    *  problem data assigned. In case the result pointer is set to SCIP_DIDNOTRUN the target SCIP will have no problem data
    *  at all.
    *
    *  The variable map and the constraint map can be used via the function SCIPgetVarCopy() and SCIPgetConsCopy(),
    *  respectively, to get for certain variables and constraints of the source SCIP the counter parts in the target
    *  SCIP. You should be very carefully in using these two methods since they could lead to an infinite loop due to
    *  recursion.
    *
    *  possible return values for *result:
    *  - SCIP_DIDNOTRUN  : the copying process was not performed 
    *  - SCIP_SUCCESS    : the copying process was successfully performed
    */
   virtual SCIP_RETCODE scip_copy(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP*              sourcescip,         /**< source SCIP main data structure */
      SCIP_HASHMAP*      varmap,             /**< a hashmap which stores the mapping of source variables to corresponding 
                                              *   target variables */
      SCIP_HASHMAP*      consmap,            /**< a hashmap which stores the mapping of source contraints to corresponding 
                                              *   target constraints */
      ObjProbData**      objprobdata,        /**< pointer to store the copied problem data object */
      SCIP_Bool          global,             /**< create a global or a local copy? */
      SCIP_RESULT*       result              /**< pointer to store the result of the call */
      )
   {  /*lint --e{715}*/
      (*objprobdata) = 0;
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
};

} /* namespace scip */



/** creates empty problem, initializes all solving data structures, and sets the user problem data to point to the
 *  given user data object
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyProbData* myprobdata = new MyProbData(...);
 *       SCIP_CALL( SCIPcreateObjProb(scip, "probname", &myprobdata, FALSE) );
 *       ... // solve the problem
 *       SCIP_CALL( SCIPfreeProb(scip) );
 *       delete myprobdata;    // delete probdata AFTER SCIPfreeProb() !
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfreeProb() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPcreateObjProb(scip, "probname", new MyProbData(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // problem is freed and destructor of MyProbData is called here
 */
EXTERN
SCIP_RETCODE SCIPcreateObjProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */
   scip::ObjProbData*    objprobdata,        /**< user problem data object */
   SCIP_Bool             deleteobject        /**< should the user problem data object be deleted when problem is freed? */
   );

/** gets user problem data object
 *  Warning! This method should only be called after a problem was created with SCIPcreateObjProb().
 *  Otherwise, a segmentation fault may arise, or an undefined pointer is returned.
 */
EXTERN
scip::ObjProbData* SCIPgetObjProbData(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
