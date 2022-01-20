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

/**@file   objprobdata.cpp
 * @brief  C++ wrapper for user problem data
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objprobdata.h"




/*
 * Data structures
 */

/** user problem data */
struct SCIP_ProbData
{
   scip::ObjProbData*    objprobdata;        /**< user problem data object */
   SCIP_Bool             deleteobject;       /**< should the user problem data object be deleted when problem is freed? */
};




/*
 * Callback methods of user problem data
 */

extern "C"
{

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probDelorigObj)
{  /*lint --e{715}*/
   assert(probdata != NULL);
   assert(*probdata != NULL);
   assert((*probdata)->objprobdata != NULL);

   /* call virtual method of probdata object */
   SCIP_CALL( (*probdata)->objprobdata->scip_delorig(scip) );

   /* free probdata object */
   if( (*probdata)->deleteobject )
      delete (*probdata)->objprobdata;

   /* free probdata data */
   delete *probdata;
   *probdata = 0; /*lint !e64*/
   
   return SCIP_OKAY;
}


/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed)
 */
static
SCIP_DECL_PROBTRANS(probTransObj)
{  /*lint --e{715}*/
   scip::ObjProbData* objprobdata; /*lint !e78 !e40 !e55 !e530 !e522*/
   SCIP_Bool deleteobject;

   assert(sourcedata != NULL);
   assert(sourcedata->objprobdata != NULL);
   assert(targetdata != NULL);
   assert(*targetdata == NULL);

   /* call virtual method of probdata object */
   SCIP_CALL( sourcedata->objprobdata->scip_trans(scip, &objprobdata, &deleteobject) ); /*lint !e40*/

   /* create transformed user problem data */
   *targetdata = new SCIP_PROBDATA;
   (*targetdata)->objprobdata = objprobdata; /*lint !e40*/
   (*targetdata)->deleteobject = deleteobject;

   return SCIP_OKAY;
}


/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probDeltransObj)
{  /*lint --e{715}*/
   assert(probdata != NULL);
   assert(*probdata != NULL);
   assert((*probdata)->objprobdata != NULL);

   /* call virtual method of probdata object */
   SCIP_CALL( (*probdata)->objprobdata->scip_deltrans(scip) );

   /* free probdata object */
   if( (*probdata)->deleteobject )
      delete (*probdata)->objprobdata;

   /* free probdata data */
   delete *probdata;
   *probdata = 0; /*lint !e64*/

   return SCIP_OKAY;
}


/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probInitsolObj)
{  /*lint --e{715}*/
   assert(probdata != NULL);
   assert(probdata->objprobdata != NULL);

   /* call virtual method of probdata object */
   SCIP_CALL( probdata->objprobdata->scip_initsol(scip) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probExitsolObj)
{  /*lint --e{715}*/
   assert(probdata != NULL);
   assert(probdata->objprobdata != NULL);

   /* call virtual method of probdata object */
   SCIP_CALL( probdata->objprobdata->scip_exitsol(scip, restart) );

   return SCIP_OKAY;
}

/** copies user data if you want to copy it to a subscip */
static
SCIP_DECL_PROBCOPY(probCopyObj)
{  /*lint --e{715}*/
   scip::ObjProbData* objprobdata; /*lint !e78 !e40 !e55 !e530 !e522*/

   assert(sourcedata != NULL);
   assert(sourcedata->objprobdata != NULL);
   assert(targetdata != NULL);
   assert(*targetdata == NULL);

   /* call virtual method of probdata object */
   SCIP_CALL( sourcedata->objprobdata->scip_copy(scip, sourcescip, varmap, consmap, &objprobdata, global, result) ); /*lint !e40*/

   if( objprobdata != 0 )
   {
      assert(*result == SCIP_SUCCESS);

      /* create trarget user problem data */
      *targetdata = new SCIP_PROBDATA;
      (*targetdata)->objprobdata = objprobdata; /*lint !e40*/
      (*targetdata)->deleteobject = TRUE; /* always delete object, because we created it */
   }
   else
   {
      assert(*result == SCIP_DIDNOTRUN);
      *targetdata = 0;
   }

   return SCIP_OKAY;
}

}



/*
 * user problem data specific interface methods
 */

/** creates empty problem, initializes all solving data structures, and sets the user problem data to point to the
 *  given user data object
 */
SCIP_RETCODE SCIPcreateObjProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */
   scip::ObjProbData*    objprobdata,        /**< user problem data object */
   SCIP_Bool             deleteobject        /**< should the user problem data object be deleted when problem is freed? */
   )
{
   SCIP_PROBDATA* probdata;

   /* create user problem data */
   probdata = new SCIP_PROBDATA;
   probdata->objprobdata = objprobdata;
   probdata->deleteobject = deleteobject;

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, name, probDelorigObj, probTransObj, probDeltransObj, 
         probInitsolObj, probExitsolObj, probCopyObj, probdata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** gets user problem data object
 *  Warning! This method should only be called after a problem was created with SCIPcreateObjProb().
 *  Otherwise, a segmentation fault may arise, or an undefined pointer is returned.
 */
scip::ObjProbData* SCIPgetObjProbData(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->objprobdata;
}

