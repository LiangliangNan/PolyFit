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

/**@file   objsepa.cpp
 * @brief  C++ wrapper for cut separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objsepa.h"




/*
 * Data structures
 */

/** cut separator data */
struct SCIP_SepaData
{
   scip::ObjSepa*        objsepa;            /**< cut separator object */
   SCIP_Bool             deleteobject;       /**< should the cut separator object be deleted when cut separator is freed? */
};



/*
 * Callback methods of cut separator
 */

extern "C"
{

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyObj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   
   assert(scip != NULL);
   
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->objsepa != NULL);
   assert(sepadata->objsepa->scip_ != scip);

   if( sepadata->objsepa->iscloneable() )
   {
      scip::ObjSepa* newobjsepa;
      newobjsepa = dynamic_cast<scip::ObjSepa*> (sepadata->objsepa->clone(scip));

      /* call include method of separator object */
      SCIP_CALL( SCIPincludeObjSepa(scip, newobjsepa, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of cut separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeObj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->objsepa != NULL);
   assert(sepadata->objsepa->scip_ == scip);

   /* call virtual method of sepa object */
   SCIP_CALL( sepadata->objsepa->scip_free(scip, sepa) );

   /* free sepa object */
   if( sepadata->deleteobject )
      delete sepadata->objsepa;

   /* free sepa data */
   delete sepadata;
   SCIPsepaSetData(sepa, NULL); /*lint !e64*/
   
   return SCIP_OKAY;
}


/** initialization method of cut separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitObj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->objsepa != NULL);
   assert(sepadata->objsepa->scip_ == scip);

   /* call virtual method of sepa object */
   SCIP_CALL( sepadata->objsepa->scip_init(scip, sepa) );

   return SCIP_OKAY;
}


/** deinitialization method of cut separator (called before transformed problem is freed) */
static
SCIP_DECL_SEPAEXIT(sepaExitObj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->objsepa != NULL);

   /* call virtual method of sepa object */
   SCIP_CALL( sepadata->objsepa->scip_exit(scip, sepa) );

   return SCIP_OKAY;
}


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
static
SCIP_DECL_SEPAINITSOL(sepaInitsolObj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->objsepa != NULL);

   /* call virtual method of sepa object */
   SCIP_CALL( sepadata->objsepa->scip_initsol(scip, sepa) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolObj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->objsepa != NULL);

   /* call virtual method of sepa object */
   SCIP_CALL( sepadata->objsepa->scip_exitsol(scip, sepa) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpObj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->objsepa != NULL);

   /* call virtual method of sepa object */
   SCIP_CALL( sepadata->objsepa->scip_execlp(scip, sepa, result, allowlocal) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolObj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->objsepa != NULL);

   /* call virtual method of sepa object */
   SCIP_CALL( sepadata->objsepa->scip_execsol(scip, sepa, sol, result, allowlocal) );

   return SCIP_OKAY;
}
}



/*
 * cut separator specific interface methods
 */

/** creates the cut separator for the given cut separator object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjSepa*        objsepa,            /**< cut separator object */
   SCIP_Bool             deleteobject        /**< should the cut separator object be deleted when cut separator is freed? */
   )
{
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);
   assert(objsepa != NULL);

   /* create cut separator data */
   sepadata = new SCIP_SEPADATA;
   sepadata->objsepa = objsepa;
   sepadata->deleteobject = deleteobject;

   /* include cut separator */
   SCIP_CALL( SCIPincludeSepa(scip, objsepa->scip_name_, objsepa->scip_desc_, objsepa->scip_priority_,
         objsepa->scip_freq_, objsepa->scip_maxbounddist_, objsepa->scip_usessubscip_, objsepa->scip_delay_,
         sepaCopyObj, sepaFreeObj, sepaInitObj, sepaExitObj, sepaInitsolObj, sepaExitsolObj,
         sepaExeclpObj, sepaExecsolObj,
         sepadata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the sepa object of the given name, or 0 if not existing */
scip::ObjSepa* SCIPfindObjSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of cut separator */
   )
{
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;

   sepa = SCIPfindSepa(scip, name);
   if( sepa == NULL )
      return 0;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   return sepadata->objsepa;
}
   
/** returns the sepa object for the given cut separator */
scip::ObjSepa* SCIPgetObjSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa                /**< cut separator */
   )
{
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   return sepadata->objsepa;
}
