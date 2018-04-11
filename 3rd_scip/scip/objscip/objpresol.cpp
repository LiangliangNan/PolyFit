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

/**@file   objpresol.cpp
 * @brief  C++ wrapper for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objpresol.h"




/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   scip::ObjPresol*      objpresol;          /**< presolver object */
   SCIP_Bool             deleteobject;       /**< should the presolver object be deleted when presolver is freed? */
};




/*
 * Callback methods of presolver
 */

extern "C"
{

/** copy method for presolver plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyObj)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   
   assert(scip != NULL);
   
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);
   assert(presoldata->objpresol->scip_ != scip);

   if( presoldata->objpresol->iscloneable() )
   {
      scip::ObjPresol* newobjpresol;
      newobjpresol = dynamic_cast<scip::ObjPresol*> (presoldata->objpresol->clone(scip));

      /* call include method of presolver object */
      SCIP_CALL( SCIPincludeObjPresol(scip, newobjpresol, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeObj)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);
   assert(presoldata->objpresol->scip_ == scip);

   /* call virtual method of presol object */
   SCIP_CALL( presoldata->objpresol->scip_free(scip, presol) );

   /* free presol object */
   if( presoldata->deleteobject )
      delete presoldata->objpresol;

   /* free presol data */
   delete presoldata;
   SCIPpresolSetData(presol, NULL); /*lint !e64*/
   
   return SCIP_OKAY;
}


/** initialization method of presolver (called after problem was transformed) */
static
SCIP_DECL_PRESOLINIT(presolInitObj)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);
   assert(presoldata->objpresol->scip_ == scip);

   /* call virtual method of presol object */
   SCIP_CALL( presoldata->objpresol->scip_init(scip, presol) );

   return SCIP_OKAY;
}


/** deinitialization method of presolver (called before transformed problem is freed) */
static
SCIP_DECL_PRESOLEXIT(presolExitObj)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);

   /* call virtual method of presol object */
   SCIP_CALL( presoldata->objpresol->scip_exit(scip, presol) );

   return SCIP_OKAY;
}


/** presolving initialization method of presolver (called when presolving is about to begin) */
static
SCIP_DECL_PRESOLINITPRE(presolInitpreObj)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);

   /* call virtual method of presol object */
   SCIP_CALL( presoldata->objpresol->scip_initpre(scip, presol) );

   return SCIP_OKAY;
}


/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreObj)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);

   /* call virtual method of presol object */
   SCIP_CALL( presoldata->objpresol->scip_exitpre(scip, presol) );

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecObj)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);

   /* call virtual method of presol object */
   SCIP_CALL( presoldata->objpresol->scip_exec(scip, presol, nrounds, presoltiming,
         nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
         nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides,
         nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes,
         ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

   return SCIP_OKAY;
}
}



/*
 * presolver specific interface methods
 */

/** creates the presolver for the given presolver object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjPresol*      objpresol,          /**< presolver object */
   SCIP_Bool             deleteobject        /**< should the presolver object be deleted when presolver is freed? */
   )
{
   SCIP_PRESOLDATA* presoldata;

   assert(scip != NULL);
   assert(objpresol != NULL);

   /* create presolver data */
   presoldata = new SCIP_PRESOLDATA;
   presoldata->objpresol = objpresol;
   presoldata->deleteobject = deleteobject;

   /* include presolver */
   SCIP_CALL( SCIPincludePresol(scip, objpresol->scip_name_, objpresol->scip_desc_,
         objpresol->scip_priority_, objpresol->scip_maxrounds_, objpresol->scip_timing_,
         presolCopyObj, presolFreeObj, presolInitObj, presolExitObj,
         presolInitpreObj, presolExitpreObj, presolExecObj,
         presoldata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the presol object of the given name, or 0 if not existing */
scip::ObjPresol* SCIPfindObjPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of presolver */
   )
{
   SCIP_PRESOL* presol;
   SCIP_PRESOLDATA* presoldata;

   presol = SCIPfindPresol(scip, name);
   if( presol == NULL )
      return 0;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   return presoldata->objpresol;
}
   
/** returns the presol object for the given presolver */
scip::ObjPresol* SCIPgetObjPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   return presoldata->objpresol;
}
