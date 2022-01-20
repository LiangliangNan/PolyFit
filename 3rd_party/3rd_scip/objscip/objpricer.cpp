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

/**@file   objpricer.cpp
 * @brief  C++ wrapper for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objpricer.h"




/*
 * Data structures
 */

/** variable pricer data */
struct SCIP_PricerData
{
   scip::ObjPricer*      objpricer;          /**< variable pricer object */
   SCIP_Bool             deleteobject;       /**< should the pricer object be deleted when pricer is freed? */
};




/*
 * Callback methods of variable pricer
 */

extern "C"
{

/** copy method for pricer plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRICERCOPY(pricerCopyObj)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;
   
   assert(scip != NULL);
   
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);
   assert(pricerdata->objpricer->scip_ != scip);

   if( pricerdata->objpricer->iscloneable() )
   {
      scip::ObjPricer* newobjpricer;
      newobjpricer = dynamic_cast<scip::ObjPricer*> (pricerdata->objpricer->clone(scip, valid));

      /* call include method of pricer object */
      SCIP_CALL( SCIPincludeObjPricer(scip, newobjpricer, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreeObj)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);
   assert(pricerdata->objpricer->scip_ == scip);

   /* call virtual method of pricer object */
   SCIP_CALL( pricerdata->objpricer->scip_free(scip, pricer) );

   /* free pricer object */
   if( pricerdata->deleteobject )
      delete pricerdata->objpricer;

   /* free pricer data */
   delete pricerdata;
   SCIPpricerSetData(pricer, NULL); /*lint !e64*/
   
   return SCIP_OKAY;
}


/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitObj)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);
   assert(pricerdata->objpricer->scip_ == scip);

   /* call virtual method of pricer object */
   SCIP_CALL( pricerdata->objpricer->scip_init(scip, pricer) );

   return SCIP_OKAY;
}


/** deinitialization method of variable pricer (called before transformed problem is freed) */
static
SCIP_DECL_PRICEREXIT(pricerExitObj)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   SCIP_CALL( pricerdata->objpricer->scip_exit(scip, pricer) );

   return SCIP_OKAY;
}


/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
static
SCIP_DECL_PRICERINITSOL(pricerInitsolObj)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   SCIP_CALL( pricerdata->objpricer->scip_initsol(scip, pricer) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolObj)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   SCIP_CALL( pricerdata->objpricer->scip_exitsol(scip, pricer) );

   return SCIP_OKAY;
}


/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostObj)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   SCIP_CALL( pricerdata->objpricer->scip_redcost(scip, pricer, lowerbound, stopearly, result) );

   return SCIP_OKAY;
}


/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasObj)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   SCIP_CALL( pricerdata->objpricer->scip_farkas(scip, pricer, result) );

   return SCIP_OKAY;
}
}



/*
 * variable pricer specific interface methods
 */

/** creates the variable pricer for the given variable pricer object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjPricer(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjPricer*      objpricer,          /**< variable pricer object */
   SCIP_Bool             deleteobject        /**< should the pricer object be deleted when pricer is freed? */
   )
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(objpricer != NULL);

   /* create variable pricer data */
   pricerdata = new SCIP_PRICERDATA;
   pricerdata->objpricer = objpricer;
   pricerdata->deleteobject = deleteobject;

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricer(scip, objpricer->scip_name_, objpricer->scip_desc_, objpricer->scip_priority_,
         objpricer->scip_delay_,
         pricerCopyObj,
         pricerFreeObj, pricerInitObj, pricerExitObj,
         pricerInitsolObj, pricerExitsolObj, pricerRedcostObj, pricerFarkasObj,
         pricerdata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the variable pricer object of the given name, or 0 if not existing */
scip::ObjPricer* SCIPfindObjPricer(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of variable pricer */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   pricer = SCIPfindPricer(scip, name);
   if( pricer == NULL )
      return 0;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->objpricer;
}
   
/** returns the variable pricer object for the given pricer */
scip::ObjPricer* SCIPgetObjPricer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer              /**< pricer */
   )
{
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->objpricer;
}
