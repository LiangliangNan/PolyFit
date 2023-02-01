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

/**@file   objtable.cpp
 * @brief  C++ wrapper for statistics tables
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objtable.h"




/*
 * Data structures
 */

/** display table data */
struct SCIP_TableData
{
   scip::ObjTable*       objtable;           /**< display statistics table object */
   SCIP_Bool             deleteobject;       /**< should the statistics table object be deleted when statistics table is freed? */
};




/*
 * Callback methods of statistics table
 */

extern "C"
{

/** copy method for statistics table plugins (called when SCIP copies plugins) */
static
SCIP_DECL_TABLECOPY(tableCopyObj)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;

   assert(scip != NULL);

   tabledata = SCIPtableGetData(table);
   assert(tabledata != NULL);
   assert(tabledata->objtable != NULL);
   assert(tabledata->objtable->scip_ != scip);

   if( tabledata->objtable->iscloneable() )
   {
      scip::ObjTable*  newobjtable;
      newobjtable = dynamic_cast<scip::ObjTable*> (tabledata->objtable->clone(scip));

      /* call include method of display column object */
      SCIP_CALL( SCIPincludeObjTable(scip, newobjtable, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of statistics table to free user data (called when SCIP is exiting) */
static
SCIP_DECL_TABLEFREE(tableFreeObj)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;

   tabledata = SCIPtableGetData(table);
   assert(tabledata != NULL);
   assert(tabledata->objtable != NULL);
   assert(tabledata->objtable->scip_ == scip);

   /* call virtual method of statistics table object */
   SCIP_CALL( tabledata->objtable->scip_free(scip, table) );

   /* free statistics table object */
   if( tabledata->deleteobject )
      delete tabledata->objtable;

   /* free statistics table data */
   delete tabledata;
   SCIPtableSetData(table, NULL); /*lint !e64*/

   return SCIP_OKAY;
}


/** initialization method of statistics table (called after problem was transformed) */
static
SCIP_DECL_TABLEINIT(tableInitObj)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;

   tabledata = SCIPtableGetData(table);
   assert(tabledata != NULL);
   assert(tabledata->objtable != NULL);
   assert(tabledata->objtable->scip_ == scip);

   /* call virtual method of statistics table object */
   SCIP_CALL( tabledata->objtable->scip_init(scip, table) );

   return SCIP_OKAY;
}


/** deinitialization method of statistics table (called before transformed problem is freed) */
static
SCIP_DECL_TABLEEXIT(tableExitObj)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;

   tabledata = SCIPtableGetData(table);
   assert(tabledata != NULL);
   assert(tabledata->objtable != NULL);

   /* call virtual method of statistics table object */
   SCIP_CALL( tabledata->objtable->scip_exit(scip, table) );

   return SCIP_OKAY;
}


/** solving process initialization method of statistics table (called when branch and bound process is about to begin) */
static
SCIP_DECL_TABLEINITSOL(tableInitsolObj)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;

   tabledata = SCIPtableGetData(table);
   assert(tabledata != NULL);
   assert(tabledata->objtable != NULL);

   /* call virtual method of statistics table object */
   SCIP_CALL( tabledata->objtable->scip_initsol(scip, table) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of statistics table (called before branch and bound process data is freed) */
static
SCIP_DECL_TABLEEXITSOL(tableExitsolObj)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;

   tabledata = SCIPtableGetData(table);
   assert(tabledata != NULL);
   assert(tabledata->objtable != NULL);

   /* call virtual method of statistics table object */
   SCIP_CALL( tabledata->objtable->scip_exitsol(scip, table) );

   return SCIP_OKAY;
}


/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputObj)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;

   tabledata = SCIPtableGetData(table);
   assert(tabledata != NULL);
   assert(tabledata->objtable != NULL);

   /* call virtual method of statistics table object */
   SCIP_CALL( tabledata->objtable->scip_output(scip, table, file) );

   return SCIP_OKAY;
}
}



/*
 * statistics table specific interface methods
 */

/** creates the statistics table for the given statistics table object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjTable(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjTable*       objtable,           /**< statistics table object */
   SCIP_Bool             deleteobject        /**< should the statistics table object be deleted when statistics table is freed? */
   )
{
   SCIP_TABLEDATA* tabledata;

   assert(scip != NULL);
   assert(objtable != NULL);

   /* create statistics table data */
   tabledata = new SCIP_TABLEDATA;
   tabledata->objtable = objtable;
   tabledata->deleteobject = deleteobject;

   /* include statistics table */
   SCIP_CALL( SCIPincludeTable(scip, objtable->scip_name_, objtable->scip_desc_, TRUE,
         tableCopyObj, tableFreeObj, tableInitObj, tableExitObj, tableInitsolObj,
         tableExitsolObj, tableOutputObj, tabledata, objtable->scip_position_, objtable->scip_earlieststage_) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the statistics table object of the given name, or 0 if not existing */
scip::ObjTable* SCIPfindObjTable(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of statistics table */
   )
{
   SCIP_TABLE* table;
   SCIP_TABLEDATA* tabledata;

   table = SCIPfindTable(scip, name);
   if( table == NULL )
      return 0;

   tabledata = SCIPtableGetData(table);
   assert(tabledata != NULL);

   return tabledata->objtable;
}

/** returns the statistics table object for the given statistics table */
scip::ObjTable* SCIPgetObjTable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TABLE*           table               /**< statistics table */
   )
{
   SCIP_TABLEDATA* tabledata;

   tabledata = SCIPtableGetData(table);
   assert(tabledata != NULL);

   return tabledata->objtable;
}
