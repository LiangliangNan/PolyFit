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
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_table.h
 * @ingroup INTERNALAPI
 * @brief  data structures for displaying statistics tables
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_TABLE_H__
#define __SCIP_STRUCT_TABLE_H__


#include "scip/def.h"
#include "scip/type_table.h"

#ifdef __cplusplus
extern "C" {
#endif

/** statistics table */
struct SCIP_Table
{
   char*                 name;               /**< name of statistics table */
   char*                 desc;               /**< description of statistics table */
   SCIP_DECL_TABLECOPY   ((*tablecopy));     /**< copy method of statistics table or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_TABLEFREE   ((*tablefree));     /**< destructor of statistics table */
   SCIP_DECL_TABLEINIT   ((*tableinit));     /**< initialize statistics table */
   SCIP_DECL_TABLEEXIT   ((*tableexit));     /**< deinitialize statistics table */
   SCIP_DECL_TABLEINITSOL ((*tableinitsol)); /**< solving process initialization method of statistics table */
   SCIP_DECL_TABLEEXITSOL ((*tableexitsol)); /**< solving process deinitialization method of statistics table */
   SCIP_DECL_TABLEOUTPUT ((*tableoutput));   /**< output method */
   SCIP_TABLEDATA*       tabledata;          /**< statistics table data */
   int                   position;           /**< relative position of statistics table */
   SCIP_STAGE            earlieststage;      /**< output of the statistics table is only printed from this stage onwards */
   SCIP_Bool             initialized;        /**< is statistics table initialized? */
   SCIP_Bool             active;             /**< should statistics table be displayed to the screen? */
};

#ifdef __cplusplus
}
#endif

#endif
