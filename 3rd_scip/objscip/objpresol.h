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

/**@file   objpresol.h
 * @brief  C++ wrapper for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJPRESOL_H__
#define __SCIP_OBJPRESOL_H__

#include <cstring>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for presolvers
 *
 *  This class defines the interface for presolvers implemented in C++. Note that there is a pure virtual
 *  function (this function has to be implemented). This function is: scip_exec().
 *
 *  - \ref PRESOL "Instructions for implementing a presolver"
 *  - \ref PRESOLVERS "List of available presolvers"
 *  - \ref type_presol.h "Corresponding C interface"
 */
class ObjPresol : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the presolver */
   char* scip_name_;

   /** description of the presolver */
   char* scip_desc_;

   /** default priority of the presolver */
   const int scip_priority_;

   /** default maximal number of presolving rounds the presolver participates in (-1: no limit) */
   const int scip_maxrounds_;

   /**< timing mask of the presolver */
   const SCIP_PRESOLTIMING scip_timing_;

   /** default constructor */
   ObjPresol(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of presolver */
      const char*        desc,               /**< description of presolver */
      int                priority,           /**< priority of the presolver */
      int                maxrounds,          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
      SCIP_PRESOLTIMING  timing              /**< timing mask of the presolver */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority),
        scip_maxrounds_(maxrounds),
        scip_timing_(timing)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** destructor */
   virtual ~ObjPresol()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** destructor of presolver to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_PRESOLFREE(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PRESOLFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of presolver (called after problem was transformed)
    *
    *  @see SCIP_DECL_PRESOLINIT(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PRESOLINIT(scip_init)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** deinitialization method of presolver (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_PRESOLEXIT(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PRESOLEXIT(scip_exit)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** presolving initialization method of presolver (called when presolving is about to begin)
    *
    *  @see SCIP_DECL_PRESOLINITPRE(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PRESOLINITPRE(scip_initpre)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** presolving deinitialization method of presolver (called after presolving has been finished)
    *
    *  @see SCIP_DECL_PRESOLEXITPRE(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PRESOLEXITPRE(scip_exitpre)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** execution method of presolver
    *
    *  @see SCIP_DECL_PRESOLEXEC(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PRESOLEXEC(scip_exec) = 0;
};

} /* namespace scip */



/** creates the presolver for the given presolver object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyPresol* mypresol = new MyPresol(...);
 *       SCIP_CALL( SCIPincludeObjPresol(scip, &mypresol, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete mypresol;    // delete presol AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjPresol(scip, new MyPresol(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyPresol is called here
 */
EXTERN
SCIP_RETCODE SCIPincludeObjPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjPresol*      objpresol,          /**< presolver object */
   SCIP_Bool             deleteobject        /**< should the presolver object be deleted when presolver is freed? */
   );

/** returns the presol object of the given name, or 0 if not existing */
EXTERN
scip::ObjPresol* SCIPfindObjPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of presolver */
   );

/** returns the presol object for the given presolver */
EXTERN
scip::ObjPresol* SCIPgetObjPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol              /**< presolver */
   );

#endif
