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

/**@file   objrelax.h
 * @brief  C++ wrapper for relaxation handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJRELAX_H__
#define __SCIP_OBJRELAX_H__

#include <cstring>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for relaxation handlers
 *
 *  This class defines the interface for relaxation handlers implemented in C++. Note that there is a pure virtual
 *  function (this function has to be implemented). This function is: scip_exec().
 *
 *  - \ref RELAX "Instructions for implementing a relaxation handler"
 *  - \ref type_relax.h "Corresponding C interface"
 */
class ObjRelax : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the relaxator */
   char* scip_name_;

   /** description of the relaxator */
   char* scip_desc_;

   /** default priority of the relaxator (negative: call after LP, non-negative: call before LP) */
   const int scip_priority_;

   /** frequency for calling relaxator */
   const int scip_freq_;

   /** does the relaxator contain all cuts in the LP? */
   const SCIP_Bool scip_includeslp_;

   /** default constructor */
   ObjRelax(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of relaxator */
      const char*        desc,               /**< description of relaxator */
      int                priority,           /**< priority of the relaxator (negative: after LP, non-negative: before LP) */
      int                freq,               /**< frequency for calling relaxator */
      SCIP_Bool          includeslp          /**< Does the relaxator contain all cuts in the LP? */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority),
        scip_freq_(freq),
        scip_includeslp_(includeslp)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** destructor */
   virtual ~ObjRelax()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** destructor of relaxator to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_RELAXFREE(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of relaxator (called after problem was transformed)
    *
    *  @see SCIP_DECL_RELAXINIT(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXINIT(scip_init)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** deinitialization method of relaxator (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_RELAXEXIT(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXEXIT(scip_exit)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of relaxator (called when branch and bound process is about to begin)
    *
    *  @see SCIP_DECL_RELAXINITSOL(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXINITSOL(scip_initsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of relaxator (called before branch and bound process data is freed)
    *
    *  @see SCIP_DECL_RELAXEXITSOL(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXEXITSOL(scip_exitsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** execution method of relaxator
    *
    *  @see SCIP_DECL_RELAXEXEC(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXEXEC(scip_exec) = 0;
};

} /* namespace scip */



/** creates the relaxator for the given relaxator object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyRelax* myrelax = new MyRelax(...);
 *       SCIP_CALL( SCIPincludeObjRelax(scip, &myrelax, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myrelax;    // delete relax AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjRelax(scip, new MyRelax(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyRelax is called here
 */
EXTERN
SCIP_RETCODE SCIPincludeObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjRelax*  objrelax,           /**< relaxator object */
   SCIP_Bool             deleteobject        /**< should the relaxator object be deleted when relaxator is freed? */
   );

/** returns the relax object of the given name, or 0 if not existing */
EXTERN
scip::ObjRelax* SCIPfindObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of relaxator */
   );

/** returns the relax object for the given relaxator */
EXTERN
scip::ObjRelax* SCIPgetObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax               /**< relaxator */
   );

#endif
