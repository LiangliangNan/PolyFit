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

/**@file   objbranchrule.h
 * @brief  C++ wrapper for branching rules
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJBRANCHRULE_H__
#define __SCIP_OBJBRANCHRULE_H__


#include <cassert>
#include <cstring>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/**
 *  @brief C++ wrapper for branching rules
 *
 *  This class defines the interface for branching rules implemented in C++.
 *
 *  - \ref BRANCH "Instructions for implementing a branching rule"
 *  - \ref BRANCHINGRULES "List of available branching rules"
 *  - \ref type_branch.h "Corresponding C interface"
 */
class ObjBranchrule : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the branching rule */
   char* scip_name_;

   /** description of the branching rule */
   char* scip_desc_;

   /** default priority of the branching rule */
   const int scip_priority_;

   /** default maximal depth for applying the branching rule */
   const int scip_maxdepth_;

   /** default maximal relative distance from current node's dual bound to primal bound
    *  compared to best node's dual bound for applying branching rule
    *  (0.0: only on current best node, 1.0: on all nodes)
    */
   const SCIP_Real scip_maxbounddist_;

   /** default constructor */
   ObjBranchrule(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of branching rule */
      const char*        desc,               /**< description of branching rule */
      int                priority,           /**< priority of the branching rule */
      int                maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
      SCIP_Real          maxbounddist        /**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying branching rule
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority),
        scip_maxdepth_(maxdepth),
        scip_maxbounddist_(maxbounddist)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** destructor */
   virtual ~ObjBranchrule()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** destructor of branching rule to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_BRANCHFREE(x) in @ref type_branch.h
    */
   virtual SCIP_DECL_BRANCHFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of branching rule (called after problem was transformed)
    *
    *  @see SCIP_DECL_BRANCHINIT(x) in @ref type_branch.h
    */
   virtual SCIP_DECL_BRANCHINIT(scip_init)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** deinitialization method of branching rule (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_BRANCHEXIT(x) in @ref type_branch.h
    */
   virtual SCIP_DECL_BRANCHEXIT(scip_exit)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of branching rule (called when branch and bound process is about to begin)
    *
    *  @see SCIP_DECL_BRANCHINITSOL(x) in @ref type_branch.h
    */
   virtual SCIP_DECL_BRANCHINITSOL(scip_initsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of branching rule (called before branch and bound process data is freed)
    *
    *  @see SCIP_DECL_BRANCHEXITSOL(x) in @ref type_branch.h
    */
   virtual SCIP_DECL_BRANCHEXITSOL(scip_exitsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** branching execution method for fractional LP solutions
    *
    *  @see SCIP_DECL_BRANCHEXECLP(x) in @ref type_branch.h
    */
   virtual SCIP_DECL_BRANCHEXECLP(scip_execlp)
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** branching execution method for external candidates
    *
    *  @see SCIP_DECL_BRANCHEXECEXT(x) in @ref type_branch.h
    */
   virtual SCIP_DECL_BRANCHEXECEXT(scip_execext)
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** branching execution method for not completely fixed pseudo solutions
    *
    *  @see SCIP_DECL_BRANCHEXECPS(x) in @ref type_branch.h
    */
   virtual SCIP_DECL_BRANCHEXECPS(scip_execps)
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
};

} /* namespace scip */


/** creates the branching rule for the given branching rule object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyBranchrule* mybranchrule = new MyBranchrule(...);
 *       SCIP_CALL( SCIPincludeObjBranchrule(scip, &mybranchrule, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete mybranchrule;    // delete branchrule AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjBranchrule(scip, new MyBranchrule(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyBranchrule is called here
 */
EXTERN
SCIP_RETCODE SCIPincludeObjBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjBranchrule*  objbranchrule,      /**< branching rule object */
   SCIP_Bool             deleteobject        /**< should the branching rule object be deleted when branching rule is freed? */
   );

/** returns the branchrule object of the given name, or 0 if not existing */
EXTERN
scip::ObjBranchrule* SCIPfindObjBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of branching rule */
   );

/** returns the branchrule object for the given branching rule */
EXTERN
scip::ObjBranchrule* SCIPgetObjBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

#endif
