/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objconshdlr.h
 * @brief  C++ wrapper for constraint handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJCONSHDLR_H__
#define __SCIP_OBJCONSHDLR_H__


#include <cassert>
#include <cstring>
#include <utility>

#include "scip/scip.h"
#include "objscip/objprobcloneable.h"

namespace scip
{

/**
 *  @brief C++ wrapper for constraint handlers
 *
 *  This class defines the interface for constraint handlers implemented in C++. Note that there are pure virtual
 *  functions (these have to be implemented). These functions are: scip_trans(), scip_enfolp(), scip_enforelax(),
 *  scip_enfops(), scip_check(), and scip_lock().
 *
 *  - \ref CONS "Instructions for implementing a constraint handler"
 *  - \ref CONSHDLRS "List of available constraint handlers"
 *  - \ref type_cons.h "Corresponding C interface"
 */
class ObjConshdlr : public ObjProbCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the constraint handler */
   char* scip_name_;

   /** description of the constraint handler */
   char* scip_desc_;

   /** default separation priority of the constraint handler */
   const int scip_sepapriority_;

   /** default enforcing priority of the constraint handler */
   const int scip_enfopriority_;

   /** default checking priority of the constraint handler */
   const int scip_checkpriority_;

   /** default separation frequency of the constraint handler */
   const int scip_sepafreq_;

   /** default propagation frequency of the constraint handler */
   const int scip_propfreq_;

   /** default frequency of the constraint handler for eager evaluations in separation, propagation and enforcement */
   const int scip_eagerfreq_;

   /** maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   const int scip_maxprerounds_;

   /** should separation method be delayed, if other separators found cuts? */
   const SCIP_Bool scip_delaysepa_;

   /** should propagation method be delayed, if other propagators found reductions? */
   const SCIP_Bool scip_delayprop_;

   /** should the constraint handler be skipped, if no constraints are available? */
   const SCIP_Bool scip_needscons_;

   /** positions in the node solving loop where propagation method of constraint handler should be executed */
   const SCIP_PROPTIMING scip_proptiming_;

   /**< timing mask of the constraint handler's presolving method */
   const SCIP_PRESOLTIMING scip_presoltiming_;

   /** default constructor */
   ObjConshdlr(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of constraint handler */
      const char*        desc,               /**< description of constraint handler */
      int                sepapriority,       /**< priority of the constraint handler for separation */
      int                enfopriority,       /**< priority of the constraint handler for constraint enforcing */
      int                checkpriority,      /**< priority of the constraint handler for checking infeasibility (and propagation) */
      int                sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
      int                propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
      int                eagerfreq,          /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
      int                maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
      SCIP_Bool          delaysepa,          /**< should separation method be delayed, if other separators found cuts? */
      SCIP_Bool          delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
      SCIP_Bool          needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
      SCIP_PROPTIMING    proptiming,         /**< positions in the node solving loop where propagation method of constraint handlers should be executed */
      SCIP_PRESOLTIMING  presoltiming        /**< timing mask of the constraint handler's presolving method */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_sepapriority_(sepapriority),
        scip_enfopriority_(enfopriority),
        scip_checkpriority_(checkpriority),
        scip_sepafreq_(sepafreq),
        scip_propfreq_(propfreq),
        scip_eagerfreq_(eagerfreq),
        scip_maxprerounds_(maxprerounds),
        scip_delaysepa_(delaysepa),
        scip_delayprop_(delayprop),
        scip_needscons_(needscons),
        scip_proptiming_(proptiming),
        scip_presoltiming_(presoltiming)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** copy constructor */
   ObjConshdlr(const ObjConshdlr& o)
       : ObjConshdlr(o.scip_, o.scip_name_, o.scip_desc_, o.scip_sepapriority_, o.scip_enfopriority_,
                     o.scip_checkpriority_, o.scip_sepafreq_, o.scip_propfreq_, o.scip_eagerfreq_, o.scip_maxprerounds_,
                     o.scip_delaysepa_, o.scip_delayprop_, o.scip_needscons_, o.scip_proptiming_, o.scip_presoltiming_)
   {
   }

   /** move constructor */
   ObjConshdlr(ObjConshdlr&& o)
       : scip_(o.scip_),
         scip_name_(0),
         scip_desc_(0),
         scip_sepapriority_(o.scip_sepapriority_),
         scip_enfopriority_(o.scip_enfopriority_),
         scip_checkpriority_(o.scip_checkpriority_),
         scip_sepafreq_(o.scip_sepafreq_),
         scip_propfreq_(o.scip_propfreq_),
         scip_eagerfreq_(o.scip_eagerfreq_),
         scip_maxprerounds_(o.scip_maxprerounds_),
         scip_delaysepa_(o.scip_delaysepa_),
         scip_delayprop_(o.scip_delayprop_),
         scip_needscons_(o.scip_needscons_),
         scip_proptiming_(o.scip_proptiming_),
         scip_presoltiming_(o.scip_presoltiming_)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
   }

   /** destructor */
   virtual ~ObjConshdlr()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjConshdlr& operator=(const ObjConshdlr& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjConshdlr& operator=(ObjConshdlr&& o) = delete;

   /** destructor of constraint handler to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_CONSFREE(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of constraint handler (called after problem has been transformed)
    *
    *  @see SCIP_DECL_CONSINIT(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSINIT(scip_init)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** deinitialization method of constraint handler (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_CONSEXIT(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSEXIT(scip_exit)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** presolving initialization method of constraint handler (called when presolving is about to begin)
    *
    *  @see SCIP_DECL_CONSINITPRE(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSINITPRE(scip_initpre)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** presolving deinitialization method of constraint handler (called after presolving has been finished)
    *
    *  @see SCIP_DECL_CONSEXITPRE(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSEXITPRE(scip_exitpre)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of constraint handler (called when branch and bound process is about to begin)
    *
    *  @see SCIP_DECL_CONSINITSOL(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSINITSOL(scip_initsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of constraint handler (called before branch and bound process data is freed)
    *
    *  @see SCIP_DECL_CONSEXITSOL(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSEXITSOL(scip_exitsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** frees specific constraint data
    *
    *  @see SCIP_DECL_CONSDELETE(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSDELETE(scip_delete)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** transforms constraint data into data belonging to the transformed problem
    *
    *  @see SCIP_DECL_CONSTRANS(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSTRANS(scip_trans) = 0;

   /** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved)
    *
    *  @see SCIP_DECL_CONSINITLP(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSINITLP(scip_initlp)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** separation method of constraint handler for LP solution
    *
    *  @see SCIP_DECL_CONSSEPALP(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSSEPALP(scip_sepalp)
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** separation method of constraint handler for arbitrary primal solution
    *
    *  @see SCIP_DECL_CONSSEPASOL(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSSEPASOL(scip_sepasol)
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** constraint enforcing method of constraint handler for LP solutions
    *
    *  @see SCIP_DECL_CONSENFOLP(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSENFOLP(scip_enfolp) = 0;

   /** constraint enforcing method of constraint handler for relaxation solutions
    *
    *  @see SCIP_DECL_CONSENFORELAX(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSENFORELAX(scip_enforelax)
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** constraint enforcing method of constraint handler for pseudo solutions
    *
    *  @see SCIP_DECL_CONSENFOPS(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSENFOPS(scip_enfops) = 0;

   /** feasibility check method of constraint handler for primal solutions
    *
    *  @see SCIP_DECL_CONSCHECK(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSCHECK(scip_check) = 0;

   /** domain propagation method of constraint handler
    *
    *  @see SCIP_DECL_CONSPROP(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSPROP(scip_prop)
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** presolving method of constraint handler
    *
    *  @see SCIP_DECL_CONSPRESOL(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSPRESOL(scip_presol)
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** propagation conflict resolving method of constraint handler
    *
    *  @see SCIP_DECL_CONSRESPROP(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSRESPROP(scip_resprop)
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   /** variable rounding lock method of constraint handler
    *
    *  @see SCIP_DECL_CONSLOCK(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSLOCK(scip_lock) = 0;

   /** constraint activation notification method of constraint handler
    *
    *  @see SCIP_DECL_CONSACTIVE(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSACTIVE(scip_active)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** constraint deactivation notification method of constraint handler
    *
    *  @see SCIP_DECL_CONSDEACTIVE(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSDEACTIVE(scip_deactive)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** constraint enabling notification method of constraint handler
    *
    *  @see SCIP_DECL_CONSENABLE(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSENABLE(scip_enable)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** constraint disabling notification method of constraint handler
    *
    *  @see SCIP_DECL_CONSDISABLE(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSDISABLE(scip_disable)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** variable deletion method of constraint handler
    *
    *  @see SCIP_DECL_CONSDELVARS(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSDELVARS(scip_delvars)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** constraint display method of constraint handler
    *
    *  @see SCIP_DECL_CONSPRINT(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSPRINT(scip_print)
   {  /*lint --e{715}*/
      if ( file == NULL )
	 fprintf(stdout, "constraint handler <%s> does not support printing constraints\n", SCIPconshdlrGetName(conshdlr));
      else
	 fprintf(file, "constraint handler <%s> does not support printing constraints\n", SCIPconshdlrGetName(conshdlr));
      return SCIP_OKAY;
   }

   /** constraint copying method of constraint handler
    *
    *  @see SCIP_DECL_CONSCOPY(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSCOPY(scip_copy)
   {  /*lint --e{715}*/
      *valid = FALSE;
      return SCIP_OKAY;
   }

   /** constraint parsing method of constraint handler
    *
    *  @see SCIP_DECL_CONSPARSE(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSPARSE(scip_parse)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** constraint method of constraint handler which returns the variables (if possible)
    *
    *  @see SCIP_DECL_CONSGETVARS(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSGETVARS(scip_getvars)
   {  /*lint --e{715}*/

      (*success) = FALSE;

      return SCIP_OKAY;
   }

   /** constraint method of constraint handler which returns the number of variables (if possible)
    *
    *  @see SCIP_DECL_CONSGETNVARS(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSGETNVARS(scip_getnvars)
   {  /*lint --e{715}*/

      (*nvars) = 0;
      (*success) = FALSE;

      return SCIP_OKAY;
   }

   /** constraint handler method to suggest dive bound changes during the generic diving algorithm
    *
    *  @see SCIP_DECL_CONSGETDIVEBDCHGS(x) in @ref type_cons.h
    */
   virtual SCIP_DECL_CONSGETDIVEBDCHGS(scip_getdivebdchgs)
   {  /*lint --e{715}*/

      (*success) = FALSE;

      return SCIP_OKAY;
   }
};

} /* namespace scip */



/** creates the constraint handler for the given constraint handler object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyConshdlr* myconshdlr = new MyConshdlr(...);
 *       SCIP_CALL( SCIPincludeObjConshdlr(scip, &myconshdlr, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myconshdlr;    // delete conshdlr AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjConshdlr(scip, new MyConshdlr(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyConshdlr is called here
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjConshdlr*    objconshdlr,        /**< constraint handler object */
   SCIP_Bool             deleteobject        /**< should the constraint handler object be deleted when conshdlr is freed? */
   );

/** returns the conshdlr object of the given name, or 0 if not existing */
SCIP_EXPORT
scip::ObjConshdlr* SCIPfindObjConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint handler */
   );

/** returns the conshdlr object for the given constraint handler */
SCIP_EXPORT
scip::ObjConshdlr* SCIPgetObjConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

#endif
