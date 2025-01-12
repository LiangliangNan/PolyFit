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

/**@file   objdisp.h
 * @brief  C++ wrapper for display columns
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJDISP_H__
#define __SCIP_OBJDISP_H__

#include <cstring>
#include <utility>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/**
 *  @brief C++ wrapper for display columns
 *
 *  This class defines the interface for display columns implemented in C++. Note that there is a pure virtual function
 *  (this function has to be implemented). This function is: scip_output().
 *
 * - \ref DISP "Instructions for implementing a display column"
 * - \ref DISPLAYS "List of available display columns"
 *  - \ref type_disp.h "Corresponding C interface"
 */
class ObjDisp : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the display column */
   char* scip_name_;

   /** description of the display column */
   char* scip_desc_;

   /** head line of the display column */
   char* scip_header_;

   /** width of the display column (no. of chars used) */
   const int scip_width_;

   /** priority of the display column */
   const int scip_priority_;

   /** relative position of the display column */
   const int scip_position_;

   /** should the column be separated with a line from its right neighbour? */
   const SCIP_Bool scip_stripline_;

   /** default constructor */
   ObjDisp(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of display column */
      const char*        desc,               /**< description of display column */
      const char*        header,             /**< head line of display column */
      int                width,              /**< width of display column (no. of chars used) */
      int                priority,           /**< priority of display column */
      int                position,           /**< relative position of display column */
      SCIP_Bool          stripline           /**< should the column be separated with a line from its right neighbour? */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_header_(0),
        scip_width_(width),
        scip_priority_(priority),
        scip_position_(position),
        scip_stripline_(stripline)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_header_, header, std::strlen(header)+1) );
   }

   /** copy constructor */
   ObjDisp(const ObjDisp& o)
       : ObjDisp(o.scip_, o.scip_name_, o.scip_desc_, o.scip_header_, o.scip_width_, o.scip_priority_, o.scip_position_,
                 o.scip_stripline_)
   {
   }

   /** move constructor */
   ObjDisp(ObjDisp&& o)
       : scip_(o.scip_),
         scip_name_(0),
         scip_desc_(0),
         scip_header_(0),
         scip_width_(o.scip_width_),
         scip_priority_(o.scip_priority_),
         scip_position_(o.scip_position_),
         scip_stripline_(o.scip_stripline_)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
      std::swap(scip_header_, o.scip_header_);
   }

   /** destructor */
   virtual ~ObjDisp()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
      SCIPfreeMemoryArray(scip_, &scip_header_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjDisp& operator=(const ObjDisp& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjDisp& operator=(ObjDisp&& o) = delete;

   /** destructor of display column to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_DISPFREE(x) in @ref type_disp.h
    */
   virtual SCIP_DECL_DISPFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of display column (called after problem was transformed)
    *
    *  @see SCIP_DECL_DISPINIT(x) in @ref type_disp.h
    */
   virtual SCIP_DECL_DISPINIT(scip_init)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** deinitialization method of display column (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_DISPEXIT(x) in @ref type_disp.h
    */
   virtual SCIP_DECL_DISPEXIT(scip_exit)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of display column (called when branch and bound process is about to begin)
    *
    *  @see SCIP_DECL_DISPINITSOL(x) in @ref type_disp.h
    */
   virtual SCIP_DECL_DISPINITSOL(scip_initsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of display column (called before branch and bound process data is freed)
    *
    *  @see SCIP_DECL_DISPEXITSOL(x) in @ref type_disp.h
    */
   virtual SCIP_DECL_DISPEXITSOL(scip_exitsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** output method of display column to output file stream 'file'
    *
    *  @see SCIP_DECL_DISPOUTPUT(x) in @ref type_disp.h
    */
   virtual SCIP_DECL_DISPOUTPUT(scip_output) = 0;
};

} /* namespace scip */



/** creates the display column for the given display column object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyDisp* mydisp = new MyDisp(...);
 *       SCIP_CALL( SCIPincludeObjDisp(scip, &mydisp, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete mydisp;    // delete disp AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjDisp(scip, new MyDisp(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyDisp is called here
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjDisp*        objdisp,            /**< display column object */
   SCIP_Bool             deleteobject        /**< should the display column object be deleted when display column is freed? */
   );

/** returns the display column object of the given name, or 0 if not existing */
SCIP_EXPORT
scip::ObjDisp* SCIPfindObjDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of display column */
   );

/** returns the display column object for the given display column */
SCIP_EXPORT
scip::ObjDisp* SCIPgetObjDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DISP*            disp                /**< display column */
   );

#endif
