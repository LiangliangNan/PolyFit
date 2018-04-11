/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  mpsinput.h
 * @brief Read MPS format files.
 */
#ifndef _MPSINPUT_H_
#define _MPSINPUT_H_

#include <iostream>

#include "spxout.h"

namespace soplex
{ 

/**@class MPSInput
   
   Reads MPS input files. A right-hand side for the objective function is
   allowed but ignored.
 */
class MPSInput
{
public:

   //-----------------------------------
   /**@name Types */
   //@{
   ///

   enum Section
   {
      NAME, OBJSEN, OBJNAME, ROWS, COLUMNS, RHS, RANGES, BOUNDS, ENDATA
   };
   ///

   /// optimization sense.
   enum Sense
   {
      MAXIMIZE = 1,
      MINIMIZE = -1
   };

   enum { MAX_LINE_LEN = 256 };

   //@}

private:

   //-----------------------------------
   /**@name Private data */
   //@{
   ///
   Section         m_section;
   /// the input stream from which the file is read
   std::istream&   m_input;
   /// line number
   int             m_lineno;
   /// objctive sense (maximization or minimization)
   Sense           m_objsense;
   /// is set to \c true upon a syntax error
   bool            m_has_error;
   /// the line buffer
   char            m_buf[MAX_LINE_LEN];
   /// first field in a line
   const char*     m_f0;
   /// second field in a line
   const char*     m_f1;
   /// third field in a line
   const char*     m_f2;
   /// fourth field in a line
   const char*     m_f3;
   /// fifth field in a line
   const char*     m_f4;
   /// sixth field in a line
   const char*     m_f5;
   /// problem name
   char            m_probname[MAX_LINE_LEN];
   /// objective name
   char            m_objname [MAX_LINE_LEN];
   ///
   bool            m_is_integer;
   /// new MPS format?
   bool            m_is_new_format;
   /// Number of already ignored entries.
   int             m_ignored;
   /// Maximal number of ignored entries for which a warning will be issued.
   static const int m_max_ignore = 1000;
   //@}

public:

   //-----------------------------------
   /**@name Construction / destruction */
   //@{
   ///
   explicit
   MPSInput( std::istream& p_input )
      : m_section       ( NAME )
      , m_input         ( p_input )
      , m_lineno        ( 0 )
      , m_objsense      ( MPSInput::MINIMIZE )
      , m_has_error     ( false )
      , m_is_integer    ( false )
      , m_is_new_format ( false )
      , m_ignored       ( 0 )
   {
      m_f0 = m_f1 = m_f2 = m_f3 = m_f4 = m_f5 = 0;

      m_buf     [0] = '\0';
      m_probname[0] = '\0';
      m_objname [0] = '\0';
   }
   ///
   virtual 
   ~MPSInput()
   {
      // only to signal to flexelint that the pointers do
      // not point to anything that has to be freed.
      m_f0 = m_f1 = m_f2 = m_f3 = m_f4 = m_f5 = 0;
   }
   //@}
   
   //-----------------------------------
   /**@name Access */
   //@{
   ///
   Section         section()   const { return m_section; }
   ///
   int             lineno()    const { return m_lineno; }
   ///
   const char*     field0()    const { return m_f0; }
   ///
   const char*     field1()    const { return m_f1; }
   ///
   const char*     field2()    const { return m_f2; }
   ///
   const char*     field3()    const { return m_f3; }
   ///
   const char*     field4()    const { return m_f4; }
   ///
   const char*     field5()    const { return m_f5; }
   ///
   const char*     probName()  const { return m_probname; }
   ///
   const char*     objName()   const { return m_objname; }
   ///
   Sense           objSense()  const { return m_objsense; }
   ///
   bool            hasError()  const { return m_has_error; }
   ///
   bool            isInteger() const { return m_is_integer; }
   //@}

   //-----------------------------------
   /**@name Modification */
   //@{
   ///
   void setSection(Section p_section)
   {
      m_section = p_section;
   }
   ///
   void setProbName(const char* p_probname)
   {
      assert(strlen(p_probname) < MAX_LINE_LEN);
      spxSnprintf(m_probname, MAX_LINE_LEN, p_probname);
   }
   ///
   void setObjName(const char* p_objname)
   {
      assert(strlen(p_objname) < MAX_LINE_LEN);
      spxSnprintf(m_objname, MAX_LINE_LEN, p_objname);
   }
   ///
   void setObjSense(Sense sense)
   {
      m_objsense = sense;
   }
   //@}

   //-----------------------------------
   /**@name Warnings and Errors */
   //@{
   ///
   void syntaxError() 
   {
      MSG_ERROR( std::cerr << "Syntax error in line " << m_lineno << std::endl; )
      m_section = ENDATA;
      m_has_error = true;
   }
   ///
   void entryIgnored(
      const char* what, const char* what_name, 
      const char* entity, const char* entity_name)
   {
      if ( m_ignored < m_max_ignore )
      {
         MSG_ERROR( std::cerr << "Warning: line " << m_lineno << ": "
                              << what << " \"" << what_name << "\""
                              << " for " << entity << " \""
                              << entity_name << "\" ignored" << std::endl; )
         ++m_ignored;

         if ( m_ignored == m_max_ignore )
            MSG_ERROR( std::cerr << "Warning: This was the " << m_max_ignore
                                 << " ignored entry. No further warnings on "
                                 << "ignored entries will be given." << std::endl; )
      }
   }
   //@}

   //-----------------------------------
   /**@name Helpers */
   //@{
   /// reads an MPS format data line and parse the fields.
   bool readLine();
   /// Inserts \p name as field 1 and shifts all other fields up.
   void insertName( const char* name, 
                    bool second = false );
   //@}
};
} // namespace soplex
#endif // _MPSINPUT_H_
