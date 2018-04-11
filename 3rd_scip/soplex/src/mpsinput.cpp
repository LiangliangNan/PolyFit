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

/**@file  mpsinput.cpp
 * @brief Read MPS format files.
 */

#include <assert.h>
#include <ctype.h>
#include <string.h>

#include "spxdefines.h"
#include "mpsinput.h"
#include "spxout.h"

#define PATCH_CHAR    '_'
#define BLANK         ' '

namespace soplex
{ 

/// fill the line from \p pos up to column 80 with blanks.
static void clear_from(char* buf, int pos)
{
   for(int i = pos; i < 80; i++)
      buf[i] = BLANK;
   buf[80] = '\0';
}

/// change all blanks inside a field to #PATCH_CHAR.
static void patch_field(char* buf, int beg, int end)
{
   while((beg <= end) && (buf[end] == BLANK))
      end--;

   while((beg <= end) && (buf[beg] == BLANK))
      beg++;

   for(int i = beg; i <= end; i++)
      if (buf[i] == BLANK)
         buf[i] = PATCH_CHAR;
}

/// read a MPS format data line and parse the fields.
bool MPSInput::readLine()
{
   int   len;
   int   space;
   char* s;
   bool  is_marker;
   bool  is_comment;

   do
   {
      m_f0 = m_f1 = m_f2 = m_f3 = m_f4 = m_f5 = 0;
      is_marker = false;
   
      // Read until we have a non-empty, non-comment line.
      do
      {
         if( !m_input.getline(m_buf, sizeof(m_buf)).good() && !m_input.eof() )
            return false;
        m_lineno++;

        MSG_DEBUG( std::cout << "DMPSIN01 Line " << m_lineno
                          << " " << m_buf << std::endl; )

        /* check if comment line */
        is_comment = true;
        if( m_buf[0] == '*' )
           continue;

        /* normalize line and check if it is empty */
        len = int(strlen(m_buf));
        for( int i = 0; i < len; i++ )
        {
           if( m_buf[i] == '\t' || m_buf[i] == '\n' || m_buf[i] == '\r' )
              m_buf[i] = BLANK;
           else if( m_buf[i] != BLANK )
              is_comment = false;
        }
      }
      while( is_comment );

      len = int(strlen(m_buf));
      if (len < 80)
         clear_from(m_buf, len);

      assert(strlen(m_buf) >= 80);

      /* Look for new section
       */
      if (*m_buf != BLANK)
      {
         m_f0 = strtok(&m_buf[0], " ");

         assert(m_f0 != 0);

         m_f1 = strtok(0, " ");

         return true;
      }

      if (!m_is_new_format)
      {
         /* Test for fixed format comments
          */
         if ((m_buf[14] == '$') && (m_buf[13] == ' '))
            clear_from(m_buf, 14);
         else if ((m_buf[39] == '$') && (m_buf[38] == ' '))
            clear_from(m_buf, 39);

         /* Test for fixed format
          */
         space = m_buf[12] | m_buf[13] 
            | m_buf[22] | m_buf[23] 
            | m_buf[36] | m_buf[37] | m_buf[38]
            | m_buf[47] | m_buf[48] 
            | m_buf[61] | m_buf[62] | m_buf[63];

         if (space == BLANK || len < 13)
         {
            /* Now we have space at the right positions.
             * But are there also the non space where they
             * should be ?
             */
            bool number = isdigit(m_buf[24]) || isdigit(m_buf[25]) 
               || isdigit(m_buf[26]) || isdigit(m_buf[27]) 
               || isdigit(m_buf[28]) || isdigit(m_buf[29]) 
               || isdigit(m_buf[30]) || isdigit(m_buf[31]) 
               || isdigit(m_buf[32]) || isdigit(m_buf[33]) 
               || isdigit(m_buf[34]) || isdigit(m_buf[35]); 
            
            /* len < 13 is handle ROW lines with embedded spaces 
             * in the names correctly
             */
            if (number || len < 13)
            {
               /* Now we assume fixed format, so we patch possible embedded spaces.
                */
               patch_field(m_buf,  4, 12);
               patch_field(m_buf, 14, 22);
               patch_field(m_buf, 39, 47);
            }
            else
            {
               if (  m_section == COLUMNS || m_section == RHS
                  || m_section == RANGES  || m_section == BOUNDS)
                  m_is_new_format = true;
            }
         }
         else
         {
            m_is_new_format = true;
         }
      }
      s = &m_buf[1];
      
      /* At this point it is not clear if we have a indicator field.
       * If there is none (e.g. empty) f1 will be the first name field.
       * If there is one, f2 will be the first name field.
       * 
       * Initially comment marks '$' ar only allowed in the beginning
       * of the 2nd and 3rd name field. We test all fields but the first.
       * This makes no difference, since if the $ is at the start of a value
       * field, the line will be errornous anyway.
       */
      do
      {
         if (0 == (m_f1 = strtok(s, " ")))
            break;
         
         if ((0 == (m_f2 = strtok(0, " "))) || (*m_f2 == '$'))
         {
            m_f2 = 0;
            break;      
         }
         if (!strcmp(m_f2, "'MARKER'"))
            is_marker = true;
            
         if ((0 == (m_f3 = strtok(0, " "))) || (*m_f3 == '$'))
         {
            m_f3 = 0;
            break;      
         }
         if (is_marker)
         {
            if (!strcmp(m_f3, "'INTORG'"))
               m_is_integer = true;
            else if (!strcmp(m_f3, "'INTEND'"))
               m_is_integer = false;
            else
               break; // unknown marker
         }

         if (!strcmp(m_f3, "'MARKER'"))
            is_marker = true;

         if ((0 == (m_f4 = strtok(0, " "))) || (*m_f4 == '$'))
         {
            m_f4 = 0;
            break;      
         }
         if (is_marker)
         {
            if (!strcmp(m_f4, "'INTORG'"))
               m_is_integer = true;
            else if (!strcmp(m_f4, "'INTEND'"))
               m_is_integer = false;
            else
               break; // unknown marker
         }

         if ((0 == (m_f5 = strtok(0, " "))) || (*m_f5 == '$'))
            m_f5 = 0;
      }
      while(false);
   }
   while(is_marker);

   MSG_DEBUG( std::cerr
             << "DMPSIN02 -----------------------------------------------"
             << std::endl
             << "DMPSIN03 f0=" << ((m_f0 == 0) ? "nil" : m_f0) << std::endl
             << "DMPSIN04 f1=" << ((m_f1 == 0) ? "nil" : m_f1) << std::endl
             << "DMPSIN05 f2=" << ((m_f2 == 0) ? "nil" : m_f2) << std::endl
             << "DMPSIN06 f3=" << ((m_f3 == 0) ? "nil" : m_f3) << std::endl
             << "DMPSIN07 f4=" << ((m_f4 == 0) ? "nil" : m_f4) << std::endl
             << "DMPSIN08 f5=" << ((m_f5 == 0) ? "nil" : m_f5) << std::endl
             << "DMPSIN09 -----------------------------------------------" 
             << std::endl;
   )

   return true;
}

/// Insert \p name as field 1 and shift all other fields up.
void MPSInput::insertName(const char* name, bool second)
{
   m_f5 = m_f4;
   m_f4 = m_f3;
   m_f3 = m_f2;

   if (second)
      m_f2 = name;
   else
   {
      m_f2 = m_f1;
      m_f1 = name;
   }
}

} // namespace soplex
