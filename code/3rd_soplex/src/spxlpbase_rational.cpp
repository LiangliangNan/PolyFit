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

/**@file  spxlpbase_rational.cpp
 * @brief Saving LPs with Rational values in a form suitable for SoPlex.
 */

#ifndef SOPLEX_LEGACY
#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include <iostream>

#include "spxdefines.h"
#include "spxlpbase.h"
#include "spxout.h"
#include "mpsinput.h"
#include "exceptions.h"

#define MAX_LINE_WRITE_LEN 65536   ///< maximum length allowed for writing lines

namespace soplex
{

template<>
void SPxLPBase<Rational>::computePrimalActivity(const VectorBase<Rational>& primal, VectorBase<Rational>& activity, const bool unscaled) const
{
   if( primal.dim() != nCols() )
   {
      throw SPxInternalCodeException("XSPXLP01 Primal vector for computing row activity has wrong dimension");
   }

   if( activity.dim() != nRows() )
   {
      throw SPxInternalCodeException("XSPXLP03 Activity vector computing row activity has wrong dimension");
   }

   int c;

   for( c = 0; c < nCols() && primal[c] == 0; c++ )
      ;

   if( c >= nCols() )
   {
      activity.clear();
      return;
   }

   activity = colVector(c);

   activity *= primal[c];
   c++;

   for( ; c < nCols(); c++ )
   {
      if( primal[c] != 0 )
      {
         activity.multAdd(primal[c], colVector(c));
      }
   }
}

template<>
void SPxLPBase<Rational>::computeDualActivity(const VectorBase<Rational>& dual, VectorBase<Rational>& activity, const bool unscaled) const
{
   if( dual.dim() != nRows() )
   {
      throw SPxInternalCodeException("XSPXLP02 Dual vector for computing dual activity has wrong dimension");
   }

   if( activity.dim() != nCols() )
   {
      throw SPxInternalCodeException("XSPXLP04 Activity vector computing dual activity has wrong dimension");
   }

   int r;

   for( r = 0; r < nRows() && dual[r] == 0; r++ )
      ;

   if( r >= nRows() )
   {
      activity.clear();
      return;
   }

   activity = rowVector(r);

   activity *= dual[r];
   r++;

   for( ; r < nRows(); r++ )
   {
      if( dual[r] != 0 )
      {
         activity.multAdd(dual[r], rowVector(r));
      }
   }
}

template<>
Rational SPxLPBase<Rational>::maxAbsNzo(bool /* unscaled */) const
{
   Rational maxi = Rational(0);

   for( int i = 0; i < nCols(); ++i )
   {
      Rational m = colVector(i).maxAbs();

      if( m > maxi )
         maxi = m;
   }

   assert(maxi >= Rational(0));

   return maxi;
}

template<>
Rational SPxLPBase<Rational>::minAbsNzo(bool /* unscaled */) const
{
   Rational mini = infinity;

   for( int i = 0; i < nCols(); ++i )
   {
      Rational m = colVector(i).minAbs();

      if( m < mini )
         mini = m;
   }

   assert(mini >= Rational(0));

   return mini;
}

// ---------------------------------------------------------------------------------------------------------------------
//  Specialization for reading LP format
// ---------------------------------------------------------------------------------------------------------------------

#define LPF_MAX_LINE_LEN  8192     ///< maximum length of a line (8190 + \\n + \\0)

/// Is \p c a \c space, \c tab, \c nl or \c cr ?
static inline bool LPFisSpace(int c)
{
   return (c == ' ') || (c == '\t') || (c == '\n') || (c == '\r');
}



/// Is there a number at the beginning of \p s ?
static bool LPFisValue(const char* s)
{
   return ((*s >= '0') && (*s <= '9')) || (*s == '+') || (*s == '-') || (*s == '.');
}



/// Is there a possible column name at the beginning of \p s ?
static bool LPFisColName(const char* s)
{
   // strchr() gives a true for the null char.
   if( *s == '\0' )
      return false;

   return ((*s >= 'A') && (*s <= 'Z'))
      || ((*s >= 'a') && (*s <= 'z'))
      || (strchr("!\"#$%&()/,;?@_'`{}|~", *s) != 0);
}



/// Is there a comparison operator at the beginning of \p s ?
static bool LPFisSense(const char* s)
{
   return (*s == '<') || (*s == '>') || (*s == '=');
}



static bool LPFisInfinity(const char* s)
{
   return ((s[0] == '-') || (s[0] == '+'))
      && (tolower(s[1]) == 'i')
      && (tolower(s[2]) == 'n')
      && (tolower(s[3]) == 'f');
}



static bool LPFisFree(const char* s)
{
   return (tolower(s[0]) == 'f')
      && ( tolower(s[1]) == 'r')
      && ( tolower(s[2]) == 'e')
      && ( tolower(s[3]) == 'e');
}



/// Read the next number and advance \p pos.
/** If only a sign is encountered, the number is assumed to be \c sign * 1.  This routine will not catch malformatted
 *  numbers like .e10 !
 */
static Rational LPFreadValue(char*& pos, SPxOut* spxout, const int lineno = -1)
{
   assert(LPFisValue(pos));

   char        tmp[LPF_MAX_LINE_LEN];
   const char* s = pos;
   char*       t;
   Rational        value = 1;
   bool        has_digits = false;
   bool        has_emptyexponent = false;
   bool        has_dot = false;
   bool        has_exponent = false;
   bool        has_emptydivisor = false;

   // 1. sign
   if( (*s == '+') || (*s == '-') )
      s++;

   // 2. Digits before the decimal dot
   while( (*s >= '0') && (*s <= '9') )
   {
      has_digits = true;
      s++;
   }

   // 3. Decimal dot
   if( *s == '.' )
   {
      has_dot = true;
      s++;

      // 4. If there was a dot, possible digit behind it
      while( (*s >= '0') && (*s <= '9') )
      {
         has_digits = true;
         s++;
      }
   }

   // 5. Exponent
   if( tolower(*s) == 'e' )
   {
      has_exponent = true;
      has_emptyexponent = true;
      s++;

      // 6. Exponent sign
      if( (*s == '+') || (*s == '-') )
         s++;

      // 7. Exponent digits
      while( (*s >= '0') && (*s <= '9') )
      {
         has_emptyexponent = false;
         s++;
      }
   }

   // 8. Division
   if( *s == '/' )
   {
      s++;
      has_emptydivisor = true;
      while( (*s >= '0') && (*s <= '9') )
      {
         has_emptydivisor = false;
         s++;
      }
      if( has_dot || has_exponent || has_emptydivisor ||
         (*s == '.') ||(*s == '+') || (*s == '-') || (tolower(*s) == 'e') )
      {
         MSG_WARNING( (*spxout), (*spxout) << "WLPFRD03 Warning: In line " << lineno << ": malformed rational value in LP file\n";)
      }
   }


   assert(s != pos);

   if( has_emptyexponent )
   {
      MSG_WARNING( (*spxout), (*spxout) << "WLPFRD01 Warning: In line " << lineno << ": found empty exponent in LP file - check for forbidden variable names with initial 'e' or 'E'\n" );
   }

   if( !has_digits )
      value = (*pos == '-') ? -1 : 1;
   else
   {
      for( t = tmp; pos != s; pos++ )
         *t++ = *pos;
      *t = '\0';

      if( !value.readString(tmp) )
      {
         MSG_WARNING( (*spxout), (*spxout) <<"WLPFRD04 Warning: In line " << lineno << ": malformed rational value in LP file\n" );
      }
   }
   pos += s - pos;

   assert(pos == s);

   MSG_DEBUG( std::cout << "DLPFRD01 LPFreadValue = " << value << std::endl; )

   if( LPFisSpace(*pos) )
      pos++;

   return value;
}



/// Read the next column name from the input.
/** The name read is looked up and if not found \p emptycol
 *  is added to \p colset. \p pos is advanced behind the name.
 *  @return The Index of the named column.
 */
static int LPFreadColName(char*& pos, NameSet* colnames, LPColSetBase<Rational>& colset, const LPColBase<Rational>* emptycol, SPxOut* spxout)
{
   assert(LPFisColName(pos));
   assert(colnames != 0);

   char        name[LPF_MAX_LINE_LEN];
   const char* s = pos;
   int         i;
   int         colidx;

   // These are the characters that are not allowed in a column name.
   while( (strchr("+-.<>= ", *s) == 0) && (*s != '\0') )
      s++;

   for( i = 0; pos != s; i++, pos++ )
      name[i] = *pos;

   name[i] = '\0';

   if( (colidx = colnames->number(name)) < 0 )
   {
      // We only add the name if we got an empty column.
      if( emptycol == 0 )
         MSG_WARNING( (*spxout), (*spxout) << "WLPFRD02 Unknown variable \"" << name << "\" "; )
      else
      {
         colidx = colnames->num();
         colnames->add(name);
         colset.add(*emptycol);
      }
   }

   MSG_DEBUG( std::cout << "DLPFRD03 LPFreadColName [" << name << "] = " << colidx << std::endl; )

   if( LPFisSpace(*pos) )
      pos++;

   return colidx;
}



/// Read the next <,>,=,==,<=,=<,>=,=> and advance \p pos.
static int LPFreadSense(char*& pos)
{
   assert(LPFisSense(pos));

   int sense = *pos++;

   if( (*pos == '<') || (*pos == '>') )
      sense = *pos++;
   else if( *pos == '=' )
      pos++;

   MSG_DEBUG( std::cout << "DLPFRD04 LPFreadSense = " << static_cast<char>(sense) << std::endl; )

   if( LPFisSpace(*pos) )
      pos++;

   return sense;
}



/// Is the \p keyword present in \p buf ? If yes, advance \p pos.
/** \p keyword should be lower case. It can contain optional sections which are enclosed in '[' ']' like "min[imize]".
 */
static bool LPFhasKeyword(char*& pos, const char* keyword)
{
   int i;
   int k;

   assert(keyword != 0);

   for( i = 0, k = 0; keyword[i] != '\0'; i++, k++ )
   {
      if( keyword[i] == '[' )
      {
         i++;

         // Here we assumed that we have a ']' for the '['.
         while( (tolower(pos[k]) == keyword[i]) && (pos[k] != '\0') )
         {
           k++;
           i++;
         }
         while( keyword[i] != ']' )
            i++;
         --k;
      }
      else
      {
         if( keyword[i] != tolower(pos[k]) )
            break;
      }
   }

   // we have to be at the end of the keyword and the word found on the line also has to end here.  Attention: The
   // LPFisSense is a kludge to allow LPFhasKeyword also to process Inf[inity] keywords in the bounds section.
   if( keyword[i] == '\0' && (pos[k] == '\0' || LPFisSpace(pos[k]) || LPFisSense(&pos[k])) )
   {
      pos += k;

      MSG_DEBUG( std::cout << "DLPFRD05 LPFhasKeyword: " << keyword << std::endl; )

      return true;
   }

   return false;
}



/// If \p buf start with "name:" store the name in \p rownames and advance \p pos.
static bool LPFhasRowName(char*& pos, NameSet* rownames)
{
   const char* s = strchr(pos, ':');

   if( s == 0 )
      return false;

   int dcolpos = int(s - pos);

   int end;
   int srt;

   // skip spaces between name and ":"
   for( end = dcolpos-1; end >= 0; end-- )
      if( pos[end] != ' ')
         break;

   // are there only spaces in front of the ":" ?
   if( end < 0 )
   {
      pos = &(pos[dcolpos+1]);
      return false;
   }

   // skip spaces in front of name
   for( srt = end-1; srt >= 0; srt-- )
      if( pos[srt] == ' ' )
         break;

   // go back to the non-space character
   srt++;

   assert(srt <= end && pos[srt] != ' ');

   char name[LPF_MAX_LINE_LEN];
   int i;
   int k = 0;

   for( i = srt; i <= end; i++ )
      name[k++] = pos[i];

   name[k] = '\0';

   if( rownames != 0 )
      rownames->add(name);

   pos = &(pos[dcolpos+1]);

   return true;
}



static Rational LPFreadInfinity(char*& pos)
{
   assert(LPFisInfinity(pos));

   Rational sense = (*pos == '-') ? -1 : 1;

   (void) LPFhasKeyword(++pos, "inf[inity]");

   sense *= Rational(infinity);
   return sense;
}



/// Read LP in "CPLEX LP File Format".
   /** The specification is taken from the ILOG CPLEX 7.0 Reference Manual, Appendix E, Page 527.
    *
    *  This routine should read (most?) valid LP format files.  What it will not do, is find all cases where a file is ill
    *  formed.  If this happens it may complain and read nothing or read "something".
    *
    *  Problem: A line ending in '+' or '-' followed by a line starting with a number, will be regarded as an error.
    *
    *  The reader will accept the keyword INT[egers] as a synonym for GEN[erals] which is an undocumented feature in CPLEX.
    *
    *  A difference to the CPLEX reader, is that no name for the objective row is required.
    *
    * The manual says the maximum allowed line length is 255 characters, but CPLEX does not complain if the lines are
    * longer.
    *
    *  @return true if the file was read correctly
    */
template <>
bool SPxLPBase<Rational>::readLPF(
   std::istream& p_input,                ///< input stream.
   NameSet*      p_rnames,               ///< row names.
   NameSet*      p_cnames,               ///< column names.
   DIdxSet*      p_intvars)              ///< integer variables.
{
   enum
   {
      START, OBJECTIVE, CONSTRAINTS, BOUNDS, INTEGERS, BINARIES
   } section = START;

   NameSet* rnames;                      ///< row names.
   NameSet* cnames;                      ///< column names.

   LPColSetBase<Rational> cset;              ///< the set of columns read.
   LPRowSetBase<Rational> rset;              ///< the set of rows read.
   LPColBase<Rational> emptycol;             ///< reusable empty column.
   LPRowBase<Rational> row;                  ///< last assembled row.
   DSVectorBase<Rational> vec;               ///< last assembled vector (from row).

   Rational val = 1;
   int colidx;
   int sense = 0;

   char buf[LPF_MAX_LINE_LEN];
   char tmp[LPF_MAX_LINE_LEN];
   char line[LPF_MAX_LINE_LEN];
   int lineno = 0;
   bool unnamed = true;
   bool finished = false;
   bool other;
   bool have_value = true;
   int i;
   int k;
   char* s;
   char* pos;
   char* pos_old = 0;

   if( p_cnames )
      cnames = p_cnames;
   else
   {
      cnames = 0;
      spx_alloc(cnames);
      cnames = new (cnames) NameSet();
   }

   cnames->clear();

   if( p_rnames )
      rnames = p_rnames;
   else
   {
      try
      {
         rnames = 0;
         spx_alloc(rnames);
         rnames = new (rnames) NameSet();
      }
      catch( const SPxMemoryException& x )
      {
         if( !p_cnames )
         {
            cnames->~NameSet();
            spx_free(cnames);
         }
         throw x;
      }
   }

   rnames->clear();

   SPxLPBase<Rational>::clear(); // clear the LP.

   //--------------------------------------------------------------------------
   //--- Main Loop
   //--------------------------------------------------------------------------
   for(;;)
   {
      // 0. Read a line from the file.
      if( !p_input.getline(buf, sizeof(buf)) )
      {
         if( strlen(buf) == LPF_MAX_LINE_LEN - 1 )
         {
            MSG_ERROR( std::cerr << "ELPFRD06 Line exceeds " << LPF_MAX_LINE_LEN - 2
                            << " characters" << std::endl; )
         }
         else
         {
            MSG_ERROR( std::cerr << "ELPFRD07 No 'End' marker found" << std::endl; )
            finished = true;
         }
         break;
      }
      lineno++;
      i   = 0;
      pos = buf;

      MSG_DEBUG( std::cout << "DLPFRD08 Reading line " << lineno
                        << " (pos=" << pos << ")" << std::endl; )

      // 1. Remove comments.
      if( 0 != (s = strchr(buf, '\\')) )
         *s = '\0';

      // 2. Look for keywords.
      if( section == START )
      {
         if( LPFhasKeyword(pos, "max[imize]") )
         {
            changeSense(SPxLPBase<Rational>::MAXIMIZE);
            section = OBJECTIVE;
         }
         else if( LPFhasKeyword(pos, "min[imize]") )
         {
            changeSense(SPxLPBase<Rational>::MINIMIZE);
            section = OBJECTIVE;
         }
      }
      else if( section == OBJECTIVE )
      {
         if( LPFhasKeyword(pos, "s[ubject][   ]t[o]")
            || LPFhasKeyword(pos, "s[uch][    ]t[hat]")
            || LPFhasKeyword(pos, "s[.][    ]t[.]")
            || LPFhasKeyword(pos, "lazy con[straints]") )
         {
            // store objective vector
            for( int j = vec.size() - 1; j >= 0; --j )
               cset.maxObj_w(vec.index(j)) = vec.value(j);
            // multiplication with -1 for minimization is done below
            vec.clear();
            have_value = true;
            val = 1;
            section = CONSTRAINTS;
         }
      }
      else if( section == CONSTRAINTS &&
              (LPFhasKeyword(pos, "s[ubject][   ]t[o]")
              || LPFhasKeyword(pos, "s[uch][    ]t[hat]")
              || LPFhasKeyword(pos, "s[.][    ]t[.]")) )
      {
         have_value = true;
         val = 1;
      }
      else
      {
         if( LPFhasKeyword(pos, "lazy con[straints]") )
            ;
         else if( LPFhasKeyword(pos, "bound[s]") )
            section = BOUNDS;
         else if( LPFhasKeyword(pos, "bin[ary]") )
            section = BINARIES;
         else if( LPFhasKeyword(pos, "bin[aries]") )
            section = BINARIES;
         else if( LPFhasKeyword(pos, "gen[erals]") )
            section = INTEGERS;
         else if( LPFhasKeyword(pos, "int[egers]") ) // this is undocumented
            section = INTEGERS;
         else if( LPFhasKeyword(pos, "end") )
         {
            finished = true;
            break;
         }
         else if( LPFhasKeyword(pos, "s[ubject][   ]t[o]") // second time
            || LPFhasKeyword(pos, "s[uch][    ]t[hat]")
            || LPFhasKeyword(pos, "s[.][    ]t[.]")
            || LPFhasKeyword(pos, "lazy con[straints]") )
         {
            // In principle this has to checked for all keywords above,
            // otherwise we just ignore any half finished constraint
            if( have_value )
               goto syntax_error;

            have_value = true;
            val = 1;
         }
      }

      // 3a. Look for row names in objective and drop it.
      if( section == OBJECTIVE )
         LPFhasRowName(pos, 0);

      // 3b. Look for row name in constraint and store it.
      if( section == CONSTRAINTS )
         if( LPFhasRowName(pos, rnames) )
            unnamed = false;

      // 4a. Remove initial spaces.
      while( LPFisSpace(pos[i]) )
         i++;

      // 4b. remove spaces if they do not appear before the name of a vaiable.
      for( k = 0; pos[i] != '\0'; i++ )
         if( !LPFisSpace(pos[i]) || LPFisColName(&pos[i + 1]) )
            tmp[k++] = pos[i];

      tmp[k] = '\0';

      // 5. Is this an empty line ?
      if( tmp[0] == '\0' )
         continue;

      // 6. Collapse sequences of '+' and '-'. e.g ++---+ => -
      for( i = 0, k = 0; tmp[i] != '\0'; i++ )
      {
         while( ((tmp[i] == '+') || (tmp[i] == '-')) && ((tmp[i + 1] == '+') || (tmp[i + 1] == '-')) )
         {
            if( tmp[i++] == '-' )
               tmp[i] = (tmp[i] == '-') ? '+' : '-';
         }
         line[k++] = tmp[i];
      }
      line[k] = '\0';

      //-----------------------------------------------------------------------
      //--- Line processing loop
      //-----------------------------------------------------------------------
      pos = line;

      MSG_DEBUG( std::cout << "DLPFRD09 pos=" << pos << std::endl; )

      // 7. We have something left to process.
      while( (pos != 0) && (*pos != '\0') )
      {
         // remember our position, so we are sure we make progress.
         pos_old = pos;

         // now process the sections
         switch( section )
         {
         case OBJECTIVE:
            if( LPFisValue(pos) )
            {
               Rational pre_sign = 1;

               /* Already having here a value could only result from being the first number in a constraint, or a sign
                * '+' or '-' as last token on the previous line.
                */
               if( have_value )
               {
                  if( spxAbs(val) != 1 )
                     goto syntax_error;

                  if( val == -1 )
                     pre_sign = val;
               }
               have_value = true;
               val = LPFreadValue(pos, spxout, lineno);
               val *= pre_sign;
            }
            if( *pos == '\0' )
               continue;

            if( !have_value || !LPFisColName(pos) )
               goto syntax_error;

            have_value = false;
            colidx = LPFreadColName(pos, cnames, cset, &emptycol, spxout);
            vec.add(colidx, val);
            break;
         case CONSTRAINTS:
            if( LPFisValue(pos) )
            {
               Rational pre_sign = 1;

               /* Already having here a value could only result from being the first number in a constraint, or a sign
                * '+' or '-' as last token on the previous line.
                */
               if( have_value )
               {
                  if( spxAbs(val) != 1 )
                     goto syntax_error;

                  if( val == -1 )
                     pre_sign = val;
               }

               have_value = true;
               val = LPFreadValue(pos, spxout, lineno);
               val *= pre_sign;

               if( sense != 0 )
               {
                  if( sense == '<' )
                  {
                     row.setLhs(-infinity);
                     row.setRhs(val);
                  }
                  else if( sense == '>' )
                  {
                     row.setLhs(val);
                     row.setRhs(infinity);
                  }
                  else
                  {
                     assert(sense == '=');

                     row.setLhs(val);
                     row.setRhs(val);
                  }
                  row.setRowVector(vec);
                  rset.add(row);
                  vec.clear();

                  if( !unnamed )
                     unnamed = true;
                  else
                  {
                     char name[16];
                     spxSnprintf(name, 16, "C%d", rset.num());
                     rnames->add(name);
                  }
                  have_value = true;
                  val = 1;
                  sense = 0;
                  pos = 0;
                  // next line
                  continue;
               }
            }
            if( *pos == '\0' )
               continue;

            if( have_value )
            {
               if( LPFisColName(pos) )
               {
                  colidx = LPFreadColName(pos, cnames, cset, &emptycol, spxout);

                  if( val != 0 )
                  {
                     // Do we have this index already in the row?
                     int n = vec.pos(colidx);

                     // if not, add it
                     if( n < 0 )
                        vec.add(colidx, val);
                     // if yes, add them up and remove the element if it amounts to zero
                     else
                     {
                        assert(vec.index(n) == colidx);

                        val += vec.value(n);

                        if( val == 0 )
                           vec.remove(n);
                        else
                           vec.value(n) = val;

                        assert(cnames->has(colidx));

                        MSG_WARNING( (*spxout), (*spxout) << "WLPFRD10 Duplicate index "
                                            << (*cnames)[colidx]
                                            << " in line " << lineno
                                            << std::endl; )
                     }
                  }
                  have_value = false;
               }
               else
               {
                  // We have a row like c1: <= 5 with no variables. We can not handle 10 <= 5; issue a syntax error.
                  if( val != 1 )
                     goto syntax_error;

                  // If the next thing is not the sense we give up also.
                  if( !LPFisSense(pos) )
                     goto syntax_error;

                  have_value = false;
               }
            }
            assert(!have_value);

            if( LPFisSense(pos) )
               sense = LPFreadSense(pos);
            break;
         case BOUNDS:
            other = false;
            sense = 0;

            if( LPFisValue(pos) )
            {
               val = LPFisInfinity(pos) ? LPFreadInfinity(pos) : LPFreadValue(pos, spxout, lineno);

               if( !LPFisSense(pos) )
                  goto syntax_error;

               sense = LPFreadSense(pos);
               other = true;
            }

            if( !LPFisColName(pos) )
               goto syntax_error;

            if( (colidx = LPFreadColName(pos, cnames, cset, 0, spxout)) < 0 )
            {
               MSG_WARNING( (*spxout), (*spxout) << "WLPFRD11 in Bounds section line "
                                   << lineno << " ignored" << std::endl; )
               *pos = '\0';
               continue;
            }

            if( sense )
            {
               if( sense == '<' )
                  cset.lower_w(colidx) = val;
               else if( sense == '>' )
                  cset.upper_w(colidx) = val;
               else
               {
                  assert(sense == '=');
                  cset.lower_w(colidx) = val;
                  cset.upper_w(colidx) = val;
               }
            }

            if( LPFisFree(pos) )
            {
               cset.lower_w(colidx) = -infinity;
               cset.upper_w(colidx) =  infinity;
               other = true;
               pos += 4;  // set position after the word "free"
            }
            else if( LPFisSense(pos) )
            {
               sense = LPFreadSense(pos);
               other = true;

               if( !LPFisValue(pos) )
                  goto syntax_error;

               val = LPFisInfinity(pos) ? LPFreadInfinity(pos) : LPFreadValue(pos,  spxout, lineno);

               if( sense == '<' )
                  cset.upper_w(colidx) = val;
               else if( sense == '>' )
                  cset.lower_w(colidx) = val;
               else
               {
                  assert(sense == '=');
                  cset.lower_w(colidx) = val;
                  cset.upper_w(colidx) = val;
               }
            }

            /* Do we have only a single column name in the input line?  We could ignore this savely, but it is probably
             * a sign of some other error.
             */
            if( !other )
               goto syntax_error;
            break;
         case BINARIES:
         case INTEGERS:
            if( (colidx = LPFreadColName(pos, cnames, cset, 0, spxout)) < 0 )
            {
               MSG_WARNING( (*spxout), (*spxout) << "WLPFRD12 in Binary/General section line " << lineno << " ignored" << std::endl; )
            }
            else
            {
               if( section == BINARIES )
               {
		  if( cset.lower(colidx) < 0 )
		  {
		     cset.lower_w(colidx) = 0;
		  }
		  if( cset.upper(colidx) > 1 )
		  {
		     cset.upper_w(colidx) = 1;
		  }
               }

               if( p_intvars != 0 )
                  p_intvars->addIdx(colidx);
            }
            break;
         case START:
            MSG_ERROR( std::cerr << "ELPFRD13 This seems to be no LP format file" << std::endl; )
            goto syntax_error;
         default:
            throw SPxInternalCodeException("XLPFRD01 This should never happen.");
         }

         if( pos == pos_old )
            goto syntax_error;
      }
   }

   assert(isConsistent());

   addCols(cset);
   assert(isConsistent());

   addRows(rset);
   assert(isConsistent());

syntax_error:
   if( finished )
   {
      MSG_INFO2( (*spxout), (*spxout) << "Finished reading " << lineno << " lines" << std::endl; )
   }
   else
      MSG_ERROR( std::cerr << "ELPFRD15 Syntax error in line " << lineno << std::endl; )

   if( p_cnames == 0 )
      spx_free(cnames);
   if( p_rnames == 0 )
      spx_free(rnames);

   return finished;
}



// ---------------------------------------------------------------------------------------------------------------------
// Specialization for reading MPS format
// ---------------------------------------------------------------------------------------------------------------------

/// Process NAME section.
static void MPSreadName(MPSInput& mps, SPxOut* spxout)
{
   do
   {
      // This has to be the Line with the NAME section.
      if( !mps.readLine() || (mps.field0() == 0) || strcmp(mps.field0(), "NAME") )
         break;

      // Sometimes the name is omitted.
      mps.setProbName((mps.field1() == 0) ? "_MPS_" : mps.field1());

      MSG_INFO2( (*spxout), (*spxout) << "IMPSRD01 Problem name   : " << mps.probName() << std::endl; )

      // This has to be a new section
      if( !mps.readLine() || (mps.field0() == 0) )
         break;

      if( !strcmp(mps.field0(), "ROWS") )
         mps.setSection(MPSInput::ROWS);
      else if( !strncmp(mps.field0(), "OBJSEN", 6) )
         mps.setSection(MPSInput::OBJSEN);
      else if( !strcmp(mps.field0(), "OBJNAME") )
         mps.setSection(MPSInput::OBJNAME);
      else
         break;

      return;
   }
   while(false);

   mps.syntaxError();
}



/// Process OBJSEN section. This Section is an ILOG extension.
static void MPSreadObjsen(MPSInput& mps)
{
   do
   {
      // This has to be the Line with MIN or MAX.
      if( !mps.readLine() || (mps.field1() == 0) )
         break;

      if( !strcmp(mps.field1(), "MIN") )
         mps.setObjSense(MPSInput::MINIMIZE);
      else if( !strcmp(mps.field1(), "MAX") )
         mps.setObjSense(MPSInput::MAXIMIZE);
      else
         break;

      // Look for ROWS or OBJNAME Section
      if( !mps.readLine() || (mps.field0() == 0) )
         break;

      if( !strcmp(mps.field0(), "ROWS") )
         mps.setSection(MPSInput::ROWS);
      else if( !strcmp(mps.field0(), "OBJNAME") )
         mps.setSection(MPSInput::OBJNAME);
      else
         break;

      return;
   }
   while(false);

   mps.syntaxError();
}



/// Process OBJNAME section. This Section is an ILOG extension.
static void MPSreadObjname(MPSInput& mps)
{
   do
   {
      // This has to be the Line with the name.
      if( !mps.readLine() || (mps.field1() == 0) )
         break;

      mps.setObjName(mps.field1());

      // Look for ROWS Section
      if( !mps.readLine() || (mps.field0() == 0) )
         break;

      if( strcmp(mps.field0(), "ROWS") )
         break;

      mps.setSection(MPSInput::ROWS);

      return;
   }
   while(false);

   mps.syntaxError();
}



/// Process ROWS section.
static void MPSreadRows(MPSInput& mps, LPRowSetBase<Rational>& rset, NameSet& rnames, SPxOut* spxout)
{
   LPRowBase<Rational> row;

   while( mps.readLine() )
   {
      if( mps.field0() != 0 )
      {
         MSG_INFO2( (*spxout), (*spxout) << "IMPSRD02 Objective name : " << mps.objName() << std::endl; )

         if( strcmp(mps.field0(), "COLUMNS") )
            break;

         mps.setSection(MPSInput::COLUMNS);

         return;
      }

      if( *mps.field1() == 'N' )
      {
         if( *mps.objName() == '\0' )
            mps.setObjName(mps.field2());
      }
      else
      {
         if( rnames.has(mps.field2()) )
            break;

         rnames.add(mps.field2());

         switch( *mps.field1() )
         {
         case 'G':
            row.setLhs(0);
            row.setRhs(infinity);
            break;
         case 'E':
            row.setLhs(0);
            row.setRhs(0);
            break;
         case 'L':
            row.setLhs(-infinity);
            row.setRhs(0);
            break;
         default:
            mps.syntaxError();
            return;
         }

         rset.add(row);
      }

      assert((*mps.field1() == 'N') || (rnames.number(mps.field2()) == rset.num() - 1));
   }

   mps.syntaxError();
}



/// Process COLUMNS section.
static void MPSreadCols(MPSInput& mps, const LPRowSetBase<Rational>& rset, const NameSet&  rnames, LPColSetBase<Rational>& cset, NameSet& cnames, DIdxSet* intvars, SPxOut* spxout)
{
   Rational val;
   int idx;
   char colname[MPSInput::MAX_LINE_LEN] = { '\0' };
   LPColBase<Rational> col(rset.num());
   DSVectorBase<Rational> vec;

   col.setObj(0);
   vec.clear();

   while( mps.readLine() )
   {
      if( mps.field0() != 0 )
      {
         if( strcmp(mps.field0(), "RHS") )
            break;

         if( colname[0] != '\0' )
         {
            col.setColVector(vec);
            cset.add(col);
         }

         mps.setSection(MPSInput::RHS);

         return;
      }

      if( (mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0) )
         break;

      // new column?
      if( strcmp(colname, mps.field1()) )
      {
         // first column?
         if( colname[0] != '\0' )
         {
            col.setColVector(vec);
            cset.add(col);
         }

         // save copy of string (make sure string ends with \0)
         spxSnprintf(colname, MPSInput::MAX_LINE_LEN-1, mps.field1());
         colname[MPSInput::MAX_LINE_LEN-1] = '\0';
         cnames.add(colname);
         vec.clear();
         col.setObj(0);
         col.setLower(0);
         col.setUpper(infinity);

         if( mps.isInteger() )
         {
            assert(cnames.number(colname) == cset.num());

            if( intvars != 0 )
               intvars->addIdx(cnames.number(colname));

            // for Integer variable the default bounds are 0/1
            col.setUpper(1);
         }
      }

      if( !val.readString(mps.field3()) )
      {
         MSG_WARNING( (*spxout), (*spxout) <<"WMPSRD01 Warning: malformed rational value in MPS file\n" );
      }


      if( !strcmp(mps.field2(), mps.objName()) )
         col.setObj(val);
      else
      {
         if( (idx = rnames.number(mps.field2())) < 0 )
            mps.entryIgnored("Column", mps.field1(), "row", mps.field2());
         else if( val != 0 )
            vec.add(idx, val);
      }

      if( mps.field5() != 0 )
      {
         assert(mps.field4() != 0);

         if( !val.readString(mps.field5()) )
         {
            MSG_WARNING( (*spxout), (*spxout) <<"WMPSRD02 Warning: malformed rational value in MPS file\n" );
         }

         if( !strcmp(mps.field4(), mps.objName()) )
            col.setObj(val);
         else
         {
            if( (idx = rnames.number(mps.field4())) < 0 )
               mps.entryIgnored("Column", mps.field1(), "row", mps.field4());
            else if( val != 0 )
               vec.add(idx, val);
         }
      }
   }

   mps.syntaxError();
}



/// Process RHS section.
static void MPSreadRhs(MPSInput& mps, LPRowSetBase<Rational>& rset, const NameSet& rnames, SPxOut* spxout)
{
   char rhsname[MPSInput::MAX_LINE_LEN] = { '\0' };
   char addname[MPSInput::MAX_LINE_LEN] = { '\0' };
   int idx;
   Rational val;

   while( mps.readLine() )
   {
      if( mps.field0() != 0 )
      {
         MSG_INFO2( (*spxout), (*spxout) << "IMPSRD03 RHS name       : " << rhsname  << std::endl; );

         if( !strcmp(mps.field0(), "RANGES") )
            mps.setSection(MPSInput::RANGES);
         else if( !strcmp(mps.field0(), "BOUNDS") )
            mps.setSection(MPSInput::BOUNDS);
         else if( !strcmp(mps.field0(), "ENDATA") )
            mps.setSection(MPSInput::ENDATA);
         else
            break;

         return;
      }

      if( ((mps.field2() != 0) && (mps.field3() == 0)) || ((mps.field4() != 0) && (mps.field5() == 0)) )
         mps.insertName("_RHS_");

      if( (mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0) )
         break;

      if( *rhsname == '\0' )
         spxSnprintf(rhsname, MPSInput::MAX_LINE_LEN , mps.field1());

      if( strcmp(rhsname, mps.field1()) )
      {
         if( strcmp(addname, mps.field1()) )
         {
            assert(strlen(mps.field1()) < MPSInput::MAX_LINE_LEN);
            spxSnprintf(addname, MPSInput::MAX_LINE_LEN, mps.field1());
            MSG_INFO3( (*spxout), (*spxout) << "IMPSRD07 RHS ignored    : " << addname << std::endl );
         }
      }
      else
      {
         if( (idx = rnames.number(mps.field2())) < 0 )
            mps.entryIgnored("RHS", mps.field1(), "row", mps.field2());
         else
         {
            if( !val.readString(mps.field3()) )
            {
               MSG_WARNING( (*spxout), (*spxout) <<"WMPSRD03 Warning: malformed rational value in MPS file\n" );
            }

            // LE or EQ
            if( double(rset.rhs(idx)) < double(infinity) )
               rset.rhs_w(idx) = val;
            // GE or EQ
            if( double(rset.lhs(idx)) > double(-infinity) )
               rset.lhs_w(idx) = val;
         }

         if( mps.field5() != 0 )
         {
            if( (idx = rnames.number(mps.field4())) < 0 )
               mps.entryIgnored("RHS", mps.field1(), "row", mps.field4());
            else
            {
               if( !val.readString(mps.field5()) )
               {
                  MSG_WARNING( (*spxout), (*spxout) <<"WMPSRD04 Warning: malformed rational value in MPS file\n" );
               }

               // LE or EQ
               if( double(rset.rhs(idx)) < double(infinity) )
                  rset.rhs_w(idx) = val;
               // GE or EQ
               if( double(rset.lhs(idx)) > double(-infinity) )
                  rset.lhs_w(idx) = val;
            }
         }
      }
   }

   mps.syntaxError();
}



/// Process RANGES section.
static void MPSreadRanges(MPSInput& mps,  LPRowSetBase<Rational>& rset, const NameSet& rnames, SPxOut* spxout)
{
   char rngname[MPSInput::MAX_LINE_LEN] = { '\0' };
   int idx;
   Rational val;

   while( mps.readLine() )
   {
      if( mps.field0() != 0 )
      {
         MSG_INFO2( (*spxout), (*spxout) << "IMPSRD04 Range name     : " << rngname << std::endl; );

         if( !strcmp(mps.field0(), "BOUNDS") )
            mps.setSection(MPSInput::BOUNDS);
         else if( !strcmp(mps.field0(), "ENDATA") )
            mps.setSection(MPSInput::ENDATA);
         else
            break;

         return;
      }

      if( ((mps.field2() != 0) && (mps.field3() == 0)) || ((mps.field4() != 0) && (mps.field5() == 0)) )
         mps.insertName("_RNG_");

      if( (mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0) )
         break;

      if( *rngname == '\0' )
      {
         assert(strlen(mps.field1()) < MPSInput::MAX_LINE_LEN);
         spxSnprintf(rngname, MPSInput::MAX_LINE_LEN, mps.field1());
      }

      /* The rules are:
       * Row Sign   LHS             RHS
       * ----------------------------------------
       *  G   +/-   rhs             rhs + |range|
       *  L   +/-   rhs - |range|   rhs
       *  E   +     rhs             rhs + range
       *  E   -     rhs + range     rhs
       * ----------------------------------------
       */
      if( !strcmp(rngname, mps.field1()) )
      {
         if( (idx = rnames.number(mps.field2())) < 0 )
            mps.entryIgnored("Range", mps.field1(), "row", mps.field2());
         else
         {
            if( !val.readString(mps.field3()) )
            {
               MSG_WARNING( (*spxout), (*spxout) <<"WMPSRD05 Warning: malformed rational value in MPS file\n" );
            }

            // EQ
            if( (double(rset.lhs(idx)) > -double(infinity)) && (double(rset.rhs_w(idx)) < double(infinity)) )
            {
               assert(rset.lhs(idx) == rset.rhs(idx));

               if( double(val) >= 0 )
                  rset.rhs_w(idx) += val;
               else
                  rset.lhs_w(idx) += val;
            }
            else
            {
               // GE
               if( double(rset.lhs(idx)) > -double(infinity) )
               {
                  rset.rhs_w(idx) = rset.lhs(idx);
                  rset.rhs_w(idx) += spxAbs(val);
               }
               // LE
               else
               {
                  rset.lhs_w(idx) = rset.rhs(idx);
                  rset.lhs_w(idx) -= spxAbs(val);
               }
            }
         }

         if( mps.field5() != 0 )
         {
            if( (idx = rnames.number(mps.field4())) < 0 )
               mps.entryIgnored("Range", mps.field1(), "row", mps.field4());
            else
            {
               if( !val.readString(mps.field5()) )
               {
                  MSG_WARNING( (*spxout), (*spxout) <<"WMPSRD06 Warning: malformed rational value in MPS file\n" );
               }

               // EQ
               if( (double(rset.lhs(idx)) > -double(infinity)) && (double(rset.rhs(idx)) <  double(infinity)) )
               {
                  assert(rset.lhs(idx) == rset.rhs(idx));

                  if( double(val) >= 0 )
                     rset.rhs_w(idx) += val;
                  else
                     rset.lhs_w(idx) += val;
               }
               else
               {
                  // GE
                  if( double(rset.lhs(idx)) > -double(infinity) )
                  {
                     rset.rhs_w(idx) = rset.lhs(idx);
                     rset.rhs_w(idx) += spxAbs(val);
                  }
                  // LE
                  else
                  {
                     rset.lhs_w(idx) = rset.rhs(idx);
                     rset.lhs_w(idx) -= spxAbs(val);
                  }
               }
            }
         }
      }
   }

   mps.syntaxError();
}



/// Process BOUNDS section.
static void MPSreadBounds(MPSInput& mps, LPColSetBase<Rational>& cset, const NameSet& cnames, DIdxSet* intvars, SPxOut* spxout)
{
   DIdxSet oldbinvars;
   char bndname[MPSInput::MAX_LINE_LEN] = { '\0' };
   int  idx;
   Rational val;

   while( mps.readLine() )
   {
      if( mps.field0() != 0 )
      {
         MSG_INFO2( (*spxout), (*spxout) << "IMPSRD05 Bound name     : " << bndname << std::endl; )

         if( strcmp(mps.field0(), "ENDATA") )
            break;

         mps.setSection(MPSInput::ENDATA);

         return;
      }

      // Is the value field used ?
      if(  (!strcmp(mps.field1(), "LO"))
         || (!strcmp(mps.field1(), "UP"))
         || (!strcmp(mps.field1(), "FX"))
         || (!strcmp(mps.field1(), "LI"))
         || (!strcmp(mps.field1(), "UI")) )
      {
         if( (mps.field3() != 0) && (mps.field4() == 0) )
            mps.insertName("_BND_", true);
      }
      else
      {
         if( (mps.field2() != 0) && (mps.field3() == 0) )
            mps.insertName("_BND_", true);
      }

      if( (mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0) )
         break;

      if( *bndname == '\0' )
      {
         assert(strlen(mps.field2()) < MPSInput::MAX_LINE_LEN);
         spxSnprintf(bndname, MPSInput::MAX_LINE_LEN, mps.field2());
      }

      // Only read the first Bound in section
      if( !strcmp(bndname, mps.field2()) )
      {
         if( (idx = cnames.number(mps.field3())) < 0 )
            mps.entryIgnored("column", mps.field3(), "bound", bndname);
         else
         {
            if( mps.field4() == 0 )
               val = 0;
            else if( !strcmp(mps.field4(), "-Inf") || !strcmp(mps.field4(), "-inf") )
               val = -infinity;
            else if( !strcmp(mps.field4(), "Inf") || !strcmp(mps.field4(), "inf") || !strcmp(mps.field4(), "+Inf") || !strcmp(mps.field4(), "+inf") )
               val = infinity;
            else if( !val.readString(mps.field4()) )
            {
               MSG_WARNING( (*spxout), (*spxout) <<"WMPSRD07 Warning: malformed rational value in MPS file line " << mps.lineno() << ": " << mps.field4() << "\n" );
               mps.syntaxError();
               return;
            }

            // ILOG extension (Integer Bound)
            if( mps.field1()[1] == 'I' )
            {
               if( intvars != 0 )
                  intvars->addIdx(idx);

               // if the variable has appeared in the MARKER section of the COLUMNS section then its default bounds were
               // set to 0,1; the first time it is declared integer we need to change to default bounds 0,infinity
               if( oldbinvars.pos(idx) < 0 )
               {
                  cset.upper_w(idx) = infinity;
                  oldbinvars.addIdx(idx);
               }
            }

            switch( *mps.field1() )
            {
            case 'L':
               cset.lower_w(idx) = val;
               break;
            case 'U':
               cset.upper_w(idx) = val;
               break;
            case 'F':
               if( mps.field1()[1] == 'X' )
               {
                  cset.lower_w(idx) = val;
                  cset.upper_w(idx) = val;
               }
               else
               {
                  cset.lower_w(idx) = -infinity;
                  cset.upper_w(idx) = infinity;
               }
               break;
            case 'M':
               cset.lower_w(idx) = -infinity;
               break;
            case 'P':
               cset.upper_w(idx) = infinity;
               break;
               // Ilog extension (Binary)
            case 'B':
               cset.lower_w(idx) = 0;
               cset.upper_w(idx) = 1;

               if( intvars != 0 )
                  intvars->addIdx(idx);
               break;
            default:
               mps.syntaxError();
               return;
            }
         }
      }
   }

   mps.syntaxError();
}



/// Read LP in MPS File Format.
/**
 *  The specification is taken from the IBM Optimization Library Guide and Reference, online available at
 *  http://www.software.ibm.com/sos/features/libuser.htm and from the ILOG CPLEX 7.0 Reference Manual, Appendix E, Page
 *  531.
 *
 *  This routine should read all valid MPS format files.  What it will not do, is find all cases where a file is ill
 *  formed.  If this happens it may complain and read nothing or read "something".
 *
 *  @return true if the file was read correctly.
 */
#define INIT_COLS 1000 ///< initialy allocated columns.
#define INIT_NZOS 5000 ///< initialy allocated non zeros.
template <>
bool SPxLPBase<Rational>::readMPS(
   std::istream& p_input,           ///< input stream.
   NameSet*      p_rnames,          ///< row names.
   NameSet*      p_cnames,          ///< column names.
   DIdxSet*      p_intvars)         ///< integer variables.
{
   LPRowSetBase<Rational>& rset = *this;
   LPColSetBase<Rational>& cset = *this;
   NameSet* rnames;
   NameSet* cnames;

   if( p_cnames )
      cnames = p_cnames;
   else
   {
      cnames = 0;
      spx_alloc(cnames);
      cnames = new (cnames) NameSet();
   }

   cnames->clear();

   if( p_rnames )
      rnames = p_rnames;
   else
   {
      try
      {
         rnames = 0;
         spx_alloc(rnames);
         rnames = new (rnames) NameSet();
      }
      catch( const SPxMemoryException& x )
      {
         if( !p_cnames )
         {
            cnames->~NameSet();
            spx_free(cnames);
         }
         throw x;
      }
   }

   rnames->clear();

   SPxLPBase<Rational>::clear(); // clear the LP.

   cset.memRemax(INIT_NZOS);
   cset.reMax(INIT_COLS);

   MPSInput mps(p_input);

   MPSreadName(mps, spxout);

   if( mps.section() == MPSInput::OBJSEN )
      MPSreadObjsen(mps);

   if( mps.section() == MPSInput::OBJNAME )
      MPSreadObjname(mps);

   if( mps.section() == MPSInput::ROWS )
      MPSreadRows(mps, rset, *rnames, spxout);

   addedRows(rset.num());

   if( mps.section() == MPSInput::COLUMNS )
      MPSreadCols(mps, rset, *rnames, cset, *cnames, p_intvars, spxout);

   if( mps.section() == MPSInput::RHS )
      MPSreadRhs(mps, rset, *rnames, spxout);

   if( mps.section() == MPSInput::RANGES )
      MPSreadRanges(mps, rset, *rnames, spxout);

   if( mps.section() == MPSInput::BOUNDS )
      MPSreadBounds(mps, cset, *cnames, p_intvars, spxout);

   if( mps.section() != MPSInput::ENDATA )
      mps.syntaxError();

   if( mps.hasError() )
      clear();
   else
   {
      changeSense(mps.objSense() == MPSInput::MINIMIZE ? SPxLPBase<Rational>::MINIMIZE : SPxLPBase<Rational>::MAXIMIZE);

      MSG_INFO2( (*spxout), (*spxout) << "IMPSRD06 Objective sense: " << ((mps.objSense() == MPSInput::MINIMIZE) ? "Minimize\n" : "Maximize\n") );

      added2Set(
         *(reinterpret_cast<SVSetBase<Rational>*>(static_cast<LPRowSetBase<Rational>*>(this))),
         *(reinterpret_cast<SVSetBase<Rational>*>(static_cast<LPColSetBase<Rational>*>(this))),
         cset.num());
      addedCols(cset.num());

      assert(isConsistent());
   }

   if( p_cnames == 0 )
   {
      cnames->~NameSet();
      spx_free(cnames);
   }

   if( p_rnames == 0 )
   {
      rnames->~NameSet();
      spx_free(rnames);
   }

   return !mps.hasError();
}



// ---------------------------------------------------------------------------------------------------------------------
// Specialization for writing LP format
// ---------------------------------------------------------------------------------------------------------------------

// get the name of a row or construct one
static const char* LPFgetRowName(
   const SPxLPBase<Rational>& p_lp,
   int                    p_idx,
   const NameSet*         p_rnames,
   char*                  p_buf,
   int                    p_num_written_rows
   )
{
   assert(p_buf != 0);
   assert(p_idx >= 0);
   assert(p_idx <  p_lp.nRows());

   if( p_rnames != 0 )
   {
      DataKey key = p_lp.rId(p_idx);

      if( p_rnames->has(key) )
         return (*p_rnames)[key];
   }

   spxSnprintf(p_buf, 16, "C%d", p_num_written_rows);

   return p_buf;
}



// get the name of a column or construct one
static const char* getColName(
   const SPxLPBase<Rational>& p_lp,
   int                    p_idx,
   const NameSet*         p_cnames,
   char*                  p_buf
   )
{
   assert(p_buf != 0);
   assert(p_idx >= 0);
   assert(p_idx <  p_lp.nCols());

   if( p_cnames != 0 )
   {
      DataKey key = p_lp.cId(p_idx);

      if( p_cnames->has(key) )
         return (*p_cnames)[key];
   }

   spxSnprintf(p_buf, 16, "x%d", p_idx);

   return p_buf;
}



// write an SVector
#define NUM_ENTRIES_PER_LINE 5
static void LPFwriteSVector(
   const SPxLPBase<Rational>&   p_lp,       ///< the LP
   std::ostream&            p_output,   ///< output stream
   const NameSet*           p_cnames,   ///< column names
   const SVectorBase<Rational>& p_svec,     ///< vector to write
   SPxOut*                  spxout      ///< out stream
   )
{

   char name[16];
   int num_coeffs = 0;
   long long pos;

   pos = p_output.tellp();

   for( int j = 0; j < p_lp.nCols(); ++j )
   {
      const Rational coeff = p_svec[j];

      if( coeff == 0 )
         continue;

      if( num_coeffs == 0 )
         p_output << coeff << " " << getColName(p_lp, j, p_cnames, name);
      else
      {
         // insert a line break every NUM_ENTRIES_PER_LINE columns or whenever max line length is nearly exceeded
         if( num_coeffs == NUM_ENTRIES_PER_LINE ||
            (long long)(p_output.tellp()) - pos + (long long)(rationalToString(coeff, false).length() + 100) > MAX_LINE_WRITE_LEN )
         {
            num_coeffs = 0;
            p_output << "\n\t";
            if( (long long)(p_output.tellp()) - pos  >  MAX_LINE_WRITE_LEN )
            {
               MSG_WARNING( (*spxout), (*spxout) << "XLPSWR01 Warning: MAX_LINE_WRITE_LEN possibly exceeded when writing LP file\n" );
            }
            pos = p_output.tellp();
         }

         if( coeff < 0 )
            p_output << " - " << -coeff;
         else
            p_output << " + " << coeff;

         p_output << " " << getColName(p_lp, j, p_cnames, name);
      }

      ++num_coeffs;
   }
}



// write the objective
static void LPFwriteObjective(
   const SPxLPBase<Rational>& p_lp,       ///< the LP
   std::ostream&          p_output,   ///< output stream
   const NameSet*         p_cnames,   ///< column names
   SPxOut*                spxout      ///< out stream
   )
{

   const int sense = p_lp.spxSense();

   p_output << ((sense == SPxLPBase<Rational>::MINIMIZE) ? "Minimize\n" : "Maximize\n");
   p_output << "  obj: ";

   const VectorBase<Rational>& obj = p_lp.maxObj();
   DSVectorBase<Rational> svec(obj.dim());
   svec.operator=(obj);
   svec *= Rational(sense);
   LPFwriteSVector(p_lp, p_output, p_cnames, svec, spxout);
   p_output << "\n";
}



// write non-ranged rows
static void LPFwriteRow(
   const SPxLPBase<Rational>&   p_lp,       ///< the LP
   std::ostream&            p_output,   ///< output stream
   const NameSet*           p_cnames,   ///< column names
   const SVectorBase<Rational>& p_svec,     ///< vector of the row
   const Rational&              p_lhs,      ///< lhs of the row
   const Rational&              p_rhs,      ///< rhs of the row
   SPxOut*                      spxout      ///< out stream
   )
{

   long long pos;
   pos = p_output.tellp();

   LPFwriteSVector(p_lp, p_output, p_cnames, p_svec, spxout);

   long long sidelen;
   sidelen = (p_lhs == p_rhs || double(p_lhs) <= double(-infinity)) ? (long long)rationalToString(p_rhs, false).length() : (long long)rationalToString(p_lhs, false).length();

   // insert a line break if max line length is in danger of being exceeded
   if( (long long)(p_output.tellp()) - pos + sidelen + (long long)100 > MAX_LINE_WRITE_LEN )
   {
      p_output << "\n\t";
      if( (long long)(p_output.tellp()) - pos  >  MAX_LINE_WRITE_LEN )
      {
         MSG_WARNING( (*spxout), (*spxout) << "XLPSWR02 Warning: MAX_LINE_WRITE_LEN possibly exceeded when writing LP file\n" );
      }
      pos = p_output.tellp();
   }

   // write bound value
   if( p_lhs == p_rhs )
      p_output << " = " << p_rhs;
   else if( double(p_lhs) <= double(-infinity) )
      p_output << " <= " << p_rhs;
   else
   {
      assert(double(p_rhs) >= double(infinity));
      p_output << " >= " << p_lhs;
   }

   p_output << "\n";

   if( (long long)(p_output.tellp()) - pos  >  MAX_LINE_WRITE_LEN )
   {
      MSG_WARNING( (*spxout), (*spxout) << "XLPSWR03 Warning: MAX_LINE_WRITE_LEN possibly exceeded when writing LP file\n" );
   }
}



// write all rows
static void LPFwriteRows(
   const SPxLPBase<Rational>& p_lp,       ///< the LP
   std::ostream&          p_output,   ///< output stream
   const NameSet*         p_rnames,   ///< row names
   const NameSet*         p_cnames,   ///< column names
   SPxOut*                spxout      ///< out stream
   )
{

   char name[16];

   p_output << "Subject To\n";

   for( int i = 0; i < p_lp.nRows(); ++i )
   {
      const Rational lhs = p_lp.lhs(i);
      const Rational rhs = p_lp.rhs(i);

      if( double(lhs) > -double(infinity) && double(rhs) < double(infinity) && lhs != rhs )
      {
         // ranged row -> write two non-ranged rows
         p_output << " " << LPFgetRowName(p_lp, i, p_rnames, name, i) << "_1 : ";
         LPFwriteRow(p_lp, p_output, p_cnames, p_lp.rowVector(i), lhs, infinity, spxout);

         p_output << " " << LPFgetRowName(p_lp, i, p_rnames, name, i) << "_2 : ";
         LPFwriteRow(p_lp, p_output, p_cnames, p_lp.rowVector(i), -infinity, rhs, spxout);
      }
      else
      {
         p_output << " " << LPFgetRowName(p_lp, i, p_rnames, name, i) << " : ";
         LPFwriteRow(p_lp, p_output, p_cnames, p_lp.rowVector(i), lhs, rhs, spxout);
      }
   }
}



// write the variable bounds
// (the default bounds 0 <= x <= infinity are not written)
static void LPFwriteBounds(
   const SPxLPBase<Rational>&   p_lp,       ///< the LP to write
   std::ostream&            p_output,   ///< output stream
   const NameSet*           p_cnames,   ///< column names
   SPxOut*                  spxout      ///< out stream
   )
{

   char name[16];
   long long pos;

   pos = p_output.tellp();

   p_output << "Bounds\n";

   for( int j = 0; j < p_lp.nCols(); ++j )
   {
      const Rational lower = p_lp.lower(j);
      const Rational upper = p_lp.upper(j);

      if( lower == upper )
      {
         p_output << "  "   << getColName(p_lp, j, p_cnames, name) << " = "  << upper << '\n';
      }
      else if( double(lower) > -double(infinity) )
      {
         if( double(upper) < double(infinity) )
         {
            // range bound
            if( lower != 0 )
               p_output << "  "   << lower << " <= "
                        << getColName(p_lp, j, p_cnames, name)
                        << " <= " << upper << '\n';
            else
               p_output << "  "   << getColName(p_lp, j, p_cnames, name)
                        << " <= " << upper << '\n';
         }
         else if( lower != 0 )
            p_output << "  " << lower << " <= "
                     << getColName(p_lp, j, p_cnames, name)
                     << '\n';
      }
      else if( double(upper) < double(infinity) )
         p_output << "   -Inf <= "
                  << getColName(p_lp, j, p_cnames, name)
                  << " <= " << upper << '\n';
      else
         p_output << "  "   << getColName(p_lp, j, p_cnames, name)
                  << " free\n";

      // check if max line length exceeded
      if( (long long)(p_output.tellp()) - pos  >  MAX_LINE_WRITE_LEN )
      {
         MSG_WARNING( (*spxout), (*spxout) << "XLPSWR04 Warning: MAX_LINE_WRITE_LEN exceeded when writing LP file\n" );
      }
      pos = p_output.tellp();
   }
}



// write the generals section
static void LPFwriteGenerals(
   const SPxLPBase<Rational>&   p_lp,         ///< the LP to write
   std::ostream&            p_output,     ///< output stream
   const NameSet*           p_cnames,     ///< column names
   const DIdxSet*           p_intvars     ///< integer variables
   )
{

   char name[16];

   if( p_intvars == NULL || p_intvars->size() <= 0 )
      return;  // no integer variables

   p_output << "Generals\n";

   for( int j = 0; j < p_lp.nCols(); ++j )
      if( p_intvars->pos(j) >= 0 )
         p_output << "  " << getColName(p_lp, j, p_cnames, name) << "\n";
}



/// Write LP in LP Format.
template <>
void SPxLPBase<Rational>::writeLPF(
   std::ostream&  p_output,          ///< output stream
   const NameSet* p_rnames,          ///< row names
   const NameSet* p_cnames,          ///< column names
   const DIdxSet* p_intvars          ///< integer variables
   ) const
{

   LPFwriteObjective(*this, p_output, p_cnames, spxout);
   LPFwriteRows(*this, p_output, p_rnames, p_cnames, spxout);
   LPFwriteBounds(*this, p_output, p_cnames, spxout);
   LPFwriteGenerals(*this, p_output, p_cnames, p_intvars);

   p_output << "End" << std::endl;
}



// ---------------------------------------------------------------------------------------------------------------------
// Specialization for writing MPS format
// ---------------------------------------------------------------------------------------------------------------------


static void MPSwriteRecord(
   std::ostream&  os,
   const char*    indicator,
   const char*    name,
   SPxOut* spxout,
   const char*    name1  = 0,
   const Rational value1 = 0,
   const char*    name2  = 0,
   const Rational value2 = 0
   )
{
   char buf[81];
   long long pos;
   pos = os.tellp();
   sprintf(buf, " %-2.2s %-8.8s", (indicator == 0) ? "" : indicator, (name == 0)      ? "" : name);
   os << buf;

   if( name1 != 0 )
   {
      spxSnprintf(buf, 81, " %-8.8s ", name1);
      os << buf << value1;
      if( name2 != 0 )
      {
         spxSnprintf(buf, 81, " %-8.8s ", name2);
         os << buf << value2;
      }
   }
   os << std::endl;

   // Warning if line is too long
   if( (long long)(os.tellp()) - pos > MAX_LINE_WRITE_LEN )
   {
      MSG_WARNING( (*spxout), (*spxout) << "XMPSWR04 Warning: MAX_LINE_WRITE_LEN exceeded when writing MPS file\n" );
   }
}



static Rational MPSgetRHS(Rational left, Rational right)
{
   Rational rhsval;

   if( double(left) > -double(infinity) ) /// This includes ranges
      rhsval = left;
   else if( double(right) <  double(infinity) )
      rhsval = right;
   else
      throw SPxInternalCodeException("XMPSWR01 This should never happen.");

   return rhsval;
}



static const char* MPSgetRowName(
   const SPxLPBase<Rational>& lp,
   int                   idx,
   const NameSet*        rnames,
   char*                 buf
   )
{
   assert(buf != 0);
   assert(idx >= 0);
   assert(idx <  lp.nRows());

   if( rnames != 0 )
   {
      DataKey key = lp.rId(idx);

      if( rnames->has(key) )
         return (*rnames)[key];
   }

   spxSnprintf(buf, 16, "C%d", idx);

   return buf;
}



/// Write LP in MPS format.
/** @note There will always be a BOUNDS section, even if there are no bounds.
 */
template <>
void SPxLPBase<Rational>::writeMPS(
   std::ostream&  p_output,          ///< output stream.
   const NameSet* p_rnames,          ///< row names.
   const NameSet* p_cnames,          ///< column names.
   const DIdxSet* p_intvars          ///< integer variables.
   ) const
{

   const char*    indicator;
   char           name [16];
   char           name1[16];
   char           name2[16];
   bool           has_ranges = false;
   int            i;
   int            k;

   // --- NAME Section ---
   p_output << "NAME          MPSDATA" << std::endl;

   // --- ROWS Section ---
   p_output << "ROWS" << std::endl;

   for( i = 0; i < nRows(); i++ )
   {
      if( lhs(i) == rhs(i) )
         indicator = "E";
      else if( (double(lhs(i)) > -double(infinity)) && (double(rhs(i)) < double(infinity)) )
      {
         indicator = "E";
         has_ranges = true;
      }
      else if( double(lhs(i)) > -double(infinity) )
         indicator = "G";
      else if( double(rhs(i)) <  double(infinity) )
         indicator = "L";
      else
         throw SPxInternalCodeException("XMPSWR02 This should never happen.");

      MPSwriteRecord(p_output, indicator, MPSgetRowName(*this, i, p_rnames, name), spxout);
   }

   MPSwriteRecord(p_output, "N", "MINIMIZE", spxout);

   // --- COLUMNS Section ---
   p_output << "COLUMNS" << std::endl;

   bool has_intvars = (p_intvars != 0) && (p_intvars->size() > 0);

   for( int j = 0; j < (has_intvars ? 2 : 1); j++ )
   {
      bool is_intrun = has_intvars && (j == 1);

      if( is_intrun )
         p_output << "    MARK0001  'MARKER'                 'INTORG'" << std::endl;

      for( i = 0; i < nCols(); i++ )
      {
         bool is_intvar = has_intvars && (p_intvars->pos(i) >= 0);

         if( ( is_intrun && !is_intvar) || (!is_intrun &&  is_intvar) )
             continue;

         const SVectorBase<Rational>& col = colVector(i);
         int colsize2 = (col.size() / 2) * 2;

         assert(colsize2 % 2 == 0);

         for( k = 0; k < colsize2; k += 2 )
            MPSwriteRecord(p_output, 0, getColName(*this, i, p_cnames, name), spxout,
               MPSgetRowName(*this, col.index(k), p_rnames, name1), col.value(k),
               MPSgetRowName(*this, col.index(k + 1), p_rnames, name2), col.value(k + 1));

         if( colsize2 != col.size() )
            MPSwriteRecord(p_output, 0, getColName(*this, i, p_cnames, name), spxout,
               MPSgetRowName(*this, col.index(k), p_rnames, name1), col.value(k));

         if( maxObj(i) != 0 )
            MPSwriteRecord(p_output, 0, getColName(*this, i, p_cnames, name), spxout, "MINIMIZE", -maxObj(i));
      }

      if( is_intrun )
         p_output << "    MARK0001  'MARKER'                 'INTEND'" << std::endl;
   }

   // --- RHS Section ---
   p_output << "RHS" << std::endl;

   i = 0;
   while( i < nRows() )
   {
      Rational rhsval1 = 0;
      Rational rhsval2 = 0;

      for( ; i < nRows(); i++ )
         if( (rhsval1 = MPSgetRHS(lhs(i), rhs(i))) != 0 )
            break;

      if( i < nRows() )
      {
         for( k = i + 1; k < nRows(); k++ )
         {
            if( (rhsval2 = MPSgetRHS(lhs(k), rhs(k))) != 0 )
               break;
         }

         if( k < nRows() )
         {
            MPSwriteRecord(p_output, 0, "RHS", spxout, MPSgetRowName(*this, i, p_rnames, name1), rhsval1,
               MPSgetRowName(*this, k, p_rnames, name2), rhsval2);
         }
         else
            MPSwriteRecord(p_output, 0, "RHS", spxout, MPSgetRowName(*this, i, p_rnames, name1), rhsval1);

         i = k + 1;
      }
   }

   // --- RANGES Section ---
   if( has_ranges )
   {
      p_output << "RANGES" << std::endl;

      for( i = 0; i < nRows(); i++ )
      {
         if( (double(lhs(i)) > -double(infinity)) && (double(rhs(i)) < double(infinity)) )
         {
            Rational range = rhs(i);
            range -= lhs(i);
            MPSwriteRecord(p_output, "", "RANGE", spxout, MPSgetRowName(*this, i, p_rnames, name1), range);
         }
      }
   }

   // --- BOUNDS Section ---
   p_output << "BOUNDS" << std::endl;

   for( i = 0; i < nCols(); i++ )
   {
      // skip variables that do not appear in the objective function or any constraint
      const SVectorBase<Rational>& col = colVector(i);

      if( col.size() == 0 && maxObj(i) == 0 )
         continue;

      if( lower(i) == upper(i) )
      {
         MPSwriteRecord(p_output, "FX", "BOUND", spxout, getColName(*this, i, p_cnames, name1), lower(i));
         continue;
      }

      if( (double(lower(i)) <= double(-infinity)) && (double(upper(i)) >= double(infinity)) )
      {
         MPSwriteRecord(p_output, "FR", "BOUND", spxout, getColName(*this, i, p_cnames, name1));
         continue;
      }

      if( lower(i) != 0 )
      {
         if( double(lower(i)) > -double(infinity) )
            MPSwriteRecord(p_output, "LO", "BOUND", spxout, getColName(*this, i, p_cnames, name1), lower(i));
         else
            MPSwriteRecord(p_output, "MI", "BOUND", spxout, getColName(*this, i, p_cnames, name1));
      }

      if( has_intvars && (p_intvars->pos(i) >= 0) )
      {
         // Integer variables have default upper bound 1, but we should write
         // it nevertheless since CPLEX seems to assume infinity otherwise.
         MPSwriteRecord(p_output, "UP", "BOUND", spxout, getColName(*this, i, p_cnames, name1), upper(i));
      }
      else
      {
         // Continous variables have default upper bound infinity
         if( double(upper(i)) < double(infinity) )
            MPSwriteRecord(p_output, "UP", "BOUND", spxout, getColName(*this, i, p_cnames, name1), upper(i));
      }
   }

   // --- ENDATA Section ---
   p_output << "ENDATA" << std::endl;

   // Output warning when writing a maximisation problem
   if( spxSense() == SPxLPBase<Rational>::MAXIMIZE )
   {
      MSG_WARNING( (*spxout), (*spxout) << "XMPSWR03 Warning: objective function inverted when writing maximization problem in MPS file format\n" );
   }
}



/// Building the dual problem from a given LP
/// @note primalRows must be as large as the number of unranged primal rows + 2 * the number of ranged primal rows.
///       dualCols must have the identical size to the primal rows.
template < >
void SPxLPBase<Rational>::buildDualProblem(SPxLPBase<Rational>& dualLP, SPxRowId primalRowIds[], SPxColId primalColIds[],
      SPxRowId dualRowIds[], SPxColId dualColIds[], int* nprimalrows, int* nprimalcols, int* ndualrows, int* ndualcols)
{
   assert(false);
   MSG_ERROR( std::cerr << "Method buildDualProblem() not implemented for Rational\n" );
}


// ---------------------------------------------------------------------------------------------------------------------
//  Explicit instantiation
// ---------------------------------------------------------------------------------------------------------------------

template class SPxLPBase < Rational >;
} // namespace soplex
#endif
