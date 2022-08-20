#ifdef SOPLEX_WITH_ZLIB

// ============================================================================
// gzstream, C++ iostream classes wrapping the zlib compression library.
// Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// File          : gzstream.h
// Revision      : $Revision: 1.8 $
// Revision_date : $Date: 2005/11/09 13:53:50 $
// Author(s)     : Deepak Bandyopadhyay, Lutz Kettner
// 
// Standard streambuf implementation following Nicolai Josuttis, "The 
// Standard C++ Library".
// ============================================================================

/**@file gzstream.h
 * @brief Utilities for handling gzipped input and output streams.
 */
#ifndef GZSTREAM_H
#define GZSTREAM_H 1

// standard C++ with new header file names and std:: namespace
#include <iostream>
#include <fstream>
#include <zlib.h>

#define GZSTREAM_NAMESPACE gzstream

#ifdef GZSTREAM_NAMESPACE
namespace GZSTREAM_NAMESPACE {
#endif

// ----------------------------------------------------------------------------
//    Internal classes to implement gzstream. See below for user classes.
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
//    class gzstreambuf
// ----------------------------------------------------------------------------

/**@class gzstreambuf
   @brief Internal class to implement gzstream.
*/
class gzstreambuf 
   : public std::streambuf 
{
private:

   //------------------------------------
   /**@name Types */
   //@{
   ///
   static const int bufferSize = 47+256;   ///< size of data buff
   // totals 512 bytes under g++ for igzstream at the end.
   //@}

   //------------------------------------
   /**@name Data */
   //@{
   gzFile           file;               ///< file handle for compressed file
   char             buffer[bufferSize]; ///< data buffer
   char             opened;             ///< open/close state of stream
   unsigned int     mode;               ///< I/O mode
   //@}

   //------------------------------------
   /**@name Internal helpers */
   //@{
   ///
   int flush_buffer();
   //@}

public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   gzstreambuf()
      : file(0)
      , opened(0)
      , mode(0)
   {
      setp( buffer, buffer + (bufferSize-1));
      setg( buffer + 4,     // beginning of putback area
            buffer + 4,     // read position
            buffer + 4);    // end position      
      // ASSERT: both input & output capabilities will not be used together
   }
   /// destructor
   ~gzstreambuf() 
   { 
      close(); 
   }
   //@}

   //------------------------------------
   /**@name Interface */
   //@{
   ///
   int is_open() 
   { 
      return opened; 
   }
   ///
   gzstreambuf* open( const char* name, int open_mode );
   ///
   gzstreambuf* close();
   ///
   virtual int     overflow( int c = EOF );
   ///
   virtual int     underflow();
   ///
   virtual int     sync();
   //@}
};

// ----------------------------------------------------------------------------
//    class gzstreambase
// ----------------------------------------------------------------------------

/**@class gzstreambase
   @brief Internal class to implement gzstream.
*/
class gzstreambase 
   : virtual public std::ios 
{
protected:

   //------------------------------------
   /**@name Data */
   //@{
   ///
   gzstreambuf buf;
   //@}

public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   gzstreambase() 
   { 
      init(&buf); 
   }
   /// full constructor
   gzstreambase( const char* _name, int _open_mode );
   /// destructor
   ~gzstreambase();
   //@}

   //------------------------------------
   /**@name Interface */
   //@{
   ///
   void open( const char* _name, int _open_mode );
   ///
   void close();
   ///
   gzstreambuf* rdbuf() 
   { 
      return &buf; 
   }
   //@}
};

// ----------------------------------------------------------------------------
// User classes. Use igzstream and ogzstream analogously to ifstream and
// ofstream respectively. They read and write files based on the gz* 
// function interface of the zlib. Files are compatible with gzip compression.
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
//    class igzstream
// ----------------------------------------------------------------------------

/**@class igzstream
   @brief Class to implement a gzip'd input stream.
*/
class igzstream 
   : public std::istream
   , public gzstreambase 
{
public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   igzstream() 
      : std::istream( &buf) 
   {} 
   /// full constructor
   igzstream( const char*  _name, 
              int          _open_mode = std::ios::in )
      : std::istream( &buf )
      , gzstreambase( _name, _open_mode )
   {}  
   //@}

   //------------------------------------
   /**@name Interface */
   //@{
   ///
   gzstreambuf* rdbuf() 
   { 
      return gzstreambase::rdbuf(); 
   }
   ///
   void open( const char*  _name, 
              int          _open_mode = std::ios::in )
   {
       gzstreambase::open( _name, _open_mode );
   }
   //@}
};

// ----------------------------------------------------------------------------
//    class ogzstream
// ----------------------------------------------------------------------------

/**@class ogzstream
   @brief Class to implement a gzip'd output stream.
*/
class ogzstream 
   : public gzstreambase
   , public std::ostream 
{
public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   ogzstream() 
      : std::ostream( &buf) 
   {}
   /// full constructor
   explicit
   ogzstream( const char* _name, 
              int         _open_mode = std::ios::out )
      : gzstreambase( _name, _open_mode )
      , std::ostream( &buf) 
   {}  
   //@}

   //------------------------------------
   /**@name Interface */
   //@{
   ///
   gzstreambuf* rdbuf() 
   { 
      return gzstreambase::rdbuf(); 
   }
   ///
   void open( const char*  _name, 
              int          _open_mode = std::ios::out ) 
   {
      gzstreambase::open( _name, _open_mode );
   }
};

#ifdef GZSTREAM_NAMESPACE
} // namespace GZSTREAM_NAMESPACE
#endif

#endif // GZSTREAM_H
// ============================================================================
// EOF //

#endif
