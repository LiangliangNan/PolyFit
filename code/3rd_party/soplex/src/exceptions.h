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

/**@file  exceptions.h
 * @brief Exception classes for SoPlex.
 */
#ifndef _EXCEPTIONS_H_
#define _EXCEPTIONS_H_

#include <string.h>

namespace soplex
{
   /**@brief   Exception base class.
    * @ingroup Elementary
    * 
    * This class implements a base class for our SoPlex exceptions
    * We provide a what() function which returns the exception message.
    */
   class SPxException
   {
   private:
      //----------------------------------------
      /**@name Private data */
      //@{
      /// Exception message.
      std::string msg;
      //@}
   public:
      //----------------------------------------
      /**@name Construction / destruction */
      //@{
      /// constructor
      /** The constructor receives an optional string as an exception message.
       */
      SPxException(const std::string& m = "") : msg(m) {}
      /// destructor
      virtual ~SPxException() {}
      //@}
      
      //----------------------------------------
      /**@name Access / modification */
      //@{
      /// returns exception message
      virtual const std::string what() const
      {
         return msg;
      }
      //@}
   };

   /**@brief   Exception class for out of memory exceptions.
    * @ingroup Elementary
    * 
    * This class is derived from the SoPlex exception base class.
    * It does not provide any new functionality.
    */
   class SPxMemoryException : public SPxException 
   {
   public:
     //----------------------------------------
      /**@name Construction / destruction */
      //@{
      /// constructor
      /** The constructor receives an optional string as an exception message.
       */
      SPxMemoryException(const std::string& m = "") : SPxException(m) {}
      //@}
   };

   /**@brief   Exception class for status exceptions during the computations
    * @ingroup Elementary
    * 
    * This class is derived from the SoPlex exception base class.
    * It does not provide any new functionality.
    */
   class SPxStatusException : public SPxException
   {
   public:
     //----------------------------------------
      /**@name Construction / destruction */
      //@{
      /// constructor
      /** The constructor receives an optional string as an exception message.
       */
      SPxStatusException(const std::string& m = "") : SPxException(m) {}
      //@}
   };

   /**@brief   Exception class for things that should NEVER happen.
    * @ingroup Elementary
    * 
    * This class is derived from the SoPlex exception base class.
    * It does not provide any new functionality. Most often it is used to replace
    * assert(false) terms in earlier code.
    */
   class SPxInternalCodeException : public SPxException
   {
   public:
     //----------------------------------------
      /**@name Construction / destruction */
      //@{
      /// constructor
      /** The constructor receives an optional string as an exception message.
       */
      SPxInternalCodeException(const std::string& m = "") : SPxException(m) {}
      //@}
   };


   /**@brief   Exception class for incorrect usage of interface methods.
    * @ingroup Elementary
    */
   class SPxInterfaceException : public SPxException
   {
   public:
     //----------------------------------------
      /**@name Construction / destruction */
      //@{
      /// constructor
      /** The constructor receives an optional string as an exception message.
       */
      SPxInterfaceException(const std::string& m = "") : SPxException(m) {}
      //@}
   };

} //namespace soplex

#endif // _EXCEPTIONS_H_
