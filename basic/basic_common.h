
#ifndef _BASIC_COMMON_H_
#define _BASIC_COMMON_H_

// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the BASIC_EXPORTS symbol.
// This symbol should not be defined on any project that uses this DLL. This way any other
// project whose source files include this file see BASIC_EXPORTS functions as being imported 
// from a DLL, whereas this DLL sees symbols defined with this macro as being exported.


// Windows DLL export macros
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
// Disable a warning message about dll. This is a temporary solution
// http://support.microsoft.com/default.aspx?scid=kb;EN-US;168958
#pragma warning( disable : 4251 )
#pragma warning( disable : 4996 )  
#pragma warning( disable : 4267 ) 
#pragma warning( disable : 4091 ) 
#pragma warning( disable : 4005 ) 
#pragma warning( disable : 4244 ) 
#pragma warning( disable : 4101 ) 
#endif


#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
#    ifdef STATIC_LIB
#        define EXPORT_LIBRARY
#        define IMPORT_LIBRARY
#    else
#        define EXPORT_LIBRARY __declspec( dllexport )
#        define IMPORT_LIBRARY __declspec( dllimport )
#    endif
#else
#    define EXPORT_LIBRARY
#    define IMPORT_LIBRARY
#endif



# ifdef BASIC_EXPORTS
#   define BASIC_API  EXPORT_LIBRARY
# else
#   define BASIC_API  IMPORT_LIBRARY
# endif


//_______________________________________________

// 32 or 64 bit environment

#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64) // Check windows
#	if defined (_WIN64) || defined (WIN64)
#		define ENV_64_BIT
#	else
#		define ENV_32_BIT
#	endif
#elif defined (__GNUC__) //	Check GCC
#	if defined (__x86_64__) || defined (__ppc64__)
#		define ENV_64_BIT
#	else
#		define ENV_32_BIT
#	endif
#else
#	error "Error: unknown environment."
#endif




#endif // _BASIC_COMMON_H_