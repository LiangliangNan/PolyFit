/* GLPK configuration file (Microsoft Visual Studio Express) */

#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
#define __WOE__ 1


#define TLS __declspec(thread)
/* thread local storage-class specifier for reentrancy */

#define ODBC_DLNAME "odbc32.dll"
/* ODBC shared library name if this feature is enabled */

#if 0
#define MYSQL_DLNAME "libmysql.dll"
/* MySQL shared library name if this feature is enabled */
#endif

#endif
/* eof */
