
#ifndef _METHOD_COMMON_H_
#define _METHOD_COMMON_H_


#include "../basic/basic_common.h"


# ifdef METHOD_EXPORTS
#   define METHOD_API  EXPORT_LIBRARY
# else
#   define METHOD_API  IMPORT_LIBRARY
# endif


#endif