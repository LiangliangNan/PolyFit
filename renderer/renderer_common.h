
#ifndef _RENDERER_COMMON_H_
#define _RENDERER_COMMON_H_


#include "../basic/basic_common.h"


# ifdef RENDERER_EXPORTS
#   define RENDERER_API  EXPORT_LIBRARY
# else
#   define RENDERER_API  IMPORT_LIBRARY
# endif


#endif