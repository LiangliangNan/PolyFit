#ifndef _VRENDER_TYPES_H
#define _VRENDER_TYPES_H

#ifdef WIN32
# include <windows.h>
#endif

#ifdef __APPLE__
# include <OpenGL/gl.h>
#else
# include <GL/gl.h>
#endif

namespace vrender
{
	typedef double FLOAT ;
	typedef GLdouble GLFLOAT ;

#ifdef A_VOIR
	typedef T_Vect3<double> DVector3 ;
	typedef T_Vect2<double> Vector2 ;
#endif

	class Primitive ;
	typedef Primitive *PtrPrimitive ;

	const float FLAT_POLYGON_EPS = 1e-5f ;
}

#endif
