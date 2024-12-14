#ifndef _VRENDER_PARSERGL_H
#define _VRENDER_PARSERGL_H

//  This class implements the conversion from OpenGL feedback buffer into more
// usable data structures such as points, segments, and polygons (See Primitive.h)

#include <vector>
#include "Primitive.h"

namespace vrender
{
	class ParserGL
	{
		public:
			void parseFeedbackBuffer(	GLfloat *,
												int size,
												std::vector<PtrPrimitive>& primitive_tab,
												VRenderParams& vparams) ;
			void printStats() const ;

			inline GLfloat xmin() const { return _xmin ; }
			inline GLfloat ymin() const { return _ymin ; }
			inline GLfloat zmin() const { return _zmin ; }
			inline GLfloat xmax() const { return _xmax ; }
			inline GLfloat ymax() const { return _ymax ; }
			inline GLfloat zmax() const { return _zmax ; }
		private:
			int nb_lines ;
			int nb_polys ;
			int nb_points ;
			int nb_degenerated_lines ;
			int nb_degenerated_polys ;
			int nb_degenerated_points ;

			GLfloat _xmin ;
			GLfloat _ymin ;
			GLfloat _zmin ;
			GLfloat _xmax ;
			GLfloat _ymax ;
			GLfloat _zmax ;
	};
}

#endif
