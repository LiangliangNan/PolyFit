#ifndef _VRENDER_H_
#define _VRENDER_H_

#include "../config.h"
#include <QTextStream>
#include <QString>

#include "../qglviewer.h"

namespace vrender
{
	class VRenderParams ;
	typedef void (*RenderCB)(void *) ;
	typedef void (*ProgressFunction)(float,const QString&) ;

	void VectorialRender(RenderCB DrawFunc, void *callback_params, VRenderParams& render_params) ;

	class VRenderParams
	{
		public:
			VRenderParams() ;
			~VRenderParams() ;

			enum VRenderSortMethod { NoSorting, BSPSort, TopologicalSort, AdvancedTopologicalSort };
			enum VRenderFormat     { EPS, PS, XFIG, SVG };

			enum VRenderOption {	CullHiddenFaces         = 0x1,
						OptimizeBackFaceCulling = 0x4,
						RenderBlackAndWhite     = 0x8,
						AddBackground           = 0x10,
						TightenBoundingBox      = 0x20 } ;

			int sortMethod()    { return _sortMethod; }
			void setSortMethod(VRenderParams::VRenderSortMethod s) { _sortMethod = s ; }

			int format()        { return _format; }
			void setFormat(VRenderFormat f) { _format = f; }

			const QString filename() { return _filename ; }
			void setFilename(const QString& filename) ;

			void setOption(VRenderOption,bool) ;
			bool isEnabled(VRenderOption) ;

			void setProgressFunction(ProgressFunction pf) { _progress_function = pf ; }

		private:
			int _error;
			VRenderSortMethod _sortMethod;
			VRenderFormat     _format ;

			ProgressFunction _progress_function ;

			unsigned int _options; // _DrawMode; _ClearBG; _TightenBB;
			QString _filename;

			friend void VectorialRender(	RenderCB render_callback,
							void *callback_params,
							VRenderParams& vparams);
			friend class ParserGL ;
			friend class Exporter ;
			friend class BSPSortMethod ;
			friend class VisibilityOptimizer ;
			friend class TopologicalSortMethod ;
			friend class TopologicalSortUtils ;

			int& error() { return _error ; }
			int& size()  { static int size=1000000; return size ; }

			void progress(float,const QString&) ;
	};
}
#endif

