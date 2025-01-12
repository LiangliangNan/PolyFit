#include <vector>
#include "VRender.h"
#include "Optimizer.h"
#include "Primitive.h"

using namespace std ;
using namespace vrender ;

// Over-simplified algorithm to check wether a polygon is front-facing or not.
// Only works for convex polygons.

void BackFaceCullingOptimizer::optimize(std::vector<PtrPrimitive>& primitives_tab,VRenderParams&)
{
	Polygone *P ;
	int nb_culled = 0 ;

	for(size_t i=0;i<primitives_tab.size();++i)
		if((P = dynamic_cast<Polygone *>(primitives_tab[i])) != nullptr)
		{
						for(unsigned int j=0;j<P->nbVertices();++j)
				if(( (P->vertex(j+2) - P->vertex(j+1))^(P->vertex(j+1) - P->vertex(j))).z() > 0.0 )
				{
					delete primitives_tab[i] ;
					primitives_tab[i] = nullptr ;
					++nb_culled ;
					break ;
				}
		}

	// Rule out gaps. This avoids testing for null primitives later.

	int j=0 ;
	for(size_t k=0;k<primitives_tab.size();++k)
		if(primitives_tab[k] != nullptr)
			primitives_tab[j++] = primitives_tab[k] ;

	primitives_tab.resize(j) ;
#ifdef DEBUG_BFC
	cout << "Backface culling: " << nb_culled << " polygons culled." << endl ;
#endif
}
