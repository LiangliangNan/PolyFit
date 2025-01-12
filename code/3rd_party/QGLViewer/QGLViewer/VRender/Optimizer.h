#ifndef _OPTIMIZER_H
#define _OPTIMIZER_H

#include "Types.h"

namespace vrender
{
	// Implements some global optimizations on the polygon sorting.

	class VRenderParams ;
	class Optimizer
	{
		public:
			virtual void optimize(std::vector<PtrPrimitive>&,VRenderParams&) = 0 ;
			virtual ~Optimizer() {} ;
	};

	//  Optimizes visibility by culling primitives which do not appear in the
	// rendered image. Computations are done analytically rather than using an item
	// buffer.

	class VisibilityOptimizer: public Optimizer
	{
		public:
			virtual void optimize(std::vector<PtrPrimitive>&,VRenderParams&) ;
			virtual ~VisibilityOptimizer() {} ;
	};

	//  Optimizes by collapsing together primitives which can be, without
	// perturbating the back to front painting algorithm.

	class PrimitiveSplitOptimizer: public Optimizer
	{
		public:
			virtual void optimize(std::vector<PtrPrimitive>&,VRenderParams&) {}
			virtual ~PrimitiveSplitOptimizer() {} ;
	};

	class BackFaceCullingOptimizer: public Optimizer
	{
		public:
			virtual void optimize(std::vector<PtrPrimitive>&,VRenderParams&) ;
			virtual ~BackFaceCullingOptimizer() {} ;
	};
}

#endif
