#ifndef _SORTMETHOD_H
#define _SORTMETHOD_H

#include <vector>
#include "Types.h"

namespace vrender
{
	// Class which implements the sorting of the primitives. An object of
	class VRenderParams ;
	class SortMethod
	{
		public:
			SortMethod() {}
			virtual ~SortMethod() {}

			virtual void sortPrimitives(std::vector<PtrPrimitive>&,VRenderParams&) = 0 ;

			void SetZDepth(FLOAT s) { zSize = s ; }
			FLOAT ZDepth() const { return zSize ; }

		protected:
			FLOAT zSize ;
	};

	class DontSortMethod: public SortMethod
	{
		public:
			DontSortMethod() {}
			virtual ~DontSortMethod() {}

			virtual void sortPrimitives(std::vector<PtrPrimitive>&,VRenderParams&) {}
	};

	class BSPSortMethod: public SortMethod
	{
		public:
			BSPSortMethod() {} ;
			virtual ~BSPSortMethod() {}

			virtual void sortPrimitives(std::vector<PtrPrimitive>&,VRenderParams&) ;
	};

	class TopologicalSortMethod: public SortMethod
	{
		public:
			TopologicalSortMethod() ;
			virtual ~TopologicalSortMethod() {}

			virtual void sortPrimitives(std::vector<PtrPrimitive>&,VRenderParams&) ;

			void setBreakCycles(bool b) { _break_cycles = b ; }
		private:
			bool _break_cycles ;
	};
}

#endif
