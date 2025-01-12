#ifndef _VRENDER_AXISALIGNEDBOX_H
#define _VRENDER_AXISALIGNEDBOX_H

namespace vrender
{
  class Vector2;
  class Vector3;

	template<class T> class AxisAlignedBox
	{
		public:
			AxisAlignedBox() ;
			AxisAlignedBox(const T& v) ;
			AxisAlignedBox(const T& v,const T& w) ;

			const T& mini() const { return _min ; }
			const T& maxi() const { return _max ; }

			void include(const T& v) ;
			void include(const AxisAlignedBox<T>& b) ;
		private:
			T _min ;
			T _max ;
	};

	typedef AxisAlignedBox< Vector2 > AxisAlignedBox_xy ;
	typedef AxisAlignedBox< Vector3 > AxisAlignedBox_xyz ;

	template<class T> AxisAlignedBox<T>::AxisAlignedBox()
	: _min(T::inf), _max(-T::inf)
	{
	}

	template<class T> AxisAlignedBox<T>::AxisAlignedBox(const T& v)
		: _min(v), _max(v)
	{
	}

	template<class T> AxisAlignedBox<T>::AxisAlignedBox(const T& v,const T& w)
		: _min(v), _max(v)
	{
		include(w) ;
	}

	template<class T> void AxisAlignedBox<T>::include(const T& v)
	{
		_min = T::mini(_min,v) ;
		_max = T::maxi(_max,v) ;
	}

	template<class T> void AxisAlignedBox<T>::include(const AxisAlignedBox<T>& b)
	{
		include(b._min) ;
		include(b._max) ;
	}
}
#endif
