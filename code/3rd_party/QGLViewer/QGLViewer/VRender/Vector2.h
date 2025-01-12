#ifndef _VRENDER_VECTOR2_H
#define _VRENDER_VECTOR2_H

#include <stdexcept>
#include <iostream>

namespace vrender
{
  class Vector3;

  class Vector2
	{
		public:
			// ---------------------------------------------------------------------------
			//! @name Constant
			//@{
			static const Vector2 inf;
			//@}

			// ---------------------------------------------------------------------------
			//! @name Constructor(s) and destructor
			//@{
			Vector2 ();
			~Vector2 ();
			Vector2 (const Vector2&);
			Vector2 (const Vector3& u);
			Vector2 (double,double);
			//@}

			// ---------------------------------------------------------------------------
			//! @name Access methods
			//@{
			inline double  x() const { return _xyz[0]; }
			inline double  y() const { return _xyz[1]; }
			inline void  setX(double r) { _xyz[0] = r; }
			inline void  setY(double r) { _xyz[1] = r; }
			inline void  setXY (double x,double y) { _xyz[0] = x; _xyz[1] = y; }
			//@}

			// ---------------------------------------------------------------------------
			//! @name Assignment
			//@{
			inline Vector2& operator= (const Vector2& u)  { _xyz[0] = u._xyz[0]; _xyz[1] = u._xyz[1]; return *this; }
			//@}

			// ---------------------------------------------------------------------------
			//! @name Comparisons
			//@{
			friend bool operator== (const Vector2&,const Vector2&);
			friend bool operator!= (const Vector2&,const Vector2&);
			//@}

			// ---------------------------------------------------------------------------
			//! @name Algebraic operations
			//@{
			inline Vector2& operator+= (const Vector2& v)
			{
				_xyz[0] += v._xyz[0];
				_xyz[1] += v._xyz[1];
				return *this;
			}

			inline Vector2& operator-= (const Vector2& v)
			{
				_xyz[0] -= v._xyz[0];
				_xyz[1] -= v._xyz[1];
				return *this;
			}

			inline Vector2& operator*= (double f) { _xyz[0] *= f; _xyz[1] *= f; return *this;}
			inline Vector2& operator/= (double f) { _xyz[0] /= f; _xyz[1] /= f; return *this;}

			friend Vector2 operator- (const Vector2&);
			static Vector2 mini(const Vector2&,const Vector2&) ;
			static Vector2 maxi(const Vector2&,const Vector2&) ;

			inline Vector2 operator+(const Vector2& u) const
			{
				return Vector2(_xyz[0]+u._xyz[0],_xyz[1]+u._xyz[1]);
			}
			inline Vector2 operator-(const Vector2& u) const
			{
				return Vector2(_xyz[0]-u._xyz[0],_xyz[1]-u._xyz[1]);
			}

			inline double    operator*(const Vector2& u) const
			{
				return _xyz[0]*u._xyz[0] + _xyz[1]*u._xyz[1] ;
			}

			inline double operator^(const Vector2& v) const
			{
				return _xyz[0]*v._xyz[1] - _xyz[1]*v._xyz[0] ;
			}

			Vector2 operator/ (double v) { return Vector2(_xyz[0]/v,_xyz[1]/v); }
			Vector2 operator* (double v) { return Vector2(_xyz[0]*v,_xyz[1]*v); }

			friend Vector2 operator* (double,const Vector2&);
			//@}

			// ---------------------------------------------------------------------------
			//! @name Metrics
			//@{
			double norm       () const;
			double squareNorm () const;
			double infNorm    () const; /// Should be used for most comparisons, for efficiency reasons.
			//@}
			// ---------------------------------------------------------------------------
			//! @name Stream overrides
			//@{
			friend std::ostream& operator<< (std::ostream&,const Vector2&);
			//@}

			double  operator[] (int i) const
			{
				if((i < 0)||(i > 1))
					throw std::runtime_error("Out of bounds in Vector2::operator[]") ;

				return _xyz[i];
			}

			double& operator[] (int i)
			{
				if((i < 0)||(i > 1))
					throw std::runtime_error("Out of bounds in Vector2::operator[]") ;

				return _xyz[i];
			}

		private:
			double _xyz[2];  //!< The 3 vector components

	}; // interface of Vector2
}

#endif // _VECTOR2_H
