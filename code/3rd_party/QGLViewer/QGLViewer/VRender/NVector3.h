#ifndef _VRENDER_NVECTOR3_H
#define _VRENDER_NVECTOR3_H

#include <iostream>
#include <stdexcept>

namespace vrender
{
  class Vector3;

  class NVector3
  {
  public:
    NVector3();
    NVector3(const NVector3& u);
    inline NVector3(double x,double y,double z,bool normalization=true)
    {
      setXYZ(x,y,z,normalization);
    }

    NVector3(const Vector3 &u,bool normalization=true);

    inline double x() const {return _n[0];}
    inline double y() const {return _n[1];}
    inline double z() const {return _n[2];}
    void setXYZ(double x,double y,double z,bool normalization=true);

    NVector3& operator=(const NVector3& u);
    /*
      inline friend bool operator==(const NVector3 &u,const Vector3  &v) {return u.isEqualTo(v);}
      inline friend bool operator==(const Vector3  &u,const NVector3 &v) {return v.isEqualTo(u);}
      inline friend bool operator==(const NVector3 &u,const NVector3 &v) {return u.isEqualTo(v);}
      inline friend bool operator!=(const NVector3 &u,const Vector3  &v) {return !(u == v);}
      inline friend bool operator!=(const Vector3  &u,const NVector3 &v) {return !(u == v);}
      inline friend bool operator!=(const NVector3 &u,const NVector3 &v) {return !(u == v);}
    */

    inline friend NVector3 operator-(const NVector3 &u) { return NVector3(-u[0],-u[1],-u[2],false); }
    //inline friend Vector3 operator+(const NVector3 &u,const Vector3  &v);
    //inline friend Vector3 operator+(const Vector3  &u,const NVector3 &v);
    //inline friend Vector3 operator+(const NVector3 &u,const NVector3 &v);
    //inline friend Vector3 operator-(const NVector3 &u,const Vector3  &v);
    //inline friend Vector3 operator-(const Vector3  &u,const NVector3 &v);
    //inline friend Vector3 operator-(const NVector3 &u,const NVector3 &v);
    friend double operator*(const NVector3 &u,const Vector3  &v);
    friend double operator*(const Vector3  &u,const NVector3 &v);
    //inline friend double operator*(const NVector3 &u,const NVector3 &v);
    //inline friend Vector3 operator*(double r,const NVector3 &u);
    //inline friend Vector3 operator/(const NVector3 &u,double r);

    //inline friend Vector3 operator^(const NVector3 &u,const Vector3  &v);
    //inline friend Vector3 operator^(const Vector3  &u,const NVector3 &v);
    //inline friend Vector3 operator^(const NVector3 &u,const NVector3 &v);

    inline double norm() const {return 1.0;}
    inline double squareNorm() const {return 1.0;}
    friend std::ostream& operator<<(std::ostream &out,const NVector3 &u);

    double operator[](int i) const
    {
      if((i < 0)||(i > 2))
	throw std::runtime_error("Out of bounds in NVector3::operator[]") ;

      return _n[i];
    }

  private:
    void normalize();

    double _n[3];  //!< normalized vector

  }; // interface of NVector3

}

#endif // _NVECTOR3_H
