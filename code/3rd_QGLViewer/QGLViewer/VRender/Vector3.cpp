#include <iostream>
#include "Vector3.h"
#include "NVector3.h"
#include <math.h>
#include <algorithm>

#ifdef WIN32
# include <windows.h>
#endif

using namespace vrender;
using namespace std;

const Vector3 Vector3::inf(FLT_MAX, FLT_MAX, FLT_MAX); 

Vector3::Vector3 ()
{
  _xyz[0] = 0.0;
  _xyz[1] = 0.0;
  _xyz[2] = 0.0;
}

// -----------------------------------------------------------------------------
//! Default destructor
Vector3::~Vector3 ()
{
}

// -----------------------------------------------------------------------------
//! Copy constructor
Vector3::Vector3 (const Vector3& u)
{
  setXYZ(u[0],u[1],u[2]);
}

// -----------------------------------------------------------------------------
//! Copy constructor from a normalized vector
Vector3::Vector3 (const NVector3& u)
{
  setXYZ(u[0],u[1],u[2]);
}

// -----------------------------------------------------------------------------
//! Create a vector from real values
Vector3::Vector3 (double x,double y,double z)
{
  setXYZ(x,y,z);
}
// -----------------------------------------------------------------------------
//! Assignment with a normalized vector
Vector3& Vector3::operator= (const NVector3& u)
{
  _xyz[0] = u[0];
  _xyz[1] = u[1];
  _xyz[2] = u[2];
  return ( *this );
}

// -----------------------------------------------------------------------------
//! Self addition with a normalized vector
Vector3& Vector3::operator+= (const NVector3& u)
{
  _xyz[0] += u[0];
  _xyz[1] += u[1];
  _xyz[2] += u[2];
  return ( *this );
}

// -----------------------------------------------------------------------------
//! Self substraction with a normalized vector
Vector3& Vector3::operator-= (const NVector3& u)
{
  _xyz[0] -= u[0];
  _xyz[1] -= u[1];
  _xyz[2] -= u[2];
  return ( *this );
}

// -----------------------------------------------------------------------------
//! Left multiplication by a real value
Vector3 vrender::operator* (double r,const Vector3& u)
{
  return ( Vector3(r*u[0], r*u[1], r*u[2]) );
}


// -----------------------------------------------------------------------------
//! Norm
double Vector3::norm () const
{
  return sqrt( _xyz[0]*_xyz[0] + _xyz[1]*_xyz[1] + _xyz[2]*_xyz[2] );
}

// -----------------------------------------------------------------------------
//! Square norm (self dot product)
double Vector3::squareNorm () const
{
  return _xyz[0]*_xyz[0] + _xyz[1]*_xyz[1] + _xyz[2]*_xyz[2];
}

// -----------------------------------------------------------------------------
//! Infinite norm
double Vector3::infNorm() const
{
  return std::max(std::max(fabs(_xyz[0]),fabs(_xyz[1])),fabs(_xyz[2])) ;
}


// -----------------------------------------------------------------------------
//! Out stream override: prints the 3 vector components
std::ostream& vrender::operator<< (std::ostream& out,const Vector3& u)
{
  out << u[0] << " " << u[1] << " " << u[2];
  return ( out );
}

Vector3 Vector3::mini(const Vector3& v1,const Vector3& v2)
{
  return Vector3(std::min(v1[0],v2[0]),std::min(v1[1],v2[1]),std::min(v1[2],v2[2])) ;
}

Vector3 Vector3::maxi(const Vector3& v1,const Vector3& v2)
{
  return Vector3(std::max(v1[0],v2[0]),std::max(v1[1],v2[1]),std::max(v1[2],v2[2])) ;
}

