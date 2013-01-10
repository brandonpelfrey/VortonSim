#ifndef Vec3_h_
#define Vec3_h_

#include <cmath>

struct Vec3;
Vec3 operator*(float r, const Vec3& v);

struct Vec3 {
 union {
  struct {
   float x,y,z;
  };
  struct {
   float r,g,b;
  };
  float D[3];
 };

	Vec3() { }
 Vec3(float _x, float _y, float _z)
  :x(_x), y(_y), z(_z)
  { }

 float& operator[](unsigned int i) {
  return D[i];
 }

 const float& operator[](unsigned int i) const {
  return D[i];
 }

 float maxComponent() const {
	float r = x;
	if(y>r) r = y;
	if(z>r) r = z;
	return r;
	}
 
	float minComponent() const {
	float r = x;
	if(y<r) r = y;
	if(z<r) r = z;
	return r;
	}

 Vec3 operator+(const Vec3& r) const {
  return Vec3(x+r.x, y+r.y, z+r.z); 
 }

 Vec3 operator-(const Vec3& r) const {
  return Vec3(x-r.x, y-r.y, z-r.z); 
 }

 Vec3 cmul(const Vec3& r) const {
  return Vec3(x*r.x, y*r.y, z*r.z);
 }

 Vec3 cdiv(const Vec3& r) const {
  return Vec3(x/r.x, y/r.y, z/r.z);
 }

 Vec3 operator*(float r) const {
  return Vec3(x*r,y*r,z*r);
 }


 Vec3 operator/(float r) const {
  return Vec3(x/r, y/r, z/r);
 }

 ////
 Vec3& operator+=(const Vec3& r) {
  x+=r.x;
  y+=r.y;
  z+=r.z;
  return *this;
 }

 Vec3& operator-=(const Vec3& r) {
  x-=r.x;
  y-=r.y;
  z-=r.z;
  return *this;
 }

 Vec3& operator*=(float r) {
  x*=r; y*=r; z*=r;
  return *this;
 }

 float operator*(const Vec3& r) const {
  return x*r.x + y*r.y + z*r.z;
 }

 float norm() const {
  return sqrtf(x*x+y*y+z*z);
 }

 float normSquared() const {
  return x*x + y*y + z*z;
 }

 Vec3 operator^(const Vec3& r) const {
  return Vec3(
   y * r.z - z * r.y, 
   z * r.x - x * r.z, 
   x * r.y - y * r.x
   );
 }

 Vec3 normalized() const {
  return *this / norm();
 }

 Vec3 reflect(const Vec3& normal) const {
  return *this - (2.f * (*this * normal)) * normal;
 }

 void tangentSpace(Vec3* a, Vec3* b) const {
  const Vec3 special(0.577350269f, 0.577350269f, 0.577350269f); // 1/sqrt(3)
  *a = special ^ normalized();
	*a = a->normalized();
  *b = *a ^ *this;
	*b = b->normalized();
 }
};
 
//*
inline Vec3 operator*(float r, const Vec3& v) {
 return Vec3(v.x*r, v.y*r, v.z*r);
}

inline Vec3 min(const Vec3& a, const Vec3& b) {
	Vec3 r = a;
	if(b.x<r.x) r.x = b.x;
	if(b.y<r.y) r.y = b.y;
	if(b.z<r.z) r.z = b.z;
	return r;
}

inline Vec3 max(const Vec3& a, const Vec3& b) {
	Vec3 r = a;
	if(b.x>r.x) r.x = b.x;
	if(b.y>r.y) r.y = b.y;
	if(b.z>r.z) r.z = b.z;
	return r;
}

#endif
