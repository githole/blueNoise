#ifndef	_VEC_3_H_
#define	_VEC_3_H_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <assert.h>

namespace hstd {

template<typename T>
struct Vec3 {
	typedef T ValueType;

	T x, y, z;
	explicit Vec3(const T x, const T y, const T z) : x(x), y(y), z(z) {}
	explicit Vec3() : x(0), y(0), z(0) {}
	
	T operator[](const std::size_t i) const {
#ifdef _DEBUG
		assert(i < 3);
#endif // _DEBUG
		return (&x)[i];
	}
	T &operator[](const std::size_t i) {
#ifdef _DEBUG
		assert(i < 3);
#endif // _DEBUG
		return (&x)[i];
	}
	
	template <class E>
	Vec3& operator=(const E& r) {
		for (std::size_t i = 0; i < 3; ++i)
			((*this)[i]) = r[i];
		return *this;
	}
	
	template <class E>
	Vec3(const E &r) {
		for (std::size_t i = 0; i < 3; ++i)
			((*this)[i]) = r[i];
	}
	
	inline Vec3& operator+=(const Vec3 &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	inline Vec3& operator-=(const Vec3 &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	
	inline Vec3& operator*=(const T b) {
		x *= b; y *= b; z *= b;
		return *this;
	}
	
	inline Vec3 operator+(const Vec3 &b) const {
		return Vec3(x + b.x, y + b.y, z + b.z);
	}
	inline Vec3 operator-(const Vec3 &b) const {
		return Vec3(x - b.x, y - b.y, z - b.z);
	}
	inline Vec3 operator*(const T b) const {
		return Vec3(x * b, y * b, z * b);
	}
	inline Vec3 operator/(const T b) const {
		return Vec3(x / b, y / b, z / b);
	}
	inline const float lengthSquared() const { 
		return x*x + y*y + z*z; 
	}
	inline const float length() const { 
		return sqrt(lengthSquared()); 
	}
	
	Vec3 operator-() const { 
		return Vec3(-x, -y, -z);
	}
};

template <typename T>
inline Vec3<T> operator*(T f, const Vec3<T> &v) { 
	return v * f; 
}

template <typename T>
inline Vec3<T> normalize(const Vec3<T> &v) {
	return v * (1 / v.length()); 
}

template <typename T>
inline const Vec3<T> times(const Vec3<T> &v1, const Vec3<T> &v2) {
	return Vec3<T>(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

template <typename T>
inline const T dot(const Vec3<T> &v1, const Vec3<T> &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
inline const Vec3<T> cross(const Vec3<T> &v1, const Vec3<T> &v2) {
	return Vec3<T>(
		(v1.y * v2.z) - (v1.z * v2.y),
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x));
}

template <typename T>
std::ostream &operator<<(std::ostream &out, const Vec3<T> &v) {
	out << "<" << v.x << "," << v.y << "," << v.z << ">";
	return out;
}

template <typename T>
inline void createOrthoNormalBasis(const Vec3<T> &normal, Vec3<T> *tangent, Vec3<T> *binormal) {
	if (abs(normal.x) > abs(normal.y))
		*tangent = normalize(cross(Vec3<T>(0, 1, 0), normal));
	else
		*tangent = normalize(cross(Vec3<T>(1, 0, 0), normal));
	*binormal = normalize(cross(normal, *tangent));
}

typedef Vec3<float> Float3;
typedef Vec3<double> Double3;
typedef Vec3<int> Int3;

};

#endif
