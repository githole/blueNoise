#ifndef _HMATH_H_
#define _HMATH_H_

#include <cmath>

#include "vec3.h"

namespace hstd {
	
template<typename T>
inline T clamp(T x, T min_value, T max_value) {
	if (x < min_value)
		return min_value;
	if (max_value < x)
		return max_value;
	return x;
}

template<typename T>
inline T saturate(T x) {
	return clamp(x, 0, 1);
}

template<typename T>
inline T lerp(T x, T y, T s) {
	return x + s * (y - x);
}

template<typename T>
inline T smoothstep (T edge0, T edge1, T x)
{
	x = saturate((x - edge0) / (edge1 - edge0)); 
	return x*x*(3-2*x);
}

/// <summary>
/// theta: [0, pi], phi: [0, 2pi]
/// </summary>
template<typename T>
inline void directionToPolarCoordinate(const Vec3<T>& dir, T *theta, T *phi) {
	*theta = acos(dir.y);
	*phi = atan2(dir.z, dir.x);
	if (*phi < 0)
		*phi += 2.0f * kPI;
}

// Y-up
template<typename T>
inline Vec3<T> polarCoordinateToDirection(T theta, T phi) {
	return Vec3<T>(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));
}

template<typename T>
inline Vec3<T> reflect(const Vec3<T> &in, const Vec3<T> &normal) {
	return normalize(in - normal * 2.0 * dot(normal, in));
}

}

#endif // _HMATH_H_