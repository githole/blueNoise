#ifndef _SAMPLING_H_
#define _SAMPLING_H_

#include "random.h"
#include "constant.h"

namespace hstd {


template <class Random, typename T>
class GenericSampling {
public:
	static void uniformCircle(Random &random, T *x, T *y) {
		const T sqrt_r= sqrt(random.next01());
		const T theta = random.next(0.0f, 2.0f * kPI);

		*x = sqrt_r * cos(theta);
		*y = sqrt_r * sin(theta);
	}
	
	static Vec3<T> uniformSphereSurface(Random &random) {
		const T tz = random.next(-1.0f, 1.0f);
		const T phi = random.next(0.0f, 2.0f * kPI);
		const T k = sqrt(1.0f - tz * tz);
		const T tx = k * cos(phi);
		const T ty = k * sin(phi);

		return Vec3<T>(tx, ty, tz);
	}
	
	static Vec3<T> uniformHemisphereSurface(Random &random, const Vec3<T> &normal, const Vec3<T> &tangent, const Vec3<T> &binormal) {
		const T tz = random.next(0.0f, 1.0f);
		const T phi = random.next(0.0f, 2.0f * kPI);
		const T k = sqrt(1.0f - tz * tz);
		const T tx = k * cos(phi);
		const T ty = k * sin(phi);

		return tz * normal + tx * tangent + ty * binormal;
	}

	static Vec3<T> cosineWeightedHemisphereSurface(Random &random, const Vec3<T> &normal, const Vec3<T> &tangent, const Vec3<T> &binormal) {
		const T phi = random.next(0.0f, 2.0f * kPI);
		const T r2 = random.next01(), r2s = sqrt(r2);

		const T tx = r2s * cos(phi);
		const T ty = r2s * sin(phi);
		const T tz = sqrt(1.0f - r2);

		return tz * normal + tx * tangent + ty * binormal;
	}
};

typedef GenericSampling<Random, float> Sampling;

}

#endif // _SAMPLING_H_