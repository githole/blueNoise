#ifndef _SPHERE_H_
#define _SPHERE_H_

#include "..\vec3.h"
#include "..\constant.h"
#include "intersection.h"
#include "bbox.h"

namespace hstd {

namespace rt {

const float kEPS = 1e-6f;

struct Sphere {
	float radius;
	Float3 position;

	Sphere(const float radius, const Float3 &position) :
	radius(radius), position(position) {}

	// 入力のrayに対する交差点までの距離を返す。交差しなかったら0を返す。
	// rayとの交差判定を行う。交差したらtrue,さもなくばfalseを返す。
	bool intersect(const Ray &ray, Hitpoint *hitpoint) const {
		const Float3 p_o = position - ray.org;
		const float b = dot(p_o, ray.dir);
		const float D4 = b * b - dot(p_o, p_o) + radius * radius;

		if (D4 < 0.0)
			return false;
		
		const float sqrt_D4 = sqrt(D4);
		const float t1 = b - sqrt_D4, t2 = b + sqrt_D4;
	
		if (t1 < kEPS && t2 < kEPS)
			return false;

		if (t1 > kEPS) {
			hitpoint->distance = t1;
		} else {
			hitpoint->distance = t2;
		}
		return true;
	}
};


} // namespace rt

} // namespace hstd

#endif