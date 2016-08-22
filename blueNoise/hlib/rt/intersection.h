#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include "..\constant.h"

namespace hstd {

namespace rt {

struct Hitpoint {
	float distance;
	int triangle_index;
	float b1, b2; // barycentric coordinate

//	Float3 normal;
//	Float3 v0, v1, v2;

	Hitpoint() : distance(kINF), triangle_index(-1), b1(-1), b2(-1) {}
};

} // namespace rt

} // namespace hstd

#endif // _INTERSECTION_H_