#ifndef _TRIANGLE_H_
#define _TRIANGLE_H_

#include "..\vec3.h"
#include "..\constant.h"
#include "bbox.h"

namespace hstd {

namespace rt {

// どこか別のメモリ上に存在する頂点データに対する参照を保持する三角形データ構造
class RefTriangle {
public:
	Float3 *p[3];
	int original_triangle_index;

	RefTriangle(Float3 *p1, Float3 *p2, Float3 *p3, int original_triangle_index) :
	 original_triangle_index(original_triangle_index) {
		 p[0] = p1;
		 p[1] = p2;
		 p[2] = p3;
	}
	BBox objectBound() const {
		return unionBBox(BBox(*p[0], *p[1]), *p[2]);
	}
};

} // namespace rt

} // namespace hstd

#endif