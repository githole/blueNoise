#ifndef _TRIANGLE_MESH_H_
#define _TRIANGLE_MESH_H_

#include <map>

#include "..\vec3.h"
#include "..\constant.h"
#include "bbox.h"
#include "intersection.h"
#include "qbvh.h"

#include "..\memfile.h"
#include "..\textutil.h"

namespace hstd {

namespace rt {

struct Material {
	Float3 diffuse;
	Float3 specular;
	float specular_coefficient;
	float metalic;

	Material(Float3& diffuse, Float3& specular, float specular_coefficient, float metalic) :
	diffuse(diffuse), specular(specular), specular_coefficient(specular_coefficient), metalic(metalic) {}
};

typedef std::map<std::string, Material> MaterialMap;

struct Triangle {
	Int3 v_index;
	Int3 vt_index;
	Int3 vn_index;
	Material* material;
};

struct MeshBody {
	std::vector<Float3> v;
	std::vector<Float3> vt;
	std::vector<Float3> vn;
	std::vector<Triangle> triangle;
	MaterialMap matmap;
};

struct TriangleElement {
	Float3* v[3];
	Float3* vt[3];
	Float3* vn[3];
	Material* material;

	TriangleElement() {
		v[0] = v[1] = v[2] = NULL;
		vn[0] = vn[1] = vn[2] = NULL;
		vt[0] = vt[1] = vt[2] = NULL;
		material = NULL;
	}
};

class TriangleMesh {
private:
	MeshBody body_;
	QBVH qbvh_;
public:
	void set(const MeshBody& body) {
		body_ = body;
	}
	MeshBody& get_body() { return body_; }

	void build(std::vector<RefTriangle>& ref_triangle) {
		qbvh_.build(ref_triangle);
	}

	TriangleElement getTriangle(const int triangle_index) {
		TriangleElement t;
		if (triangle_index >= body_.triangle.size())
			return t;

		Triangle& now_t = body_.triangle[triangle_index];
		for (int i = 0; i < 3; ++i) {
			t.v[i] = &body_.v[now_t.v_index[i]];
			if (0 <= now_t.vn_index[i] && now_t.vn_index[i] < body_.vn.size())
				t.vn[i] = &body_.vn[now_t.vn_index[i]];
			if (0 <= now_t.vt_index[i] && now_t.vt_index[i] < body_.vt.size())
				t.vt[i] = &body_.vt[now_t.vt_index[i]];
		}
		t.material = now_t.material;

		return t;
	}

	virtual bool intersect(const Ray &ray, Hitpoint* hitpoint) {
		*hitpoint = Hitpoint();
		return qbvh_.intersect(ray, hitpoint);
	}
};


} // namespace rt

} // namespace hstd

#endif // _TRIANGLE_MESH_H_