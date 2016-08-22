#ifndef _OBJ_MESH_
#define _OBJ_MESH_

#include "..\vec3.h"
#include "..\constant.h"
#include "..\memfile.h"
#include "..\textutil.h"

#include "bbox.h"
#include "intersection.h"
#include "triangleMesh.h"
#include "qbvh.h"

namespace hstd {

namespace rt {
class OBJOperator {
private:
	static Int3 face_to_data(const std::string &input) {
		std::string tmp;
		int data[3] = {-1, -1, -1};
		int now_data_index = 0;
		int now_index = 0;
		for (;now_index < input.size(); ++now_index) {
			if (input[now_index] == '/') {
				
				if (tmp != "") {
					data[now_data_index] = atoi(tmp.c_str());
					tmp = "";
				}
				++now_data_index;

				continue;
			}

			tmp += input[now_index];
		}

		if (tmp != "") {
			data[now_data_index] = atoi(tmp.c_str());
			++now_data_index;
		}

		return Int3(data[0], data[1], data[2]);
	}
public:
	static bool load_material(const char *filename, MaterialMap *matmap) {
		FileManager manager;
		if (!manager.load(filename))
			return false;

		std::string now_material_name;
		Float3 diffuse;
		Float3 specular;
		float specular_coefficient = 0;
		float metalic = 0;
		
		for (;;) {
			std::string now_line;
			if (!manager.gets(&now_line))
				break;
			
			if (now_line.size() == 0 || now_line[0] == '#') {
				continue;
			}

			std::vector<std::string> ret = split(now_line, ' ');

			if (ret[0] == "newmtl") {
				now_material_name = ret[1];
			} else if (ret[0] == "endmtl") {
				matmap->insert(std::pair<std::string, Material>(now_material_name, Material(diffuse, specular, specular_coefficient, metalic)));
				now_material_name = "";
			} else if (ret[0] == "diffuse" && ret.size() >= 4) {
				const float x = atof(ret[1].c_str());
				const float y = atof(ret[2].c_str());
				const float z = atof(ret[3].c_str());
				diffuse = Float3(x, y, z);
			} else if (ret[0] == "specular" && ret.size() >= 4) {
				const float x = atof(ret[1].c_str());
				const float y = atof(ret[2].c_str());
				const float z = atof(ret[3].c_str());
				specular = Float3(x, y, z);
			} else if (ret[0] == "metalic" && ret.size() >= 2) {
				const float x = atof(ret[1].c_str());
				metalic = x;
			} else if (ret[0] == "specular_coefficient" && ret.size() >= 2) {
				const float x = atof(ret[1].c_str());
				specular_coefficient = x;
			}
		}

		return true;
	}

	static bool load(const char *filename, TriangleMesh *mesh) {
		// 一行ずつ処理すんぞ
		FileManager manager;
		if (!manager.load(filename))
			return false;

		MeshBody mesh_body;
		Material *now_material = NULL;

		for (;;) {
			std::string now_line;
			if (!manager.gets(&now_line))
				break;


			if (now_line.size() == 0 || now_line[0] == '#') {
				continue;
			}

			std::vector<std::string> ret = split(now_line, ' ');

			if (ret[0] == "mtllib") {
				// マテリアルの読み込み
				if (!load_material(ret[1].c_str(), &mesh_body.matmap))
					return false;

				/*
				MaterialMap::iterator result = mesh_body.matmap.find("Material.001");
				std::cout << &(result->second) << std::endl;
				*/

			} else if (ret[0] == "usemtl") {

				MaterialMap::iterator result = mesh_body.matmap.find(ret[1]);
				if (result == mesh_body.matmap.end()) {
					now_material = NULL;
				} else {
					now_material = &(result->second);

//					std::cout << result->second.specular << std::endl;

//					now_material = new Material(result->second);
				}
			} else if (ret[0] == "o") {
				// めんどいからとりあえずoはなしね
			} else if (ret[0] == "v" && ret.size() >= 4) {
				const float x = atof(ret[1].c_str());
				const float y = atof(ret[2].c_str());
				const float z = atof(ret[3].c_str());
				mesh_body.v.push_back(Float3(x, y, z));
			} else if (ret[0] == "vn" && ret.size() >= 4) {
				const float x = atof(ret[1].c_str());
				const float y = atof(ret[2].c_str());
				const float z = atof(ret[3].c_str());
				mesh_body.vn.push_back(Float3(x, y, z));
			} else if (ret[0] == "vt" && ret.size() >= 4) {
				const float x = atof(ret[1].c_str());
				const float y = atof(ret[2].c_str());
				const float z = atof(ret[3].c_str());
				mesh_body.vt.push_back(Float3(x, y, z));
			} else if (ret[0] == "f" && ret.size() >= 4) {
				const Int3 f1 = face_to_data(ret[1]);
				const Int3 f2 = face_to_data(ret[2]);
				const Int3 f3 = face_to_data(ret[3]);
				Triangle t;

				t.v_index  = Int3(f1.x, f2.x, f3.x) - Int3(1, 1, 1);
				t.vt_index = Int3(f1.y, f2.y, f3.y) - Int3(1, 1, 1);
				t.vn_index = Int3(f1.z, f2.z, f3.z) - Int3(1, 1, 1);
				t.material = now_material;

				mesh_body.triangle.push_back(t); 
			}
		}
		/*
		MaterialMap::iterator result = mesh_body.matmap.find("Material.001");
		std::cout << &(result->second) << std::endl;
		std::cout << result->second.specular << std::endl;
		*/
		

		mesh->set(mesh_body);

		// BVH作る
		std::vector<RefTriangle> ref_triangle;

		for (int i = 0; i < mesh_body.triangle.size(); ++i) {
			const Int3 v_index = mesh_body.triangle[i].v_index;
			const Float3 v0 = mesh_body.v[v_index.x];
			const Float3 v1 = mesh_body.v[v_index.y];
			const Float3 v2 = mesh_body.v[v_index.z];

			RefTriangle tr(&mesh_body.v[v_index.x], &mesh_body.v[v_index.y], &mesh_body.v[v_index.z], i);
			ref_triangle.push_back(tr);
		}

		mesh->build(ref_triangle);

		return true;
	}
};



} // namespace rt

} // namespace hstd


#endif // _OBJ_MESH_