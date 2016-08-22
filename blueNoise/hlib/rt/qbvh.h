#ifndef _QBVH_H_
#define _QBVH_H_

#include <vector>
#include <algorithm>
#include <queue>

/*
#ifdef _MSC_VER
#ifndef _M_IX86
#error "No SSE!"
#endif
#else
#ifndef __SSE__
#error "No SSE!"
#endif
#endif
*/

#include <xmmintrin.h>

#include "..\vec3.h"
#include "..\constant.h"

#include "bbox.h"
#include "triangle.h"
#include "triangleMesh.h"

namespace hstd {

namespace rt {
	
static const int OrderTable[] = {
	0x44444,0x44444,0x44444,0x44444,0x44444,0x44444,0x44444,0x44444,
	0x44440,0x44440,0x44440,0x44440,0x44440,0x44440,0x44440,0x44440,
	0x44441,0x44441,0x44441,0x44441,0x44441,0x44441,0x44441,0x44441,
	0x44401,0x44401,0x44410,0x44410,0x44401,0x44401,0x44410,0x44410,
	0x44442,0x44442,0x44442,0x44442,0x44442,0x44442,0x44442,0x44442,
	0x44402,0x44402,0x44402,0x44402,0x44420,0x44420,0x44420,0x44420,
	0x44412,0x44412,0x44412,0x44412,0x44421,0x44421,0x44421,0x44421,
	0x44012,0x44012,0x44102,0x44102,0x44201,0x44201,0x44210,0x44210,
	0x44443,0x44443,0x44443,0x44443,0x44443,0x44443,0x44443,0x44443,
	0x44403,0x44403,0x44403,0x44403,0x44430,0x44430,0x44430,0x44430,
	0x44413,0x44413,0x44413,0x44413,0x44431,0x44431,0x44431,0x44431,
	0x44013,0x44013,0x44103,0x44103,0x44301,0x44301,0x44310,0x44310,
	0x44423,0x44432,0x44423,0x44432,0x44423,0x44432,0x44423,0x44432,
	0x44023,0x44032,0x44023,0x44032,0x44230,0x44320,0x44230,0x44320,
	0x44123,0x44132,0x44123,0x44132,0x44231,0x44321,0x44231,0x44321,
	0x40123,0x40132,0x41023,0x41032,0x42301,0x43201,0x42310,0x43210,
};


const float kEPS = 1e-5;

class QBVH {
private:
	int maxPrimsInNode_;

	
	struct BVHPrimitiveInfo {
		int primitiveNumber;
		Float3 centroid;
		BBox bounds;

		BVHPrimitiveInfo(const int pn, const BBox& b) :
			primitiveNumber(pn), bounds(b){
				centroid = 0.5f * b.pmin + 0.5f * b.pmax;
		}
	};


	// 四つの三角形をパックする
	struct SIMDTrianglePack {
		__m128 x[3];
		__m128 y[3];
		__m128 z[3];
		int idx[4];


	//	char padding_[12];
	};

	struct BVHBuildNode  {
		BBox bounds;
		BVHBuildNode* children[2];
		int splitAxis, firstPrimOffset, nPrimitives;

	//	SIMDTrianglePack *simd_[2];
		int simdTrisIdx;

		BVHBuildNode() {
			children[0] = children[1] = NULL;
		//	children[0] = children[1] = BVHBuildNodePtr(NULL);
		}

		void InitLeaf(int first, int n, const BBox& b, const int asimdTrisIdx){
			firstPrimOffset = first;
			nPrimitives = n;
			bounds = b;
			simdTrisIdx = asimdTrisIdx;
			splitAxis = 0;
		}

		void InitInterior(int axis, BVHBuildNode* c0, BVHBuildNode* c1) {
			children[0] = c0;
			children[1] = c1;
			bounds = unionBBox(c0->bounds, c1->bounds);
			splitAxis = axis;
			firstPrimOffset = -1;
			nPrimitives = 0;
		}
	};

	struct ComparePoints {
		int dim;
		ComparePoints(const int d) : dim(d){}
		bool operator()(const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b) const {
			return a.centroid[dim] < b.centroid[dim];
		}
	};

	struct CompareToBucket {
		int splitBucket, nBuckets, dim;
		const BBox &centroidBounds;

		CompareToBucket(int split, int num, int d, const BBox &b)
			: centroidBounds(b) {
				splitBucket = split;
				nBuckets = num;
				dim = d;
		}

		bool operator()(const BVHPrimitiveInfo &p) const {
			int b = (int)(nBuckets * ((p.centroid[dim] - centroidBounds.pmin[dim]) / 
								(centroidBounds.pmax[dim] - centroidBounds.pmin[dim])));
			if (b == nBuckets) 
				b = nBuckets - 1;
			return b <= splitBucket;
		};
	};

	struct Children {
		union {
			struct Node {
				unsigned flag : 1;
				unsigned index: 31;
			} node;
			struct Leaf {
				unsigned flag : 1;
				unsigned nPrimitives: 3;
				unsigned index: 28;
			} leaf;

			unsigned int raw;
		};
	};

	struct SIMDBVHNode{
		__m128 bboxes[2][3];//4 float min-max xyz
		Children children[4]; //4 children
		int axis_top;       //top axis
		int axis_left;      //left axis
		int axis_right;     //right axis
		int reserved;       //padding
	};

public:
	std::vector<RefTriangle*> orderedTris;
	std::vector<RefTriangle> tris;
	std::vector<SIMDTrianglePack*> simdTris;
	std::vector<SIMDBVHNode*> simdNodes;
//	std::vector<SIMDBVHNode, AlignmentAllocator<SIMDBVHNode, 16> > simdNodes_;
	
	int stats[16];

	__m128 zero;
	__m128 one;
	__m128 inf;
	__m128 keps;

	QBVH() {
		for (int i = 0; i < 16; ++i)
			stats[i] = 0;
		//std::cout << sizeof(SIMDTrianglePack) << std::endl;
		__declspec(align(32)) float one_f[4] = {1.0f, 1.0f, 1.0f, 1.0f};
		__declspec(align(32)) float inf_f[4] = {kINF, kINF, kINF, kINF};
		__declspec(align(32)) float zero_f[4] = {0.0f, 0.0f, 0.0f, 0.0f};
		__declspec(align(32)) float keps_f[4] = {kEPS, kEPS, kEPS, kEPS};
		
		zero = _mm_load_ps(zero_f);
		one = _mm_load_ps(one_f);
		inf = _mm_load_ps(inf_f);
		keps = _mm_load_ps(keps_f);
	};

	BVHBuildNode* recursiveBuild(std::vector<BVHPrimitiveInfo> &buildData, int start, int end, int *totalNodes, std::vector<RefTriangle*> &orderedPrims) {
		(*totalNodes) ++;
		BVHBuildNode* node = new BVHBuildNode;

		BBox bbox;
		for (int i = start; i < end; ++i)
			bbox = unionBBox(bbox, buildData[i].bounds);

		int nPrimitives = end - start;
		if (nPrimitives <= 4) {
			//std::cout << nPrimitives << " ";
			// 葉
			int firstPrimOffset = orderedPrims.size();

			SIMDTrianglePack *simdt = (SIMDTrianglePack*)_aligned_malloc(sizeof(SIMDTrianglePack), 16);
			
			__declspec(align(16)) float x[4 * 3] = {0};
			__declspec(align(16)) float y[4 * 3] = {0};
			__declspec(align(16)) float z[4 * 3] = {0};

			int cnt = 0;
			for (int i = start; i < end; ++i , ++cnt){
				const int idx = buildData[i].primitiveNumber;
				orderedPrims.push_back(&tris[idx]);
				
				int t = cnt % 4;

				simdt->idx[t] = firstPrimOffset + cnt;
				x[t] = tris[idx].p[0]->x;
				x[4 + t] = tris[idx].p[1]->x;
				x[8 + t] = tris[idx].p[2]->x;
				y[t] = tris[idx].p[0]->y;
				y[4 + t] = tris[idx].p[1]->y;
				y[8 + t] = tris[idx].p[2]->y;
				z[t] = tris[idx].p[0]->z;
				z[4 + t] = tris[idx].p[1]->z;
				z[8 + t] = tris[idx].p[2]->z;
			}
			for (; cnt < 4; ++cnt) {
				simdt->idx[cnt%4] = -1;
			}

			for (int i = 0; i < 3; ++i) {
				simdt->x[i] = _mm_load_ps(x + 4 * i);
				simdt->y[i] = _mm_load_ps(y + 4 * i);
				simdt->z[i] = _mm_load_ps(z + 4 * i);
			}

			simdTris.push_back(simdt);
			node->InitLeaf(firstPrimOffset, nPrimitives, bbox, simdTris.size() - 1);
//			node->InitLeaf(firstPrimOffset, nPrimitives, bbox, NULL);
		} else {
			// 中間ノード
			BBox centroidBounds;
			for (int i = start; i < end; ++i)
				centroidBounds = unionBBox(centroidBounds, buildData[i].centroid);
			int dim = centroidBounds.maximumExtent();
			int mid = (start + end) / 2;
			const float denom = (centroidBounds.pmax[dim] - centroidBounds.pmin[dim]);

			if (nPrimitives <= 16 || denom <= 0.0f) {
				// Equally-sized subsets
				std::nth_element(&buildData[start], &buildData[mid], &buildData[end-1] + 1, ComparePoints(dim));
			} else {
				// SAH
				const int nBuckets = 12;

				struct BucketInfo {
					BucketInfo() { count = 0; }
					int count;
					BBox bounds;
				};
				BucketInfo buckets[nBuckets];
				for (int i = start; i < end;  ++i ) {
					int b = (int)(nBuckets * ((buildData[i].centroid[dim] - centroidBounds.pmin[dim]) / denom));
					if (b == nBuckets) b = nBuckets - 1;
					
					buckets[b].count ++;
					buckets[b].bounds = unionBBox(buckets[b].bounds, buildData[i].bounds);
				}
				float cost[nBuckets - 1];
				for (int i = 0; i < nBuckets - 1; i ++) {
					BBox b0, b1;
					int count0 = 0, count1  = 0;
					for (int j = 0; j <= i; j ++) {
						b0 = unionBBox(b0, buckets[j].bounds);
						count0 += buckets[j].count;
					}
					for (int j = i + 1; j < nBuckets; j ++) {
						b1 = unionBBox(b1, buckets[j].bounds);
						count1 += buckets[j].count;
					}
					cost[i] += 0.125f + (count0 * b0.surfaceArea() + count1 * b1.surfaceArea()) / bbox.surfaceArea();
				}

				float minCost = cost[0];
				int minCostSplit = 0;
				for (int i = 1; i < nBuckets - 1; ++ i) {
					if (cost[i] < minCost) {
						minCost = cost[i];
						minCostSplit = i;
					}
				}

				if (nPrimitives > 0 || minCost < nPrimitives) {
					BVHPrimitiveInfo* pmid = std::partition(&buildData[start], &buildData[end-1] + 1, CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
					mid = pmid - &buildData[0];
				}
			}
			node->InitInterior(dim,
				recursiveBuild(buildData, start, mid, totalNodes, orderedPrims),
				recursiveBuild(buildData, mid, end, totalNodes, orderedPrims));
		}

		return node;
	}
	void collapse2QBVH(BVHBuildNode* node) {
		BVHBuildNode *lc = node->children[0];
		BVHBuildNode *rc = node->children[1];

		BVHBuildNode *c[4] = {0};
		
		SIMDBVHNode *n;
		n = (SIMDBVHNode*)_aligned_malloc(sizeof(SIMDBVHNode), 16);
		simdNodes.push_back(n);
		n->axis_top = node->splitAxis;
		n->axis_left = n->axis_right = 0;

		if (lc != NULL) {
			n->axis_left = lc->splitAxis;
			if (lc->nPrimitives == 0) {
				c[0] = lc->children[0];
				c[1] = lc->children[1];
			} else {
				c[0] = lc;
			}
		}
		if (rc != NULL) {
			n->axis_right = rc->splitAxis;
			if (rc->nPrimitives == 0) {
				c[2] = rc->children[0];
				c[3] = rc->children[1];
			} else {
				c[2] = rc;
			}
		}
		__declspec(align(16)) float bboxes[2][3][4];
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 4; ++k) {
				if (c[k] != NULL) {
					bboxes[0][j][k] = c[k]->bounds.pmin[j];
					bboxes[1][j][k] = c[k]->bounds.pmax[j];
				}
			}
		}
		for(int m = 0; m < 2; m++){//minmax
			for(int j = 0; j < 3; j++){//xyz
				n->bboxes[m][j] = _mm_load_ps(bboxes[m][j]);
			}
		}

		for (int i = 0; i < 4; ++i) {
			if (c[i] == NULL) {
				n->children[i].leaf.flag = 1;
				n->children[i].leaf.nPrimitives = 0;
				n->children[i].leaf.index = 0;
			} else {
				if (c[i]->nPrimitives == 0) {
					n->children[i].node.flag = 0;
					n->children[i].node.index= simdNodes.size();
					collapse2QBVH(c[i]);
				} else {
					// std::cout << c[i]->nPrimitives_ << " ";
					n->children[i].leaf.flag = 1;
					n->children[i].leaf.nPrimitives = c[i]->nPrimitives;
					n->children[i].leaf.index = c[i]->simdTrisIdx;
				}
			}
		}

		return;
	}
	BVHBuildNode *rootNode;
	void build(std::vector<RefTriangle>& atris) {
//		tris = atris;
		for (int i = 0; i < atris.size(); ++i) {
			tris.push_back(atris[i]);
		}

		orderedTris.clear();
		maxPrimsInNode_ = 32;

		std::vector<BVHPrimitiveInfo> buildData;
		for (unsigned int i = 0; i < tris.size(); i ++) {
			BBox b = tris[i].objectBound();

			buildData.push_back(BVHPrimitiveInfo(i, b));
		}
		int totalNodes = 0;
		orderedTris.reserve(tris.size());

		rootNode = recursiveBuild(buildData, 0, tris.size(), &totalNodes, orderedTris);

		// collapse
		collapse2QBVH(rootNode);

		SIMDBVHNode *root = simdNodes[0];
		// std::cerr << "QBVH: " << simdNodes.size() << std::endl;
		// print(root, 0);
	}

	inline int test_AABB(
		const __m128 bboxes[2][3],  //4boxes : min-max[2] of xyz[3] of boxes[4]

		const __m128 org[3],        //ray origin
		const __m128 idir[3],       //ray inveresed direction
		const int sign[3],          //ray xyz direction -> +:0,-:1
		__m128 tmin, __m128 tmax    //ray range tmin-tmax 
		)
	{
		// x coordinate
		tmin = _mm_max_ps(
			tmin,
			_mm_mul_ps(_mm_sub_ps(bboxes[sign[0]][0],org[0]), idir[0])
			);
		tmax = _mm_min_ps(
			tmax,
			_mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[0]][0], org[0]), idir[0])
			);

		// y coordinate
		tmin = _mm_max_ps(
			tmin,
			_mm_mul_ps(_mm_sub_ps(bboxes[sign[1]][1],org[1]), idir[1])
			);
		tmax = _mm_min_ps(
			tmax,
			_mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[1]][1], org[1]), idir[1])
			);

		// z coordinate
		tmin = _mm_max_ps(
			tmin,
			_mm_mul_ps(_mm_sub_ps(bboxes[sign[2]][2],org[2]), idir[2])
			);
		tmax = _mm_min_ps(
			tmax,
			_mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[2]][2], org[2]), idir[2])
			);
		return _mm_movemask_ps(_mm_cmpge_ps(tmax, tmin));//tmin<tmaxとなれば交差
	}
	
	bool intersect(const Ray &ray, Hitpoint* hitpoint) {

		__m128 sseOrg[3];
		__m128 sseiDir[3];
		int sign[3];
		Float3 idir(1.0f / ray.dir.x, 1.0f / ray.dir.y, 1.0f / ray.dir.z);
		__declspec(align(16)) float r_idir_x[4] = {idir.x, idir.x, idir.x, idir.x};
		__declspec(align(16)) float r_idir_y[4] = {idir.y, idir.y, idir.y, idir.y};
		__declspec(align(16)) float r_idir_z[4] = {idir.z, idir.z, idir.z, idir.z};
					
		__declspec(align(16)) float r_org_x[4] = {ray.org.x, ray.org.x, ray.org.x, ray.org.x};
		__declspec(align(16)) float r_org_y[4] = {ray.org.y, ray.org.y, ray.org.y, ray.org.y};
		__declspec(align(16)) float r_org_z[4] = {ray.org.z, ray.org.z, ray.org.z, ray.org.z};

		__declspec(align(16)) float r_dir_x[4] = {ray.dir.x, ray.dir.x, ray.dir.x, ray.dir.x};
		__declspec(align(16)) float r_dir_y[4] = {ray.dir.y, ray.dir.y, ray.dir.y, ray.dir.y};
		__declspec(align(16)) float r_dir_z[4] = {ray.dir.z, ray.dir.z, ray.dir.z, ray.dir.z};
					
		__m128 dir_x = _mm_load_ps(r_dir_x);
		__m128 dir_y = _mm_load_ps(r_dir_y);
		__m128 dir_z = _mm_load_ps(r_dir_z);

		sseOrg[0] = _mm_load_ps(r_org_x);
		sseOrg[1] = _mm_load_ps(r_org_y);
		sseOrg[2] = _mm_load_ps(r_org_z);

		sseiDir[0] = _mm_load_ps(r_idir_x);
		sseiDir[1] = _mm_load_ps(r_idir_y);
		sseiDir[2] = _mm_load_ps(r_idir_z);
		
		sign[0] = idir[0] < 0;
		sign[1] = idir[1] < 0;
		sign[2] = idir[2] < 0;

		const SIMDBVHNode* nodes = simdNodes[0];

		Children nodeStack[40];
		int todoNode = 0;

		nodeStack[0].raw = 0;

		bool hit = false;
		int triangle_index = -1;

		int cnt = 0;

		while (todoNode >= 0) {
			Children item = nodeStack[todoNode];
			todoNode--;//pop stack

			if(item.node.flag == 0){
				const SIMDBVHNode& node = *(simdNodes[item.node.index]);
				__declspec(align(16)) float now_distance_f[4] = {hitpoint->distance, hitpoint->distance, hitpoint->distance, hitpoint->distance};
				__m128 now_distance = _mm_load_ps(now_distance_f);
				const int HitMask = test_AABB(node.bboxes, sseOrg, sseiDir, sign, zero, now_distance);

				if (HitMask) {
					const int nodeIdx = (sign[node.axis_top] << 2) | (sign[node.axis_left] << 1) | (sign[node.axis_right]);
					int bboxOrder = OrderTable[HitMask * 8 + nodeIdx];
					

					for (int i = 0; i < 4; ++i) {
						if (bboxOrder & 0x4)
							break;
						++todoNode;
						nodeStack[todoNode] = node.children[bboxOrder & 0x3];
						bboxOrder >>= 4;
					}
				}

			} else {

				// __declspec(align(16)) float no_hit_f[4];
				__declspec(align(16)) float t_f[4];
				__declspec(align(16)) float b1_f[4];
				__declspec(align(16)) float b2_f[4];
				int nohitmask;
				SIMDTrianglePack *s = simdTris[item.leaf.index];
				
				float eps = 1e-4f;
				__declspec(align(32)) float t0_f[4] = {0.0f - eps, 0.0f - eps, 0.0f - eps, 0.0f - eps};
				__declspec(align(32)) float t1_f[4] = {1.0f + eps, 1.0f + eps, 1.0f + eps, 1.0f + eps};
		
				__m128 t0 = _mm_load_ps(t0_f);
				__m128 t1 = _mm_load_ps(t1_f);
				
				__m128 e1_x = _mm_sub_ps(s->x[1], s->x[0]);
				__m128 e1_y = _mm_sub_ps(s->y[1], s->y[0]);
				__m128 e1_z = _mm_sub_ps(s->z[1], s->z[0]);
					
				__m128 e2_x = _mm_sub_ps(s->x[2], s->x[0]);
				__m128 e2_y = _mm_sub_ps(s->y[2], s->y[0]);
				__m128 e2_z = _mm_sub_ps(s->z[2], s->z[0]);

				__m128 s1_x = _mm_sub_ps(_mm_mul_ps(dir_y, e2_z), _mm_mul_ps(dir_z, e2_y));
				__m128 s1_y = _mm_sub_ps(_mm_mul_ps(dir_z, e2_x), _mm_mul_ps(dir_x, e2_z));
				__m128 s1_z = _mm_sub_ps(_mm_mul_ps(dir_x, e2_y), _mm_mul_ps(dir_y, e2_x));

				__m128 divisor = _mm_add_ps(_mm_add_ps(_mm_mul_ps(s1_x, e1_x), _mm_mul_ps(s1_y, e1_y)), _mm_mul_ps(s1_z, e1_z));
				__m128 no_hit  = _mm_cmpeq_ps(divisor, zero);	

				__m128 invDivisor = _mm_rcp_ps(divisor);

				__m128 d_x = _mm_sub_ps(sseOrg[0], s->x[0]);
				__m128 d_y = _mm_sub_ps(sseOrg[1], s->y[0]);
				__m128 d_z = _mm_sub_ps(sseOrg[2], s->z[0]);

				__m128 b1 = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(d_x, s1_x), _mm_mul_ps(d_y, s1_y)), _mm_mul_ps(d_z, s1_z)),
										invDivisor);
				no_hit = _mm_or_ps(no_hit, _mm_or_ps(_mm_cmplt_ps(b1, t0), _mm_cmpgt_ps(b1, t1)));
					
				__m128 s2_x = _mm_sub_ps(_mm_mul_ps(d_y, e1_z), _mm_mul_ps(d_z, e1_y));
				__m128 s2_y = _mm_sub_ps(_mm_mul_ps(d_z, e1_x), _mm_mul_ps(d_x, e1_z));
				__m128 s2_z = _mm_sub_ps(_mm_mul_ps(d_x, e1_y), _mm_mul_ps(d_y, e1_x));

				__m128 b2 = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(dir_x, s2_x), _mm_mul_ps(dir_y, s2_y)), _mm_mul_ps(dir_z, s2_z)),
										invDivisor);
				no_hit = _mm_or_ps(no_hit, _mm_or_ps(_mm_cmplt_ps(b2, t0), _mm_cmpgt_ps(_mm_add_ps(b1, b2), t1)));

				__m128 t = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(e2_x, s2_x), _mm_mul_ps(e2_y, s2_y)), _mm_mul_ps(e2_z, s2_z)),
										invDivisor);
				
				no_hit = _mm_or_ps(no_hit, _mm_cmplt_ps(t, keps));
				
				nohitmask = _mm_movemask_ps(no_hit);
				_mm_store_ps(t_f, t);
				
				for (int i = 0; i < 4; ++i) {
					if ((nohitmask & (1 << i)) == 0 && hitpoint->distance > t_f[i]) {
						hit = true;
						triangle_index = s->idx[i];
						hitpoint->distance = t_f[i];
						
						_mm_store_ps(b1_f, b1);
						_mm_store_ps(b2_f, b2);
						hitpoint->b1 = b1_f[i];
						hitpoint->b2 = b2_f[i];
					}
				}
			}
		}
		if (hit) {
			RefTriangle* t = orderedTris[triangle_index];
			hitpoint->triangle_index = t->original_triangle_index;
			
			
			/*
			hitpoint->v0 = *t->p[0];
			hitpoint->v1 = *t->p[1];
			hitpoint->v2 = *t->p[2];
			hitpoint->normal   = normalize(cross((*t->p[2]) - (*t->p[0]), (*t->p[1]) - (*t->p[0])));
			*/
		}

		return hit;
	}

	virtual ~QBVH() {
	}

};


} // namespace rt

} // namespace hstd

#endif // _QBVH_H_
