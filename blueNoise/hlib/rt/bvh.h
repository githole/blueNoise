#ifndef _BVH_H_
#define _BVH_H_

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

#include "vec3.h"
#include "bbox.h"
#include "triangle.h"
#include "triangleMesh.h"
#include "constant.h"

namespace hstd {

namespace rt {

#if 0
const float kEPS = 1e-5;

struct BVHPrimitiveInfo {
	int primitiveNumber;
	Double3 centroid;
	BBox bounds;

	BVHPrimitiveInfo(int pn, const BBox& b) :
		primitiveNumber(pn), bounds(b){
			centroid = 0.5f * b.pmin + 0.5f * b.pmax;
	}
};

struct BVHBuildNode;
typedef BVHBuildNode* BVHBuildNodePtr;
struct BVHBuildNode  {
	BBox bounds;
	BVHBuildNodePtr children[2];
	int splitAxis, firstPrimOffset, nPrimitives;

	BVHBuildNode() {
	//	children[0] = children[1] = BVHBuildNodePtr(NULL);
	}

	void InitLeaf(int first, int n, const BBox& b) {
		firstPrimOffset = first;
		nPrimitives = n;
		bounds = b;
	}

	void InitInterior(int axis, BVHBuildNodePtr c0, BVHBuildNodePtr c1) {
		children[0] = c0;
		children[1] = c1;
		bounds = unionBBox(c0->bounds, c1->bounds);
		splitAxis = axis;
		nPrimitives = 0;
	}
};

struct ComparePoints {
	int dim;
	ComparePoints(int d) { dim = d; }
	bool operator()(const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b) const {
		return a.centroid[dim] < b.centroid[dim];
	}
};

struct CompareToBucket {
	CompareToBucket(int split, int num, int d, const BBox &b)
		: centroidBounds(b) {
			splitBucket = split;
			nBuckets = num;
			dim = d;
	}

	bool operator()(const BVHPrimitiveInfo &p) const {
		int b = nBuckets * ((p.centroid[dim] - centroidBounds.pmin[dim]) / 
							(centroidBounds.pmax[dim]- centroidBounds.pmin[dim]));
		if (b == nBuckets) b = nBuckets - 1;
		return b <= splitBucket;
	};

	int splitBucket, nBuckets, dim;
	const BBox &centroidBounds;
};

struct LinearBVHNode {
	BBox bounds;
	union {
		int primitiveOffset;
		int secondChildOffset;
	};
	char nPrimitives;
	char axis;
	char pad[2];
};


class Triangle2BVH {
private:
	int maxPrimsInNode;
public:
	std::vector<RefTriangle*> orderedTris;
	std::vector<RefTriangle*> tris;
	
	LinearBVHNode* nodes;

	Triangle2BVH() : nodes(NULL) {
	};

	BVHBuildNodePtr recursiveBuild(std::vector<BVHPrimitiveInfo> &buildData, int start, int end, int *totalNodes, std::vector<RefTriangle*> &orderedPrims) {
		(*totalNodes) ++;
		BVHBuildNodePtr node = BVHBuildNodePtr(new BVHBuildNode);

		BBox bbox;
		for (int i = start; i < end; ++i)
			bbox = unionBBox(bbox, buildData[i].bounds);

		int nPrimitives = end - start;
		if (nPrimitives == 1) {
			// 葉
			int firstPrimOffset = orderedPrims.size();
			for (int i = start; i < end; ++i ){
				int primNum = buildData[i].primitiveNumber;
				orderedPrims.push_back(tris[primNum]);
			}
			node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
		} else {
			// 中間ノード
			BBox centroidBounds;
			for (int i = start; i < end; ++i)
				centroidBounds = unionBBox(centroidBounds, buildData[i].centroid);
			int dim = centroidBounds.maximumExtent();

			int mid = (start + end) / 2;
			if (centroidBounds.pmax[dim] == centroidBounds.pmin[dim]) {
				// 葉
				int firstPrimOffset = orderedPrims.size();
				for (int i = start; i < end; ++i ){
					int primNum = buildData[i].primitiveNumber;
					orderedPrims.push_back(tris[primNum]);
				}
				node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
				return node;
			}

			// 分割
			if (nPrimitives <= 4) {
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
					int b = nBuckets * ((buildData[i].centroid[dim] - centroidBounds.pmin[dim]) /
										(centroidBounds.pmax[dim] - centroidBounds.pmin[dim]));
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

				if (nPrimitives > maxPrimsInNode || minCost < nPrimitives) {
					BVHPrimitiveInfo* pmid = std::partition(&buildData[start], &buildData[end-1] + 1, CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
					mid = pmid - &buildData[0];
				} else {
					// 葉
					int firstPrimOffset = orderedPrims.size();
					for (int i = start; i < end; ++i ){
						int primNum = buildData[i].primitiveNumber;
						orderedPrims.push_back(tris[primNum]);
					}
					node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
				}
			}

			node->InitInterior(dim,
				recursiveBuild(buildData, start, mid, totalNodes, orderedPrims),
								recursiveBuild(buildData, mid, end, totalNodes, orderedPrims));
		}

		return node;
	}

	int flattenBVHTree(BVHBuildNodePtr node, int *offset) {
		LinearBVHNode* linearNode = &nodes[*offset];
		linearNode->bounds = node->bounds;
		int myOffset = (*offset) ++;
		if (node->nPrimitives > 0) {
			linearNode->primitiveOffset = node->firstPrimOffset;
			linearNode->nPrimitives = node->nPrimitives;
		} else {
			linearNode->axis = node->splitAxis;
			linearNode->nPrimitives = 0;
			flattenBVHTree(node->children[0], offset);
			linearNode->secondChildOffset = flattenBVHTree(node->children[1], offset);
		}

		return myOffset;
	}

	void CreateBVHFromTriangle2s(const std::vector<RefTriangle>& tris_) {
		//tris = tris_;
		orderedTris.clear();
		maxPrimsInNode = 32;

		std::vector<BVHPrimitiveInfo> buildData;
		for (int i = 0; i < tris.size(); i ++) {
			BBox b = tris[i]->objectBound();

			buildData.push_back(BVHPrimitiveInfo(i, b));
		}
		int totalNodes = 0;
		orderedTris.reserve(tris.size());

		BVHBuildNodePtr root = recursiveBuild(buildData, 0, tris.size(), &totalNodes, orderedTris);

		tris.swap(orderedTris);

		nodes = new LinearBVHNode[totalNodes];
		int offset = 0;
		flattenBVHTree(root, &offset);
	}
	
	inline bool IntersectP(const BBox &bounds, const Ray &ray, const Double3& invDir, const int dirIsNeg[3]) {
		float tmin =  (bounds[  dirIsNeg[0]].x - ray.org.x) * invDir.x;
		float tmax =  (bounds[1-dirIsNeg[0]].x - ray.org.x) * invDir.x;
		float tymin = (bounds[  dirIsNeg[1]].y - ray.org.y) * invDir.y;
		float tymax = (bounds[1-dirIsNeg[1]].y - ray.org.y) * invDir.y;
		if ((tmin > tymax) || (tymin > tmax))
			return false;
		if (tymin > tmin) tmin = tymin;
		if (tymax < tmax) tmax = tymax;

		float tzmin = (bounds[  dirIsNeg[2]].z - ray.org.z) * invDir.z;
		float tzmax = (bounds[1-dirIsNeg[2]].z - ray.org.z) * invDir.z;
		if ((tmin > tzmax) || (tzmin > tmax))
			return false;
		if (tzmin > tmin)
			tmin = tzmin;
		if (tzmax < tmax)
			tmax = tzmax;
		return (tmin < kINF) && (tmax > 0.0);
	}
	
	virtual bool Intersect(const Ray &ray, Hitpoint* hitpoint) {
		// BVHで探索
		Double3 invDir(1.0f / ray.dir.x, 1.0f / ray.dir.y, 1.0f / ray.dir.z);
		int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0 };

		int todoOffset = 0, nodeNum = 0;
		int todo[256];

		bool hit = false;

		while (true) {
			const LinearBVHNode *node = &nodes[nodeNum];

			if (IntersectP(node->bounds, ray, invDir, dirIsNeg)) {

				if (node->nPrimitives > 0) {
					hit += node->nPrimitives;
					for (int j = 0; j < node->nPrimitives; ++j) {
						/*
						Intersection in;
						Triangle2* tri = bvh.tris[node->primitiveOffset + j];
						// Triangle2 と交差判定
						if (tri->Intersect(ray, &in, reverse)) {
							hit = true;
							if (intersection->thit > in.thit) {
								*intersection = in;
							}
						}
						*/
					}

					if (todoOffset == 0) break;
					nodeNum = todo[--todoOffset];
				} else {
					if (dirIsNeg[node->axis]) {
						todo[todoOffset ++] = nodeNum + 1;
						nodeNum = node->secondChildOffset;
					} else {
						todo[todoOffset ++] = node->secondChildOffset;
						nodeNum = nodeNum + 1;
					}
				}
			} else {
				if (todoOffset == 0) break;
				nodeNum = todo[--todoOffset];
			}
		}
		return hit;
	}

	virtual ~Triangle2BVH() {
		Delete();
	}

	void Delete() {
		if (nodes != NULL)
			delete[] nodes;
	}

};

#endif


} // namespace rt

} // namespace hstd

#endif // _BVH_H_
