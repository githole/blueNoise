#ifndef _KD_TREE_H_
#define _KD_TREE_H_

#include <algorithm>
#include <vector>
#include <queue>

#include "vec3.h"

namespace hstd {

namespace rt {

// KD-tree
template<typename T>
class KDTree {
public:
	// k-NN searchのクエリ
	struct Query {
		double max_distance2; // 探索の最大半径
		size_t max_search_num; // 最大探索点数
		Float3 search_position; // 探索中心
		Float3 normal; // 探索中心における法線
		Query(const Float3 &search_position_, const Float3 &normal_, const double max_distance2_, const size_t max_search_num_) :
		max_distance2(max_distance2_), normal(normal_), max_search_num(max_search_num_), search_position(search_position_) {}
	};
	// 結果のQueueに乗せるためのデータ構造。
	struct ElementForQueue {
		const T *point;
		double distance2;
		ElementForQueue(const T *point_, const double distance2_) : point(point_), distance2(distance2_) {}
		bool operator<(const ElementForQueue &b) const {
			return distance2 < b.distance2;
		}
	};
	// KNNの結果を格納するキュー
	typedef std::priority_queue<ElementForQueue, std::vector<ElementForQueue> > ResultQueue;
private:
	std::vector<T> points;
	struct KDTreeNode {
		T* point;
		KDTreeNode* left;
		KDTreeNode* right;
		int axis;
	};
	KDTreeNode* root;
	void delete_kdtree(KDTreeNode* node) {
		if (node == NULL)
			return;
		delete_kdtree(node->left);
		delete_kdtree(node->right);
		delete node;
	}

	// フツーのk-NN search。
	void locate_points(typename KDTree<T>::ResultQueue* pqueue, KDTreeNode* node, typename KDTree<T>::Query &query) {
		if (node == NULL)
				return;
		const int axis = node->axis;

		double delta;
		switch (axis) {
		case 0: delta = query.search_position.x - node->point->position.x; break;
		case 1: delta = query.search_position.y - node->point->position.y; break;
		case 2: delta = query.search_position.z - node->point->position.z; break;
		}

		// 対象点<->探索中心の距離が設定半径以下　かつ　対象点<->探索中心の法線方向の距離が一定以下　という条件ならその対象点格納
		const Float3 dir = node->point->position - query.search_position;
		const double distance2 = dir.LengthSquared();
		const double dt = Dot(query.normal, dir / sqrt(distance2));
		if (distance2 < query.max_distance2 && fabs(dt) <= query.max_distance2 * 0.01) {
			pqueue->push(ElementForQueue(node->point, distance2));
			if (pqueue->size() > query.max_search_num) {
				pqueue->pop();
				query.max_distance2 = pqueue->top().distance2;
			}
		}
		if (delta > 0.0) { // みぎ
			locate_points(pqueue,node->right, query);
			if (delta * delta < query.max_distance2) {
				locate_points(pqueue, node->left, query);
			}
		} else { // ひだり
			locate_points(pqueue,node->left, query);
			if (delta * delta < query.max_distance2) {
				locate_points(pqueue, node->right, query);
			}
		}

	}

	static bool kdtree_less_operator_x(const T& left, const T& right) {
		return left.position.x < right.position.x;
	}
	static bool kdtree_less_operator_y(const T& left, const T& right) {
		return left.position.y < right.position.y;
	}
	static bool kdtree_less_operator_z(const T& left, const T& right) {
		return left.position.z < right.position.z;
	}

	KDTreeNode* create_kdtree_sub(typename std::vector<T>::iterator begin,typename std::vector<T>::iterator end, int depth) {
		if (end - begin <= 0) {
			return NULL;
		}
		const int axis = depth % 3;
		// 中央値
		switch (axis) {
		case 0: std::sort(begin, end, kdtree_less_operator_x); break;
		case 1: std::sort(begin, end, kdtree_less_operator_y); break;
		case 2: std::sort(begin, end, kdtree_less_operator_z); break;
		}
		const int median = (end - begin) / 2;
		KDTreeNode* node = new KDTreeNode;
		node->axis = axis;
		node->point = &(*(begin + median));
		// 子供
		node->left = create_kdtree_sub(begin, begin + median, depth + 1);
		node->right = create_kdtree_sub(begin + median + 1, end, depth + 1);
		return node;
	}
public:
	KDTree() {
		root = NULL;
	}
	virtual ~KDTree() {
		delete_kdtree(root);
	}
	size_t Size() {
		return points.size();
	}
	void SearchKNN(typename KDTree::ResultQueue* pqueue, typename KDTree<T>::Query &query) {
		locate_points(pqueue, root, query);
	}
	void AddPoint(const T &point) {
		points.push_back(point);
	}
	void CreateKDtree() {
		root = create_kdtree_sub(points.begin(), points.end(), 0);
	}
};


} // namespace rt

} // namespace hstd


#endif // _KD_TREE_H_