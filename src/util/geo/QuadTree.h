// Copyright 2020, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_QUADTREE_H_
#define UTIL_GEO_QUADTREE_H_

#include <map>
#include <set>
#include <vector>
#include "util/geo/Geo.h"

namespace util {
namespace geo {

template <typename V, typename T>
struct QuadValue {
  V val;           // the actual value of this entry
  Point<T> point;  // the value's position

  int64_t nextValue;  // index of the next quad value, -1 means no next value
};

template <typename T>
struct QuadNode {
  int64_t numEls;  // number of elements, -1 if this is not a leaf node
  int64_t childs;  // for leafs, points to the first value contained. for
                   // other nodes, points to the array block containing the
                   // 4 childs
  Box<T> bbox;
};

// QuadTree for point data (and only point data)
template <typename V, typename T>
class QuadTree {
 public:
  // initialization of a quad tree with maximum depth d
  QuadTree(size_t d, const Box<T>& bbox);

  // insert into the tree
  void insert(const V& val, const Point<T>& point);

  // insert into a specific node
  void insert(int64_t vid, int64_t nid);

  bool shouldSplit(const QuadNode<T>& nd, const QuadValue<V, T>& newVal) const;

  size_t size() const;

  const std::vector<QuadNode<T>>& getNds() const;
  const QuadNode<T>& getNd(size_t nid) const;

 private:
  size_t _maxDepth;
  std::vector<QuadValue<V, T>> _vals;
  std::vector<QuadNode<T>> _nds;

  // split a node
  void split(size_t nid);

};

#include "util/geo/QuadTree.tpp"

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_QUADTREE_H_
