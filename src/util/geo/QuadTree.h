// Copyright 2020, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_QUADTREE_H_
#define UTIL_GEO_QUADTREE_H_

#include <map>
#include <set>
#include <vector>
#include "util/geo/Geo.h"
#include "util/geo/output/GeoJsonOutput.h"

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

template <typename V, typename T>
struct SplitFunc {
  virtual ~SplitFunc() = default;
  virtual bool operator()(const QuadNode<T>& nd,
                       const QuadValue<V, T>& newVal) const = 0;
};

template <typename V, typename T>
struct CapaSplitFunc : SplitFunc<V, T> {
  CapaSplitFunc(size_t c) : _c(c) {}
  virtual bool operator()(const QuadNode<T>& nd,
                       const QuadValue<V, T>& newVal) const {
    UNUSED(newVal);
    return static_cast<size_t>(nd.numEls) + 1 > _c;
  }
  size_t _c;
};

// QuadTree for point data (and only point data)
template <typename V, typename T>
class QuadTree {
 public:
  // initialization of a quad tree with maximum depth d and maximum node // capacity c
  QuadTree(size_t d, size_t c, const Box<T>& bbox);

  QuadTree(size_t d, const SplitFunc<V, T>& splitF, const Box<T>& bbox);

  // insert into the tree
  void insert(const V& val, const Point<T>& point);

  // insert into a specific node
  void insert(int64_t vid, int64_t nid, size_t d);

  size_t size() const;

  const std::vector<QuadNode<T>>& getNds() const;
  const QuadNode<T>& getNd(size_t nid) const;

  // GeoJSON output
  void print(std::ostream& o) const;

 private:
  size_t _maxDepth;
  std::vector<QuadValue<V, T>> _vals;
  std::vector<QuadNode<T>> _nds;

  CapaSplitFunc<V, T> _capaFunc;

  const SplitFunc<V, T>& _splFunc;


  // split a node
  void split(size_t nid, size_t d);
};

#include "util/geo/QuadTree.tpp"

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_QUADTREE_H_
