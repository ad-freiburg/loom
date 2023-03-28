// Copyright 2023, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_RTREE_H_
#define UTIL_GEO_RTREE_H_

#include <unordered_map>
#include <set>
#include <vector>
#include "util/3rdparty/RTree.h"

namespace util {
namespace geo {

class RTreeException : public std::runtime_error {
 public:
  RTreeException(std::string const& msg) : std::runtime_error(msg) {}
};

template <typename V, template <typename> class G, typename T>
class RTree {
 public:
  // RTree(const RTree<V, G, T>&) = delete;
  // RTree(RTree<V, G, T>&& o) = delete;

  // the empty RTree
  RTree() {};

  // add object t to this RTree
  void add(G<T> geom, V val);

  void get(const Box<T>& btbox, std::set<V>* s) const;
  void get(const G<T>& geom, double d, std::set<V>* s) const;
  void remove(V val);

  void getNeighbors(const V& val, double d, std::set<V>* s) const;
 private:
  // make mutable, 3rd party implementation is not really const-correct
  mutable DA::RTree<V, T, 2, T> _rtree;

  std::unordered_map<V, util::geo::Box<T>> _valIdx;

  static bool searchCb(V val, void* arg);
};

#include "util/geo/RTree.tpp"

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_RTREE_H_
