// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRIDGRAPH_NODECOST_H_
#define OCTI_GRIDGRAPH_NODECOST_H_


#include <sstream>


namespace octi {
namespace gridgraph {

struct NodeCost {
  double _cost[8] = {0, 0, 0, 0, 0, 0, 0, 0};

  void normalize() {
    double smallest = std::numeric_limits<double>::max();
    for (int i = 0; i < 8; i++) {
      if (_cost[i] > -1 && _cost[i] < smallest) {
        smallest = _cost[i];
      }
    }

    for (int i = 0; i < 8; i++) {
      if (_cost[i] > -1) {
        _cost[i] -= smallest;
      }
    }
  }

  std::string toString() {
    std::stringstream ret;
    ret << "{";
    for (size_t i = 0; i < 8; i++) {
      if (_cost[i] < -1)  ret << "<B>" << (i != 7 ? ", " : "}");
      else  ret << _cost[i] << (i != 7 ? ", " : "}");
    }
    return ret.str();
  }

  double operator[](size_t i) const {
    return _cost[i];
  }

  double& operator[](size_t i) {
    return _cost[i];
  }

  NodeCost operator+(const NodeCost& other) const {
    NodeCost ret = (*this);
    for (size_t i = 0; i < 8; i++) {
      if (ret[i] < - 1) continue;
      ret[i] += other[i];
    }
    return ret;
  }

  void operator+=(const NodeCost& other) {
    (*this) = (*this) + other;
  }
};

}  // gridgraph
}  // octi



#endif  // OCTI_GRIDGRAPH_NODECOST_H_
