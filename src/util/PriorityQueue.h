// Copyright 2019, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_PRIORITYQUEUE_H_
#define UTIL_PRIORITYQUEUE_H_

#include<iomanip>
#include<queue>
#include<iostream>

namespace util {
template <typename K, typename V>
class PriorityQueue {
  struct _ByFirst {
    bool operator()(const std::pair<K, V>& a, const std::pair<K, V>& b) {
      return a.first > b.first;
    }
  };


 public:
  PriorityQueue() : _last(std::numeric_limits<K>::lowest()) {}
  void push(K k, const V& v);
  const K topKey() ;
  const V& topVal() ;
  void pop();

  bool empty() const;

 private:
  K _last;
  std::priority_queue<std::pair<K, V>, std::vector<std::pair<K, V>>, _ByFirst>
      _pq;
};
#include "util/PriorityQueue.tpp"
}  // namespace util

#endif
