// Copyright 2019, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename K, typename V>
void PriorityQueue<K, V>::push(K key, const V& val) {
  // if (_last - key > 0) key = _last; 
  if (key < _last) {
    std::cout << std::setprecision(10) << _last << " vs " << key << std::endl;
    key = _last;
  }
  // assert(key >= _last);
  _pq.emplace(std::pair<K, V>(key, val));
}

// _____________________________________________________________________________
template <typename K, typename V>
const K PriorityQueue<K, V>::topKey() {
  _last = _pq.top().first;
  return _pq.top().first;
}

// _____________________________________________________________________________
template <typename K, typename V>
const V& PriorityQueue<K, V>::topVal() {
  _last = _pq.top().first;
  return _pq.top().second;
}

// _____________________________________________________________________________
template <typename K, typename V>
bool PriorityQueue<K, V>::empty() const {
  return _pq.empty();
}

// _____________________________________________________________________________
template <typename K, typename V>
void PriorityQueue<K, V>::pop() {
  _pq.pop();
}
