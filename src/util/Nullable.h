// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdexcept>

#ifndef UTIL_NULLABLE_H_
#define UTIL_NULLABLE_H_

namespace util {

template<typename T>
class Nullable {
 public:
  Nullable()
   : val(), null(true) {}
  Nullable(T* valPointer)
   : val(), null(true) {
   if (valPointer) {
     assign(*valPointer);
   }
  }
  Nullable(const T& value)
   : val(value), null(false) {}
  Nullable(const Nullable& other)
   : val(other.val), null(other.isNull()) {}

  Nullable& operator=(const Nullable& other) {
    if (!other.isNull()) val = other.get();
    null = other.isNull();
    return *this;
  }

  T operator=(const T& other) {
    assign(other);
    return val;
  }

  /**
   * Passing through comparision operators
   */

  bool operator==(const Nullable& other) const {
    return (other.isNull() && isNull()) || other.get() == get();
  }

  bool operator!=(const Nullable& other) const {
    return !(*this == other);
  }

  bool operator<(const Nullable& other) const {
    return !other.isNull() && !isNull() && get() < other.get();
  }

  bool operator>(const Nullable& other) const {
    return !(*this < other || *this == other);
  }

  bool operator<=(const Nullable& other) const {
    return *this < other || *this == other;
  }

  bool operator>=(const Nullable& other) const {
    return *this > other || *this == other;
  }

  bool operator==(const T& other) const {
    return !isNull() && other == get();
  }

  bool operator!=(const T& other) const {
    return !(*this == other);
  }

  bool operator<(const T& other) const {
    return !isNull() && get() < other;
  }

  bool operator>(const T& other) const {
    return !(*this < other || *this == other);
  }

  bool operator<=(const T& other) const {
    return *this < other || *this == other;
  }

  bool operator>=(const T& other) const {
    return *this > other || *this == other;
  }

  operator T() const {
    return get();
  }

  bool isNull() const {
    return null;
  }

  T get() const {
    if (!isNull()) return val;
    else throw std::runtime_error("Trying to retrieve value of NULL object.");
  }

private:
  void assign(T v) {
    val = v;
    null = false;
  }

  T val;
  bool null;
};

}

#endif  // UTIL_NULLABLE_H_
