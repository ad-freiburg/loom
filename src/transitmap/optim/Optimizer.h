// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_OPTIMIZER_H_
#define TRANSITMAP_OPTIM_OPTIMIZER_H_

namespace transitmapper {
namespace optim {

class Optimizer {
 public:
  virtual int optimize() const = 0;
};
}  // namespace optim
}  // namespace transitmapper

#endif  // TRANSITMAP_OPTIM_OPTIMIZER_H_
