// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_SCORER_H_
#define TRANSITMAP_OPTIM_SCORER_H_

#include <string>
#include <vector>
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/graph/Node.h"
#include "transitmap/graph/Penalties.h"

using transitmapper::graph::TransitGraph;
using transitmapper::graph::Node;
using transitmapper::graph::Penalties;
using transitmapper::graph::OrderingConfig;

namespace transitmapper {
namespace optim {

class Scorer {
 public:
  Scorer(const TransitGraph* g) : _g(g) {}

  double getScore(const Penalties& pens) const;
  double getScore(const Penalties& pens, const OrderingConfig& c) const;

  double getCrossScore(const Penalties& pens) const;
  double getCrossScore(const Penalties& pens, const OrderingConfig& c) const;

  double getSeparationScore(const Penalties& pens) const;
  double getSeparationScore(const Penalties& pens, const OrderingConfig& c) const;

  size_t getNumCrossings() const;
  size_t getNumCrossings(const OrderingConfig& c) const;

  size_t getNumSeparations() const;
  size_t getNumSeparations(const OrderingConfig& c) const;

  double getNumPossSolutions() const;
  size_t getMaxCrossPenalty() const;
  size_t getMaxSplitPenalty() const;

  double getScore(const Node* n, const Penalties& pens, bool crossAdjPen, bool splitAdjPen, const graph::OrderingConfig& cfg) const;
  size_t getNumCrossings(const Node* n, const graph::OrderingConfig& c) const;
  size_t getNumSeparations(const Node* n, const graph::OrderingConfig& c) const;
  size_t getSeparationScore(const Node* n, const OrderingConfig& c, const Penalties& pens, bool adjpen) const;
  double getCrossingScore(const Node* n, const OrderingConfig& c, const Penalties& pens, bool adjpen) const;

  int getCrossingPenaltySameSeg(const Node* n, const Penalties& pens, bool adjpen) const;
  int getCrossingPenaltyDiffSeg(const Node* n, const Penalties& pens, bool adjpen) const;
  int getSplittingPenalty(const Node* n, const Penalties& pens, bool adjpen) const;

 private:
 	const TransitGraph* _g;
};

}}

#endif  // TRANSITMAP_GRAPH_ROUTE_H_
