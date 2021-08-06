// Copyright 2016
// Author: Patrick Brosi
//

#include "shared/rendergraph/RenderGraph.h"
#include "loom/optim/CombOptimizer.h"
#include "loom/config/LoomConfig.h"

// _____________________________________________________________________________
int main(int argc, char** argv) {
  UNUSED(argc);
  UNUSED(argv);

  loom::config::Config cfg;
  cfg.untangleGraph = false;
  cfg.collapseLinePartners = false;
  cfg.separationOpt = false;
  cfg.createCoreOptimGraph = false;
  cfg.optimRuns = 1;

  shared::rendergraph::Penalties pens{1, 1, 1, 1, 1, 1, 1, 1, false, false};

  // TODO
  // TEST(g.numEdgs() == 18);
  // TEST(g.numNds() == 18);

  loom::optim::ExhaustiveOptimizer exhausOptim(&cfg, pens);
  loom::optim::ILPOptimizer ilpOptim(&cfg, pens);
  loom::optim::ILPEdgeOrderOptimizer ilpImprOptim(&cfg, pens);
  loom::optim::CombOptimizer combOptim(&cfg, pens);

  std::vector<loom::optim::Optimizer*> optimizers;
  optimizers.push_back(&exhausOptim);
  optimizers.push_back(&ilpOptim);
  optimizers.push_back(&ilpImprOptim);
  optimizers.push_back(&combOptim);

  for (auto optim : optimizers) {
    shared::rendergraph::RenderGraph g(5, 5);

    std::ifstream input;
    input.open("/home/patrick/repos/loom/src/loom/tests/datasets/y-splitting.json");
    g.readFromJson(&input, 3);

    auto res = optim->optimize(&g);

    TEST(res.sameSegCrossings, ==, 0);
    TEST(res.diffSegCrossings, ==, 0);
    TEST(res.separations, ==, 0);
  }
}
