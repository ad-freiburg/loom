// Copyright 2016
// Author: Patrick Brosi
//

#include <vector>
#include "loom/config/LoomConfig.h"
#include "loom/optim/CombOptimizer.h"
#include "shared/rendergraph/RenderGraph.h"

struct FileTest {
  std::string fname;
  size_t sameSegCrossings, diffSegCrossings, separations;
};

static const std::vector<FileTest> fileTests(
    {{"/home/patrick/repos/loom/src/loom/tests/datasets/y-splitting.json", 0, 0,
      0}});

// _____________________________________________________________________________
int main(int argc, char** argv) {
  UNUSED(argc);
  UNUSED(argv);

  loom::config::Config baseCfg;
  baseCfg.untangleGraph = false;
  baseCfg.collapseLinePartners = false;
  baseCfg.separationOpt = false;
  baseCfg.createCoreOptimGraph = false;
  baseCfg.optimRuns = 1;

  shared::rendergraph::Penalties pens{1, 1, 1, 1, 1, 1, 1, 1, false, false};

  std::vector<loom::config::Config> configs;
  configs.push_back(baseCfg);

  baseCfg.collapseLinePartners = true;
  baseCfg.createCoreOptimGraph = true;
  configs.push_back(baseCfg);

  baseCfg.untangleGraph = true;
  configs.push_back(baseCfg);

  // TODO
  // TEST(g.numEdgs() == 18);
  // TEST(g.numNds() == 18);

  for (const auto& cfg : configs) {
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
      for (const auto& test : fileTests) {
        shared::rendergraph::RenderGraph g(5, 5);

        std::ifstream input;
        input.open(test.fname);
        g.readFromJson(&input, 3);

        auto res = optim->optimize(&g);

        TEST(res.sameSegCrossings, ==, test.sameSegCrossings);
        TEST(res.diffSegCrossings, ==, test.diffSegCrossings);
        TEST(res.separations, ==, test.separations);
      }
    }
  }
}
