// Copyright 2016
// Author: Patrick Brosi
//

#include <vector>
#include "loom/config/LoomConfig.h"
#include "loom/optim/CombOptimizer.h"
#include "shared/rendergraph/RenderGraph.h"

struct FileTest {
  std::string fname;
  size_t sameSegCrossings, diffSegCrossings, separations, numNds, numTopoNds,
      numEdgs;
};

static const std::vector<FileTest> fileTests(
    {{
         "/home/patrick/repos/loom/src/loom/tests/datasets/single.json",
         0,  // number same segment crossings
         0,  // number diff segment crossings
         0,  // number separations
         2,  // number nodes
         0,  // number topological nodes
         1   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/simplify.json",
         0,  // number same segment crossings
         0,  // number diff segment crossings
         0,  // number separations
         4,  // number nodes
         0,  // number topological nodes
         3   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/simplify-2.json",
         0,  // number same segment crossings
         0,  // number diff segment crossings
         0,  // number separations
         6,  // number nodes
         0,  // number topological nodes
         5   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/simplify-3.json",
         0,  // number same segment crossings
         0,  // number diff segment crossings
         0,  // number separations
         6,  // number nodes
         0,  // number topological nodes
         5   // number edges
     },
     // Y splitting
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/y-splitting.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         28,  // number nodes
         0,   // number topological nodes
         22   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/y-splitting-2.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         23,  // number nodes
         0,   // number topological nodes
         18   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/y-splitting-3.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         38,  // number nodes
         0,   // number topological nodes
         30   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/y-splitting-4.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         45,  // number nodes
         0,   // number topological nodes
         37   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/y-splitting-5.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         24,  // number nodes
         0,   // number topological nodes
         20   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/y-splitting-6.json",
         0,  // number same segment crossings
         0,  // number diff segment crossings
         0,  // number separations
         4,  // number nodes
         0,  // number topological nodes
         3   // number edges
     },

     // partial Y splitting
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "y-splitting-partial.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         28,  // number nodes
         0,   // number topological nodes
         22   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "y-splitting-partial-2.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         23,  // number nodes
         0,   // number topological nodes
         18   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "y-splitting-partial-3.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         38,  // number nodes
         0,   // number topological nodes
         30   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "y-splitting-partial-4.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         47,  // number nodes
         8,   // number topological nodes
         39   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "y-splitting-rec.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         24,  // number nodes
         0,   // number topological nodes
         20   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "y-splitting-rec-2.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         24,  // number nodes
         0,   // number topological nodes
         20   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "y-splitting-rec-3.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         15,  // number nodes
         0,   // number topological nodes
         12   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "y-splitting-rec-4.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         16,  // number nodes
         0,   // number topological nodes
         13   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "dog-bone-splitting.json",
         0,   // number same segment crossings
         2,   // number diff segment crossings
         0,   // number separations
         12,  // number nodes
         0,   // number topological nodes
         10   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "dog-bone-splitting-2.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         33,  // number nodes
         0,   // number topological nodes
         28   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "dog-bone-splitting-3.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         66,  // number nodes
         0,   // number topological nodes
         59   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "dog-bone-splitting-4.json",
         0,   // number same segment crossings
         18,  // number diff segment crossings
         0,   // number separations
         66,  // number nodes
         5,   // number topological nodes
         59   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "dog-bone-splitting-partial.json",
         0,   // number same segment crossings
         2,   // number diff segment crossings
         0,   // number separations
         12,  // number nodes
         0,   // number topological nodes
         10   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "dog-bone-splitting-partial-2.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         33,  // number nodes
         0,   // number topological nodes
         28   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "dog-bone-splitting-partial-3.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         66,  // number nodes
         0,   // number topological nodes
         59   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "dog-bone-splitting-partial-4.json",
         0,   // number same segment crossings
         21,  // number diff segment crossings
         0,   // number separations
         86,  // number nodes
         7,   // number topological nodes
         77   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "dog-bone-splitting-partial-5.json",
         0,  // number same segment crossings
         0,  // number diff segment crossings
         0,  // number separations
         8,  // number nodes
         1,  // number topological nodes
         6   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/"
         "dog-bone-splitting-partial-6.json",
         0,   // number same segment crossings
         2,   // number diff segment crossings
         0,   // number separations
         12,  // number nodes
         0,   // number topological nodes
         10   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/full-cross.json",
         0,   // number same segment crossings
         0,   // number diff segment crossings
         0,   // number separations
         17,  // number nodes
         0,   // number topological nodes
         14   // number edges
     },
     {
         "/home/patrick/repos/loom/src/loom/tests/datasets/full-cross-2.json",
         1,   // number same segment crossings
         3,   // number diff segment crossings
         0,   // number separations
         19,  // number nodes
         0,   // number topological nodes
         16   // number edges
     }

    });

// _____________________________________________________________________________
int main(int argc, char** argv) {
  UNUSED(argc);
  UNUSED(argv);

  loom::config::Config baseCfg;
  baseCfg.untangleGraph = false;
  baseCfg.collapseLinePartners = false;
  baseCfg.createCoreOptimGraph = false;
  baseCfg.optimRuns = 1;

  shared::rendergraph::Penalties pens{1, 0, 1, 1, 0, 1, 1, 0, false, false};

  // without separation penalty

  std::vector<loom::config::Config> configs;
  configs.push_back(baseCfg);

  baseCfg.collapseLinePartners = true;
  baseCfg.createCoreOptimGraph = true;
  configs.push_back(baseCfg);

  baseCfg.untangleGraph = true;
  configs.push_back(baseCfg);

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

        TEST(g.numEdgs(), ==, test.numEdgs);
        TEST(g.numNds(), ==, test.numNds);
        TEST(g.numNds(true), ==, test.numTopoNds);

        if (optim == &exhausOptim && g.searchSpaceSize() > 50000) continue;
        if (optim == &ilpOptim && g.searchSpaceSize() > 500000) continue;

        auto res = optim->optimize(&g);

        TEST(res.sameSegCrossings, ==, test.sameSegCrossings);
        TEST(res.diffSegCrossings, ==, test.diffSegCrossings);
      }
    }
  }

  // miscellaneous
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
      shared::rendergraph::RenderGraph g(5, 5);

      std::ifstream input;
      input.open(
          "/home/patrick/repos/loom/src/loom/tests/datasets/"
          "freiburg-tram.json");
      g.readFromJson(&input, 3);

      TEST(g.numEdgs(), ==, 78);
      TEST(g.numNds(), ==, 77);
      TEST(g.numNds(true), ==, 5);

      if (optim == &exhausOptim && g.searchSpaceSize() > 50000) continue;
      if (optim == &ilpOptim && g.searchSpaceSize() > 500000) continue;

      auto res = optim->optimize(&g);

      TEST(res.sameSegCrossings + res.diffSegCrossings, ==, 4);
      TEST(res.score, ==, 4);
    }
  }

  {
    shared::rendergraph::Penalties pensLoc = pens;
    pensLoc.diffSegCrossPen = 100;
    pensLoc.inStatCrossPenDiffSeg = 200;

    for (const auto& cfg : configs) {
      loom::optim::ExhaustiveOptimizer exhausOptim(&cfg, pensLoc);
      loom::optim::ILPOptimizer ilpOptim(&cfg, pensLoc);
      loom::optim::ILPEdgeOrderOptimizer ilpImprOptim(&cfg, pensLoc);
      loom::optim::CombOptimizer combOptim(&cfg, pensLoc);

      std::vector<loom::optim::Optimizer*> optimizers;
      optimizers.push_back(&exhausOptim);
      optimizers.push_back(&ilpOptim);
      optimizers.push_back(&ilpImprOptim);
      optimizers.push_back(&combOptim);

      for (auto optim : optimizers) {
        shared::rendergraph::RenderGraph g(5, 5);

        std::ifstream input;
        input.open(
            "/home/patrick/repos/loom/src/loom/tests/datasets/"
            "freiburg-tram.json");
        g.readFromJson(&input, 3);

        TEST(g.numEdgs(), ==, 78);
        TEST(g.numNds(), ==, 77);
        TEST(g.numNds(true), ==, 5);

        if (optim == &exhausOptim && g.searchSpaceSize() > 50000) continue;
        if (optim == &ilpOptim && g.searchSpaceSize() > 500000) continue;

        auto res = optim->optimize(&g);

        TEST(res.sameSegCrossings + res.diffSegCrossings, ==, 4);
        TEST(res.diffSegCrossings, ==, 0);
        TEST(res.score, ==, 4);
      }
    }
  }

  {
    shared::rendergraph::Penalties pensLoc = pens;
    pensLoc.diffSegCrossPen = 100;
    pensLoc.inStatCrossPenDiffSeg = 200;
    pensLoc.inStatCrossPenSameSeg = 5;
    pensLoc.inStatCrossPenDegTwo = 5;
    pensLoc.crossAdjPen = true;

    for (const auto& cfg : configs) {
      loom::optim::ExhaustiveOptimizer exhausOptim(&cfg, pensLoc);
      loom::optim::ILPOptimizer ilpOptim(&cfg, pensLoc);
      loom::optim::ILPEdgeOrderOptimizer ilpImprOptim(&cfg, pensLoc);
      loom::optim::CombOptimizer combOptim(&cfg, pensLoc);

      std::vector<loom::optim::Optimizer*> optimizers;
      optimizers.push_back(&exhausOptim);
      optimizers.push_back(&ilpOptim);
      optimizers.push_back(&ilpImprOptim);
      optimizers.push_back(&combOptim);

      for (auto optim : optimizers) {
        shared::rendergraph::RenderGraph g(5, 5);

        std::ifstream input;
        input.open(
            "/home/patrick/repos/loom/src/loom/tests/datasets/"
            "freiburg-tram.json");
        g.readFromJson(&input, 3);

        TEST(g.numEdgs(), ==, 78);
        TEST(g.numNds(), ==, 77);
        TEST(g.numNds(true), ==, 5);

        if (optim == &exhausOptim && g.searchSpaceSize() > 50000) continue;
        if (optim == &ilpOptim && g.searchSpaceSize() > 500000) continue;

        auto res = optim->optimize(&g);

        TEST(res.sameSegCrossings + res.diffSegCrossings, ==, 4);
        TEST(res.diffSegCrossings, ==, 0);
        TEST(res.score, ==, 18);
      }
    }
  }

  // with separation penalty

  pens.inStatSplitPenDegTwo = 1;
  pens.inStatSplitPen = 1;
  pens.splitPen = 1;

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

        TEST(g.numEdgs(), ==, test.numEdgs);
        TEST(g.numNds(), ==, test.numNds);
        TEST(g.numNds(true), ==, test.numTopoNds);

        if (optim == &exhausOptim && g.searchSpaceSize() > 50000) continue;
        if (optim == &ilpOptim && g.searchSpaceSize() > 500000) continue;

        auto res = optim->optimize(&g);

        TEST(res.sameSegCrossings, ==, test.sameSegCrossings);
        TEST(res.diffSegCrossings, ==, test.diffSegCrossings);
        TEST(res.separations, ==, test.separations);
      }
    }
  }

  {
    shared::rendergraph::Penalties pensLoc = pens;
    pensLoc.diffSegCrossPen = 100;
    pensLoc.inStatCrossPenDiffSeg = 200;
    pensLoc.inStatCrossPenSameSeg = 5;
    pensLoc.inStatCrossPenDegTwo = 5;
    pensLoc.inStatSplitPenDegTwo = 300;
    pensLoc.inStatSplitPen = 300;
    pensLoc.splitPen = 500;
    pensLoc.crossAdjPen = true;
    pensLoc.splitAdjPen = true;

    for (const auto& cfg : configs) {
      loom::optim::ExhaustiveOptimizer exhausOptim(&cfg, pensLoc);
      loom::optim::ILPOptimizer ilpOptim(&cfg, pensLoc);
      loom::optim::ILPEdgeOrderOptimizer ilpImprOptim(&cfg, pensLoc);
      loom::optim::CombOptimizer combOptim(&cfg, pensLoc);

      std::vector<loom::optim::Optimizer*> optimizers;
      optimizers.push_back(&exhausOptim);
      optimizers.push_back(&ilpOptim);
      optimizers.push_back(&ilpImprOptim);
      optimizers.push_back(&combOptim);

      for (auto optim : optimizers) {
        shared::rendergraph::RenderGraph g(5, 5);

        std::ifstream input;
        input.open(
            "/home/patrick/repos/loom/src/loom/tests/datasets/"
            "freiburg-tram.json");
        g.readFromJson(&input, 3);

        TEST(g.numEdgs(), ==, 78);
        TEST(g.numNds(), ==, 77);
        TEST(g.numNds(true), ==, 5);

        if (optim == &exhausOptim && g.searchSpaceSize() > 50000) continue;
        if (optim == &ilpOptim && g.searchSpaceSize() > 500000) continue;

        auto res = optim->optimize(&g);

        TEST(res.sameSegCrossings + res.diffSegCrossings, ==, 4);
        TEST(res.diffSegCrossings, ==, 0);
        TEST(res.separations, ==, 2);
        TEST(res.score, ==, 620);
      }
    }
  }
}
