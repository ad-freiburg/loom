// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_OPTIM_ILPSOLVPROV_H_
#define SHARED_OPTIM_ILPSOLVPROV_H_

#include "shared/optim/COINSolver.h"
#include "shared/optim/GLPKSolver.h"
#include "shared/optim/GurobiSolver.h"
#include "shared/optim/ILPSolver.h"
#include "util/log/Log.h"

namespace shared {
namespace optim {

class ILPProviderErr : public std::exception {
 public:
  ILPProviderErr(const std::string& msg) : _msg(msg) {}
  const char* what() { return _msg.c_str(); }

 private:
  std::string _msg;
};

inline ILPSolver* getSolver(std::string wish, shared::optim::DirType dir) {
  UNUSED(dir);  // prevent warning if no ILP solver is present
  bool force = (wish.back() == '!');

  if (force) {
    wish.pop_back();
  }
  ILPSolver* lp = 0;

  // first try to consider wish

  // aliases
  if (wish == "cbc") wish = "coin";

  try {
#ifdef GUROBI_FOUND
    if (wish == "gurobi") lp = new shared::optim::GurobiSolver(dir);
#endif

#if COIN_FOUND
    if (wish == "coin") lp = new shared::optim::COINSolver(dir);
#endif

#if GLPK_FOUND
    if (wish == "glpk") lp = new shared::optim::GLPKSolver(dir);
#endif
  } catch (std::exception& e) {
    LOG(ERROR) << e.what();
  }

  if (!lp && force) {
    throw ILPProviderErr(
        "The configured ILP solver could not be provided");
  }

  // fallbacks

#ifdef GUROBI_FOUND
  try {
    if (!lp && wish != "gurobi") lp = new shared::optim::GurobiSolver(dir);
  } catch (std::exception& e) {
    LOG(ERROR) << e.what();
  }
#endif

#ifdef COIN_FOUND
  try {
    if (!lp && wish != "coin") lp = new shared::optim::COINSolver(dir);
  } catch (std::exception& e) {
    LOG(ERROR) << e.what();
  }
#endif

#if GLPK_FOUND
  try {
    if (!lp && wish != "glpk") lp = new shared::optim::GLPKSolver(dir);
  } catch (std::exception& e) {
    LOG(ERROR) << e.what();
  }
#endif

  if (!lp) {
    throw ILPProviderErr("No ILP solver found.");
  }

  return lp;
}

}  // namespace optim
}  // namespace shared

#endif
