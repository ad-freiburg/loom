// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include "dot/Parser.h"

using namespace dot;
using namespace parser;

// _____________________________________________________________________________
Parser::Parser(std::istream* is) : _is(is), _s(NONE), _has(true) {}

// _____________________________________________________________________________
bool Parser::has() { return _has; }

// _____________________________________________________________________________
const Entity& Parser::get() {
  _ret.ids.clear();
  _ret.attrs.clear();
  _ret.type = EMPTY;
  std::string tmp;

  while (_is->get(_c)) {
    std::cout << "state: " << _s << " char: " << _c << std::endl;
    switch (_s) {
      case NONE:
        if (std::isspace(_c)) continue;
        tmp.clear();
        if (_c == 'd' || _c == 'D') {
          _s = KW_DIRGRAPH;
          continue;
        }
        if (_c == 'g' || _c == 'G') {
          _s = KW_GRAPH;
          continue;
        }
        if (_c == 's' || _c == 'S') {
          _s = KW_STRICT;
          continue;
        }

        std::cerr << "Expected keywords 'strict', 'graph' or 'digraph'" << std::endl;
        exit(1);
      case KW_STRICT:
        tmp += _c;
        while (_is->get(_c)) {
          if (std::isalpha(_c)) {
            tmp += std::tolower(_c);
          } else break;
        }

        if (tmp == "trict" && std::isspace(_c)) {
          _ret.graphType = STRICT_GRAPH;
          tmp.clear();
          _s = AW_KW_GRAPH;
          continue;
        }

        std::cerr << "Expected keyword 'strict'" << std::endl;
        exit(1);
      case KW_GRAPH:
        tmp += _c;
        while (_is->get(_c)) {
          if (std::isalpha(_c)) {
            tmp += std::tolower(_c);
            if (tmp == "raph") break;
          } else break;
        }

        if (tmp == "raph") {
          if (_ret.graphType == STRICT_GRAPH) {
            _ret.graphType = STRICT_GRAPH;
          } else {
            _ret.graphType = GRAPH;
          }
          tmp.clear();
          _s = AW_GRAPH_ID;
          continue;
        }

        std::cerr << "Expected keyword 'graph'" << std::endl;
        exit(1);
      case KW_DIRGRAPH:
        tmp += _c;
        while (_is->get(_c)) {
          if (std::isalpha(_c)) {
            tmp += std::tolower(_c);
            if (tmp == "igraph") break;
          } else break;
        }

        if (tmp == "igraph") {
          if (_ret.graphType == STRICT_GRAPH) {
            _ret.graphType = STRICT_DIGRAPH;
          } else {
            _ret.graphType = DIGRAPH;
          }
          tmp.clear();
          _s = AW_GRAPH_ID;
          continue;
        }

        std::cerr << tmp << std::endl;

        std::cerr << "Expected keyword 'digraph'" << std::endl;
        exit(1);
      case AW_KW_GRAPH:
        if (std::isspace(_c)) continue;
        if (_c == 'g' || _c == 'G') {
          _s = KW_GRAPH;
          continue;
        }
        if (_c == 'd' || _c == 'D') {
          _s = KW_DIRGRAPH;
          continue;
        }
      case AW_GRAPH_ID:
        if (std::isspace(_c)) continue;
        if (isIDChar(_c)) {
          tmp += _c;
          _s = IN_GRAPH_ID;
          continue;
        } else if (_c == '{') {
          _s = AW_LL_ID;
          continue;
        }

        std::cerr << "Expected graph id or opening {" << std::endl;
        exit(1);
      case IN_GRAPH_ID:
        if (isIDChar(_c)) {
          tmp += _c;
          continue;
        }

        if (_c == '{') {
          _ret.graphName = tmp;
          tmp.clear();
          _s = AW_LL_ID;
          continue;
        }

        if (std::isspace(_c)) {
          _ret.graphName = tmp;
          tmp.clear();
          _s = AW_OPEN;
          continue;
        }

        std::cerr << "Expected ID characters" << std::endl;
        exit(1);
      case AW_OPEN:
        if (std::isspace(_c)) continue;
        if (_c == '{') {
          tmp.clear();
          _s = AW_LL_ID;
          continue;
        }
        std::cerr << "Expected opening {" << std::endl;
        exit(1);
      case AW_LL_ID:
        if (std::isspace(_c)) continue;
        if (isIDChar(_c)) {
          tmp += _c;
          _s = IN_ID;
          continue;
        }
        if (_c == '"') {
          _s = IN_QUOTED_ID;
          continue;
        }
        if (_c == '}') {
          _has = false;
          return _ret;
        }

        std::cerr << "Expected statement" << std::endl;
        exit(1);

      case AW_ID:
        if (std::isspace(_c)) continue;
        if (isIDChar(_c)) {
          tmp += _c;
          _s = IN_ID;
          continue;
        }
        if (_c == '"') {
          _s = IN_QUOTED_ID;
          continue;
        }

        std::cerr << "Expected ID" << std::endl;
        exit(1);

      case IN_ID:
        if (isIDChar(_c)) {
          tmp += _c;
          continue;
        }
        _ret.ids.push_back(tmp);
        tmp.clear();

        _s = AW_STMT_DEC;
        // slip down

      case AW_STMT_DEC:
        if (std::isspace(_c)) continue;
        if (_c == '=') {
          if (_ret.type != EMPTY) {
            std::cerr << "Syntax error." << std::endl;
            exit(1);
          }
          tmp.clear();
          _ret.type = ATTR;
          _s = AW_ID;
          continue;
        }

        if (_c == '-') {
          if (_ret.type != EMPTY && _ret.type != EDGE) {
            std::cerr << "Syntax error." << std::endl;
            exit(1);
          }
          _s = AW_EDGE_OP;
          continue;
        }

        if (_c == '[') {
          _s = AW_ATTR_KEY;
          continue;
        }

        if (_c == '"') {
          _s = IN_QUOTED_ID;
          return _ret;
        }

        if (_c == ';' || _c == ',') {
          _s = AW_LL_ID;
          return _ret;
        }

        if (_c == '}') {
          _has = false;
          return _ret;
        }

        std::cerr << "Expected edge operator, = or [" << std::endl;
        exit(1);

      case IN_QUOTED_ID:
        if (_c != '"' || tmp.back() == '\\') {
          tmp += _c;
          continue;
        } else {
          _ret.ids.push_back(tmp);
          tmp.clear();
          _s = AW_STMT_DEC;
          continue;
        }

      case AW_EDGE_OP:
        if (_c == '>') {
          if (_ret.graphType == STRICT_GRAPH || _ret.graphType == GRAPH) {
            std::cerr << "No directed edges allowed in undirected graph." << std::endl;
            exit(1);
          }
          _ret.type = EDGE;
          _s = AW_ID;
          continue;
        }

        if (_c == '-') {
          if (_ret.graphType == STRICT_DIGRAPH || _ret.graphType == DIGRAPH) {
            std::cerr << "No undirected edges allowed in directed graph." << std::endl;
            exit(1);
          }
          _ret.type = EDGE;
          _s = AW_ID;
          continue;
        }

        std::cerr << "Expected edge operator" << std::endl;
        exit(1);

    }
  }
}

// _____________________________________________________________________________
bool Parser::isIDChar(const char& c) const { return std::isalnum(c) || c == '_'; }
