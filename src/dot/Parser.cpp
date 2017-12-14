// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
#include <iostream>
#include <cstring>
#include <map>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "dot/Parser.h"

using namespace dot;
using namespace parser;

// _____________________________________________________________________________
File::File(const std::string& path) : _s(NONE), _c(0), _lastBytes(0), _which(0) {
}

// _____________________________________________________________________________
bool File::has() {
  return _tagStack.size();
}

// _____________________________________________________________________________
const Tag& File::get() {
  _ret.name = 0;
  _ret.attrs.clear();
  while (_lastBytes) {
    for (; _c - _buffer[_which] < _lastBytes; ++_c) {
      switch (_s) {
        case NONE:
          if (*_c == '<') {
            _s = IN_TAG_TENTATIVE;
            continue;
          }

          if (std::isspace(*_c)) continue;
        case IN_TAG_TENTATIVE:
          if (*_c == '/') {
            _s = IN_TAG_NAME_CLOSE;
            _tmp = _c + 1;
            continue;
          }
          if (*_c == '?') {
            _s = IN_TAG_NAME_META;
            continue;
          }
          if (std::isalnum(*_c)) {
            _s = IN_TAG_NAME;
            _ret.name = _c;
            continue;
          }

        case IN_TAG:
          if (std::isspace(*_c)) continue;

          if (std::isalnum(*_c)) {
            _s = IN_ATTRKEY;
            _tmp = _c;
            continue;
          }
          if (*_c == '/') {
            _s = AW_CLOSING;
            continue;
          }
          if (*_c == '>') {
            _tagStack.push(_ret.name);
            _s = WS_SKIP;
            continue;
          }

        case IN_ATTRVAL_SQ:
          if (*_c == '\'') {
            _s = IN_TAG;
            *_c = 0;
            _ret.attrs[_tmp] = _tmp2;
            continue;
          }

          continue;

        case IN_ATTRVAL_DQ:
          if (*_c == '"') {
            _s = IN_TAG;
            *_c = 0;
            _ret.attrs[_tmp] = _tmp2;
            continue;
          }

          continue;

        case AW_IN_ATTRVAL:
          if (std::isspace(*_c)) continue;

          if (*_c == '\'') {
            _s = IN_ATTRVAL_SQ;
            _tmp2 = _c + 1;
            continue;
          }

          if (*_c == '"') {
            _s = IN_ATTRVAL_DQ;
            _tmp2 = _c + 1;
            continue;
          }

          // TODO: error!
          std::cout << "ERROR 3!" << std::endl;
          std::cout << *_c << std::endl;
          exit(1);

        case IN_ATTRKEY:
          if (std::isspace(*_c)) {
            *_c = 0;
            _s = AFTER_ATTRKEY;
            continue;
          }
          if (std::isalnum(*_c)) {
            continue;
          }
          if (*_c == '=') {
            *_c = 0;
            _s = AW_IN_ATTRVAL;
            continue;
          }

          std::cout << "ERROR 5" << std::endl;
          exit(0);

        case AFTER_ATTRKEY:
          if (std::isspace(*_c)) continue;

          if (*_c == '=') {
            _s = AW_IN_ATTRVAL;
            continue;
          }

          // TODO: error
          continue;

        case IN_TAG_NAME:
          if (std::isspace(*_c)) {
            *_c = 0;
            _s = IN_TAG;
            continue;
          }
          if (*_c == '>') {
            *_c = 0;
            _tagStack.push(_ret.name);
            _s = WS_SKIP;
            continue;
          }
          if (*_c == '/') {
            *_c = 0;
            _s = AW_CLOSING;
            continue;
          }
          if (std::isalnum(*_c)) {
            continue;
          }

        case IN_TAG_NAME_META:
          // TODO: read meta tags!
          if (*_c == '>') {
            _s = NONE;
            continue;
          }

          continue;

        case IN_TAG_NAME_CLOSE:
          if (std::isspace(*_c)) {
            *_c = 0;
            _s = IN_TAG_CLOSE;
            continue;
          }

          if (std::isalnum(*_c)) {
            continue;
          }

          if (*_c == '>') {
            *_c = 0;
            if (_tmp != _tagStack.top()) {
              // TODO: error!
              std::cout << "ERROR 6!" << std::endl;
              exit(1);
            }
            _tagStack.pop();
            _s = NONE;
            continue;
          }

        case IN_TAG_CLOSE:
          if (std::isspace(*_c)) continue;

          if (*_c == '>') {
            if (_tmp != _tagStack.top()) {
              // TODO: error!
              std::cout << "ERROR 1!" << std::endl;
              exit(1);
            }
            _tagStack.pop();
            _s = NONE;
            continue;
          }

          // TODO: error
          std::cout << "ERROR 4" << std::endl;
          exit(1);

        case AW_CLOSING:
          if (*_c == '>') {
            _s = WS_SKIP;
            continue;
          }

        case WS_SKIP:
          if (std::isspace(*_c)) continue;
          else {
            _s = NONE;
            return _ret;
          }
      }
    }

 }
}
