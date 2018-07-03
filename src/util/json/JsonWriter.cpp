// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <iomanip>
#include "JsonWriter.h"
#include "util/String.h"
using namespace util;
using namespace json;

using std::ostream;
using std::string;
using std::map;

// _____________________________________________________________________________
JsonWriter::JsonWriter(std::ostream* out)
    : _out(out), _pretty(false), _indent(2) {
  *_out << std::setprecision(10);
}

// _____________________________________________________________________________
JsonWriter::JsonWriter(std::ostream* out, size_t prec)
    : _out(out), _pretty(false), _indent(2) {
  *_out << std::setprecision(prec);
}

// _____________________________________________________________________________
JsonWriter::JsonWriter(std::ostream* out, size_t prec, bool pret)
    : _out(out), _pretty(pret), _indent(2) {
  *_out << std::setprecision(prec);
}

// _____________________________________________________________________________
JsonWriter::JsonWriter(std::ostream* out, size_t prec, bool pret, size_t indent)
    : _out(out), _pretty(pret), _indent(indent) {
  *_out << std::setprecision(prec);
}

// _____________________________________________________________________________
void JsonWriter::obj() {
  if (!_stack.empty() && _stack.top().type == OBJ)
    throw JsonWriterException("Object not allowed as key");
  if (!_stack.empty() && _stack.top().type == KEY) _stack.pop();
  if (!_stack.empty() && _stack.top().type == ARR) valCheck();
  if (_stack.size() && _stack.top().type == ARR) prettor();
  (*_out) << "{";
  _stack.push({OBJ, 1});
}

// _____________________________________________________________________________
void JsonWriter::obj(const std::map<std::string, std::string> kvs) {
  obj();
  for (const auto& kv : kvs) {
    key(kv.first);
    val(kv.second);
  }
  close();
}

// _____________________________________________________________________________
void JsonWriter::key(const std::string& k) {
  if (_stack.empty() || _stack.top().type != OBJ)
    throw JsonWriterException("Keys only allowed in objects.");
  if (!_stack.top().empty) (*_out) << "," << (_pretty ? " " : "");
  _stack.top().empty = 0;
  prettor();
  (*_out) << "\"" << k << "\""
          << ":" << (_pretty ? " " : "");
  _stack.push({KEY, 1});
}

// _____________________________________________________________________________
void JsonWriter::keyVal(const std::string& k, const std::string& v) {
  key(k);
  val(v);
}

// _____________________________________________________________________________
void JsonWriter::valCheck() {
  if (_stack.empty() || (_stack.top().type != KEY && _stack.top().type != ARR))
    throw JsonWriterException("Value not allowed here.");
  if (!_stack.empty() && _stack.top().type == KEY) _stack.pop();
  if (!_stack.empty() && _stack.top().type == ARR) {
    if (!_stack.top().empty) (*_out) << "," << (_pretty ? " " : "");
    _stack.top().empty = 0;
  }
}

// _____________________________________________________________________________
void JsonWriter::val(const std::string& v) {
  valCheck();
  (*_out) << "\"" << util::jsonStringEscape(v) << "\"";
}

// _____________________________________________________________________________
void JsonWriter::val(int v) {
  valCheck();
  (*_out) << v;
}

// _____________________________________________________________________________
void JsonWriter::val(double v) {
  valCheck();
  (*_out) << v;
}

// _____________________________________________________________________________
void JsonWriter::arr() {
  if (!_stack.empty() && _stack.top().type == OBJ)
    throw JsonWriterException("Array not allowed as key");
  if (!_stack.empty() && _stack.top().type == KEY) _stack.pop();
  if (!_stack.empty() && _stack.top().type == ARR) valCheck();
  (*_out) << "[";
  _stack.push({ARR, 1});
}

// _____________________________________________________________________________
void JsonWriter::prettor() {
  if (_pretty) {
    (*_out) << "\n";
    for (size_t i = 0; i < _indent * _stack.size(); i++) (*_out) << " ";
  }
}

// _____________________________________________________________________________
void JsonWriter::closeAll() {
  while (!_stack.empty()) close();
}

// _____________________________________________________________________________
void JsonWriter::close() {
  if (_stack.empty()) return;
  switch (_stack.top().type) {
    case OBJ:
      _stack.pop();
      prettor();
      (*_out) << "}";
      break;
    case ARR:
      _stack.pop();
      (*_out) << "]";
      break;
    case KEY:
      throw JsonWriterException("Missing value.");
  }
}
