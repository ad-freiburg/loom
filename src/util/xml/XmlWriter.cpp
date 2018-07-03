// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <map>
#include <ostream>
#include <stack>
#include <string>
#include "XmlWriter.h"

using namespace util;
using namespace xml;

using std::ostream;
using std::string;
using std::map;

// _____________________________________________________________________________
XmlWriter::XmlWriter(std::ostream* out)
    : _out(out), _pretty(false), _indent(4) {}

// _____________________________________________________________________________
XmlWriter::XmlWriter(std::ostream* out, bool pret)
    : _out(out), _pretty(pret), _indent(4) {}

// _____________________________________________________________________________
XmlWriter::XmlWriter(std::ostream* out, bool pret, size_t indent)
    : _out(out), _pretty(pret), _indent(indent) {}

// _____________________________________________________________________________
void XmlWriter::openTag(const string& tag, const map<string, string>& attrs) {
  if (!_nstack.empty() && _nstack.top().t == COMMENT) {
    throw XmlWriterException("Opening tags not allowed while inside comment.");
  }

  checkTagName(tag);
  closeHanging();
  doIndent();

  *_out << "<" << tag;

  for (auto kv : attrs) {
    *_out << " ";
    putEsced(_out, kv.first, '"');
    *_out << "=\"";
    putEsced(_out, kv.second, '"');
    *_out << "\"";
  }

  _nstack.push(XmlNode(TAG, tag, true));
}

// _____________________________________________________________________________
void XmlWriter::openTag(const string& tag) {
  openTag(tag, map<string, string>());
}

// _____________________________________________________________________________
void XmlWriter::openTag(const string& tag, const string& k, const string& v) {
  map<string, string> kv;
  kv[k] = v;
  openTag(tag, kv);
}

// _____________________________________________________________________________
void XmlWriter::openComment() {
  // don't allow nested comments
  if (!_nstack.empty() && _nstack.top().t == COMMENT) return;

  closeHanging();
  doIndent();

  *_out << "<!-- ";

  _nstack.push(XmlNode(COMMENT, "", false));
}

// _____________________________________________________________________________
void XmlWriter::writeText(const string& text) {
  if (_nstack.empty()) {
    throw XmlWriterException("Text content not allowed in prolog / trailing.");
  }
  closeHanging();
  doIndent();
  putEsced(_out, text, ' ');
}

// _____________________________________________________________________________
void XmlWriter::closeTag() {
  while (!_nstack.empty() && _nstack.top().t == TEXT) _nstack.pop();

  if (_nstack.empty()) return;

  if (_nstack.top().t == COMMENT) {
    _nstack.pop();
    doIndent();
    *_out << " -->";
  } else if (_nstack.top().t == TAG) {
    if (_nstack.top().hanging) {
      *_out << " />";
      _nstack.pop();
    } else {
      string tag = _nstack.top().pload;
      _nstack.pop();
      doIndent();
      *_out << "</" << tag << ">";
    }
  }
}

// _____________________________________________________________________________
void XmlWriter::closeTags() {
  while (!_nstack.empty()) closeTag();
}

// _____________________________________________________________________________
void XmlWriter::doIndent() {
  if (_pretty) {
    *_out << std::endl;
    for (size_t i = 0; i < _nstack.size() * _indent; i++) *_out << " ";
  }
}

// _____________________________________________________________________________
void XmlWriter::closeHanging() {
  if (_nstack.empty()) return;

  if (_nstack.top().hanging) {
    *_out << ">";
    _nstack.top().hanging = false;
  } else if (_nstack.top().t == TEXT) {
    _nstack.pop();
  }
}

// _____________________________________________________________________________
void XmlWriter::putEsced(ostream* out, const string& str, char quot) {
  if (!_nstack.empty() && _nstack.top().t == COMMENT) {
    *out << str;
    return;
  }

  for (const char& c : str) {
    if (quot == '"' && c == '"')
      *out << "&quot;";
    else if (quot == '\'' && c == '\'')
      *out << "&apos;";
    else if (c == '<')
      *out << "&lt;";
    else if (c == '>')
      *out << "&gt;";
    else if (c == '&')
      *out << "&amp;";
    else
      *out << c;
  }
}

// _____________________________________________________________________________
void XmlWriter::checkTagName(const string& str) const {
  if (!isalpha(str[0]) && str[0] != '_')
    throw XmlWriterException(
        "XML elements must start with either a letter "
        "or an underscore");

  string begin = str.substr(0, 3);
  std::transform(begin.begin(), begin.end(), begin.begin(), ::tolower);
  if (begin == "xml")
    throw XmlWriterException(
        "XML elements cannot start with"
        " XML, xml, Xml etc.");

  for (const char& c : str) {
    // we allow colons in tag names for primitive namespace support
    if (!isalpha(c) && !isdigit(c) && c != '-' && c != '_' && c != '.' &&
        c != ':')
      throw XmlWriterException(
          "XML elements can only contain letters, "
          "digits, hyphens, underscores and periods.");
  }
}
