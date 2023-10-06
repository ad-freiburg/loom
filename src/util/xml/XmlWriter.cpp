// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <map>
#include <cstring>
#include <ostream>
#include <stack>
#include <string>
#ifdef BZLIB_FOUND
#include <bzlib.h>
#endif

#include "XmlWriter.h"

using namespace util;
using namespace xml;

using std::map;
using std::ostream;
using std::string;

// _____________________________________________________________________________
XmlWriter::XmlWriter(std::ostream* out)
    : _out(out), _pretty(false), _indent(4), _gzfile(0), _bzfile(0) {}

// _____________________________________________________________________________
XmlWriter::XmlWriter(std::ostream* out, bool pret)
    : _out(out), _pretty(pret), _indent(4), _gzfile(0), _bzfile(0) {}

// _____________________________________________________________________________
XmlWriter::XmlWriter(std::ostream* out, bool pret, size_t indent)
    : _out(out), _pretty(pret), _indent(indent), _gzfile(0), _bzfile(0) {}

// _____________________________________________________________________________
XmlWriter::XmlWriter(const std::string& file)
    : XmlWriter::XmlWriter(file, false, 4) {}

// _____________________________________________________________________________
XmlWriter::XmlWriter(const std::string& file, bool pret)
    : XmlWriter::XmlWriter(file, pret, 4) {}

// _____________________________________________________________________________
XmlWriter::XmlWriter(const std::string& file, bool pret, size_t indent)
    : _out(0), _pretty(pret), _indent(indent), _gzfile(0), _bzfile(0) {
  if (file.size() > 2 && file[file.size() - 1] == 'z' &&
      file[file.size() - 2] == 'g' && file[file.size() - 3] == '.') {
#ifdef ZLIB_FOUND
    _gzfile = gzopen(file.c_str(), "w");
    if (_gzfile == Z_NULL) {
      throw std::runtime_error("Could not open file for writing.");
    }
#else
    throw std::runtime_error(
        "Could not open gzip file for writing, was compiled without gzip "
        "support");
#endif
  } else if (file.size() > 3 && file[file.size() - 1] == '2' &&
             file[file.size() - 2] == 'z' && file[file.size() - 3] == 'b' &&
             file[file.size() - 4] == '.') {
#ifdef BZLIB_FOUND
    _bzbuf = new char[BUFFER_S];

    FILE* f = fopen(file.c_str(), "w");
    int err;
    if (!f) throw std::runtime_error("Could not open file for writing.");

    _bzfile = BZ2_bzWriteOpen(&err, f, 9, 0, 30);

    if (err != BZ_OK) {
      throw std::runtime_error("Could not open bzip file for writing.");
    }
#else
    throw std::runtime_error(
        "Could not open bzip file for writing, was compiled without bzip "
        "support");
#endif
  } else {
    _outs.open(file);
    if (_outs.fail()) {
      throw std::runtime_error("Could not open file for writing.");
    }
    _out = &_outs;
  }
}

// _____________________________________________________________________________
void XmlWriter::openTag(const string& tag, const map<string, string>& attrs) {
  if (!_nstack.empty() && _nstack.top().t == COMMENT) {
    throw XmlWriterException("Opening tags not allowed while inside comment.");
  }

  checkTagName(tag);
  closeHanging();
  doIndent();

  put("<");
  put(tag);

  for (auto kv : attrs) {
    put(" ");
    putEsced(kv.first, '"');
    put("=\"");
    putEsced(kv.second, '"');
    put("\"");
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

  put("<!-- ");

  _nstack.push(XmlNode(COMMENT, "", false));
}

// _____________________________________________________________________________
void XmlWriter::writeText(const string& text) {
  if (_nstack.empty()) {
    throw XmlWriterException("Text content not allowed in prolog / trailing.");
  }
  closeHanging();
  doIndent();
  putEsced(text, ' ');
}

// _____________________________________________________________________________
void XmlWriter::closeTag() {
  while (!_nstack.empty() && _nstack.top().t == TEXT) _nstack.pop();

  if (_nstack.empty()) return;

  if (_nstack.top().t == COMMENT) {
    _nstack.pop();
    doIndent();
    put(" -->");
  } else if (_nstack.top().t == TAG) {
    if (_nstack.top().hanging) {
      put(" />");
      _nstack.pop();
    } else {
      string tag = _nstack.top().pload;
      _nstack.pop();
      doIndent();
      put("</");
      put(tag);
      put(">");
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
    put("\n");
    for (size_t i = 0; i < _nstack.size() * _indent; i++) put(" ");
  }
}

// _____________________________________________________________________________
void XmlWriter::closeHanging() {
  if (_nstack.empty()) return;

  if (_nstack.top().hanging) {
    put(">");
    _nstack.top().hanging = false;
  } else if (_nstack.top().t == TEXT) {
    _nstack.pop();
  }
}

// _____________________________________________________________________________
void XmlWriter::put(const string& str) {
  if (_gzfile) {
#ifdef ZLIB_FOUND
    gzwrite(_gzfile, str.c_str(), str.size());
#endif
  } else if (_bzfile) {
#ifdef BZLIB_FOUND
    if (_bzbufpos == BUFFER_S || _bzbufpos + str.size() > BUFFER_S) flushBzip();
    memcpy( _bzbuf + _bzbufpos, str.c_str(), str.size());
    _bzbufpos += str.size();
#endif
  } else {
    _out->write(str.c_str(), str.size());
  }
}

// _____________________________________________________________________________
void XmlWriter::flushBzip() {
#ifdef BZLIB_FOUND
    int err = 0;
    BZ2_bzWrite(&err, _bzfile, _bzbuf, _bzbufpos);
    if (err == BZ_IO_ERROR) {
      BZ2_bzWriteClose(&err, _bzfile, 0, 0, 0);
      throw std::runtime_error("Could not write to file.");
    }

    _bzbufpos = 0;
  #endif
}

// _____________________________________________________________________________
void XmlWriter::put(const char c) {
  if (_gzfile) {
#ifdef ZLIB_FOUND
    gzputc(_gzfile, c);

#endif
  } else if (_bzfile) {
#ifdef BZLIB_FOUND
    _bzbuf[_bzbufpos++] = c;
    if (_bzbufpos == BUFFER_S) flushBzip();
#endif
  } else {
    _out->put(c);
  }
}

// _____________________________________________________________________________
void XmlWriter::putEsced(const string& str, char quot) {
  if (!_nstack.empty() && _nstack.top().t == COMMENT) {
    put(str);
    return;
  }

  for (const char& c : str) {
    if (quot == '"' && c == '"')
      put("&quot;");
    else if (quot == '\'' && c == '\'')
      put("&apos;");
    else if (c == '<')
      put("&lt;");
    else if (c == '>')
      put("&gt;");
    else if (c == '&')
      put("&amp;");
    else
      put(c);
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
