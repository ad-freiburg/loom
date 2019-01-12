// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef ZLIB_CONST
#define ZLIB_CONST
#endif

#include <fcntl.h>
#include <netdb.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <algorithm>
#include <memory>
#include <sstream>
#include <thread>
#include <unordered_map>
#ifdef ZLIB_FOUND
#include <zlib.h>
#endif
#include <vector>
#include "Server.h"
#include "util/String.h"

using util::http::Socket;
using util::http::Queue;
using util::http::Req;
using util::http::HttpErr;
using util::http::HttpServer;
using util::http::HeaderState;

// _____________________________________________________________________________
Socket::Socket(int port) {
  int y = 1;
  _sock = socket(PF_INET, SOCK_STREAM, 0);
  if (_sock < 0)
    throw std::runtime_error(std::string("Could not create socket (") +
                             std::strerror(errno) + ")");

  struct sockaddr_in addr;

  addr.sin_family = AF_INET;
  addr.sin_port = htons(port);
  addr.sin_addr.s_addr = INADDR_ANY;
  memset(&(addr.sin_zero), '\0', 8);

  setsockopt(_sock, SOL_SOCKET, SO_REUSEADDR, &y, sizeof(y));

  if (bind(_sock, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) < 0) {
    throw std::runtime_error(std::string("Could not bind to port ") +
                             std::to_string(port) + " (" +
                             std::strerror(errno) + ")");
  }
}

// _____________________________________________________________________________
Socket::~Socket() { close(_sock); }

// _____________________________________________________________________________
int Socket::wait() {
  if (listen(_sock, BLOG) < 0)
    throw std::runtime_error(std::string("Cannot listen to socket (") +
                             std::strerror(errno) + ")");
  sockaddr_in cli_addr;
  socklen_t clilen = sizeof(cli_addr);
  int sock = accept(_sock, reinterpret_cast<sockaddr*>(&cli_addr), &clilen);
  return sock;
}

// _____________________________________________________________________________
void HttpServer::send(int sock, Answer* aw) {
  std::string enc = "identity";
  if (aw->gzip) aw->pl = compress(aw->pl, &enc);

  aw->params["Content-Encoding"] = enc;
  aw->params["Content-Length"] = std::to_string(aw->pl.size());

  std::stringstream ss;
  ss << "HTTP/1.1 " << aw->status << "\r\n";
  for (const auto& kv : aw->params)
    ss << kv.first << ": " << kv.second << "\r\n";
  ss << "\r\n" << aw->pl;
  std::string buff = ss.str();

  size_t writes = 0;

  while (writes != buff.size()) {
    int64_t out = write(sock, buff.c_str() + writes, buff.size() - writes);
    if (out < 0) {
      if (errno == EWOULDBLOCK || errno == EAGAIN || errno == EINTR) continue;
      throw std::runtime_error("Failed to write to socket");
    }
    writes += out;
  }
}

// _____________________________________________________________________________
void HttpServer::handle() {
  int connection = -1;
  while ((connection = _jobs.get()) != -1) {
    Answer answ;

    try {
      Req req = getReq(connection);
      answ = _handler->handle(req, connection);
      answ.gzip = gzipSupport(req);
    } catch (HttpErr err) {
      answ = Answer(err.what(), err.what());
    } catch (...) {
      // catch everything to make sure the server continues running
      answ = Answer(
          "500 Internal Server Error", "500 Internal Server Error");
    }

    send(connection, &answ);
    close(connection);
  }
}

// _____________________________________________________________________________
bool HttpServer::gzipSupport(const Req& req) {
  bool accepts = false;
  // decide according to
  // http://www.w3.org/Protocols/rfc2616/rfc2616-sec14.html
  for (const auto& kv : req.params) {
    if (kv.first == "Accept-Encoding") {
      for (const auto& encoding : split(kv.second, ',')) {
        std::vector<std::string> parts = split(encoding, ';');
        for (size_t i = 0; i < parts.size(); i++) {
          parts[i] = trim(parts[i]);
        }
        if (parts[0] == "*" && ((parts.size() == 1) || parts[1] != "q=0"))
          accepts = true;
        if (parts[0] == "gzip") accepts = true;
        if (parts.size() > 1 && parts[1] == "q=0") accepts = false;
      }
    }
  }
  return accepts;
}

// _____________________________________________________________________________
Req HttpServer::getReq(int connection) {
  char buf[BSIZE + 1];
  size_t rcvd = 0;
  int64_t curRcvd = 0;
  HeaderState state = NONE;
  Req ret;
  char *tmp = 0;
  char *tmp2 = 0;
  char* brk = 0;

  while ((curRcvd = read(connection, buf + rcvd, BSIZE - rcvd))) {
    if (curRcvd < 0) {
      if (errno == EAGAIN || errno == EINTR) continue;
      throw HttpErr("500 Internal Server Error");
    }

    // parse request
    for (int i = 0; i < curRcvd; i++) {
      if (brk) break;
      char* c = buf + rcvd + i;
      switch (state) {
        case NONE:
          state = I_COM;
          tmp = c;
          continue;
        case I_VER:
          if (*c == '\n') {
            *c = 0;
            ret.ver = trim(tmp);
            state = A_KEY;
          }
          continue;
        case I_URL:
          if (*c == ' ') {
            *c = 0, ret.url = trim(tmp);
            tmp = c + 1;
            state = I_VER;
          } else if (*c == '\n') {
            *c = 0, ret.url = trim(tmp);
            state = A_KEY;
          }
          continue;
        case I_COM:
          if (*c == ' ') {
            *c = 0, ret.cmd = trim(tmp);
            tmp = c + 1;
            state = I_URL;
          } else if (*c == '\n') {
            *c = 0, ret.cmd = trim(tmp);
            state = A_KEY;
          }
          continue;
        case A_KEY:
          if (*c == '\r') *c = ' ';
          if (*c == '\n')
            brk = c + 1;
          else if (*c != ' ') {
            state = I_KEY;
            tmp = c;
          }
          continue;
        case I_KEY:
          if (*c == ':') {
            *c = 0;
            state = A_VAL;
          }
          continue;
        case A_VAL:
          if (*c != ' ') {
            state = I_VAL;
            tmp2 = c;
          }
          continue;
        case I_VAL:
          if (*c == '\r') *c = ' ';
          if (*c == '\n') {
            *c = 0;
            ret.params[tmp] = trim(tmp2);
            state = A_KEY;
          }
          continue;
      }
    }

    rcvd += curRcvd;

    // buffer is full
    if (rcvd == BSIZE) throw HttpErr("431 Request Header Fields Too Large");
    if (brk) break;
  }

  // POST payload
  if (ret.cmd == "POST") {
    size_t size = 0;
    if (ret.params.count("Content-Length"))
      size = atoi(ret.params["Content-Length"].c_str());
    if (size) {
      char* postBuf = new char[size + 1];
      postBuf[size] = 0;
      size_t rem = 0;

      // copy existing to new buffer
      if ((int)rcvd > brk - buf) {
        rem = std::min(size, rcvd - (brk - buf));
        memcpy(postBuf, brk, rem);
      }

      rcvd = 0;

      if (rem < size) {
        while ((curRcvd = read(connection, postBuf + rcvd + rem, size - rem))) {
          if (curRcvd == -1 && (errno == EAGAIN || errno == EINTR)) continue;
          if (curRcvd == -1) {
            postBuf[rcvd + 1] = 0;
            break;
          }
          rcvd += curRcvd;
          if (rcvd == size - rem) break;
        }
      }

      ret.payload = postBuf;
      delete[] postBuf;
    }
  }

  return ret;
}

// _____________________________________________________________________________
std::string HttpServer::compress(const std::string& str, std::string* enc) {
#ifdef ZLIB_FOUND
  // do not compress small payloads
  if (str.size() < 500) return str;

  std::string ret;

  // based on http://www.zlib.net/zlib_how.html
  z_stream defStr;
  defStr.zalloc = Z_NULL;
  defStr.zfree = Z_NULL;
  defStr.opaque = Z_NULL;
  defStr.avail_in = 0;
  defStr.next_in = Z_NULL;

  // fail silently with no compression at all
  if (deflateInit2(&defStr, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8,
                   Z_DEFAULT_STRATEGY) != Z_OK)
    return str;

  defStr.next_in = reinterpret_cast<z_const Bytef*>(str.c_str());
  defStr.avail_in = static_cast<unsigned int>(str.size());

  size_t cSize = 0;
  do {
    if (ret.size() < (cSize + BSIZE_C)) ret.resize(cSize + BSIZE_C);
    defStr.avail_out = BSIZE_C;
    defStr.next_out = reinterpret_cast<Bytef*>(&ret[0] + cSize);
    deflate(&defStr, Z_FINISH);
    cSize += BSIZE_C - defStr.avail_out;
  } while (defStr.avail_out == 0);

  deflateEnd(&defStr);
  ret.resize(cSize);

  if (ret.size() > str.size()) return str;
  *enc = "gzip";
  return ret;
#else
  return str;
#endif
}

// _____________________________________________________________________________
void HttpServer::run() {
  Socket socket(_port);

  std::vector<std::thread> thrds(_threads);
  for (auto& thr : thrds) thr = std::thread(&HttpServer::handle, this);

  while (1) _jobs.add(socket.wait());
}

// _____________________________________________________________________________
void Queue::add(int c) {
  if (c < 0) return;
  {
    std::unique_lock<std::mutex> lock(_mut);
    _jobs.push(c);
  }
  _hasNew.notify_one();
}

// _____________________________________________________________________________
int Queue::get() {
  std::unique_lock<std::mutex> lock(_mut);
  while (_jobs.empty()) _hasNew.wait(lock);
  int next = _jobs.front();
  _jobs.pop();
  return next;
}
