// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <string>
#include <fstream>
#include <boost/filesystem.hpp>
#include "parser.h"
#include "csvparser.h"
#include "gtfs/agency.h"
#include "gtfs/stop.h"
#include "gtfs/route.h"

using namespace gtfsparser;
using namespace gtfs;
using namespace boost;

// ____________________________________________________________________________
bool Parser::parse(gtfs::Feed* targetFeed, std::string path) const {
  parseAgency(targetFeed, path);
  parseStops(targetFeed, path);
  parseRoutes(targetFeed, path);

  return true;
}

// ____________________________________________________________________________
void Parser::parseAgency(gtfs::Feed* targetFeed, std::string path) const {
  std::ifstream file_stream;
  filesystem::path gtfs_path(path);
  filesystem::path file_path("agency.txt");
  filesystem::path full_path = gtfs_path / file_path;

  file_stream.open(full_path.c_str());

  if (!file_stream.good()) fileNotFound(full_path);

  CsvParser csvp(&file_stream);
  Agency* a = 0;

  while (csvp.readNextLine()) {
    a = new Agency(
      getString(full_path, csvp, "agency_id", ""),
      getString(full_path, csvp, "agency_name"),
      getString(full_path, csvp, "agency_url"),
      getString(full_path, csvp, "agency_timezone"),
      getString(full_path, csvp, "agency_lang", ""),
      getString(full_path, csvp, "agency_phone", ""),
      getString(full_path, csvp, "agency_fare_url", ""),
      getString(full_path, csvp, "agency_email", "")
    );

    if (!targetFeed->addAgency(a)) {
      std::stringstream msg;
      msg << "'agency_id' must be dataset unique. Collision with id '"
        << a->getId() << "')";
      throw ParserException(msg.str(), "agency_id",
        csvp.getCurLine(), std::string(full_path.c_str()));
    }
  }

  if (!a) {
    throw ParserException("The feed has no agency defined."
      " This is a required field.", "", 1, std::string(full_path.c_str()));
  }
}

// ____________________________________________________________________________
void Parser::parseStops(gtfs::Feed* targetFeed, std::string path) const {
  std::ifstream file_stream;
  filesystem::path gtfs_path(path);
  filesystem::path file_path("stops.txt");
  filesystem::path full_path = gtfs_path / file_path;

  file_stream.open(full_path.c_str());

  if (!file_stream.good()) fileNotFound(full_path);

  CsvParser csvp(&file_stream);

  while (csvp.readNextLine()) {
    Stop* s = new Stop(
      getString(full_path, csvp, "stop_id", ""),
      getString(full_path, csvp, "stop_code"),
      getString(full_path, csvp, "stop_name"),
      getString(full_path, csvp, "stop_desc"),
      getDouble(full_path, csvp, "stop_lat"),
      getDouble(full_path, csvp, "stop_lon"),
      getString(full_path, csvp, "zone_id", ""),
      getString(full_path, csvp, "stop_url", ""),
      static_cast<Stop::LOCATION_TYPE>(
        getRangeInteger(full_path, csvp, "location_type", 0, 1, 0)
      ),
      getString(full_path, csvp, "parent_station", ""),
      getString(full_path, csvp, "stop_timezone", ""),
      static_cast<Stop::WHEELCHAIR_BOARDING>(
        getRangeInteger(full_path, csvp, "wheelchair_boarding", 0, 2, 0)
      )
    );

    if (!targetFeed->addStop(s)) {
      std::stringstream msg;
      msg << "'stop_id' must be dataset unique. Collision with id '"
        << s->getId() << "')";
      throw ParserException(msg.str(), "stop_id",
        csvp.getCurLine(), std::string(full_path.c_str()));
    }
  }
}

// ____________________________________________________________________________
void Parser::parseRoutes(gtfs::Feed* targetFeed, std::string path) const {
  std::ifstream file_stream;
  filesystem::path gtfs_path(path);
  filesystem::path file_path("routes.txt");
  filesystem::path full_path = gtfs_path / file_path;

  file_stream.open(full_path.c_str());

  if (!file_stream.good()) fileNotFound(full_path);

  CsvParser csvp(&file_stream);

  while (csvp.readNextLine()) {
    std::string agencyId = getString(full_path, csvp, "agency_id", "");
    Agency* routeAgency = 0;

    if (!agencyId.empty()) {
      routeAgency = targetFeed->getAgencyById(agencyId);
      if (!routeAgency) {
        std::stringstream msg;
        msg << "No agency with id '" << agencyId << "' defined, cannot "
          << "reference here.";
        throw ParserException(msg.str(), "agency_id", csvp.getCurLine(),
          full_path.c_str());
      }
    }

    Route* r = new Route(
      getString(full_path, csvp, "route_id"),
      routeAgency,
      getString(full_path, csvp, "route_short_name", ""),
      getString(full_path, csvp, "route_long_name"),
      getString(full_path, csvp, "route_desc", ""),
      static_cast<Route::TYPE> (
        getRangeInteger(full_path, csvp, "route_type", 0, 7)
      ),
      getString(full_path, csvp, "route_url", ""),
      getColorFromHexString(full_path, csvp, "route_color", "FFFFFF"),
      getColorFromHexString(full_path, csvp, "route_text_color", "000000")
    );

    if (!targetFeed->addRoute(r)) {
      std::stringstream msg;
      msg << "'route_id' must be dataset unique. Collision with id '"
        << r->getId() << "')";
      throw ParserException(msg.str(), "route_id",
        csvp.getCurLine(), std::string(full_path.c_str()));
    }
  }
}

// ___________________________________________________________________________
void Parser::fileNotFound(filesystem::path file) const {
  throw ParserException("File not found",
     "", 1, std::string(file.c_str()));
}

// ___________________________________________________________________________
std::string Parser::getString(const filesystem::path& fn, const CsvParser& csv,
  const std::string& field) const {
  return csv.getTString(field.c_str());
}

// ___________________________________________________________________________
std::string Parser::getString(const filesystem::path& fn, const CsvParser& csv,
  const std::string& fld, const std::string& def) const {
  std::string ret = def;

  if (csv.hasItem(fld.c_str()) && !csv.fieldIsEmpty(fld.c_str())) {
    ret = csv.getTString(fld.c_str());
  }

  return ret;
}


// ___________________________________________________________________________
double Parser::getDouble(const filesystem::path& fn, const CsvParser& csv,
  const std::string& field) const {
  return csv.getDouble(field.c_str());
}

// ___________________________________________________________________________
double Parser::getDouble(const filesystem::path& fn, const CsvParser& csv,
  const std::string& field, double ret) const {
  if (csv.hasItem(field.c_str()) && !csv.fieldIsEmpty(field.c_str())) {
    ret = csv.getDouble(field.c_str());
  }

  return ret;
}


// ___________________________________________________________________________
int64_t Parser::getRangeInteger(const filesystem::path& fn, const CsvParser& csv,
  const std::string& field, int64_t minv, int64_t maxv) const {
  int64_t ret = csv.getLong(field.c_str());

  if (ret < minv || ret > maxv) {
    std::stringstream msg;
    msg << "Expected integer in range [" << minv << "," << maxv << "]";
    throw ParserException(msg.str(), field, csv.getCurLine(), fn.c_str());
  }

  return ret;
}

// ___________________________________________________________________________
int64_t Parser::getRangeInteger(const filesystem::path& fn, const CsvParser& csv,
  const std::string& field, int64_t minv, int64_t maxv, int64_t def) const {
  int64_t ret = def;

  if (csv.hasItem(field.c_str()) && !csv.fieldIsEmpty(field.c_str())) {
    ret = csv.getLong(field.c_str());
  }

  if (ret < minv || ret > maxv) {
    std::stringstream msg;
    msg << "Expected integer in range [" << minv << "," << maxv << "]";
    throw ParserException(msg.str(), field, csv.getCurLine(), fn.c_str());
  }

  return ret;
}

// ___________________________________________________________________________
uint32_t Parser::getColorFromHexString(const filesystem::path& fn,
  const CsvParser& csv, const std::string& field, const std::string& def)
  const {
  std::string color_string;

  if (csv.hasItem(field.c_str())) {
    color_string= csv.getTString(field.c_str());
  }

  if (color_string.empty()) color_string = def;

  size_t chars_processed = 0;
  uint32_t ret = 0;

  try {
    ret = std::stoul("0x"+color_string, &chars_processed, 16);
  } catch (std::exception e) {
   std::stringstream msg;
   msg << "Expected a 6-character hexadecimal color string, found '"
      << color_string << "' instead. (Error while parsing was: "
      << e.what() << ")";
   throw ParserException(msg.str(), field, csv.getCurLine(), fn.c_str());
  }

  if (color_string.size() != 6 || chars_processed != 8) {
    std::stringstream msg;
    msg << "Expected a 6-character hexadecimal color string, found '"
      << color_string << "' instead.";
    throw ParserException(msg.str(), field, csv.getCurLine(), fn.c_str());
  }

  return ret;
}


