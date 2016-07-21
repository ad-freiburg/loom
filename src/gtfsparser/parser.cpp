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
  std::ifstream fs;
  filesystem::path gtfsPath(path);
  filesystem::path curFile;

  try {
    curFile = gtfsPath / "agency.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) fileNotFound(curFile);
    parseAgency(targetFeed, &fs);
    fs.close();

    curFile = gtfsPath / "stops.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) fileNotFound(curFile);
    parseStops(targetFeed, &fs);
    fs.close();

    curFile = gtfsPath / "routes.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) fileNotFound(curFile);
    parseRoutes(targetFeed, &fs);
    fs.close();

    curFile = gtfsPath / "calendar.txt";
    fs.open(curFile.c_str());
    if (fs.good()) {
      parseCalendar(targetFeed, &fs);
      fs.close();
    }

    curFile = gtfsPath / "calendar_dates.txt";
    fs.open(curFile.c_str());
    if (fs.good()) {
      parseCalendarDates(targetFeed, &fs);
      fs.close();
    }

    curFile = gtfsPath / "trips.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) fileNotFound(curFile);
    parseTrips(targetFeed, &fs);
    fs.close();

    curFile = gtfsPath / "stop_times.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) fileNotFound(curFile);
    parseStopTimes(targetFeed, &fs);
    fs.close();

 } catch (const CsvParserException& e) {
    throw ParserException(e.getMsg(), e.getFieldName(), e.getLine(),
      curFile.c_str());
  } catch (const ParserException& e) {
    // augment with file namoe
    ParserException fe = e;
    fe.setFileName(curFile.c_str());
    throw fe;
  }

  return true;
}

// ____________________________________________________________________________
void Parser::parseAgency(gtfs::Feed* targetFeed, std::istream* s) const {
  CsvParser csvp(s);
  Agency* a = 0;

  while (csvp.readNextLine()) {
    a = new Agency(
      getString(csvp, "agency_id", ""),
      getString(csvp, "agency_name"),
      getString(csvp, "agency_url"),
      getString(csvp, "agency_timezone"),
      getString(csvp, "agency_lang", ""),
      getString(csvp, "agency_phone", ""),
      getString(csvp, "agency_fare_url", ""),
      getString(csvp, "agency_email", "")
    );

    if (!targetFeed->addAgency(a)) {
      std::stringstream msg;
      msg << "'agency_id' must be dataset unique. Collision with id '"
        << a->getId() << "')";
      throw ParserException(msg.str(), "agency_id",
        csvp.getCurLine());
    }
  }

  if (!a) {
    throw ParserException("the feed has no agency defined."
      " This is a required field.", "", 1);
  }
}

// ____________________________________________________________________________
void Parser::parseShapes(gtfs::Feed* targetFeed, std::istream* s) const {
  CsvParser csvp(s);

  while (csvp.readNextLine()) {

  }
}

// ____________________________________________________________________________
void Parser::parseStops(gtfs::Feed* targetFeed, std::istream* s) const {
  CsvParser csvp(s);

  while (csvp.readNextLine()) {
    Stop* s = new Stop(
      getString(csvp, "stop_id", ""),
      getString(csvp, "stop_code", ""),
      getString(csvp, "stop_name"),
      getString(csvp, "stop_desc", ""),
      getDouble(csvp, "stop_lat"),
      getDouble(csvp, "stop_lon"),
      getString(csvp, "zone_id", ""),
      getString(csvp, "stop_url", ""),
      static_cast<Stop::LOCATION_TYPE>(
        getRangeInteger(csvp, "location_type", 0, 1, 0)
      ),
      getString(csvp, "parent_station", ""),
      getString(csvp, "stop_timezone", ""),
      static_cast<Stop::WHEELCHAIR_BOARDING>(
        getRangeInteger(csvp, "wheelchair_boarding", 0, 2, 0)
      )
    );

    if (!targetFeed->addStop(s)) {
      std::stringstream msg;
      msg << "'stop_id' must be dataset unique. Collision with id '"
        << s->getId() << "')";
      throw ParserException(msg.str(), "stop_id",
        csvp.getCurLine());
    }
  }
}

// ____________________________________________________________________________
void Parser::parseRoutes(gtfs::Feed* targetFeed, std::istream* s) const {
  CsvParser csvp(s);

  while (csvp.readNextLine()) {
    std::string agencyId = getString(csvp, "agency_id", "");
    Agency* routeAgency = 0;

    if (!agencyId.empty()) {
      routeAgency = targetFeed->getAgencyById(agencyId);
      if (!routeAgency) {
        std::stringstream msg;
        msg << "no agency with id '" << agencyId << "' defined, cannot "
          << "reference here.";
        throw ParserException(msg.str(), "agency_id", csvp.getCurLine());
      }
    }

    Route* r = new Route(
      getString(csvp, "route_id"),
      routeAgency,
      getString(csvp, "route_short_name", ""),
      getString(csvp, "route_long_name"),
      getString(csvp, "route_desc", ""),
      static_cast<Route::TYPE> (
        getRangeInteger(csvp, "route_type", 0, 7)
      ),
      getString(csvp, "route_url", ""),
      getColorFromHexString(csvp, "route_color", "FFFFFF"),
      getColorFromHexString(csvp, "route_text_color", "000000")
    );

    if (!targetFeed->addRoute(r)) {
      std::stringstream msg;
      msg << "'route_id' must be dataset unique. Collision with id '"
        << r->getId() << "')";
      throw ParserException(msg.str(), "route_id",
        csvp.getCurLine());
    }
  }
}

// ____________________________________________________________________________
void Parser::parseCalendar(gtfs::Feed* targetFeed, std::istream* s) const {
  CsvParser csvp(s);

  while (csvp.readNextLine()) {
    std::string serviceId = getString(csvp, "service_id");

    Service* s = new Service(
      serviceId,
      (0 << getRangeInteger(csvp,  "monday", 0, 1)) |
      (0 << getRangeInteger(csvp,  "tuesday", 0, 1) * 2) |
      (0 << getRangeInteger(csvp,  "wednesday", 0, 1) * 3) |
      (0 << getRangeInteger(csvp,  "thursday", 0, 1) * 4) |
      (0 << getRangeInteger(csvp,  "friday", 0, 1) * 5) |
      (0 << getRangeInteger(csvp,  "saturday", 0, 1) * 6) |
      (0 << getRangeInteger(csvp,  "sunday", 0, 1) * 7),
      getServiceDate(csvp, "start_date"),
      getServiceDate(csvp, "end_date")
    );

    if (!targetFeed->addService(s)) {
      std::stringstream msg;
      msg << "'service_id' must be unique in calendars.txt. Collision with id '"
        << s->getId() << "')";
      throw ParserException(msg.str(), "service_id",
        csvp.getCurLine());
    }
  }

}

// ____________________________________________________________________________
void Parser::parseCalendarDates(gtfs::Feed* targetFeed, std::istream* s) const {
  CsvParser csvp(s);

  while (csvp.readNextLine()) {
    std::string serviceId = getString(csvp, "service_id");
    ServiceDate d = getServiceDate(csvp, "date");
    Service::EXCEPTION_TYPE t = static_cast<Service::EXCEPTION_TYPE>(
      getRangeInteger(csvp, "exception_type", 1, 2)
    );

    Service* e = targetFeed->getServiceById(serviceId);

    if (!e) {
      e = new Service(serviceId);
      targetFeed->addService(e);
    }

    e->addException(d, t);
  }
}

// ____________________________________________________________________________
void Parser::parseTrips(gtfs::Feed* targetFeed, std::istream* s) const {
  CsvParser csvp(s);

  while (csvp.readNextLine()) {
    std::string routeId = getString(csvp, "route_id");
    Route* tripRoute = 0;

    tripRoute = targetFeed->getRouteById(routeId);
    if (!tripRoute) {
      std::stringstream msg;
      msg << "no route with id '" << routeId << "' defined, cannot "
        << "reference here.";
      throw ParserException(msg.str(), "route_id", csvp.getCurLine());
    }

    std::string shapeId = getString(csvp, "shape_id", "");
    Shape* tripShape = 0;

    if (false && !shapeId.empty()) {
      tripShape = targetFeed->getShapeById(shapeId);
      if (!tripShape) {
        std::stringstream msg;
        msg << "no shape with id '" << shapeId << "' defined, cannot "
          << "reference here.";
        throw ParserException(msg.str(), "shape_id", csvp.getCurLine());
      }
    }

    std::string serviceId = getString(csvp, "service_id");
    Service* tripService = targetFeed->getServiceById(serviceId);
    if (!tripService) {
      std::stringstream msg;
      msg << "no service with id '" << serviceId << "' defined, cannot "
        << "reference here.";
      throw ParserException(msg.str(), "service_id", csvp.getCurLine());
    }

    Trip* t = new Trip(
      getString(csvp, "trip_id"),
      tripRoute,
      tripService,
      getString(csvp, "trip_headsign", ""),
      getString(csvp, "trip_short_name", ""),
      static_cast<Trip::DIRECTION>(
        getRangeInteger(csvp, "direction_id", 0, 1, 2)
      ),
      getString(csvp, "block_id", ""),
      tripShape,
      static_cast<Trip::WC_BIKE_ACCESSIBLE>(
        getRangeInteger(csvp, "wheelchair_accessible", 0, 2, 0)
      ),
      static_cast<Trip::WC_BIKE_ACCESSIBLE>(
        getRangeInteger(csvp, "bikes_allowed", 0, 2, 0)
      )
    );

    if (!targetFeed->addTrip(t)) {
      std::stringstream msg;
      msg << "'trip_id' must be dataset unique. Collision with id '"
        << t->getId() << "')";
      throw ParserException(msg.str(), "trip_id",
        csvp.getCurLine());
    }
  }
}

// ____________________________________________________________________________
void Parser::parseStopTimes(gtfs::Feed* targetFeed, std::istream* s) const {
  CsvParser csvp(s);

  while (csvp.readNextLine()) {
    Stop* stop = 0;
    Trip* trip = 0;

    const std::string& stopId = getString(csvp, "stop_id");
    const std::string& tripId = getString(csvp, "trip_id");

    stop = targetFeed->getStopById(stopId);
    trip = targetFeed->getTripById(tripId);

    if (!stop) {
      std::stringstream msg;
      msg << "no stop with id '" << stopId << "' defined in stops.txt, cannot "
        << "reference here.";
      throw ParserException(msg.str(), "stop_id", csvp.getCurLine());
    }

    if (!trip) {
      std::stringstream msg;
      msg << "no trip with id '" << tripId << "' defined in trips.txt, cannot "
        << "reference here.";
      throw ParserException(msg.str(), "trip_id", csvp.getCurLine());
    }

    std::string rawDist =  getString(csvp, "shape_dist_traveled", "");

    double dist = -1;  // using -1 as a NULL value here

    if (!rawDist.empty()) {
      dist = getDouble(csvp, "shape_dist_traveled");
      if (dist < -0.01) { // TODO: better double comp
         throw ParserException("negative values not supported for distances"
           " (value was: " + std::to_string(dist),
          "trip_id",
         csvp.getCurLine());
      }
    }

    StopTime st(
      getTime(csvp, "arrival_time"),
      getTime(csvp, "departure_time"),
      stop,
      getRangeInteger(csvp, "stop_sequence", 0, UINT16_MAX),
      getString(csvp, "stop_headsign", ""),
      static_cast<StopTime::PU_DO_TYPE>(
        getRangeInteger(csvp, "drop_off_type", 0, 3, 0)
      ),
      static_cast<StopTime::PU_DO_TYPE>(
        getRangeInteger(csvp, "pick_up_type", 0, 3, 0)
      ),
      dist,
      getRangeInteger(csvp, "timepoint", 0, 1, 1)
    );

    if (st.getArrivalTime() > st.getDepartureTime()) {
     throw ParserException("arrival time '" + st.getArrivalTime().toString() +
      "' is later than departure time '" + st.getDepartureTime().toString() +
      "'. You cannot depart earlier than you arrive.",
      "departure_time",
      csvp.getCurLine());
    }

    if (!trip->addStopTime(st)) {
     throw ParserException("stop_sequence collision, stop_sequence has "
      "to be increasing for a single trip.",
      "stop_sequence",
      csvp.getCurLine());
    }

  }
}

// ___________________________________________________________________________
void Parser::fileNotFound(filesystem::path file) const {
  throw ParserException("File not found",
     "", -1, std::string(file.c_str()));
}

// ___________________________________________________________________________
std::string Parser::getString(const CsvParser& csv,
  const std::string& field) const {
  return csv.getTString(field.c_str());
}

// ___________________________________________________________________________
std::string Parser::getString(const CsvParser& csv,
  const std::string& fld, const std::string& def) const {
  std::string ret = def;

  if (csv.hasItem(fld.c_str()) && !csv.fieldIsEmpty(fld.c_str())) {
    ret = csv.getTString(fld.c_str());
  }

  return ret;
}

// ___________________________________________________________________________
double Parser::getDouble(const CsvParser& csv,
  const std::string& field) const {
  return csv.getDouble(field.c_str());
}

// ___________________________________________________________________________
double Parser::getDouble(const CsvParser& csv,
  const std::string& field, double ret) const {
  if (csv.hasItem(field.c_str()) && !csv.fieldIsEmpty(field.c_str())) {
    ret = csv.getDouble(field.c_str());
  }

  return ret;
}

// ___________________________________________________________________________
int64_t Parser::getRangeInteger(const CsvParser& csv,
  const std::string& field, int64_t minv, int64_t maxv) const {
  int64_t ret = csv.getLong(field.c_str());

  if (ret < minv || ret > maxv) {
    std::stringstream msg;
    msg << "expected integer in range [" << minv << "," << maxv << "]";
    throw ParserException(msg.str(), field, csv.getCurLine());
  }

  return ret;
}

// ___________________________________________________________________________
int64_t Parser::getRangeInteger(const CsvParser& csv,
  const std::string& field, int64_t minv, int64_t maxv, int64_t def) const {
  int64_t ret;

  if (csv.hasItem(field.c_str()) && !csv.fieldIsEmpty(field.c_str())) {
    ret = csv.getLong(field.c_str());

    if (ret < minv || ret > maxv) {
      std::stringstream msg;
      msg << "expected integer in range [" << minv << "," << maxv << "]";
      throw ParserException(msg.str(), field, csv.getCurLine());
    }

    return ret;
  }

  return def;
}

// ___________________________________________________________________________
uint32_t Parser::getColorFromHexString(
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
   msg << "expected a 6-character hexadecimal color string, found '"
      << color_string << "' instead. (Error while parsing was: "
      << e.what() << ")";
   throw ParserException(msg.str(), field, csv.getCurLine());
  }

  if (color_string.size() != 6 || chars_processed != 8) {
    std::stringstream msg;
    msg << "expected a 6-character hexadecimal color string, found '"
      << color_string << "' instead.";
    throw ParserException(msg.str(), field, csv.getCurLine());
  }

  return ret;
}

// ____________________________________________________________________________
ServiceDate Parser::getServiceDate(
  const CsvParser& csv, const std::string& field) const {
  size_t p;
  std::string val(csv.getTString(field.c_str()));

  try {
    int32_t yyyymmdd = std::stoul(val, &p, 10);
    if (p != val.length() || yyyymmdd > 99999999) {
      std::stringstream msg;
      msg << "expected a date in the YYYYMMDD format, found '"
        << val << "' instead.";
      throw ParserException(msg.str(), field, csv.getCurLine());
    }
    return ServiceDate(yyyymmdd);
  } catch (const std::out_of_range& e) {
     std::stringstream msg;
     msg << "expected a date in the YYYYMMDD format, found '"
        << val << "' instead. (Integer out of range).";
     throw ParserException(msg.str(), field, csv.getCurLine());
  } catch (const std::invalid_argument& e) {
     std::stringstream msg;
     msg << "expected a date in the YYYYMMDD format, found '"
        << val << "' instead.";
     throw ParserException(msg.str(), field, csv.getCurLine());
  }
}

// ____________________________________________________________________________
gtfs::Time Parser::getTime(
  const CsvParser& csv, const std::string& field) const {
  size_t p;
  std::string val(csv.getTString(field.c_str()));

  try {
    uint64_t h = std::stoul(val, &p, 10);
    if (h > 255) throw std::out_of_range(
      "only hour-values up to 255 are "
      "supported. (read " + std::to_string(h) + ")");
    val.erase(0, p+1);

    uint64_t m = std::stoul(val, &p, 10);
    if (p == 1) throw std::invalid_argument(
      "one-digit minute values are not allowed.");
    // allow values of 60, although standard forbids it
    if (m > 60) throw std::out_of_range(
      "only minute-values up to 60 are "
      "allowed. (read " + std::to_string(m) + ")");
    val.erase(0, p+1);

    uint64_t s = std::stoul(val, &p, 10);
    if (p == 0) s = 0; // support HH:MM format (although standard forbids it)
    if (p == 1) throw std::invalid_argument(
      "one-digit second values are not allowed.");
     // allow values of 60, although standard forbids it
    if (s > 60) throw std::out_of_range(
      "only second-values up to 60 are "
      "allowed. (read " + std::to_string(s) + ")");

    return gtfs::Time(h, m % 60, s % 60);
  } catch (const std::exception& e) {
    std::stringstream msg;
    msg << "expected a time in HH:MM:SS (or H:MM:SS) format, found '"
      << val << "' instead. (" <<  e.what() << ")";
    throw ParserException(msg.str(), field, csv.getCurLine());
  }
}

