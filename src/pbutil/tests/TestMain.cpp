// Copyright 2016
// Author: Patrick Brosi
//

#include <string>
#include "./lest.h"
#include "./../Nullable.h"
#include "./../String.h"

using lest::approx;

// define LEST cases
const lest::test specification[] = {

// ___________________________________________________________________________
CASE("url decode") {
  EXPECT("zürich" == pbutil::urlDecode("z%C3%BCrich"));
  EXPECT("!@$%^*()" == pbutil::urlDecode("!%40%24%25%5E*()"));
  EXPECT("Løkken" == pbutil::urlDecode("L%C3%B8kken"));
  EXPECT("á é" == pbutil::urlDecode("%C3%A1%20%C3%A9"));
  EXPECT("á é" == pbutil::urlDecode("%C3%A1+%C3%A9"));
},

// ___________________________________________________________________________
CASE("json escape") {
  EXPECT("Hello\\\\Goodbye!" == pbutil::jsonStringEscape("Hello\\Goodbye!"));
  EXPECT("\\\"Hello\\\"" == pbutil::jsonStringEscape("\"Hello\""));
},

// ___________________________________________________________________________
CASE("toString") {
  EXPECT(pbutil::toString(34) == "34");
  EXPECT(pbutil::toString("34") == "34");
},

// ___________________________________________________________________________
CASE("replace") {
  std::string a("lorem ipsum ipsum lorem");

  EXPECT(pbutil::replace(a, "ips", "aa"));
  EXPECT(a == "lorem aaum ipsum lorem");

  EXPECT(!pbutil::replace(a, "blablabla", ""));
  EXPECT(a == "lorem aaum ipsum lorem");

  EXPECT(pbutil::replace(a, "m", ""));
  EXPECT(a == "lore aaum ipsum lorem");

  EXPECT(!pbutil::replace(a, "", ""));
  EXPECT(a == "lore aaum ipsum lorem");

  std::string b("lorem ipsum ipsum lorem");
  EXPECT(pbutil::replaceAll(b, "ips", "aa"));
  EXPECT(b == "lorem aaum aaum lorem");

  EXPECT(pbutil::replaceAll(b, "m", ""));
  EXPECT(b == "lore aau aau lore");

  EXPECT(pbutil::replaceAll(b, "a", "aa"));
  EXPECT(b == "lore aaaau aaaau lore");

  EXPECT(pbutil::replaceAll(b, "e", "e"));
  EXPECT(b == "lore aaaau aaaau lore");

  EXPECT(pbutil::replaceAll(b, "e", "ee"));
  EXPECT(b == "loree aaaau aaaau loree");

  EXPECT(!pbutil::replaceAll(b, "", "ee"));
  EXPECT(b == "loree aaaau aaaau loree");
},

// ___________________________________________________________________________
CASE("nullable") {
  {
    pbutil::Nullable<std::string> nullable;
    EXPECT(nullable.isNull());
  }

  {
    pbutil::Nullable<std::string> nullable(0);
    EXPECT(nullable.isNull());
  }

  {
    std::string str = "aa";
    pbutil::Nullable<std::string> nullable(&str);
    EXPECT(!nullable.isNull());

    EXPECT(nullable == "aa");
    EXPECT(!(nullable == "aaa"));
    EXPECT(!(nullable != "aa"));
    EXPECT(nullable == "aa");

    EXPECT(nullable.get() == "aa");
    EXPECT(std::string(nullable) == "aa");
  }

  {
    int a = 23;
    pbutil::Nullable<int> nullable(a);
    pbutil::Nullable<int> nullable2(24);
    EXPECT(!nullable.isNull());

    EXPECT(nullable == 23);
    EXPECT(nullable >= 23);
    EXPECT(nullable <= 23);
    EXPECT(nullable < 24);
    EXPECT(nullable < 24);
    EXPECT(!(nullable < 22));
    EXPECT(nullable != nullable2);
    EXPECT(nullable < nullable2);
    EXPECT(nullable2 > nullable);

    pbutil::Nullable<int> nullable3(nullable);
    EXPECT(nullable == nullable3);

    nullable3 = nullable2;
    EXPECT(nullable2 == nullable3);
    EXPECT(nullable3 == 24);
    EXPECT(nullable2 == 24);
    EXPECT(nullable2 == nullable2.get());
    EXPECT(int(nullable2) == nullable2.get());
    EXPECT(!nullable3.isNull());
    EXPECT(!nullable2.isNull());

    pbutil::Nullable<int> voidnull;
    EXPECT(voidnull.isNull());

    EXPECT_THROWS(nullable == voidnull);
  }
}

};

// _____________________________________________________________________________
int main(int argc, char** argv) {
  return(lest::run(specification, argc, argv));
}
