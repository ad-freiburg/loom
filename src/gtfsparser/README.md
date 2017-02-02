# Usage

```
#include "gtfsparser/Parser.h"

[...]

gtfsparser::Parser parser;
gtfsparser::gtfs::Feed feed;

parser.parse(&feed, cfg.inputFeedPath);
```
