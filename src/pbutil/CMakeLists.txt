file(GLOB_RECURSE pbutil_SRC *.cpp)
list(REMOVE_ITEM pbutil_SRC TestMain.cpp)
add_library(pbutil ${pbutil_SRC})

add_subdirectory(tests)
