file(GLOB_RECURSE util_SRC *.cpp)
list(REMOVE_ITEM util_SRC TestMain.cpp)
add_library(util ${util_SRC})

add_subdirectory(tests)
