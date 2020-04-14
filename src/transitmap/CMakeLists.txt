file(GLOB_RECURSE transitmap_SRC *.cpp)

set(transitmap_main TransitMapMain.cpp)

list(REMOVE_ITEM transitmap_SRC ${transitmap_main})
list(REMOVE_ITEM transitmap_SRC TestMain.cpp)

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
	${GUROBI_INCLUDE_DIR}
	${GLPK_INCLUDE_DIR}
	${COIN_INCLUDE_DIR}
)

add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(transitmap ${transitmap_main})
add_library(transitmap_dep ${transitmap_SRC})

target_link_libraries(transitmap transitmap_dep shared_dep dot_dep util ${Boost_LIBRARIES} ${GLPK_LIBRARY} ${GUROBI_LIBRARY} ${COIN_LIBRARIES} -lproj -lpthread)
