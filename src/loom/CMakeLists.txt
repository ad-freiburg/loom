file(GLOB_RECURSE loom_SRC *.cpp)

set(loom_main TransitMapMain.cpp)

list(REMOVE_ITEM loom_SRC ${loom_main})
list(REMOVE_ITEM loom_SRC TestMain.cpp)

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

add_executable(loom ${loom_main})
add_library(loom_dep ${loom_SRC})

target_link_libraries(loom loom_dep shared_dep dot_dep util ${Boost_LIBRARIES} ${GLPK_LIBRARY} ${GUROBI_LIBRARY} ${COIN_LIBRARIES} -lpthread)
