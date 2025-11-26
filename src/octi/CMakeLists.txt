file(GLOB_RECURSE octi_SRC *.cpp)

set(octi_main OctiMain.cpp)

list(REMOVE_ITEM octi_SRC ${octi_main})
list(REMOVE_ITEM octi_SRC TestMain.cpp)

include_directories(
	SYSTEM ${GUROBI_INCLUDE_DIR}
	SYSTEM ${GLPK_INCLUDE_DIR}
	SYSTEM ${COIN_INCLUDE_DIR}
)

add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(octi ${octi_main})
add_library(octi_dep ${octi_SRC})

target_link_libraries(octi octi_dep shared_dep pb_util_geo pb_util_graph pb_util_json pb_util dot_dep ${GLPK_LIBRARY} ${GUROBI_LIBRARY} ${COIN_LIBRARIES} -lpthread)
