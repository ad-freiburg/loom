file(GLOB_RECURSE shared_SRC *.cpp)

include_directories(
	${LOOM_INCLUDE_DIR}
	SYSTEM ${GUROBI_INCLUDE_DIR}
	SYSTEM ${GLPK_INCLUDE_DIR}
	SYSTEM ${COIN_INCLUDE_DIR}
)

add_subdirectory(tests)

add_library(shared_dep ${shared_SRC})

target_link_libraries(shared_dep pb_util_graph pb_util_geo pb_util_json pb_util_xml pb_util)
