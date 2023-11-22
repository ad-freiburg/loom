file(GLOB_RECURSE shared_SRC *.cpp)

include_directories(
	${LOOM_INCLUDE_DIR}
	SYSTEM ${GUROBI_INCLUDE_DIR}
	SYSTEM ${GLPK_INCLUDE_DIR}
	SYSTEM ${COIN_INCLUDE_DIR}
)

add_subdirectory(tests)

add_library(shared_dep ${shared_SRC})
