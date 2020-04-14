file(GLOB_RECURSE shared_SRC *.cpp)

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
	${GUROBI_INCLUDE_DIR}
	${GLPK_INCLUDE_DIR}
	${COIN_INCLUDE_DIR}
)

add_subdirectory(tests)

add_library(shared_dep ${shared_SRC})
